//! Adjudication records for the corpus's two REFUSAL classes: shapes the
//! recommendations prohibit, and alleles whose members claim one territory.
//!
//! # What these are, and what they are not
//!
//! `spec_conformance_axis` pins three counters whose meaning is "ferro did not
//! refuse something": `prohibited_absolute_accepted` (32),
//! `prohibited_conditional_accepted` (40) and `conflicts_accepted` (72), plus a
//! fourth, `prohibition_violating_outputs` (32), which counts outputs that
//! *violate* a prohibition. Those are totals. This file decomposes them by
//! clause and by geometry, states what ferro **should** do per clause with the
//! clause quoted, and pins both halves.
//!
//! Per `CLAUDE.md`, a test that merely pins today's output is a change detector,
//! not an adjudication record. Every expectation below therefore carries the
//! `file.md:line` it is argued from, and every one states whether it is
//! **adjudicated-correct** (ferro already does the right thing — the assertion
//! is a guard) or a **pinned defect** (ferro does not — the assertion records
//! the target).
//!
//! # The single most important correction this file makes
//!
//! **Three of those four counters are LENIENT-MODE figures, and nothing said
//! so.** The axis census normalizes with `NormalizeConfig::default()`, which is
//! `ErrorConfig::lenient()` — and that is worth checking rather than assuming,
//! because `ErrorMode`'s own `#[default]` is `Strict` and `ErrorConfig::default()`
//! is `strict()`. Only `NormalizeConfig`'s `Default` impl (`src/normalize/config.rs:91`)
//! substitutes the lenient config, so a reader who checks either of the other two
//! reaches the opposite conclusion. Measured here:
//!
//! | census counter | in lenient mode | in strict mode |
//! |---|---|---|
//! | `conflicts_accepted` 72 | accepted, re-emitted verbatim | **all 72 refused**, `W5002` / `OverlapConflictingEdits` |
//! | `prohibited_conditional_accepted` 40 | accepted | 16 of 40 refused (`W4007`, bare-transcript intronic); 24 accepted (`X`) |
//! | `prohibited_absolute_accepted` 32 | accepted | **all 32 still accepted** |
//!
//! So "72 conflicting alleles accepted instead of refused" is not the finding it
//! reads as: ferro has already adjudicated every one of those geometries as
//! un-denotable and refuses them under strict input hygiene. What remains open
//! is the narrower question of whether *lenient* mode may hand back a
//! description that denotes no sequence, which `#1406` answered deliberately
//! (preserve the conflict verbatim so strict re-reads it as one) and which this
//! file does not reopen.
//!
//! The 32 **absolute** acceptances are a different matter: they are refused in
//! neither mode, and they are the residue that a fix has to reach.
//!
//! # The clause-level decomposition, which the totals hide
//!
//! Measured over the corpus's 164 `RowKind::Prohibited` rows (see
//! [`the_absolute_prohibitions_ferro_accepts_are_three_clauses_and_no_others`]):
//!
//! ```text
//! ACCEPTED, absolute (24)   checklist.md:33  ins6            24   <- rank-1 defect (#1627)
//! ACCEPTED, conditional(16) checklist.md:20  bare NM_ intron 16   <- strict refuses; correct
//! REFUSED at parse (92)     checklist.md:31  c.10insT        24
//!                           checklist.md:49  del3            24
//!                           general.md:96    internal space  24
//!                           checklist.md:26  c.20+2_+5del    16
//!                           checklist.md:16  g. `*` offset    4   <- the GRAMMAR, not a rule
//! REFUSED, mode-gated (32)  standards.md:39  `X` base        24   <- closed by #1627
//!                           checklist.md:16  g. `+` offset    4   <- closed by #1628
//!                           checklist.md:45  g. hyphen range  4   <- closed by #1628
//! ```
//!
//! The last block is the shape the enforcement-stage ruling asks for: strict
//! refuses at PARSE, lenient and silent accept the input and then fail at
//! NORMALIZE. Those 32 rows are absent from **both** acceptance maps in
//! [`the_absolute_prohibitions_ferro_accepts_are_three_clauses_and_no_others`],
//! which is how a mode-gated refusal reads from here.
//!
//! **The corpus's rule ids name the BULLET, the citations here name the
//! PROHIBITIVE LINE, and the two differ by one or two lines.** The corpus calls
//! the `ins6` rule `checklist.md:32-insertion-states-its-sequence` because line
//! 32 opens the bullet ("**do you provide the inserted sequence?**"); the words
//! "is not allowed" are on line 33. Same for
//! `checklist.md:30-insertion-needs-two-anchors`, whose prohibition is on line
//! 31. Every citation in this file is the line carrying the quoted words,
//! verified against spec checkout `6f85311`; where a rule id appears verbatim
//! (in the decomposition assertions) it is the corpus's id and is not a
//! citation.
//!
//! Two things fall out that no total could show.
//!
//! **One numbered checklist item is half-enforced.** Item 3, "Insertions", has
//! two bullets: `:31` rejects `c.52insT` and `:33` rejects `c.5439_5430ins6`.
//! Ferro refuses the first at parse, with a citation in the error message, and
//! accepts the second in both modes.
//!
//! **One clause naming three symbols was enforced for one of them, and the one
//! it enforced was an accident.** `checklist.md:16` says a `g.` reference "can
//! not have nucleotides with additions like a `+`, `-`, or `*`". `g.*10del` is
//! refused because no production matches it — a *grammar* refusal, in every
//! mode, carrying no clause citation — while `g.266+2del` and `g.266-268del`
//! reached the AST and were accepted and re-emitted. #1628 closed the second
//! half with a real rule (`W4009`), keyed on `GenomePos::offset` and citing
//! `background/numbering.md:6`, which states the prohibition in the definitional
//! file rather than the checklist (and states it again at `:8` for `o.` and
//! `:11` for `m.`, which is why the rule covers all three axes).
//!
//! # Where the 32 prohibition-violating OUTPUTS come from
//!
//! All of them were re-emissions, and **there are none left**.
//! `prohibition_violating_outputs` decomposed as 24 × `standards.md:39` (the `X`
//! rows, closed by #1627) and 8 × `checklist.md:16` (the `+`-offset and
//! hyphen-range rows, closed by #1628) — exactly the set of accepted rows whose
//! *input* already violated one of the three textual checks. Ferro never
//! manufactured a prohibition violation from a legal input; it propagated one,
//! which is why refusing the inputs drove the counter to **0** rather than
//! needing an output-side filter.
//!
//! # Rulings recorded alongside this file
//!
//! Four clause conflicts turned up that a test cannot carry on its own, so they
//! are recorded in `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`
//! and pinned in `hgvs_spec_normalization_tests::RULING_STATUSES`:
//!
//! - `alignment-only-symbol-in-a-description` — `standards.md:36` lists `X` in
//!   the DNA symbol table; `:39` daggers it as alignment-only. **Decided** for
//!   `:39`.
//! - `bare-transcript-intronic-position` — `checklist.md:20` says a bare `NM_`
//!   cannot express an intron; `checklist.md:45` glosses one anyway.
//!   **Decided** for `:20`, as a conditional clause.
//! - `conflicting-member-geometry-refusal-scope` — whether `general.md:58`
//!   reaches the class or only the `del`+`dup` pair it names. **Decided** for
//!   `DNA/alleles.md:5`, with `general.md:56` cited to record that it does
//!   **not** reach a multi-member allele.
//! - `absolute-prohibition-enforcement-stage` — whether an absolute prohibition
//!   is refused at parse (unconditional) or at strict-mode normalize
//!   (opt-outable). The spec addresses descriptions, not stages, so this is a
//!   project decision. **Decided**: the stage is **mode-dependent** — strict
//!   fails at parse, lenient does not validate input conformance and fails only
//!   when it cannot normalize, silent is lenient without messages. See
//!   [`the_decided_target_is_a_mode_gated_refusal`], which is still `#[ignore]`d
//!   — but no longer for this clause's sake. #1630 moved `checklist.md:20`'s
//!   strict arm to parse, so that row passes; what keeps the `#[ignore]` is the
//!   one row still open, `checklist.md:33`'s `ins6` (#1627). See that test's own
//!   doc comment for the per-row state.
//!
//! # One measurement in this file was taken through the wrong door
//!
//! Every `parse_hgvs` call below — and every one in `conformance::spec_corpus`
//! — is the **config-less** entry point. Its own doc comment claims "Uses
//! strict error handling mode by default", and it does not: it applies no
//! `ErrorConfig` at all, which is #1632. So `parse_hgvs(x).is_err()` answers
//! "the grammar refuses `x`", not "every mode refuses `x`", and the assertions
//! here that read as the latter are the former.
//!
//! Re-measured through `parse_hgvs_with_config`, two of the five prohibitions
//! [`the_absolute_prohibitions_ferro_refuses_at_parse`] groups together are in
//! fact already mode-gated exactly as the ruling asks — `g.266del3` is refused
//! in strict and **repaired** to `g.266_268del` in lenient and silent, and
//! `c.10_12 del` is refused in strict and repaired to `c.10_12del`. The other
//! three (`c.10insT`, `c.20+2_+5del`, `g.*10del`) are refused in every mode
//! because no production exists for the shape. That distinction is what the
//! ruling turns on, so it is stated here rather than left in a PR body.

use std::collections::BTreeMap;

use ferro_hgvs::conformance::spec_corpus::{
    corpus, corpus_cores, denotation_of, CorpusBounds, Denotation, Frame, RefShape, RowKind,
    SpecCorpus, DENSE_CORE_LEN,
};
use ferro_hgvs::error_handling::{ErrorConfig, ErrorMode, ErrorType};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

// ---------------------------------------------------------------------------
// Harness
// ---------------------------------------------------------------------------

/// The `AT`-alphabet core the corpus draws first, shared with
/// `spec_corpus_regressions.rs` so the two files talk about one reference.
fn at_core() -> String {
    corpus_cores(1, DENSE_CORE_LEN).remove(0)
}

/// Both shuffle directions. Every adjudication below is a validity question, and
/// validity cannot depend on which way an ambiguous edit is shifted — so each is
/// checked in both, which is coverage `spec_corpus_regressions.rs` does not have.
const DIRECTIONS: [ShuffleDirection; 2] =
    [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime];

/// Lenient (default) normalization — the mode the axis census measures in.
fn lenient(frame: &Frame, input: &str, direction: ShuffleDirection) -> Result<String, String> {
    let parsed = parse_hgvs(input).map_err(|e| e.to_string())?;
    Normalizer::with_config(
        frame.provider().clone(),
        NormalizeConfig::default().with_direction(direction),
    )
    .normalize(&parsed)
    .map(|v| v.to_string())
    .map_err(|e| e.to_string())
}

/// Strict input-hygiene normalization.
fn strict(frame: &Frame, input: &str, direction: ShuffleDirection) -> Result<String, String> {
    let parsed = parse_hgvs(input).map_err(|e| e.to_string())?;
    Normalizer::with_config(
        frame.provider().clone(),
        NormalizeConfig::default()
            .with_direction(direction)
            .with_error_mode(ErrorMode::Strict),
    )
    .normalize(&parsed)
    .map(|v| v.to_string())
    .map_err(|e| e.to_string())
}

/// Silent-mode normalization — the third arm of the enforcement-stage ruling.
///
/// Used only by [`the_decided_target_is_a_mode_gated_refusal`], because the
/// three modes are byte-identical on every row this file pins: measured over
/// all 164 prohibited rows in both directions, strict, lenient and silent
/// produce the same string with zero disagreements. That is worth having a
/// helper for rather than a comment, so the ruling's third arm is asserted
/// rather than assumed.
fn silent(frame: &Frame, input: &str, direction: ShuffleDirection) -> Result<String, String> {
    let parsed = parse_hgvs(input).map_err(|e| e.to_string())?;
    Normalizer::with_config(
        frame.provider().clone(),
        NormalizeConfig::default()
            .with_direction(direction)
            .with_error_mode(ErrorMode::Silent),
    )
    .normalize(&parsed)
    .map(|v| v.to_string())
    .map_err(|e| e.to_string())
}

/// Whether `output` violates a prohibition the spec states in as many words.
///
/// A deliberate duplicate of `spec_conformance_axis::violated_prohibition`, kept
/// literal rather than shared: the axis test pins the TOTAL, this file pins the
/// per-clause decomposition, and a shared helper would let one silently redefine
/// what the other measures. Both are a closed three-entry list, so drift is a
/// three-line diff.
fn violated_prohibition(output: &str) -> Option<&'static str> {
    if output.contains(' ') {
        return Some("general.md:96");
    }
    if output.contains('X') || output.contains('x') {
        return Some("standards.md:39");
    }
    if let Some(body) = output.split_once(":g.").map(|(_, body)| body) {
        if body.contains('+') || body.contains('*') || body.contains('-') {
            return Some("checklist.md:16");
        }
    }
    None
}

/// The corpus, built once per test that needs it. ~0.3 s.
fn built() -> SpecCorpus {
    corpus(&CorpusBounds::default())
}

// ===========================================================================
// TASK 1 — the prohibitions
// ===========================================================================

// ---------------------------------------------------------------------------
// Adjudicated CORRECT: five absolute prohibitions refused at parse
// ---------------------------------------------------------------------------

/// **Question.** Which of the checklist's absolute prohibitions does ferro
/// already enforce, and where?
///
/// **Five, all at the `parse_hgvs` entry, and inside an allele as well as
/// alone.** This is the adjudicated-correct half, pinned as a guard: each of
/// these is rank-1 validity under the operator's precedence order
/// (`rulings[adjudication-precedence-order]`, "the spec where it speaks"), so a
/// later relaxation to accept-and-warn would be a rank-1 regression, not a
/// representation choice.
///
/// **This assertion does NOT say "in every mode", and it used to.** `parse_hgvs`
/// applies no `ErrorConfig` (#1632), so it answers for the grammar and not for a
/// mode. Measured through `parse_hgvs_with_config`, two of the five are already
/// mode-gated as `rulings[absolute-prohibition-enforcement-stage]` decided:
/// `NC_TEST.1:g.266del3` is refused in strict and repaired to
/// `NC_TEST.1:g.266_268del` in lenient and silent, and `NM_TEST.1:c.10_12 del`
/// is refused in strict and repaired to `NM_TEST.1:c.10_12del`. See the module
/// docs and [`the_decided_target_is_a_mode_gated_refusal`].
///
/// The clauses, verbatim from spec checkout `6f85311`:
///
/// - `checklist.md:31` — "The format `c.52insT` is **ambiguous**, and not
///   allowed."
/// - `checklist.md:49` — "Descriptions like `g.123del3` are not allowed, correct
///   is `g.123_125del`".
/// - `general.md:96` — "**NOTE:** spaces are *not* permitted in any HGVS
///   description".
/// - `checklist.md:26` — "The format `c.123-65_123-50` is correct, the format
///   `c.123-65_-50` is not, it is **incomplete**."
/// - `checklist.md:16` — "genomic (`g.`) reference sequences start with
///   nucleotide 1 and can not have nucleotides with additions like a `+`, `-`,
///   or `*`." Only the `*` half is refused *here*, and by the grammar rather
///   than by a rule; the `+`/`-` half is a mode-gated conformance refusal, so it
///   is not an every-mode parse refusal and does not belong in this list. See
///   [`the_clauses_three_offset_symbols_are_refused_alike`].
///
/// The allele form is checked for each because a hygiene check that runs once
/// per description rather than once per member stops firing exactly there, and
/// the corpus generates both forms of every prohibition to tell the two apart.
#[test]
fn the_absolute_prohibitions_ferro_refuses_at_parse() {
    for (input, in_allele, clause) in [
        (
            "NM_TEST.1:c.10insT",
            "NM_TEST.1:c.[10insT;5del]",
            "checklist.md:31 — `c.52insT` is ambiguous, and not allowed",
        ),
        (
            "NM_TEST.1:c.10del3",
            "NM_TEST.1:c.[10del3;5del]",
            "checklist.md:49 — descriptions like `g.123del3` are not allowed",
        ),
        (
            "NM_TEST.1:c.10_12 del",
            "NM_TEST.1:c.[10_12 del;5del]",
            "general.md:96 — spaces are not permitted in any HGVS description",
        ),
        (
            "NM_TEST.1:c.20+2_+5del",
            "NM_TEST.1:c.[20+2_+5del;24del]",
            "checklist.md:26 — `c.123-65_-50` is not correct, it is incomplete",
        ),
        (
            "NC_TEST.1:g.*10del",
            "NC_TEST.1:g.[*10del;265del]",
            "checklist.md:16 — a `g.` reference admits no `+`/`-`/`*` addition",
        ),
    ] {
        for form in [input, in_allele] {
            assert!(
                parse_hgvs(form).is_err(),
                "ADJUDICATED CORRECT, REGRESSED: `{form}` must be refused ({clause}). \
                 Rank-1 validity is never traded for a representation."
            );
        }
    }
}

/// **Question.** Does the refusal say *why*, or does it just fail?
///
/// **Two of the five name their clause in the message**, which is what makes
/// them adjudication records at the API boundary rather than incidental parser
/// limits — a user gets the citation and the corrected spelling. Pinned so a
/// refactor cannot quietly degrade them into a bare syntax error.
#[test]
fn the_refusals_that_cite_their_clause_keep_citing_it() {
    let deletion = parse_hgvs("NC_TEST.1:g.266del3").expect_err("del3 is refused");
    let message = deletion.to_string();
    assert!(
        message.contains("checklist.md:49"),
        "the `del3` refusal must cite checklist.md:49; got: {message}"
    );
    assert!(
        message.contains("NC_TEST.1:g.266_268del"),
        "the `del3` refusal must offer the corrected spelling; got: {message}"
    );

    let insertion = parse_hgvs("NM_TEST.1:c.10insT").expect_err("a single-anchor ins is refused");
    let message = insertion.to_string();
    assert!(
        message.contains("DNA/insertion.md"),
        "the single-anchor refusal must cite the insertion recommendation; got: {message}"
    );
}

// ---------------------------------------------------------------------------
// PINNED DEFECT: one checklist item, one bullet enforced and one not
// ---------------------------------------------------------------------------

/// **Question.** `checklist.md` item 3 ("Insertions") states two prohibitions,
/// four lines apart. Are they enforced alike?
///
/// **No.** `:31` is refused at parse with a citation; `:33` is accepted in both
/// modes and re-emitted unchanged. The two bullets are the same numbered item,
/// the same edit type and the same prohibitive phrasing:
///
/// - `:31` — "The format `c.52insT` is **ambiguous**, and not allowed."
/// - `:33` — "Describing a variant as `c.5439_5430ins6` is not allowed, the
///   inserted sequence (for `ins6`, e.g., `TGCCAT`) should be specified."
///
/// **Ruling: refuse, in both modes.** `:33` is absolute in as many words, so it
/// is rank-1 under the precedence order and not tradeable against anything. It
/// is also rank-1 a second way, independently of the clause: an insertion that
/// states a length rather than a sequence **denotes no sequence**, which
/// [`an_insertion_stating_a_length_denotes_no_sequence`] asserts rather than
/// argues. There is no lenient reading available — lenient mode exists for
/// inputs that can be understood and repaired, and `ins6` cannot be repaired
/// without inventing six bases.
///
/// 24 of the census's 32 absolute acceptances are this one shape.
#[test]
fn one_checklist_insertion_bullet_is_enforced_and_its_neighbour_is_not() {
    let core = at_core();
    let coding = Frame::build(RefShape::CodingSingleExon, &core);

    // `:31`, enforced — and enforced BECAUSE the anchor is one position, which is
    // what the clause is about. A bare `is_err()` would not say that: it passes
    // for any refusal at all, including one caused by the accession, the `c.`
    // prefix or the inserted base, so it could not distinguish "the clause is
    // enforced" from "this string happens not to parse". The negative control is
    // the attribution — same accession, same coordinate zone, same inserted base,
    // two-position anchor — and it must PARSE.
    assert!(parse_hgvs("NM_TEST.1:c.10insT").is_err());
    assert!(
        parse_hgvs("NM_TEST.1:c.10_11insT").is_ok(),
        "the control must parse, or the refusal above is not attributable to the \
         single-position anchor `:31` names"
    );

    // `:33`, not enforced — in either mode, either direction, alone or in an
    // allele.
    for direction in DIRECTIONS {
        assert_eq!(
            lenient(&coding, "NM_TEST.1:c.10_11ins6", direction).as_deref(),
            Ok("NM_TEST.1:c.10_11ins6"),
            "PINNED DEFECT — checklist.md:33 says an insertion sized by a number `is not \
             allowed, the inserted sequence … should be specified`. Correct behaviour: refuse."
        );
        assert_eq!(
            strict(&coding, "NM_TEST.1:c.10_11ins6", direction).as_deref(),
            Ok("NM_TEST.1:c.10_11ins6"),
            "PINNED DEFECT — and STRICT mode does not catch it either, so this is not a \
             hygiene-mode question: it is refused nowhere."
        );
    }

    // Inside an allele, both directions. The member order differs by direction
    // because the legal sibling shifts; the offending member survives either way,
    // which is the point.
    assert_eq!(
        lenient(
            &coding,
            "NM_TEST.1:c.[10_11ins6;9del]",
            ShuffleDirection::ThreePrime
        )
        .as_deref(),
        Ok("NM_TEST.1:c.[9del;10_11ins6]"),
        "PINNED DEFECT — a per-description hygiene check would have to fire here too"
    );
    assert_eq!(
        lenient(
            &coding,
            "NM_TEST.1:c.[10_11ins6;9del]",
            ShuffleDirection::FivePrime
        )
        .as_deref(),
        Ok("NM_TEST.1:c.[1del;10_11ins6]"),
        "PINNED DEFECT — the 5' direction is not covered by spec_corpus_regressions.rs"
    );
}

/// **Question.** Is `c.10_11ins6` merely *prohibited*, or is it not a
/// description of a variant at all?
///
/// **Not a description of a variant.** The corpus's own applier — which does not
/// consult the normalizer — cannot produce a resulting sequence for it. That
/// makes the ruling above independent of how strongly `checklist.md:33` is
/// worded: even reading `:33` as a style preference, a description that denotes
/// nothing fails rank-1 validity by the same test the axis census applies to
/// outputs (`outputs_denoting_no_sequence`, asserted to be a rank-1 defect).
///
/// Recorded because the tempting counter-argument — "`ins6` is legacy shorthand,
/// accept it leniently and let the consumer decide" — dies on this measurement:
/// there is nothing to be lenient *toward*.
#[test]
fn an_insertion_stating_a_length_denotes_no_sequence() {
    let core = at_core();
    let coding = Frame::build(RefShape::CodingSingleExon, &core);
    let denoted = denotation_of(coding.provider(), coding.served(), "NM_TEST.1:c.10_11ins6");
    assert!(
        !matches!(denoted, Denotation::Sequence(_)),
        "`ins6` must denote no sequence for the rank-1 argument to hold; got {denoted:?}"
    );

    // The control, so the assertion above is about `ins6` and not about the
    // applier being unable to handle insertions at all.
    assert!(
        matches!(
            denotation_of(
                coding.provider(),
                coding.served(),
                "NM_TEST.1:c.10_11insTGCCAT"
            ),
            Denotation::Sequence(_)
        ),
        "the same insertion WITH its sequence stated must denote one"
    );
}

// ---------------------------------------------------------------------------
// PINNED DEFECT: one clause, three symbols, two behaviours
// ---------------------------------------------------------------------------

/// **Question.** `checklist.md:16` names three symbols in one breath — "genomic
/// (`g.`) reference sequences start with nucleotide 1 and can not have
/// nucleotides with additions like a `+`, `-`, or `*`." Are they enforced alike?
///
/// **They are now. They were not.** Before #1628's parse/normalize half, `*` was
/// refused and `+`/`-` were accepted *and re-emitted*, in all three modes — the
/// 8 `checklist.md:16` rows that were the whole of
/// `prohibition_violating_outputs`.
///
/// **This test was a pinned defect and is now an adjudicated-correct guard.**
/// The stage is the one `rulings[absolute-prohibition-enforcement-stage]` names:
/// strict refuses at PARSE (`W4009`), lenient and silent accept the input and
/// then fail at NORMALIZE, because a genomic accession carries no exon boundary
/// for an offset to be measured from. `background/numbering.md:6` is the primary
/// authority — it states the prohibition for `g.` in the definitional file, and
/// `:8`/`:11` state it again for `o.` and `m.` — with `checklist.md:16` saying
/// the same thing in the checklist's register.
///
/// **The three symbols are still not symmetric, and that is worth keeping.**
/// `*` and a leading `-` reach no production of the genomic position grammar, so
/// they are refused as *syntax*, in every mode, and never carry `W4009`. Only
/// `+` and a trailing `-` reach the AST and need a conformance rule. Reading the
/// `*` refusal as enforcement is precisely what made this clause look
/// one-third-done for as long as it did.
#[test]
fn the_clauses_three_offset_symbols_are_refused_alike() {
    let core = at_core();
    let genomic = Frame::build(RefShape::Genomic, &core);

    // `*` — refused by the GRAMMAR, alone and in an allele. The allele case
    // needs the control more than the lone one does: `g.[*10del;265del]` has a
    // second member and allele brackets, so a bare `is_err()` there passes just
    // as happily if the allele syntax or the sibling member is what failed, and
    // the claim being made is specifically that the `*` is what does it.
    // Dropping `*` from the same allele must therefore PARSE.
    //
    // The messages are deliberately not pinned. Unlike the `:31` insertion
    // refusal, which names its clause (`DNA/insertion.md:95-101`), these two
    // surface a raw combinator failure — `Error(Error { input: "*10del", code:
    // Tag })` — so pinning the text would pin an internal that carries no
    // meaning and would churn on any parser refactor. The control pins the
    // attribution, which is the part the claim rests on.
    assert!(parse_hgvs("NC_TEST.1:g.*10del").is_err());
    assert!(parse_hgvs("NC_TEST.1:g.[*10del;265del]").is_err());
    assert!(
        parse_hgvs("NC_TEST.1:g.[10del;265del]").is_ok(),
        "the control must parse, or the allele refusal above is not attributable \
         to the `*` — it could be the brackets or the sibling member"
    );

    // `+` and a trailing `-` — refused by the CONFORMANCE RULE, in every mode,
    // both directions, alone and in an allele.
    for (input, allele, clause) in [
        (
            "NC_TEST.1:g.266+2del",
            "NC_TEST.1:g.[266+2del;265del]",
            "numbering.md:6 / checklist.md:16 — a `g.` reference admits no `+` addition",
        ),
        (
            "NC_TEST.1:g.266-268del",
            "NC_TEST.1:g.[266-268del;265del]",
            "checklist.md:45 — `c.12-14del` is `Not correct`; on `g.` it has no \
             intronic reading either, so numbering.md:6 also applies",
        ),
    ] {
        for form in [input, allele] {
            for direction in DIRECTIONS {
                for (mode, produced) in [
                    ("lenient", lenient(&genomic, form, direction)),
                    ("strict", strict(&genomic, form, direction)),
                    ("silent", silent(&genomic, form, direction)),
                ] {
                    let refusal = produced.unwrap_or_else(|e| {
                        assert!(
                            e.contains("W4009") && e.contains("numbering.md"),
                            "{mode}'s refusal of `{form}` must cite the clause ({clause}); \
                             got: {e}"
                        );
                        String::new()
                    });
                    assert!(
                        refusal.is_empty(),
                        "ADJUDICATED CORRECT, REGRESSED: {mode} mode emitted `{refusal}` for \
                         `{form}` instead of refusing it ({clause}). Rule 1 of the README \
                         ruleset is about OUTPUT and has no mode escape."
                    );
                }
            }
        }
    }
}

/// **Question.** `checklist.md:45` calls `c.12-14del` "Not correct". Must ferro
/// refuse a hyphen between two positions wherever it appears?
///
/// **No — and this is the one place the answer is axis-dependent.** Read the
/// clause in full: "Not correct is `c.12-14del`, this describes a deletion of
/// nucleotide -14 in the intron directly 5' of nucleotide `c.12`". The spec is
/// not saying the string is unparseable; it is saying the string **parses to
/// something else than the author meant**. On a `c.` axis it is a well-formed
/// intronic description and a normalizer cannot see the author's intent, so
/// there is nothing to refuse — refusing it would refuse every legitimate
/// intronic deletion of the same shape.
///
/// On a `g.` axis there is no intron to read it into, so it has no legal reading
/// at all and `checklist.md:16` refuses it independently. That is the ruling:
/// **`:45` is a prohibition on the author, not on the description; `:16` is the
/// one a normalizer can enforce.**
///
/// The evidence that the axis really is the discriminator, rather than a
/// convenient story: on the `c.` axis ferro's strict mode refuses it, but as an
/// **intronic-on-bare-transcript** finding (`W4007`), which is
/// `checklist.md:20`'s clause — i.e. it is being read exactly as `:45` says it
/// reads.
#[test]
fn the_hyphen_range_is_a_legal_intronic_description_on_a_coding_axis() {
    let core = at_core();
    let coding = Frame::build(RefShape::CodingSingleExon, &core);

    for direction in DIRECTIONS {
        // ADJUDICATED CORRECT: lenient mode leaves it alone, because on `c.` it
        // is a legal description of an intronic deletion.
        assert_eq!(
            lenient(&coding, "NM_TEST.1:c.12-14del", direction).as_deref(),
            Ok("NM_TEST.1:c.12-14del"),
            "checklist.md:45 says `c.12-14del` DESCRIBES a deletion in the intron 5' of \
             c.12 — a normalizer cannot refuse it for meaning what the spec says it means"
        );

        // And strict mode refuses it under checklist.md:20, not under :45.
        let refusal = strict(&coding, "NM_TEST.1:c.12-14del", direction)
            .expect_err("strict refuses a bare-transcript intronic position");
        assert!(
            refusal.contains("W4007") || refusal.contains("IntronicOnBareTranscript"),
            "the strict refusal must be the bare-transcript finding (checklist.md:20), not a \
             hyphen-range one; got: {refusal}"
        );
    }
}

// ---------------------------------------------------------------------------
// PINNED DEFECT: the alignment-only symbols
// ---------------------------------------------------------------------------

/// **Question.** `background/standards.md` lists both `X` and `-` in the DNA
/// symbol table and daggers both, footnoting at `:39` "† used in alignment
/// only." Are they treated alike?
///
/// **Ruling: refuse both — and the STAGE is mode-dependent.** Two decided
/// records settle the two halves. `rulings[alignment-only-symbol-in-a-description]`
/// settles *what*: neither symbol may appear in a description. The strength
/// label is moot, because the alphabet half of rank-1 validity settles it
/// without the prohibition half — `general.md:48` requires "nucleotides in
/// CAPITALS using IUPAC-IUBMB assigned nucleotide symbols", and `X` is not one
/// of the fifteen (the corpus's own `DNA_SYMBOLS` excludes it for exactly this
/// reason). `rulings[absolute-prohibition-enforcement-stage]` settles *where*:
/// strict fails at PARSE, lenient does not validate input conformance and fails
/// only when it cannot NORMALIZE, silent is lenient without messages.
///
/// **This test was a pinned defect and is now an adjudicated-correct guard**
/// (#1627). Before the fix all three modes accepted every shape below and
/// re-emitted it verbatim with an *empty* warning vector — normalization was
/// not impossible, it was vacuous, which is worse, because the output looked
/// normalized while carrying a spelling the recommendations prohibit.
///
/// **The two symbols are NOT symmetric, and that is the finding worth keeping.**
/// `X` reaches the AST (as `InsertedSequence::Named`) and needed a rule. `-`
/// reaches it in no position whatsoever — leading, trailing or interior, `ins`
/// or `delins` — because it matches no arm of `parse_inserted_sequence` and the
/// variant-level trailing-character check refuses what is left unconsumed. That
/// is a *grammar* accident rather than a conformance rule, so it is asserted
/// below as a standing fact and not as the template the `X` fix copied.
#[test]
fn an_alignment_only_symbol_is_refused_in_every_mode_for_both_x_and_dash() {
    let core = at_core();
    let coding = Frame::build(RefShape::CodingSingleExon, &core);
    let genomic = Frame::build(RefShape::Genomic, &core);
    let multi = Frame::build(RefShape::CodingMultiExon(Strand::Plus), &core);

    // `-` — the grammar refuses it in every position, so in every mode. No
    // conformance rule is reached and none is needed.
    for input in [
        "NM_TEST.1:c.10delins-",
        "NM_TEST.1:c.10delins-ACGT",
        "NM_TEST.1:c.10delinsACGT-",
        "NM_TEST.1:c.10delinsAC-GT",
        "NM_TEST.1:c.10_11ins-",
    ] {
        assert!(
            parse_hgvs(input).is_err(),
            "{input}: `-` is footnoted identically to `X` at standards.md:39"
        );
        for config in [
            ErrorConfig::strict(),
            ErrorConfig::lenient(),
            ErrorConfig::silent(),
        ] {
            assert!(
                parse_hgvs_with_config(input, config).is_err(),
                "{input}: `-` must be refused in every mode"
            );
        }
    }

    // `X` — refused everywhere too, now, but at the stage the ruling names.
    // The list carries the EMBEDDED shapes as well as the lone symbol: one
    // stray `X` reclassifies an otherwise-literal run wholesale, so a guard
    // covering only `delinsX` would reproduce the undercount the census had.
    for (frame, input, allele) in [
        (
            &coding,
            "NM_TEST.1:c.10delinsX",
            "NM_TEST.1:c.[10delinsX;9del]",
        ),
        (
            &coding,
            "NM_TEST.1:c.10delinsACGTX",
            "NM_TEST.1:c.[10delinsACGTX;9del]",
        ),
        (
            &coding,
            "NM_TEST.1:c.10delinsXACGT",
            "NM_TEST.1:c.[10delinsXACGT;9del]",
        ),
        (
            &coding,
            "NM_TEST.1:c.10delinsACXGT",
            "NM_TEST.1:c.[10delinsACXGT;9del]",
        ),
        (
            &genomic,
            "NC_TEST.1:g.266delinsX",
            "NC_TEST.1:g.[266delinsX;265del]",
        ),
        (
            &genomic,
            "NC_TEST.1:g.266delinsACGTX",
            "NC_TEST.1:g.[266delinsACGTX;265del]",
        ),
        (
            &multi,
            "NM_TEST.1:c.10delinsX",
            "NM_TEST.1:c.[10delinsX;9del]",
        ),
        (
            &multi,
            "NM_TEST.1:c.10delinsACGTX",
            "NM_TEST.1:c.[10delinsACGTX;9del]",
        ),
    ] {
        // STRICT — refused at PARSE. Strict validates input conformance, not
        // merely parseability, so the refusal is reached before the AST is
        // handed to any consumer.
        for probe in [input, allele] {
            let err = parse_hgvs_with_config(probe, ErrorConfig::strict())
                .expect_err("strict refuses an alignment-only symbol at parse")
                .to_string();
            assert!(
                err.contains("W3028"),
                "{probe}: strict must refuse with W3028; got: {err}"
            );
        }

        // LENIENT and SILENT — accepted at parse, because neither validates
        // input conformance, and both then fail at NORMALIZE. That is the
        // ruling's own ground for lenient: it fails only when it cannot
        // normalize, and a masked base names no nucleotide to normalize.
        for config in [ErrorConfig::lenient(), ErrorConfig::silent()] {
            for probe in [input, allele] {
                assert!(
                    parse_hgvs_with_config(probe, config.clone()).is_ok(),
                    "{probe}: lenient/silent do not validate input conformance"
                );
            }
        }

        // The parse-stage difference between lenient and silent is the message,
        // not the verdict: lenient carries W3028, silent is quiet.
        let lenient_parse = parse_hgvs_with_config(input, ErrorConfig::lenient())
            .expect("lenient accepts at parse");
        assert!(
            lenient_parse
                .warnings
                .iter()
                .any(|w| w.error_type.code() == "W3028"),
            "{input}: lenient must say why it will fail later"
        );
        let silent_parse =
            parse_hgvs_with_config(input, ErrorConfig::silent()).expect("silent accepts at parse");
        assert!(
            !silent_parse
                .warnings
                .iter()
                .any(|w| w.error_type.code() == "W3028"),
            "{input}: silent is lenient without messages"
        );

        for direction in DIRECTIONS {
            for probe in [input, allele] {
                // All three modes refuse at normalize. Output conformance is
                // rule 1 of the README ruleset and has no mode escape, so this
                // rung is deliberately NOT mode-gated.
                for (label, outcome) in [
                    ("strict", strict(frame, probe, direction)),
                    ("lenient", lenient(frame, probe, direction)),
                    ("silent", silent(frame, probe, direction)),
                ] {
                    let err = outcome.expect_err(&format!(
                        "{probe}: {label} must refuse to normalize an alignment-only symbol"
                    ));
                    assert!(
                        err.contains("W3028"),
                        "{probe}: {label} refusal must name W3028; got: {err}"
                    );
                }
            }
        }
    }

    // The constraint on the narrowing: the named-element arm exists for these,
    // none carries a daggered symbol, and all three modes must leave them
    // exactly as they were.
    for input in [
        "NC_TEST.1:g.266delinsAluYb8",
        "NC_TEST.1:g.266delinsLINE1",
        "NC_TEST.1:g.266delinsL1",
        "NC_TEST.1:g.266delinsAlu",
    ] {
        for config in [
            ErrorConfig::strict(),
            ErrorConfig::lenient(),
            ErrorConfig::silent(),
        ] {
            assert!(
                parse_hgvs_with_config(input, config).is_ok(),
                "{input} is a legitimate mobile-element name"
            );
        }
        // Across all three modes, like the positive block above. Checking only
        // lenient here would leave the negative control weaker than the thing
        // it is controlling for: the normalize-stage refusal is deliberately
        // NOT mode-gated, so a predicate that over-matched would redden strict
        // and silent while this assertion stayed green.
        for direction in DIRECTIONS {
            for (label, outcome) in [
                ("strict", strict(&genomic, input, direction)),
                ("lenient", lenient(&genomic, input, direction)),
                ("silent", silent(&genomic, input, direction)),
            ] {
                assert_eq!(
                    outcome.as_deref(),
                    Ok(input),
                    "{input} must still round-trip in {label} mode"
                );
            }
        }
    }

    // The rank-1 argument, asserted rather than argued: `X` is not a base, so
    // the description denotes nothing.
    assert!(
        !matches!(
            denotation_of(coding.provider(), coding.served(), "NM_TEST.1:c.10delinsX"),
            Denotation::Sequence(_)
        ),
        "a `delins` inserting `X` must denote no sequence"
    );
}

// ---------------------------------------------------------------------------
// ADJUDICATED CORRECT: the conditional clause, and why the mode split is right
// ---------------------------------------------------------------------------

/// **Question.** `checklist.md:20` says an `NM_` "can only be used to describe
/// variants in introns using a `c.` prefix when a genomic reference sequence is
/// given on which the coding DNA reference sequence is annotated". Must ferro
/// refuse a bare `NM_…:c.20+2del`?
///
/// **In strict mode yes, and it does (`W4007`). In lenient mode no.** That
/// split is the correct answer rather than an accident, and it is recorded as
/// `rulings[bare-transcript-intronic-position]`, decided for `checklist.md:20`:
///
/// - The clause is **conditional in form** — "can only be used … when" states a
///   condition on the reference sequence, not a prohibition on a shape. That is
///   precisely the class input-hygiene mode exists for, and it is why the corpus
///   classifies it `Strength::Conditional` rather than `Absolute`.
/// - The spec itself reads a bare-`c.` intronic description four items later:
///   `checklist.md:45` glosses `c.12-14del` as "a deletion of nucleotide -14 in
///   the intron directly 5' of nucleotide `c.12`", with no genomic wrapper in
///   sight. A clause the spec does not apply to its own worked gloss cannot be
///   read as an absolute bar.
/// - The wrapper form the clause asks for is accepted, so the difference is a
///   refusal of the bare form and not missing support.
///
/// 16 of the census's 40 conditional acceptances are this shape, and all 16 are
/// refused by strict mode — so that counter is a lenient-mode figure, not a
/// statement that ferro does not enforce the clause.
///
/// **What this does NOT settle**, and the record says so: ferro's normalizer
/// *emits* bare-`NM_` intronic descriptions (`outputs_leaving_the_transcript`,
/// 371 at 3'), so lenient mode manufactures a description strict mode refuses.
/// That is the reverse of the laundering `#1406` legislated against, and it is a
/// defect in the shift, not in this clause's enforcement.
#[test]
fn a_bare_transcript_intronic_position_is_refused_in_strict_only() {
    let core = at_core();
    let multi = Frame::build(RefShape::CodingMultiExon(Strand::Plus), &core);

    for direction in DIRECTIONS {
        for input in ["NM_TEST.1:c.20+2del", "NM_TEST.1:c.[20+2del;24del]"] {
            assert_eq!(
                lenient(&multi, input, direction).as_deref(),
                Ok(input),
                "lenient mode accepts the bare form, per the ruling — checklist.md:20 is \
                 conditional and the spec glosses a bare-`c.` intronic description at :45"
            );
            let refusal = strict(&multi, input, direction)
                .expect_err("ADJUDICATED CORRECT, REGRESSED: strict must refuse the bare form");
            assert!(
                refusal.contains("W4007") || refusal.contains("IntronicOnBareTranscript"),
                "the refusal must name checklist.md:20's finding; got: {refusal}"
            );
        }

        // The wrapper form `checklist.md:20` asks for is accepted in BOTH modes,
        // which is what makes the refusal above a refusal of the bare spelling
        // rather than of intronic positions.
        assert_eq!(
            strict(&multi, "NC_SYNTH.1(NM_TEST.1):c.20+2del", direction).as_deref(),
            Ok("NC_SYNTH.1(NM_TEST.1):c.20+2del"),
            "the NC_(NM_) form is the one the clause names as correct"
        );
    }
}

/// **Question.** `rulings[absolute-prohibition-enforcement-stage]` puts an
/// input-conformance check at **parse**, gated by mode. `checklist.md:20` is
/// enforced with the right mode split but at the wrong **stage** — strict
/// refuses at *normalize*. Does it refuse at parse?
///
/// **It does now (#1630).** This is the third of the three defects the record's
/// own census names, and the last one whose fix is a stage move rather than a
/// relaxation:
///
/// > **`checklist.md:20`'s refusal is at normalize, not at parse.** It is a
/// > conditional clause and its own record (`bare-transcript-intronic-position`)
/// > ratifies the mode split; moving the stage is cosmetic there, but it should
/// > move for uniformity.
///
/// **"Cosmetic" is about the VERDICT, not about what a caller sees**, and the
/// distinction is the whole reason this is a guard rather than a refactor. A
/// caller that parses strictly and never normalizes — `parse_hgvs_with_config`
/// is public API and validating an identifier is its most obvious use — was
/// told a bare `NM_…:c.20+2del` conforms. It does not, and `checklist.md:20`
/// is the clause that says so. The record's ground is exactly this: "whether
/// the INPUT conforms is answered before the input is accepted, not part-way
/// through normalizing it."
///
/// **The mode split is unchanged, and that is deliberate.** Only the strict
/// arm's stage moves. Lenient still accepts, now with the `W4007` it always
/// owed a parse-only caller; silent still accepts in silence, which is what
/// "lenient with no error messages" means. `rulings[bare-transcript-intronic-
/// position]` decided the split itself and nothing here reopens it.
///
/// **The normalize rung stays**, and is still reachable: this file's own
/// [`strict`] helper parses through the config-less `parse_hgvs` and then
/// normalizes strictly, which is the shape
/// [`a_bare_transcript_intronic_position_is_refused_in_strict_only`] pins. A
/// caller that parses leniently and normalizes strictly takes the same path.
/// So the two rungs are not redundant — they answer for two different callers,
/// and removing either would drop one of them.
///
/// Both controls are load-bearing. Without the **wrapper** control a bare
/// `is_err()` is satisfied by refusing intronic positions outright; without the
/// **exonic** control it is satisfied by refusing bare `NM_` accessions
/// outright. The claim is narrower than either: it is the *conjunction* that
/// `checklist.md:20` names.
#[test]
fn a_bare_transcript_intronic_position_is_refused_at_parse_in_strict_mode() {
    // Alone and inside an allele, because a hygiene check that runs once per
    // description rather than once per member stops firing exactly there.
    for input in [
        "NM_TEST.1:c.20+2del",
        "NM_TEST.1:c.[20+2del;24del]",
        // The `n.` arm of the same rule, on a bare non-coding transcript.
        "NR_TEST.1:n.20+2del",
    ] {
        // STRICT — refused at PARSE, because strict validates input conformance.
        let refusal = parse_hgvs_with_config(input, ErrorConfig::strict())
            .expect_err(
                "strict mode must refuse a bare-transcript intronic position at PARSE \
                 (checklist.md:20, rulings[absolute-prohibition-enforcement-stage])",
            )
            .to_string();
        assert!(
            refusal.contains("W4007"),
            "the parse refusal of `{input}` must name checklist.md:20's finding by code; \
             got: {refusal}"
        );

        // LENIENT — accepts, and says so. Silence here is what left a
        // parse-only caller with no signal at all.
        let lenient_parse = parse_hgvs_with_config(input, ErrorConfig::lenient())
            .expect("lenient mode does not validate input conformance, so it must accept");
        assert!(
            lenient_parse
                .warnings
                .iter()
                .any(|w| w.error_type == ErrorType::IntronicOnBareTranscript),
            "lenient mode must accept `{input}` WITH a W4007 warning; got warnings: {:?}",
            lenient_parse.warnings
        );

        // SILENT — lenient with no error messages, per the ruling's third arm.
        let silent_parse = parse_hgvs_with_config(input, ErrorConfig::silent())
            .expect("silent mode is lenient, so it must accept");
        assert!(
            !silent_parse
                .warnings
                .iter()
                .any(|w| w.error_type == ErrorType::IntronicOnBareTranscript),
            "silent mode must accept `{input}` WITHOUT a warning; got: {:?}",
            silent_parse.warnings
        );
    }

    // CONTROL 1 — the wrapper form `checklist.md:20` names as correct must
    // still parse in strict mode. Without this the refusal above is
    // attributable to intronic positions rather than to the bare reference.
    assert!(
        parse_hgvs_with_config("NC_SYNTH.1(NM_TEST.1):c.20+2del", ErrorConfig::strict()).is_ok(),
        "the NC_(NM_) form is the one the clause names as correct; refusing it would \
         make this rule a ban on introns"
    );

    // CONTROL 2 — an EXONIC position on the same bare accession must still
    // parse in strict mode. Without this the refusal is attributable to the
    // bare `NM_` accession rather than to the intronic offset on it.
    assert!(
        parse_hgvs_with_config("NM_TEST.1:c.20del", ErrorConfig::strict()).is_ok(),
        "a bare NM_ accession is fine; it is the INTRONIC offset on one that \
         checklist.md:20 conditions"
    );
}

// ---------------------------------------------------------------------------
// The decomposition, measured over the corpus
// ---------------------------------------------------------------------------

/// **Question.** The census pins 32 absolute and 16 conditional acceptances.
/// *Which clauses*, and does either figure change under strict mode?
///
/// **Three clauses hold the whole absolute figure, and strict mode does not move
/// it at all.** The conditional figure drops to zero under strict, because
/// `checklist.md:20` is enforced there.
///
/// This is the test that turns two totals into an adjudication: a fix that
/// refuses `ins6` moves 24 of the 32, and one that refuses `g.` offsets moves
/// the remaining 8. Nothing else contributes, so nothing else needs a ruling —
/// and if a new clause ever appears in this map, the corpus has grown a shape
/// nobody has adjudicated.
///
/// **The conditional figure was 40 and is now 16 (#1627).**
/// `standards.md:39`'s 24 rows have left BOTH maps: strict refuses them at
/// parse, lenient and silent at normalize. That is the shape of a
/// mode-dependent enforcement stage as seen from here — a clause that is gone
/// from the strict map but present in the lenient one is refused *only* by
/// strict input hygiene (`checklist.md:20`); a clause gone from both is refused
/// in every mode, either by the grammar or, as here, because normalization
/// cannot proceed.
///
/// Measured over the 164 `RowKind::Prohibited` rows only, so it costs ~0.3 s
/// rather than the full census's 37 s.
#[test]
fn the_absolute_prohibitions_ferro_accepts_are_three_clauses_and_no_others() {
    let built = built();
    for direction in DIRECTIONS {
        let mut lenient_accepted: BTreeMap<(&str, &str), usize> = BTreeMap::new();
        let mut strict_accepted: BTreeMap<(&str, &str), usize> = BTreeMap::new();
        let mut rows = 0usize;
        for row in &built.rows {
            if row.kind != RowKind::Prohibited {
                continue;
            }
            rows += 1;
            let frame = row.frame();
            let (clause, strength) = row.prohibition.expect("a prohibited row names its clause");
            let key = (clause, strength.label());
            for spelling in &row.spellings {
                if lenient(&frame, spelling, direction).is_ok() {
                    *lenient_accepted.entry(key).or_default() += 1;
                }
                if strict(&frame, spelling, direction).is_ok() {
                    *strict_accepted.entry(key).or_default() += 1;
                }
            }
        }
        assert_eq!(
            rows, 164,
            "VACUOUS if this drops: the prohibited stratum is the whole denominator"
        );

        assert_eq!(
            lenient_accepted,
            BTreeMap::from([
                (
                    (
                        "checklist.md:20-intron-needs-a-genomic-wrapper",
                        "conditional"
                    ),
                    16
                ),
                (
                    ("checklist.md:32-insertion-states-its-sequence", "absolute"),
                    24
                ),
                // `checklist.md:16-genomic-has-no-offsets` (4) and
                // `checklist.md:45-range-is-underscore` (4) used to sit here and
                // no longer do (#1628): every mode now refuses them. Strict at
                // parse (`W4009`), lenient and silent at normalize, because a
                // genomic accession carries no exon boundary for the offset to
                // be measured from. See
                // `the_clauses_three_offset_symbols_are_refused_alike`.
                // `standards.md:39`'s 24 rows used to sit here and no longer do
                // (#1627): every mode now refuses them. Strict at parse, lenient
                // and silent at normalize. See
                // `an_alignment_only_symbol_is_refused_in_every_mode_for_both_x_and_dash`.
            ]),
            "the per-clause decomposition of the census's 32 absolute + 16 conditional \
             acceptances moved ({direction:?}). Every entry is adjudicated in this file; a NEW \
             clause here is a shape nobody has ruled on."
        );

        // And the decomposition must still SUM to the census it decomposes. The
        // map above pins the parts; without this, a re-bless of the axis pins
        // leaves the parts behind with nothing failing, which is the same stale-
        // literal trap `CENSUS_TOTAL` avoids.
        let pinned = match direction {
            ShuffleDirection::ThreePrime => &crate::spec_conformance_axis::THREE_PRIME,
            _ => &crate::spec_conformance_axis::FIVE_PRIME,
        };
        for (strength, total) in [
            ("absolute", pinned.prohibited_absolute_accepted),
            ("conditional", pinned.prohibited_conditional_accepted),
        ] {
            let summed: usize = lenient_accepted
                .iter()
                .filter(|((_, s), _)| *s == strength)
                .map(|(_, count)| count)
                .sum();
            assert_eq!(
                summed, total,
                "the per-clause {strength} counts sum to {summed}, but \
                 `spec_conformance_axis` pins {total} ({direction:?})"
            );
        }

        // Strict mode: `checklist.md:20` drops out, nothing else does.
        assert_eq!(
            strict_accepted,
            BTreeMap::from([
                (
                    ("checklist.md:32-insertion-states-its-sequence", "absolute"),
                    24
                ),
                // The two genomic-offset clauses used to survive strict mode
                // too, and no longer do (#1628): strict refuses them at PARSE,
                // which is where the enforcement-stage ruling puts an
                // input-conformance check. They are now absent from BOTH maps
                // rather than only this one.
                // `standards.md:39` used to survive strict mode too, and no
                // longer does (#1627): strict refuses it at PARSE, which is
                // where the enforcement-stage ruling puts an input-conformance
                // check. It is now absent from BOTH maps rather than only this
                // one, because lenient and silent refuse it at normalize.
            ]),
            "strict mode must refuse checklist.md:20 and standards.md:39 and nothing else \
             ({direction:?}). The 32 ABSOLUTE acceptances survive strict mode, which is why \
             they are the residue a fix has to reach."
        );
    }
}

/// **Question.** Where do the 32 `prohibition_violating_outputs` come from —
/// does normalization *manufacture* a prohibited shape, or only propagate one?
///
/// **Only propagate.** Every violating output is the re-emission of a row whose
/// input already violated the same clause: 8 × `checklist.md:16` and
/// 24 × `standards.md:39`, and the two figures sum to the census's pinned 32 in
/// both directions. Since the census counts violations over **all** 58,552
/// outputs and the prohibited stratum alone accounts for all 32, no legal input
/// produces a prohibited output.
///
/// That entailment is the reason this test can be cheap, and it is stated
/// explicitly because it is load-bearing: it holds only while
/// `spec_conformance_axis`'s pin is 32. If that pin moves, this test's sum must
/// move with it, and if it does not, a legal input has started producing a
/// prohibited output — which is a strictly worse defect than any pinned here.
///
/// **Consequence for the fix.** An output-side filter would be treating a
/// symptom. Refusing the three input shapes adjudicated above drives this
/// counter to zero as a side effect.
#[test]
fn every_prohibition_violating_output_is_a_re_emitted_prohibited_input() {
    /// Taken from the axis pin rather than re-typed, so the compiler carries
    /// the entailment the doc above calls load-bearing. A literal copy would
    /// survive a re-bless of the census and then report a stale constant as if
    /// a legal input had started producing a prohibited output.
    const CENSUS_TOTAL: usize =
        crate::spec_conformance_axis::THREE_PRIME.prohibition_violating_outputs;
    const _: () = assert!(
        CENSUS_TOTAL == crate::spec_conformance_axis::FIVE_PRIME.prohibition_violating_outputs,
        "both directions must pin the same prohibition_violating_outputs for this test's \
         direction-independent sum to hold"
    );

    let built = built();
    for direction in DIRECTIONS {
        let mut by_clause: BTreeMap<&str, usize> = BTreeMap::new();
        for row in &built.rows {
            if row.kind != RowKind::Prohibited {
                continue;
            }
            let frame = row.frame();
            for spelling in &row.spellings {
                let Ok(output) = lenient(&frame, spelling, direction) else {
                    continue;
                };
                if let Some(clause) = violated_prohibition(&output) {
                    *by_clause.entry(clause).or_default() += 1;
                }
            }
        }
        assert_eq!(
            by_clause,
            // `standards.md:39`'s 24 dropped out in #1627 and
            // `checklist.md:16`'s 8 in #1628 — none of those rows produces an
            // output at all any more, in any mode, so there is nothing left to
            // violate. The counter is now EMPTY, which is the end state this
            // test was written to reach: ferro emits no description that
            // violates a prohibition the spec states in as many words.
            BTreeMap::new(),
            "the clause decomposition of the violating outputs moved ({direction:?})"
        );
        assert_eq!(
            by_clause.values().sum::<usize>(),
            CENSUS_TOTAL,
            "the prohibited stratum no longer accounts for the census's {CENSUS_TOTAL} \
             prohibition-violating outputs ({direction:?}). Either the census pin moved — \
             re-bless both together — or a LEGAL input has started producing a prohibited \
             output, which is a new and worse defect."
        );
        assert!(
            !by_clause.contains_key("general.md:96"),
            "no output should carry a space; general.md:96 is refused at parse"
        );
    }
}

// ---------------------------------------------------------------------------
// The DECIDED TARGET, asserted directly
// ---------------------------------------------------------------------------

/// Whether `output` still carries the prohibited token of `clause`.
///
/// Deliberately separate from [`violated_prohibition`], which is a literal
/// duplicate of the axis test's three-entry list and must stay that way. This
/// one is keyed per clause and covers `checklist.md:33`'s `ins<n>`, which the
/// axis list does not — the axis counts *output* violations it can recognise
/// textually, and a numeric insertion payload is a shape its list predates.
fn output_still_violates(clause: &str, output: &str) -> bool {
    match clause {
        "checklist.md:33" => output.contains("ins6"),
        "standards.md:39" => output.contains('X'),
        "checklist.md:16" | "checklist.md:45" => output
            .split_once(":g.")
            .is_some_and(|(_, body)| body.contains('+') || body.contains('-')),
        other => panic!("no output check for {other}"),
    }
}

/// **The decided target**, asserted directly: an absolute prohibition is
/// refused at **parse in strict mode**, and no mode may ever *emit* the
/// prohibited spelling.
///
/// **Authority.** The `OPERATOR RULING, 2026-08-10` paragraph of
/// `rulings[absolute-prohibition-enforcement-stage]`. The stage is
/// mode-dependent: strict fails at parse because strict validates input
/// conformance rather than merely parseability; lenient does not validate input
/// conformance and fails only when it cannot normalize; silent is lenient
/// without messages.
///
/// **The second half is not a mode question and is asserted in all three.**
/// Rule 1 of the README ruleset — "Output follows the HGVS recommendations.
/// Absolute — never traded." — is about OUTPUT, so it has no mode escape. That
/// is precisely why the ruling's mode gate costs nothing: accepting a
/// non-conformant input and *normalizing it to a conformant output* trades
/// nothing away, while re-emitting the prohibited spelling is a rule-1
/// violation whatever mode produced it.
///
/// **`#[ignore]`d because ONE of its four rows is still open.** Three are now
/// closed — `standards.md:39`'s `delinsX` by #1627, and both `g.`-offset rows by
/// #1628 — leaving `checklist.md:33`'s `ins6`, which is #1627's remaining half.
/// Delete the `#[ignore]` once that lands; the other three rows already pass and
/// are guarded individually by
/// [`an_alignment_only_symbol_is_refused_in_every_mode_for_both_x_and_dash`] and
/// [`the_clauses_three_offset_symbols_are_refused_alike`].
///
/// The `g.` offset rows were the ones no un-normalizability argument reached.
/// `ins6` and `delinsX` denote nothing (`Denotation::Inexpressible`, and
/// `to_spdi` refuses both by name), so a lenient mode that failed when it could
/// not normalize would refuse them for that reason alone. `g.266+2del` appeared
/// to denote a sequence only because `to_spdi` silently **discarded** the offset
/// and answered for `g.266del`; that half is closed by #1641 and #1734, and the
/// normalize half by the `W4009` rule, so both now refuse rather than disagree.
#[test]
#[ignore = "one row still open: checklist.md:33's `ins6` — see #1627. The \
            other three rows pass (#1627's `delinsX`, #1628's two `g.` offsets), \
            as does checklist.md:20 (#1630), which this test does not yet carry \
            a row for"]
fn the_decided_target_is_a_mode_gated_refusal() {
    let core = at_core();
    let coding = Frame::build(RefShape::CodingSingleExon, &core);
    let genomic = Frame::build(RefShape::Genomic, &core);

    for (frame, input, clause) in [
        (
            &coding,
            "NM_TEST.1:c.10_11ins6",
            "checklist.md:33", // "is not allowed, the inserted sequence"
        ),
        (
            &coding,
            "NM_TEST.1:c.10delinsX",
            "standards.md:39", // "used in alignment only"
        ),
        (
            &genomic,
            "NC_TEST.1:g.266+2del",
            "checklist.md:16", // "can not have nucleotides with additions"
        ),
        (
            &genomic,
            "NC_TEST.1:g.266-268del",
            "checklist.md:45", // and checklist.md:16 on a `g.` axis
        ),
    ] {
        // STRICT — refused at parse, because strict validates input conformance.
        assert!(
            parse_hgvs_with_config(input, ErrorConfig::strict()).is_err(),
            "`{input}` must be refused at parse in STRICT mode ({clause})"
        );

        // LENIENT and SILENT — free not to judge the input, never free to emit
        // the prohibited spelling. Either outcome below is conformant; what is
        // not is an `Ok` still carrying the token.
        for direction in DIRECTIONS {
            for (mode, produced) in [
                ("lenient", lenient(frame, input, direction)),
                ("silent", silent(frame, input, direction)),
            ] {
                if let Ok(output) = produced {
                    assert!(
                        !output_still_violates(clause, &output),
                        "{mode} mode emitted `{output}`, which still violates {clause}. Rule 1 \
                         of the README ruleset is about OUTPUT and has no mode escape: lenient \
                         may decline to validate the input, but it must then either normalize \
                         it to a conformant description or fail because it cannot."
                    );
                }
            }
        }
    }
}

// ===========================================================================
// TASK 2 — conflicting member geometries
// ===========================================================================

// ---------------------------------------------------------------------------
// The correction: `conflicts_accepted` is a lenient-mode figure
// ---------------------------------------------------------------------------

/// **Question.** The census reports "72 conflicting alleles accepted instead of
/// refused". Does ferro have no refusal for these?
///
/// **It has one, and it fires on all 72.** Every conflicting row the census
/// counts as accepted is refused under `ErrorMode::Strict` as
/// `OverlapConflictingEdits` / `W5002` — "HGVS spec defines no canonical form for
/// this case". The census normalizes with `NormalizeConfig::default()`, which is
/// lenient, and nothing in its module docs says so.
///
/// **This does not make the counter wrong; it makes it narrower.** The open
/// question is not "are these geometries adjudicated?" — they are, as
/// un-denotable — but "may the DEFAULT mode return a description that denotes no
/// sequence?". `#1406` answered that deliberately: preserve the conflict
/// verbatim, so that strict mode re-reading the *output* reaches the same verdict
/// it reached on the input, rather than picking a winner among orderings the spec
/// does not rank. That contract is asserted here as the second half, because it
/// is what makes lenient acceptance safe rather than lossy.
///
/// Keyed on the diagnostic, never on `is_err()`: a bare error check would count
/// an out-of-range coordinate or a stated-reference mismatch as a successful
/// refusal.
#[test]
fn every_conflicting_allele_the_census_calls_accepted_is_refused_in_strict_mode() {
    let built = built();
    for direction in DIRECTIONS {
        let mut rows = 0usize;
        let mut lenient_accepted = 0usize;
        let mut refused_at_parse = 0usize;
        let mut refused_by_strict = 0usize;
        let mut by_geometry: BTreeMap<&str, usize> = BTreeMap::new();
        for row in &built.rows {
            if row.kind != RowKind::Conflict {
                continue;
            }
            rows += 1;
            let frame = row.frame();
            for spelling in &row.spellings {
                if parse_hgvs(spelling).is_err() {
                    refused_at_parse += 1;
                    continue;
                }
                let output = match lenient(&frame, spelling, direction) {
                    Ok(output) => output,
                    Err(_) => continue,
                };
                lenient_accepted += 1;
                *by_geometry.entry(row.geometry.label()).or_default() += 1;

                let refusal = strict(&frame, spelling, direction).expect_err(
                    "ADJUDICATED: a description whose members claim one territory denotes no \
                     sequence and must be refused under strict input hygiene",
                );
                assert!(
                    refusal.contains("W5002") || refusal.contains("OverlapConflictingEdits"),
                    "`{spelling}` was refused, but not as an overlap conflict: {refusal}"
                );
                refused_by_strict += 1;

                // #1406's contract: the lenient OUTPUT must still read as a
                // conflict, or lenient mode has laundered a description strict
                // mode refuses into one it admits.
                //
                // Keyed on the diagnostic, exactly like the refusal above, and
                // for the reason this module's doc comment already gives: a bare
                // `is_err()` counts an out-of-range coordinate or a
                // stated-reference mismatch as a successful refusal. Here that
                // would be worse than elsewhere, because the string being
                // re-read is one ferro just EMITTED — so an `is_err()` that
                // fired on a malformed output would report the laundering
                // contract as upheld by the very defect that broke it.
                let laundered = strict(&frame, &output, direction).expect_err(
                    "lenient mode turned `{spelling}` into an output strict mode ACCEPTS — \
                     the conflict was laundered out of the description (#1406)",
                );
                assert!(
                    laundered.contains("W5002") || laundered.contains("OverlapConflictingEdits"),
                    "lenient mode turned `{spelling}` into `{output}`, which strict mode \
                     refuses — but NOT as an overlap conflict: {laundered}. The #1406 contract \
                     is that the output still reads as the conflict it came from, not merely \
                     that it fails to round-trip"
                );
            }
        }
        // The axis pin, not a literal copy of it — same reason as `CENSUS_TOTAL`
        // above. Both directions pin the same figure.
        let census_conflicts = crate::spec_conformance_axis::THREE_PRIME.conflicts_accepted;
        assert_eq!(
            census_conflicts,
            crate::spec_conformance_axis::FIVE_PRIME.conflicts_accepted,
            "the two directions' `conflicts_accepted` pins have diverged, so this test can no \
             longer compare against one of them"
        );
        assert_eq!(rows, 96, "VACUOUS if the conflict stratum shrinks");
        assert_eq!(
            lenient_accepted, census_conflicts,
            "this must equal the census's `conflicts_accepted` ({direction:?})"
        );
        assert_eq!(
            refused_by_strict, census_conflicts,
            "all {census_conflicts} must be strict refusals"
        );
        assert_eq!(
            refused_at_parse, 24,
            "the self-replacement rows are refused earlier still — see \
             `the_self_replacement_geometry_the_clause_names_is_refused_at_parse`"
        );
        assert_eq!(
            by_geometry,
            BTreeMap::from([
                ("coincident-insertions", 24),
                ("nested", 24),
                ("overlapping", 24),
            ]),
            "the geometry decomposition of the 72 moved ({direction:?})"
        );
    }
}

// ---------------------------------------------------------------------------
// The geometries, adjudicated one at a time
// ---------------------------------------------------------------------------

/// **Question.** `general.md:58` — "descriptions removing part of a reference
/// sequence and replacing it with part of the same sequence are not allowed
/// (e.g., `NM_004006.2:c.[762_768del;767_774dup]`)" — is the one absolute
/// prohibition the spec states *between two members of one allele*. Is it
/// enforced as a clause, or as a spelling?
///
/// **As a spelling.** The exact edit-type pair the clause uses as its example —
/// a `del` and an overlapping `dup` — is refused at **parse**, in both modes, on
/// every axis ("Self-cancelling allele"). Swap either member for another edit
/// type with the same geometry and nothing fires until strict-mode normalize.
///
/// **Ruling: correct outcome, under-general rule.** The refusal is
/// adjudicated-correct and pinned as a guard; what is pinned as a defect is that
/// it is keyed on `(del, dup)` rather than on the geometry. `:58`'s stated
/// ground is that the description removes reference and replaces it with the
/// same reference — a property of the *footprints*, which `c.[9_12del;11_14del]`
/// has just as much as the spec's example does.
///
/// **Do not cite `general.md:56` here.** `:56` — "when a description is possible
/// according to several types, the preferred description is: (1) substitution,
/// (2) deletion, …" — ranks competing descriptions of **one** span. It says
/// nothing about two members of one allele, and citing it against a multi-member
/// allele is this repository's recorded cautionary error. `:58` is the
/// member-vs-member clause; `:56` is its neighbour and is not interchangeable
/// with it.
#[test]
fn the_self_replacement_geometry_the_clause_names_is_refused_at_parse() {
    for input in [
        // The spec's own example shape, on each axis the corpus builds.
        "NM_TEST.1:c.[9_12del;10_13dup]",
        "NC_TEST.1:g.[265_268del;266_269dup]",
        // …and the clause's literal example.
        "NM_004006.2:c.[762_768del;767_774dup]",
    ] {
        let refusal = parse_hgvs(input).expect_err(
            "ADJUDICATED CORRECT, REGRESSED: general.md:58 says a description removing part of \
             a reference and replacing it with part of the same sequence is not allowed",
        );
        assert!(
            refusal.to_string().contains("Self-cancelling"),
            "`{input}` must be refused AS a self-replacement; got: {refusal}"
        );
    }

    // The under-generality, pinned: the same geometry with two deletions parses
    // fine. That is what makes the refusal above a spelling rule rather than a
    // reading of `:58`.
    assert!(
        parse_hgvs("NM_TEST.1:c.[9_12del;11_14del]").is_ok(),
        "PINNED DEFECT — general.md:58's ground is that the footprints overlap, not that the \
         edit types are `del` and `dup`. Correct behaviour: refuse this too."
    );
}

/// **Question.** Are the three geometries the census counts genuinely
/// un-denotable, or merely redundant / order-dependent? Adjudicated one at a
/// time, because the answers are argued differently.
///
/// **All three must be refused, and all three are — in strict mode.** The
/// assertions below are the adjudicated-correct guard; the reasoning per
/// geometry is what makes this a record:
///
/// - **Overlapping** (`c.[9_12del;11_14del]`). Directly `general.md:58`'s shape.
///   The clause's own example `NM_004006.2:c.[762_768del;767_774dup]` is a
///   partial intersection of two footprints — `767_774` reaches into `762_768` —
///   so the geometry, not the edit types, is what the clause describes. The
///   second member also *states* reference bases (`11_14`) that the first has
///   already removed, so the description asserts two incompatible things about
///   one interval.
///
/// - **Coincident insertions** (`c.[8_9insAT;8_9insTC]`). Genuinely
///   un-denotable, and for a reason `:58` does not cover: two zero-width members
///   at one interbase have **no order**, and the two orders give different
///   sequences (`…ATTC…` and `…TCAT…`). `DNA/alleles.md:5` defines an allele as
///   "a series of variants on one chromosome" — a series is ordered, and this
///   pair supplies no order. `background/basics.md:38` puts **unequivocal**
///   among the four values the recommendations are designed for. Refusal is the
///   only answer that does not invent an ordering the spec does not rank, which
///   is `#1406`'s standing reasoning applied to a new geometry.
///
/// - **Nested** (`c.[9_12del;10_11del]`). The weakest of the three, and the one
///   with a real counter-argument: the inner deletion is *entailed* by the outer,
///   so one could normalize the pair to `c.9_12del` and call it redundant rather
///   than invalid. **That reading loses**, and
///   [`a_nested_pair_is_refused_on_geometry_not_on_redundancy`] is why: the
///   redundancy is a property of this particular pair, not of the geometry, and
///   no clause distinguishes the redundant nesting from the contradictory one.
///   A rule that absorbed `[9_12del;10_11del]` would have to refuse
///   `[9_12del;10_11delinsGG]`, giving two spellings one position apart opposite
///   treatments with nothing to key on. The rule has to be geometric.
///
/// Checked on all three reference shapes and in both directions, which is
/// coverage `spec_corpus_regressions::conflicting_member_geometries_are_normalized_instead_of_refused`
/// does not have — it fixes `CodingSingleExon` and the 3' direction.
#[test]
fn the_three_conflicting_geometries_are_refused_in_strict_mode_on_every_axis() {
    let core = at_core();
    let cases: [(RefShape, [(&str, &str); 3]); 3] = [
        (
            RefShape::CodingSingleExon,
            [
                ("nested", "NM_TEST.1:c.[9_12del;10_11del]"),
                ("overlapping", "NM_TEST.1:c.[9_12del;11_14del]"),
                ("coincident insertions", "NM_TEST.1:c.[8_9insAT;8_9insTC]"),
            ],
        ),
        (
            RefShape::Genomic,
            [
                ("nested", "NC_TEST.1:g.[265_268del;266_267del]"),
                ("overlapping", "NC_TEST.1:g.[265_268del;267_270del]"),
                (
                    "coincident insertions",
                    "NC_TEST.1:g.[264_265insAT;264_265insTC]",
                ),
            ],
        ),
        (
            RefShape::CodingMultiExon(Strand::Plus),
            [
                ("nested", "NM_TEST.1:c.[9_12del;10_11del]"),
                ("overlapping", "NM_TEST.1:c.[9_12del;11_14del]"),
                ("coincident insertions", "NM_TEST.1:c.[8_9insAT;8_9insTC]"),
            ],
        ),
    ];

    for (shape, rows) in cases {
        let frame = Frame::build(shape, &core);
        for (geometry, input) in rows {
            // The premise, so the test asserts a defect rather than a preference.
            assert_eq!(
                denotation_of(frame.provider(), frame.served(), input),
                Denotation::NoSequence,
                "{geometry} on {}: {input} must denote no sequence",
                shape.label()
            );

            for direction in DIRECTIONS {
                let refusal = strict(&frame, input, direction).expect_err(
                    "ADJUDICATED CORRECT, REGRESSED: a description whose members claim one \
                     territory denotes no sequence — general.md:58 for the overlapping case, \
                     DNA/alleles.md:5 and background/basics.md:38 for the rest",
                );
                assert!(
                    refusal.contains("W5002") || refusal.contains("OverlapConflictingEdits"),
                    "{geometry} on {}: refused, but not as a conflict: {refusal}",
                    shape.label()
                );

                // Lenient mode preserves it verbatim — #1406's contract, and the
                // reason the census counts it as "accepted".
                assert_eq!(
                    lenient(&frame, input, direction).as_deref(),
                    Ok(input),
                    "{geometry} on {}: lenient mode must preserve the conflict verbatim, not \
                     repair it into a description that denotes something",
                    shape.label()
                );
            }
        }
    }
}

/// **Question.** Could the nested geometry be *normalized* instead of refused,
/// on the ground that the inner member is redundant?
///
/// **No, and this is the evidence.** Take the same nesting and change the inner
/// member from `del` to `delinsGG`: `c.[9_12del;10_11delinsGG]`. Nothing is
/// redundant now — the inner member asserts two bases the outer member deleted —
/// and the pair is un-denotable outright. The two inputs differ by four
/// characters, have identical geometry, and no clause distinguishes them.
///
/// So a rule that absorbed the redundant nesting would have to inspect payloads
/// to decide whether to absorb or refuse, and would give neighbouring spellings
/// opposite verdicts. **The rule has to be geometric: two members of one allele
/// may not claim intersecting reference territory, whatever they do with it.**
/// That is also how ferro already behaves — both are `W5002` — so this test
/// pins the *reason*, which is the part a future simplification could lose.
///
/// Recorded per `CLAUDE.md`'s "record what was refuted, not only what was
/// decided": "nested is just redundant, absorb it" is the plausible belief this
/// kills.
#[test]
fn a_nested_pair_is_refused_on_geometry_not_on_redundancy() {
    let core = at_core();
    let frame = Frame::build(RefShape::CodingSingleExon, &core);

    for (why, input) in [
        (
            "redundant nesting — the inner del is entailed",
            "NM_TEST.1:c.[9_12del;10_11del]",
        ),
        (
            "contradictory nesting — the inner delins asserts bases the outer deleted",
            "NM_TEST.1:c.[9_12del;10_11delinsGG]",
        ),
    ] {
        assert_eq!(
            denotation_of(frame.provider(), frame.served(), input),
            Denotation::NoSequence,
            "{why}: {input} must denote no sequence"
        );
        for direction in DIRECTIONS {
            let refusal = strict(&frame, input, direction)
                .expect_err("both nestings must reach the same verdict");
            assert!(
                refusal.contains("W5002") || refusal.contains("OverlapConflictingEdits"),
                "{why}: {refusal}"
            );
        }
    }
}

/// **Question.** Coincident insertions with *different* payloads are equivocal
/// because the two orders give different sequences. What about identical
/// payloads — `c.[8_9insAT;8_9insAT]`, where both orders agree?
///
/// **Still equivocal, and still refused.** The ambiguity moves from order to
/// multiplicity: the description can be read as two insertions of `AT` (giving
/// `ATAT`) or as one insertion stated twice (giving `AT`). `DNA/alleles.md:5`'s
/// "a series of variants" does not say whether a repeated member is two variants
/// or one, and `basics.md:38`'s **unequivocal** is the value at stake either way.
///
/// **Ruling: refuse — the same verdict, reached by a different argument**, which
/// is exactly why it is worth pinning: a fix that keyed refusal on "the two
/// orders disagree" would let this one through, and it would be wrong.
///
/// The exactly-duplicated non-insertion member (`c.[9_12del;9_12del]`) is
/// included as the sibling case, where the two readings agree on the sequence
/// and the description is merely redundant — and ferro refuses it too, which is
/// the consistent answer given the nested ruling above.
#[test]
fn coincident_insertions_with_identical_payloads_are_still_equivocal() {
    let core = at_core();
    let frame = Frame::build(RefShape::CodingSingleExon, &core);

    for (why, input) in [
        (
            "two orders, two sequences",
            "NM_TEST.1:c.[8_9insAT;8_9insTC]",
        ),
        (
            "one order, two multiplicities",
            "NM_TEST.1:c.[8_9insAT;8_9insAT]",
        ),
        (
            "an exactly repeated deletion",
            "NM_TEST.1:c.[9_12del;9_12del]",
        ),
    ] {
        for direction in DIRECTIONS {
            let refusal = strict(&frame, input, direction).expect_err(
                "ADJUDICATED: a repeated or coincident member is equivocal, whether the \
                 ambiguity is in the ORDER or in the MULTIPLICITY",
            );
            assert!(
                refusal.contains("W5002") || refusal.contains("OverlapConflictingEdits"),
                "{why}: {refusal}"
            );
        }
    }
}

/// **Question.** Does the verdict depend on the order the members are written
/// in?
///
/// **No, and it must not** — a conflict is a property of the description, not of
/// how its author happened to sort it (`#1406`'s title, applied to geometry).
/// Pinned in both modes: strict refuses either order, lenient preserves either
/// order verbatim without re-sorting one into the other.
///
/// The lenient half is the load-bearing one. Ferro *does* sort members by
/// position in the ordinary case — `c.[10_11ins6;9del]` comes back as
/// `c.[9del;10_11ins6]` — so "preserved verbatim" here is a deliberate exception
/// for conflicting alleles and not an absence of behaviour.
#[test]
fn member_order_does_not_change_the_conflict_verdict() {
    let core = at_core();
    let frame = Frame::build(RefShape::CodingSingleExon, &core);

    for (a, b) in [
        (
            "NM_TEST.1:c.[9_12del;10_11del]",
            "NM_TEST.1:c.[10_11del;9_12del]",
        ),
        (
            "NM_TEST.1:c.[9_12del;11_14del]",
            "NM_TEST.1:c.[11_14del;9_12del]",
        ),
        (
            "NM_TEST.1:c.[8_9insAT;8_9insTC]",
            "NM_TEST.1:c.[8_9insTC;8_9insAT]",
        ),
    ] {
        for direction in DIRECTIONS {
            for input in [a, b] {
                assert!(
                    strict(&frame, input, direction).is_err(),
                    "`{input}` must be refused whichever order its members are written in"
                );
                assert_eq!(
                    lenient(&frame, input, direction).as_deref(),
                    Ok(input),
                    "`{input}` must be preserved verbatim, not re-sorted into its sibling"
                );
            }
        }
    }
}
