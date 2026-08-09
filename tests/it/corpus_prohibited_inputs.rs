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
//! ACCEPTED, absolute (32)   checklist.md:33  ins6            24   <- rank-1 defect
//!                           checklist.md:16  g. `+` offset    4   <- rank-1 defect
//!                           checklist.md:45  g. hyphen range  4   <- rank-1 defect
//! ACCEPTED, conditional(40) standards.md:39  `X` base        24   <- rank-1 defect (alphabet)
//!                           checklist.md:20  bare NM_ intron 16   <- strict refuses; correct
//! REFUSED at parse (92)     checklist.md:31  c.10insT        24
//!                           checklist.md:49  del3            24
//!                           general.md:96    internal space  24
//!                           checklist.md:26  c.20+2_+5del    16
//!                           checklist.md:16  g. `*` offset    4
//! ```
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
//! **One clause naming three symbols is enforced for one of them.**
//! `checklist.md:16` says a `g.` reference "can not have nucleotides with
//! additions like a `+`, `-`, or `*`". `g.*10del` is refused at parse;
//! `g.266+2del` and `g.266-268del` are accepted and re-emitted.
//!
//! # Where the 32 prohibition-violating OUTPUTS come from
//!
//! All of them are re-emissions. `prohibition_violating_outputs` decomposes as
//! 8 × `checklist.md:16` (the `+`-offset and hyphen-range rows, whose output
//! body still carries the offset) and 24 × `standards.md:39` (the `X` rows) —
//! which is exactly the set of accepted rows whose *input* already violated one
//! of the three textual checks. Ferro never manufactures a prohibition
//! violation from a legal input; it propagates one. That matters for where a fix
//! goes: an output-side filter would be treating the symptom, and the input gate
//! is the whole of the defect.
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
//!   (opt-outable). The spec addresses descriptions, not stages.
//!   **Undecided.**

use std::collections::BTreeMap;

use ferro_hgvs::conformance::spec_corpus::{
    corpus, corpus_cores, denotation_of, CorpusBounds, Denotation, Frame, RefShape, RowKind,
    SpecCorpus, DENSE_CORE_LEN,
};
use ferro_hgvs::error_handling::ErrorMode;
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
/// **Five, all at parse — so in both input-hygiene modes, and inside an allele
/// as well as alone.** This is the adjudicated-correct half, pinned as a guard:
/// each of these is rank-1 validity under the operator's precedence order
/// (`rulings[adjudication-precedence-order]`, "the spec where it speaks"), so a
/// later relaxation to accept-and-warn would be a rank-1 regression, not a
/// representation choice.
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
///   or `*`." — enforced for `*` only; see
///   [`one_clause_names_three_offset_symbols_and_ferro_refuses_one`].
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
/// **No: `*` is refused at parse, `+` and `-` are accepted and re-emitted.**
///
/// **Ruling: refuse all three, at parse, in both modes.** The clause is absolute
/// ("can not have") and it is rank-1 twice over, because the grammar half of
/// rank-1 validity says the same thing — `docs/syntax.yaml`'s genomic position
/// admits no offset, which is why `*` is already rejected. There is no reading on
/// which `+` and `-` are legal here and `*` is not; the difference is an
/// artifact of which spellings the parser happens to have a production for.
///
/// This is also the whole of the output-side defect for this clause: the 8
/// `checklist.md:16` violations in `prohibition_violating_outputs` are these two
/// shapes, alone and in an allele, on both cores.
#[test]
fn one_clause_names_three_offset_symbols_and_ferro_refuses_one() {
    let core = at_core();
    let genomic = Frame::build(RefShape::Genomic, &core);

    // `*` — refused, alone and in an allele. The allele case needs the control
    // more than the lone one does: `g.[*10del;265del]` has a second member and
    // allele brackets, so a bare `is_err()` there passes just as happily if the
    // allele syntax or the sibling member is what failed, and the claim being
    // made is specifically that the `*` is what does it. Dropping `*` from the
    // same allele must therefore PARSE.
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

    // `+` and `-` — accepted, both modes, both directions, alone and in an
    // allele, and the offset survives into the output.
    for (input, allele, three_prime_allele, five_prime_allele, clause) in [
        (
            "NC_TEST.1:g.266+2del",
            "NC_TEST.1:g.[266+2del;265del]",
            "NC_TEST.1:g.[265del;266+2del]",
            "NC_TEST.1:g.[256del;266+2del]",
            "checklist.md:16 — a `g.` reference admits no `+` addition",
        ),
        (
            "NC_TEST.1:g.266-268del",
            "NC_TEST.1:g.[266-268del;265del]",
            "NC_TEST.1:g.[265del;266-268del]",
            "NC_TEST.1:g.[256del;266-268del]",
            "checklist.md:45 — `c.12-14del` is `Not correct`; on `g.` it has no \
             intronic reading either, so checklist.md:16 also applies",
        ),
    ] {
        for direction in DIRECTIONS {
            assert_eq!(
                lenient(&genomic, input, direction).as_deref(),
                Ok(input),
                "PINNED DEFECT — {clause}. Correct behaviour: refuse at parse."
            );
            assert_eq!(
                strict(&genomic, input, direction).as_deref(),
                Ok(input),
                "PINNED DEFECT — strict mode does not catch it either ({clause})"
            );
        }
        assert_eq!(
            lenient(&genomic, allele, ShuffleDirection::ThreePrime).as_deref(),
            Ok(three_prime_allele),
            "PINNED DEFECT — the allele form is accepted too ({clause})"
        );
        assert_eq!(
            lenient(&genomic, allele, ShuffleDirection::FivePrime).as_deref(),
            Ok(five_prime_allele),
            "PINNED DEFECT — and in the 5' direction ({clause})"
        );

        // The output-side half: the emitted description still carries the offset,
        // so the prohibition is enforced on neither side.
        assert_eq!(
            violated_prohibition(three_prime_allele),
            Some("checklist.md:16"),
            "the re-emitted allele must still violate the clause, or this row does not \
             belong to the 8 counted in `prohibition_violating_outputs`"
        );
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
/// **No: `-` is refused at parse on every axis, `X` is accepted in both modes on
/// every axis.** One clause, two symbols, two behaviours — and the repository
/// itself holds two incompatible readings of that clause: the corpus classifies
/// an `X` input `Strength::Conditional`, while `spec_conformance_axis`'s
/// `violated_prohibition` counts an `X` output as an **absolute** violation.
///
/// **Ruling: refuse, as `-` already is.** Recorded as
/// `rulings[alignment-only-symbol-in-a-description]`, decided for
/// `standards.md:39`. The strength label is moot, because the alphabet half of
/// rank-1 validity settles it without the prohibition half: `general.md:48`
/// requires "nucleotides in CAPITALS using IUPAC-IUBMB assigned nucleotide
/// symbols", and `X` is not one of the fifteen (the corpus's own `DNA_SYMBOLS`
/// excludes it for exactly this reason). A `delins` whose inserted sequence is
/// `X` therefore denotes no sequence, asserted below.
///
/// 24 of the census's 40 conditional acceptances, and 24 of its 32
/// prohibition-violating outputs, are this shape.
#[test]
fn an_alignment_only_symbol_is_refused_for_dash_and_accepted_for_x() {
    let core = at_core();
    let coding = Frame::build(RefShape::CodingSingleExon, &core);
    let genomic = Frame::build(RefShape::Genomic, &core);
    let multi = Frame::build(RefShape::CodingMultiExon(Strand::Plus), &core);

    // `-` — refused at parse, so in both modes.
    assert!(
        parse_hgvs("NM_TEST.1:c.10delins-").is_err(),
        "`-` is footnoted identically to `X` at standards.md:39"
    );

    // `X` — accepted everywhere.
    for (frame, input, allele) in [
        (
            &coding,
            "NM_TEST.1:c.10delinsX",
            "NM_TEST.1:c.[10delinsX;9del]",
        ),
        (
            &genomic,
            "NC_TEST.1:g.266delinsX",
            "NC_TEST.1:g.[266delinsX;265del]",
        ),
        (
            &multi,
            "NM_TEST.1:c.10delinsX",
            "NM_TEST.1:c.[10delinsX;9del]",
        ),
    ] {
        for direction in DIRECTIONS {
            assert_eq!(
                lenient(frame, input, direction).as_deref(),
                Ok(input),
                "PINNED DEFECT — standards.md:39 footnotes `X` as `used in alignment only`, \
                 and general.md:48 admits only IUPAC-IUBMB symbols. Correct behaviour: refuse, \
                 as `-` already is."
            );
            assert_eq!(
                strict(frame, input, direction).as_deref(),
                Ok(input),
                "PINNED DEFECT — strict mode does not catch it either"
            );
            // The output-side half: `X` survives into the emitted description.
            let out = lenient(frame, allele, direction).expect("the allele form is accepted");
            assert_eq!(
                violated_prohibition(&out),
                Some("standards.md:39"),
                "the re-emitted allele must still carry `X`"
            );
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

// ---------------------------------------------------------------------------
// The decomposition, measured over the corpus
// ---------------------------------------------------------------------------

/// **Question.** The census pins 32 absolute and 40 conditional acceptances.
/// *Which clauses*, and does either figure change under strict mode?
///
/// **Three clauses hold the whole absolute figure, and strict mode does not move
/// it at all.** The conditional figure halves under strict, because
/// `checklist.md:20` is enforced there.
///
/// This is the test that turns two totals into an adjudication: a fix that
/// refuses `ins6` moves 24 of the 32, and one that refuses `g.` offsets moves
/// the remaining 8. Nothing else contributes, so nothing else needs a ruling —
/// and if a new clause ever appears in this map, the corpus has grown a shape
/// nobody has adjudicated.
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
                (("checklist.md:16-genomic-has-no-offsets", "absolute"), 4),
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
                (("checklist.md:45-range-is-underscore", "absolute"), 4),
                (
                    ("standards.md:39-alignment-only-symbols", "conditional"),
                    24
                ),
            ]),
            "the per-clause decomposition of the census's 32 absolute + 40 conditional \
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
                (("checklist.md:16-genomic-has-no-offsets", "absolute"), 4),
                (
                    ("checklist.md:32-insertion-states-its-sequence", "absolute"),
                    24
                ),
                (("checklist.md:45-range-is-underscore", "absolute"), 4),
                (
                    ("standards.md:39-alignment-only-symbols", "conditional"),
                    24
                ),
            ]),
            "strict mode must refuse checklist.md:20 and only checklist.md:20 ({direction:?}). \
             The 32 ABSOLUTE acceptances survive strict mode, which is why they are the residue \
             a fix has to reach."
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
            BTreeMap::from([("checklist.md:16", 8), ("standards.md:39", 24)]),
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
