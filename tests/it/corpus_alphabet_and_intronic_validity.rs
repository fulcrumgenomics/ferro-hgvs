//! Allele-level characterization tests for two rank-1 validity defects the
//! exhaustive spec corpus found: #51 (an alignment-only `X` accepted as
//! nucleotide sequence content) and #52 (an intronic position accepted on a
//! bare transcript accession). No production fix here — the fix locus
//! (`src/normalize/`, `src/hgvs/`) is being edited elsewhere; this file is
//! diagnosis plus committed tests, in the house style of
//! `tests/it/spec_corpus_regressions.rs`.
//!
//! # Why a separate file from `spec_corpus_regressions.rs`
//!
//! That file already pins the **single-member** shape of both defects:
//! `an_alignment_only_x_base_is_accepted_as_a_nucleotide` (standalone
//! `c.10delinsX`) and `a_bare_transcript_accession_accepts_an_intronic_position`
//! (standalone `c.20+2del`, plus the genomic-wrapped positive control). The
//! corpus's actual exemplars for #51/#52 are **allele** forms —
//! `c.[10delinsX;9del]` and `c.[20+2del;24del]` — and the corpus counts (24 and
//! 16 outputs respectively) are allele-level counts. This file pins those
//! allele forms, both shuffle directions, axis parity (g./c./n./r., and the
//! protein axis's *different* history), and the hygiene-mode / projector
//! questions the investigation was asked to settle.
//!
//! # #51 — root cause is NOT nucleotide-alphabet inclusion
//!
//! `background/standards.md:19-37` lists the DNA symbol table; the 15 real
//! symbols (`A,C,G,T,B,D,H,K,M,N,R,S,V,W,Y`) sit at rows 21-35, and `X` and `-`
//! sit at rows 36-37 marked `†`. `standards.md:39` footnotes: "† used in
//! alignment only." So neither `X` nor `-` is a base a description may state.
//!
//! Reading `src/hgvs/parser/edit.rs` shows `X` is not admitted through the
//! nucleotide class at all: `is_iupac_base` (`:27-62`), the function
//! `parse_inserted_sequence` (`:317-473`) uses for plain sequence content,
//! does not include `X` (or `-`). `X` "survives" only because that function
//! also has a fallback arm for **named/mobile-element** spellings (e.g.
//! `AluYb8`, `LINE1`) at `:450-468`: a lone uppercase letter with nothing
//! alphanumeric following matches it, and becomes
//! `InsertedSequence::Named("X")`, which `Display`s verbatim. This is not
//! `X`-specific — any lone uppercase letter that is not otherwise consumed
//! hits the same arm — which one test below pins directly.
//!
//! **That describes one of the two arms that produce `Named`, and it
//! understates the reach — corrected 2026-08-11.** `parse_inserted_sequence`
//! dispatches on the *first* byte, so `X…` takes the arm above, but an
//! **IUPAC** first byte takes the `is_iupac_base` arm, which walks the whole
//! run and sets `has_non_iupac` on any uppercase byte that is not an IUPAC
//! base. So a run that is literal apart from one stray symbol is reclassified
//! *wholesale*: `g.10delinsACGTX` is `Named("ACGTX")`, not `Literal(ACGT)` plus
//! a rejected `X`, and the same holds with the stray at the start or in the
//! middle. No corpus row exercises that shape, so the prohibited-row counts
//! quoted for `standards.md:39` are lone-`X` counts. Measured and pinned in
//! `tests/it/issue_1627_named_element_alphabet_reach.rs`.
//!
//! `-` fails every
//! arm (not `[`, `(`, digit/`*`, IUPAC, or uppercase) and is left unconsumed,
//! which is why `c.10delins-` is rejected: not because `-` is recognized as a
//! gap symbol and refused, but because it is simply never consumed.
//!
//! No axis-specific alphabet exists: `parse_na_edit` (c./n./r.) and
//! `parse_genome_na_edit` (g./m./o.) both dispatch to one shared
//! `parse_na_edit_inner` (`edit.rs:1716-1761`), so the finding is identical on
//! every DNA/RNA axis. Protein is different code entirely: single-letter `X`
//! is the parser-native `Xaa` ("unknown residue", `location.rs:684`), and the
//! *historical stop-codon* spelling `Trp26X` is caught by a separate hygiene
//! pass (`ErrorType::DeprecatedStopCodonX`, W3008) that strict mode rejects —
//! see `checklist.md:63` and `standards.md:213`. Protein therefore already
//! has the guard the DNA/RNA axes lack; a fix must not accidentally treat
//! protein `X` as the same defect.
//!
//! # #52 — root cause is an absent check, not a lenient one
//!
//! `checklist.md:20`: "NM reference sequences cover mature transcripts and
//! **do not contain** intron and gene flanking sequences, and can only be
//! used to describe variants in introns using a `c.` prefix when a genomic
//! reference sequence is given on which the coding DNA reference sequence is
//! annotated, e.g., `NC_000023.10(NM_004006.2):c.94-2A>G` or
//! `LRG_199t1:c.94-2A>G`." A bare `NM_...:c.123+65del` is therefore an
//! **invalid input**, not merely an unusual one.
//!
//! **Scope this claim to the PARSE stage — normalize is a different answer.**
//! `parse_hgvs_with_config` accepts the bare form in all three `ErrorMode`s,
//! which is what
//! `no_error_mode_rejects_the_bare_intronic_allele_at_parse_time` below pins.
//! Strict-mode **normalize** does refuse it, reporting
//! `IntronicOnBareTranscript` / `W4007` — see
//! `corpus_prohibited_inputs.rs::a_bare_transcript_intronic_position_is_refused_in_strict_only`
//! and the `bare-transcript-intronic-position` ruling in
//! `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`, which
//! records lenient-accepts / strict-refuses as ferro's current behaviour. An
//! earlier form of this paragraph generalized the parse-stage result into a
//! claim about the whole catalog and concluded "there is nothing to turn on",
//! which contradicts both of those artifacts.
//!
//! No `ErrorMode` variant makes the PARSER reject it (checked the full
//! `ErrorType` catalog in `src/error_handling/types.rs`); a dead, uncalled
//! module
//! (`src/hgvs/validation.rs::validate_variant_consistency`) gestures at an
//! "intronic positions in exon-only contexts" check but has zero callers
//! anywhere in `src/` or `tests/` — it is not wired into parse, normalize, or
//! project. `n.` (`NR_`) shares the identical position-parsing path as `c.`,
//! so the gap is not `c.`-specific either.
//!
//! `VariantProjector::project_to_genomic` *does* decline on a bare accession
//! ("no parent reference (genomic_context) ... see #327") — but that decline
//! is a **general** missing-genomic-parent constraint, not an intronic-offset
//! check: the already-committed test
//! `tests/it/issue_851_compound_utr_projection.rs::allele_with_unparented_member_declines_whole`
//! gets the identical decline message for a purely **exonic** bare-`NM_`
//! allele (`c.[1A>T;4C>A]`, no offset at all). So the projector's refusal is
//! orthogonal to this defect and is not duplicated here as a new test — it
//! would prove nothing not already proven by that committed test, and would
//! not touch the intronic-ness question at all.

use ferro_hgvs::conformance::spec_corpus::{corpus_cores, Frame, RefShape, DENSE_CORE_LEN};
use ferro_hgvs::error_handling::{ErrorConfig, ErrorType};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

/// The `AT`-alphabet core the corpus draws first, matching the choice
/// `spec_corpus_regressions.rs` makes for the same reason: a two-letter
/// alphabet makes long runs, which is where junction-crossing shifts live.
fn at_core() -> String {
    corpus_cores(1, DENSE_CORE_LEN).remove(0)
}

/// Normalize `input` against `frame`, in the given direction.
fn normalize(frame: &Frame, input: &str, direction: ShuffleDirection) -> String {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
    let normalizer = Normalizer::with_config(
        frame.provider().clone(),
        NormalizeConfig::default().with_direction(direction),
    );
    normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("{input} must normalize: {e}"))
        .to_string()
}

/// Normalize in the default (3') direction.
fn normalize_3prime(frame: &Frame, input: &str) -> String {
    normalize(frame, input, ShuffleDirection::ThreePrime)
}

/// Normalize in the 5' direction.
fn normalize_5prime(frame: &Frame, input: &str) -> String {
    normalize(frame, input, ShuffleDirection::FivePrime)
}

/// The message normalization refuses `input` with, for a shape that must be
/// refused rather than emitted.
///
/// **Panics if `input` is accepted**, naming the string that was emitted. A
/// helper that returned the output on success would let a regression back to
/// acceptance fail as a bare `contains(..)` mismatch, which reads as a changed
/// message rather than as a lost refusal.
fn refusal(frame: &Frame, input: &str, direction: ShuffleDirection) -> String {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
    match Normalizer::with_config(
        frame.provider().clone(),
        NormalizeConfig::default().with_direction(direction),
    )
    .normalize(&variant)
    {
        Ok(output) => panic!("{input} must be REFUSED, but normalized to {output}"),
        Err(e) => e.to_string(),
    }
}

// ---------------------------------------------------------------------------
// #51 — `X` inside an allele member
// ---------------------------------------------------------------------------

/// **Question.** Does the alignment-only `X` finding hold when `X` sits in one
/// member of a multi-member allele, after the 3' rule has reordered the members?
///
/// **It used to. It is now refused, in every mode (#1627).** These two tests
/// were PINNED FINDINGS recording that `NM_TEST.1:c.[10delinsX;9del]`
/// normalized to `NM_TEST.1:c.[9del;10delinsX]` — the members reordered by
/// position (the 3' rule doing its job, not the defect) with `X` surviving
/// inside the moved member. `standards.md:39` footnotes `X` (row 36) "used in
/// alignment only", so it is not a base a description may state, in an allele
/// member or otherwise; the decided
/// `rulings[alignment-only-symbol-in-a-description]` settles that, and
/// `rulings[absolute-prohibition-enforcement-stage]` puts the refusal at strict
/// parse and at lenient/silent normalize.
///
/// **The allele form is still pinned separately, for the reason it always was:**
/// a hygiene check that runs once per description rather than once per member
/// has its own signature, and a fix that reached the standalone form only would
/// pass the sibling test in `spec_corpus_regressions.rs` while failing here.
/// The direction is kept as a parameter because validity cannot depend on which
/// way an ambiguous edit is shifted — and because the 5' row used to have its
/// own recorded output (`[1del;10delinsX]`, the `del` member sliding left on
/// this repetitive `AT` core), so the two directions genuinely differed in
/// what they emitted while agreeing that they emitted something.
#[test]
fn an_allele_member_x_is_refused_after_a_3prime_reorder() {
    let core = at_core();
    let coding = Frame::build(RefShape::CodingSingleExon, &core);
    assert!(
        refusal(
            &coding,
            "NM_TEST.1:c.[10delinsX;9del]",
            ShuffleDirection::ThreePrime
        )
        .contains("W3028"),
        "standards.md:39 — the allele form must be refused like the standalone one"
    );
}

/// The 5' half of the pair directly above. See its doc comment.
#[test]
fn an_allele_member_x_is_refused_after_a_5prime_reorder() {
    let core = at_core();
    let coding = Frame::build(RefShape::CodingSingleExon, &core);
    assert!(
        refusal(
            &coding,
            "NM_TEST.1:c.[10delinsX;9del]",
            ShuffleDirection::FivePrime
        )
        .contains("W3028"),
        "standards.md:39 — validity does not depend on the shuffle direction"
    );
}

/// **Question.** Is `X`'s acceptance an `X`-specific carve-out, or does *any*
/// lone uppercase letter hit the same code path?
///
/// **Any lone uppercase letter reaches the code path**, and this test still
/// passes unchanged after #1627 — which is the point of keeping it.
/// `parse_inserted_sequence`'s named/mobile-element fallback arm matches any
/// bare uppercase letter with nothing alphanumeric following; it has no
/// knowledge of `X`. `Z` is not in the IUPAC table either
/// (`standards.md:19-37` lists none of `E,F,I,J,L,O,P,Q,U,X,Z`).
///
/// **CORRECTED, 2026-08-11.** This doc used to tell a future fixer that "the
/// fix cannot be 'special-case reject the letter X'". The decided
/// `rulings[alignment-only-symbol-in-a-description]` says the opposite: it is
/// scoped to `standards.md:39`'s **two daggered rows**, `X` and `-`, and #1627
/// implements exactly that. `Z` is a different question and **nobody has ruled
/// on it** — refusing it would mean refusing every `InsertedSequence::Named`
/// that is not a recognised mobile element, which reaches `AluYb8` and `LINE1`
/// and is the un-adjudicated territory `DNA/complex.md:169` gestures at ("No,
/// not really, it is not exact", about `insL1.603bp`). So this test is now an
/// **adjudicated-scope guard**: it fails if a later change widens the `X` rule
/// into an alphabet rule without a ruling to stand on.
///
/// **If you are here because this test went red, do not re-bless it.** Going
/// red means the refusal has widened past `standards.md:39`'s two daggered
/// rows, and the prerequisite is a **new `decided` record in
/// `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`** — not a
/// judgement call in a PR. The cost the record has to price is concrete and is
/// the reason this is not a free tightening: the only predicate that catches
/// `Z` is "the insert is an `InsertedSequence::Named` that is not a recognised
/// mobile element", and ferro has no such recogniser, so the same widening
/// **refuses `AluYb8`, `LINE1`, `L1` and `Alu`** — the shapes the named-element
/// arm exists to serve, and which
/// `genuine_mobile_element_names_state_no_alignment_symbol` pins as untouched.
#[test]
fn a_lone_non_daggered_uppercase_letter_is_still_accepted() {
    let core = at_core();
    let coding = Frame::build(RefShape::CodingSingleExon, &core);
    assert_eq!(
        normalize_3prime(&coding, "NM_TEST.1:c.[10delinsZ;9del]"),
        "NM_TEST.1:c.[9del;10delinsZ]",
        "SCOPE GUARD — `Z` is not an IUPAC DNA symbol (standards.md:19-37) but it is not one \
         of :39's two DAGGERED symbols either, and no ruling covers it. If this fails, the \
         #1627 refusal has widened from `standards.md:39` into a general alphabet rule; that \
         needs an adjudication first, because the same widening reaches `AluYb8` and `LINE1`."
    );
}

/// **Question.** Does the allele form relax the `-` rejection that the
/// standalone form already pins?
///
/// **No.** `standards.md:39` footnotes `-` (row 37) identically to `X` — "†
/// used in alignment only" — but ferro's parser rejects it (the byte is never
/// consumed by any arm of `parse_inserted_sequence`, so the parser's
/// require-full-consumption check fails) both standalone and inside an
/// allele. Kept here as the allele-level positive control: the allele bracket
/// syntax does not bypass a rejection that already holds per-member.
#[test]
fn the_alignment_only_gap_symbol_is_still_rejected_inside_an_allele() {
    assert!(
        parse_hgvs("NM_TEST.1:c.[10delins-;9del]").is_err(),
        "`-` is rejected inside an allele exactly as it is standalone, though standards.md:39 \
         footnotes both `X` and `-` identically"
    );
}

/// **Question.** Does `X`'s treatment differ by axis — genomic (`g.`), coding
/// (`c.`), non-coding (`n.`), or RNA (`r.`)?
///
/// **No — identical on every DNA/RNA axis, before and after #1627.**
/// `parse_na_edit` (c./n./r.) and `parse_genome_na_edit` (g./m./o.) both
/// dispatch to the single shared `parse_na_edit_inner`, so there is no
/// axis-specific nucleotide alphabet; the refusal is applied on the AST and is
/// likewise axis-independent. The parity claim is what this test is for, and it
/// is worth just as much now that the parity is "all four refuse" as it was
/// when it was "all four accept" — an axis that quietly kept accepting would be
/// a partial fix that no single-axis test could see.
#[test]
fn dna_and_rna_axes_refuse_x_identically_inside_an_allele() {
    let core = at_core();

    let genomic = Frame::build(RefShape::Genomic, &core);
    let coding = Frame::build(RefShape::CodingSingleExon, &core);
    let noncoding = Frame::build(RefShape::NonCodingMultiExon(Strand::Plus), &core);

    for (frame, input, axis) in [
        (&genomic, "NC_TEST.1:g.[266delinsX;264del]", "g."),
        (&coding, "NM_TEST.1:c.[10delinsX;9del]", "c."),
        (&noncoding, "NR_TEST.1:n.[10delinsX;9del]", "n."),
        (&coding, "NM_TEST.1:r.[10delinsX;9del]", "r."),
    ] {
        assert!(
            refusal(frame, input, ShuffleDirection::ThreePrime).contains("W3028"),
            "the {axis} axis must refuse X exactly as the others do"
        );
    }

    // Lowercase `x` is refused on every axis too, but by the GRAMMAR and for an
    // unrelated reason — the fallback arm requires an ASCII-uppercase first
    // byte. Pinned so the two routes are not conflated: `standards.md:39`'s
    // table is uppercase, and #1627 deliberately does not rule on lowercase.
    assert!(
        parse_hgvs("NM_TEST.1:r.[10delinsx;9del]").is_err(),
        "lowercase `x` (the RNA-case-convention spelling) is refused by the grammar"
    );
}

/// **Question.** Does protein's `X` behave the same way — silently accepted
/// as sequence content — when it appears inside a `p.[...]` allele?
///
/// **No — protein already has a guard the DNA/RNA axes lack.** Single-letter
/// `X` is the parser-native `Xaa` ("unknown amino acid", `location.rs:684`),
/// a different symbol than a stop codon. The *historical stop-codon* spelling
/// (`Trp26X`) is caught by a dedicated hygiene pass
/// (`ErrorType::DeprecatedStopCodonX`, W3008): `checklist.md:63` — "**`Ter`
/// or `*` should be used** to indicate a translation stop codon; the `X`
/// should not be used" — and strict mode rejects it, inside an allele exactly
/// as standalone. This is deliberately NOT filed as the same defect as #51:
/// it is cited here only to keep the axes separated, per the task's
/// instruction not to conflate DNA-axis `X` (unguarded) with protein-axis `X`
/// (guarded, and guarded for an unrelated historical reason).
#[test]
fn protein_x_inside_an_allele_is_a_different_and_already_guarded_symbol() {
    let strict_result =
        parse_hgvs_with_config("NP_TEST.1:p.[Trp26X;Lys52Asn]", ErrorConfig::strict());
    let err = strict_result
        .expect_err("strict mode must reject the deprecated stop-codon `X` inside an allele");
    assert!(
        err.to_string().contains("Deprecated protein notation 'X'")
            || err.to_string().to_lowercase().contains("deprecated"),
        "expected a deprecated-notation rejection, got: {err}"
    );

    let lenient_result =
        parse_hgvs_with_config("NP_TEST.1:p.[Trp26X;Lys52Asn]", ErrorConfig::lenient())
            .unwrap_or_else(|e| panic!("lenient mode must accept and rewrite: {e}"));
    assert_eq!(
        lenient_result.result.to_string(),
        "NP_TEST.1:p.[Trp26Ter;Lys52Asn]",
        "lenient mode rewrites the deprecated `X` to `Ter` inside the allele, member order \
         and the sibling member untouched"
    );
    assert!(
        lenient_result
            .warnings
            .iter()
            .any(|w| w.error_type == ErrorType::DeprecatedStopCodonX),
        "expected a W3008 DeprecatedStopCodonX warning; got: {:?}",
        lenient_result.warnings
    );
}

// ---------------------------------------------------------------------------
// #52 — an intronic position inside an allele, on a bare transcript accession
// ---------------------------------------------------------------------------

/// **Question.** Does the bare-accession intronic finding
/// (`a_bare_transcript_accession_accepts_an_intronic_position` in
/// `spec_corpus_regressions.rs`) hold when the intronic member sits inside a
/// multi-member allele?
///
/// **Yes, unchanged.** `NM_TEST.1:c.[20+2del;24del]` normalizes (3') to
/// itself: both members are already in position order (20+2 < 24 on this
/// plus-strand frame, once the intron length is accounted for), so nothing
/// about this test is testing normalization *drift* into an intron (that is
/// the different, already-pinned `general.md:44` shape in
/// `a_minus_strand_junction_shift_leaves_the_transcript`) — the intronic
/// position is exactly what the INPUT named. `checklist.md:20`: "NM reference
/// sequences cover mature transcripts and **do not contain** intron and gene
/// flanking sequences, and can only be used to describe variants in introns
/// using a `c.` prefix when a genomic reference sequence is given ...". A
/// bare `NM_...:c.20+2del` is an invalid input, accepted anyway.
#[test]
fn an_allele_member_keeps_an_intronic_offset_on_a_bare_nm_accession() {
    let core = at_core();
    let frame = Frame::build(RefShape::CodingMultiExon(Strand::Plus), &core);
    assert_eq!(
        normalize_3prime(&frame, "NM_TEST.1:c.[20+2del;24del]"),
        "NM_TEST.1:c.[20+2del;24del]",
        "PINNED FINDING — checklist.md:20 permits an intronic `c.` position only with a genomic \
         wrapper; ferro accepts the bare form both standalone and, as pinned here, inside an \
         allele. Correct behaviour: refuse, or demand the NC_(NM_) wrapper form."
    );
}

/// **Question.** Same allele, 5' direction — does the corpus's second count
/// (both directions, per the #51 precedent) hold here too?
///
/// **Yes.** No junction- or CDS-boundary clamp applies to either member here,
/// so the 5' direction leaves the allele exactly as the 3' direction does,
/// intronic offset and all.
#[test]
fn an_allele_member_keeps_an_intronic_offset_in_the_5prime_direction_too() {
    let core = at_core();
    let frame = Frame::build(RefShape::CodingMultiExon(Strand::Plus), &core);
    assert_eq!(
        normalize_5prime(&frame, "NM_TEST.1:c.[20+2del;24del]"),
        "NM_TEST.1:c.[20+2del;24del]",
        "PINNED FINDING — same defect, 5' direction."
    );
}

/// **Question.** Is the gap `n.`-specific, or does `NR_` (non-coding
/// transcript) admit the identical bare-accession intronic allele?
///
/// **Identical.** `n.` position parsing shares the same offset machinery as
/// `c.` (`edit.rs:1716-1761`); there is no separate `n.`-only intronic gate.
#[test]
fn an_allele_member_keeps_an_intronic_offset_on_a_bare_nr_accession_too() {
    let core = at_core();
    let frame = Frame::build(RefShape::NonCodingMultiExon(Strand::Plus), &core);
    assert_eq!(
        normalize_3prime(&frame, "NR_TEST.1:n.[20+2del;24del]"),
        "NR_TEST.1:n.[20+2del;24del]",
        "PINNED FINDING — checklist.md:20 talks about NM_ specifically, but NR_ (non-coding \
         mature transcript) is exactly as unable to carry intron sequence, and ferro's `n.` \
         offset parser is the same code path as `c.`'s. Correct behaviour: refuse, same as NM_."
    );
}

/// **Question.** Does the genomic-wrapped form the checklist actually asks
/// for (`NC_...(NM_...):c...`) work at the allele level, so the eventual fix
/// has a positive control that is not just the single-member case?
///
/// **Yes.** `NC_SYNTH.1(NM_TEST.1):c.[20+2del;24del]` parses successfully —
/// the wrapper syntax (`try_compound_suffix`, `accession.rs:477-489`) is
/// generic over every accession family and is unaffected by how many members
/// the allele has. This is the positive control: whatever fix lands for #52
/// must keep this accepted while making the bare form (above) refuse.
#[test]
fn the_genomic_wrapped_allele_form_is_the_positive_control() {
    assert!(
        parse_hgvs("NC_SYNTH.1(NM_TEST.1):c.[20+2del;24del]").is_ok(),
        "the checklist.md:20-sanctioned wrapper form must remain accepted at the allele level; \
         a fix for the bare form must not collaterally break this"
    );
}

/// **Question.** Does any `ErrorMode` (strict/lenient/silent) reject the
/// bare-accession intronic allele **at parse time**?
///
/// **Strict does, and only strict (#1630). This was a PINNED FINDING and is now
/// an adjudicated-correct guard.** It used to assert that all three modes parsed
/// it identically, on the ground that the `ErrorType` catalog had no variant for
/// "intronic offset without a genomic wrapper" and so the parser had no check to
/// turn on. `ErrorType::IntronicOnBareTranscript` / `W4007` is now that variant,
/// applied at the parse seam by `apply_bare_transcript_intron_rule`.
///
/// **The mode split is the ruling's, and it is unchanged.**
/// `rulings[absolute-prohibition-enforcement-stage]` puts an input-conformance
/// check at parse, gated by mode: strict validates input conformance and so
/// fails there; lenient does not validate it and accepts; silent is lenient
/// without messages. `rulings[bare-transcript-intronic-position]` decided the
/// split itself — `checklist.md:20` is conditional in form — and #1630 moved
/// only the *stage* of the strict arm.
///
/// **The allele form is the point of this test**, and it is not redundant with
/// the single-member guard: a hygiene check that runs once per description
/// rather than once per member stops firing exactly here. Its positive control
/// is [`the_genomic_wrapped_allele_form_is_the_positive_control`] directly
/// above, which must stay accepted.
///
/// Strict-mode *normalize* still refuses the same input too, with the same
/// `W4007` finding carrying the `EINTRONIC` tag — pinned by
/// `corpus_prohibited_inputs.rs::a_bare_transcript_intronic_position_is_refused_in_strict_only`.
/// The two rungs answer for different callers; see
/// `crate::hgvs::bare_transcript_introns`.
#[test]
fn only_strict_mode_rejects_the_bare_intronic_allele_at_parse_time() {
    let input = "NM_TEST.1:c.[20+2del;24del]";

    let refusal = parse_hgvs_with_config(input, ErrorConfig::strict())
        .expect_err(
            "ADJUDICATED CORRECT, REGRESSED: strict validates input conformance, so \
             checklist.md:20 is refused at PARSE",
        )
        .to_string();
    assert!(
        refusal.contains("W4007"),
        "the strict parse refusal must name checklist.md:20's finding; got: {refusal}"
    );

    for (mode_name, config) in [
        ("lenient", ErrorConfig::lenient()),
        ("silent", ErrorConfig::silent()),
    ] {
        let result = parse_hgvs_with_config(input, config);
        assert!(
            result.is_ok(),
            "{mode_name} mode does not validate input conformance, so it must accept \
             the bare-accession intronic allele; got: {result:?}",
        );
    }
}
