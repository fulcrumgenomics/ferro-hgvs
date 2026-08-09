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
//! hits the same arm — which one test below pins directly. `-` fails every
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

// ---------------------------------------------------------------------------
// #51 — `X` inside an allele member
// ---------------------------------------------------------------------------

/// **Question.** Does the alignment-only `X` finding
/// (`an_alignment_only_x_base_is_accepted_as_a_nucleotide` in
/// `spec_corpus_regressions.rs`) hold when `X` sits in one member of a
/// multi-member allele, after the 3' rule has reordered the members?
///
/// **Yes.** `NM_TEST.1:c.[10delinsX;9del]` normalizes (3') to
/// `NM_TEST.1:c.[9del;10delinsX]` — the members are reordered by position (the
/// 3' rule doing its job, not the defect), and `X` survives inside the moved
/// member exactly as it did standalone. `standards.md:19-37` lists the DNA
/// symbol table; `standards.md:39` footnotes `X` (row 36) "used in alignment
/// only" — so it is not a base a description may state, in an allele member
/// or otherwise.
#[test]
fn an_allele_member_x_survives_a_3prime_reorder() {
    let core = at_core();
    let coding = Frame::build(RefShape::CodingSingleExon, &core);
    assert_eq!(
        normalize_3prime(&coding, "NM_TEST.1:c.[10delinsX;9del]"),
        "NM_TEST.1:c.[9del;10delinsX]",
        "PINNED FINDING — standards.md:39 footnotes `X` as 'used in alignment only'; ferro \
         accepts it as sequence content on the DNA axis (root cause: the mobile-element-name \
         parser fallback, not nucleotide-alphabet inclusion — see module doc comment), and the \
         allele form is no different from the standalone form."
    );
}

/// **Question.** Does the same allele survive `X` in the 5' direction too?
///
/// **Yes**, though the `del` member additionally slides left on this
/// repetitive `AT`-alphabet core (5' shifting is free to walk a plain
/// deletion across identical bases, which is unrelated to the defect): the
/// output is `NM_TEST.1:c.[1del;10delinsX]` rather than the 3' direction's
/// `[9del;10delinsX]`. What matters for this test is the part the 5' shift
/// cannot touch — `X` still survives, unchanged, in the `delinsX` member.
#[test]
fn an_allele_member_x_survives_a_5prime_reorder() {
    let core = at_core();
    let coding = Frame::build(RefShape::CodingSingleExon, &core);
    assert_eq!(
        normalize_5prime(&coding, "NM_TEST.1:c.[10delinsX;9del]"),
        "NM_TEST.1:c.[1del;10delinsX]",
        "PINNED FINDING (measured) — same defect, 5' direction; the corpus counts this class in \
         both directions (24 outputs total). The `del` member also slides to `1del` on this \
         repetitive core (an unrelated 5'-shift effect), but `X` survives in `delinsX` either way."
    );
}

/// **Question.** Is `X`'s acceptance an `X`-specific carve-out, or does *any*
/// lone uppercase letter hit the same code path?
///
/// **Any lone uppercase letter.** `parse_inserted_sequence`'s named/mobile-
/// element fallback arm (`edit.rs:450-468`) matches any bare uppercase letter
/// with nothing alphanumeric following — it has no knowledge of `X`
/// specifically. `Z` is not in the IUPAC table either (`standards.md:19-37`
/// lists none of `E,F,I,J,L,O,P,Q,U,X,Z`), and it survives identically. This
/// matters for whoever fixes #51: the fix cannot be "special-case reject the
/// letter X", because the same fallback would still accept every other
/// non-IUPAC uppercase letter as a fake "named element" of length one.
#[test]
fn any_lone_uppercase_non_iupac_letter_survives_the_same_way_as_x() {
    let core = at_core();
    let coding = Frame::build(RefShape::CodingSingleExon, &core);
    assert_eq!(
        normalize_3prime(&coding, "NM_TEST.1:c.[10delinsZ;9del]"),
        "NM_TEST.1:c.[9del;10delinsZ]",
        "PINNED FINDING — `Z` is not an IUPAC DNA symbol either (standards.md:19-37) and hits \
         the identical named-mobile-element fallback arm as `X`; the defect is 'any lone \
         uppercase letter is accepted as sequence content', not 'X is special-cased'."
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

/// **Question.** Does `X`'s acceptance-as-sequence-content differ by axis —
/// genomic (`g.`), coding (`c.`), non-coding (`n.`), or RNA (`r.`)?
///
/// **No — identical on every DNA/RNA axis.** `parse_na_edit` (c./n./r.) and
/// `parse_genome_na_edit` (g./m./o.) both dispatch to the single shared
/// `parse_na_edit_inner` (`edit.rs:1716-1761`); there is no axis-specific
/// nucleotide alphabet. Each case below is an allele with the reorder the 3'
/// rule performs, so this also confirms axis-parity holds at the allele level.
#[test]
fn dna_and_rna_axes_accept_x_identically_inside_an_allele() {
    let core = at_core();

    let genomic = Frame::build(RefShape::Genomic, &core);
    let coding = Frame::build(RefShape::CodingSingleExon, &core);
    let noncoding = Frame::build(RefShape::NonCodingMultiExon(Strand::Plus), &core);

    // g. — genomic axis.
    assert!(
        normalize_3prime(&genomic, "NC_TEST.1:g.[266delinsX;264del]").contains('X'),
        "genomic axis must retain X exactly as the coding axis does"
    );

    // c. — coding axis (the exemplar axis, re-asserted here for the loop's sake).
    assert!(
        normalize_3prime(&coding, "NM_TEST.1:c.[10delinsX;9del]").contains('X'),
        "coding axis must retain X"
    );

    // n. — non-coding axis.
    let noncoding_output = normalize_3prime(&noncoding, "NR_TEST.1:n.[10delinsX;9del]");
    assert!(
        noncoding_output.contains('X'),
        "non-coding axis must retain X, got {noncoding_output}"
    );

    // r. — RNA axis, same NM_ accession as coding but the `r.` prefix.
    assert!(
        parse_hgvs("NM_TEST.1:r.[10delinsx;9del]").is_err(),
        "lowercase `x` (the RNA-case-convention spelling) is rejected on every axis — the \
         mobile-element fallback arm requires an ASCII-uppercase letter"
    );
    let rna_variant = parse_hgvs("NM_TEST.1:r.[10delinsX;9del]")
        .unwrap_or_else(|e| panic!("uppercase X must parse on the r. axis too: {e}"));
    assert!(
        rna_variant.to_string().contains('X'),
        "uppercase X survives on the r. axis exactly as on g./c./n., because the parser has no \
         axis-specific nucleotide alphabet at all"
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

/// **Question.** Does any `ErrorMode` (strict/lenient/silent) already reject
/// the bare-accession intronic allele **at parse time**, making this a
/// hygiene-mode gap in the parser rather than an absent check there?
///
/// **No — all three modes parse it identically.** The full `ErrorType` catalog
/// (`src/error_handling/types.rs`) has no variant for "intronic offset without
/// a genomic wrapper" (the nearest relative, `DeprecatedIvsNotation`, covers
/// only the retracted `IVS` spelling). This pins that the parser has no check
/// to turn on.
///
/// **The scope is the parser, and only the parser.** Strict-mode *normalize*
/// does refuse the same input, with `IntronicOnBareTranscript` / `W4007`; that
/// is pinned by
/// `corpus_prohibited_inputs.rs::a_bare_transcript_intronic_position_is_refused_in_strict_only`
/// and recorded in the `bare-transcript-intronic-position` ruling. The test is
/// named for the stage so the two cannot be read as contradicting each other.
#[test]
fn no_error_mode_rejects_the_bare_intronic_allele_at_parse_time() {
    for (mode_name, config) in [
        ("strict", ErrorConfig::strict()),
        ("lenient", ErrorConfig::lenient()),
        ("silent", ErrorConfig::silent()),
    ] {
        let result = parse_hgvs_with_config("NM_TEST.1:c.[20+2del;24del]", config);
        assert!(
            result.is_ok(),
            "PINNED FINDING — {mode_name} mode accepts the bare-accession intronic allele \
             exactly as the other two do; no hygiene mode implements checklist.md:20's \
             wrapper requirement, got: {result:?}",
        );
    }
}
