//! HGVS parser using nom
//!
//! This module provides a complete parser for HGVS variant nomenclature.
//!
//! # Error Handling Modes
//!
//! The parser supports three error handling modes:
//!
//! - **Strict**: Reject all non-standard input
//! - **Lenient**: Auto-correct common errors with warnings
//! - **Silent**: Auto-correct common errors without warnings
//!
//! A mode is applied only by [`parse_hgvs_with_config`]; strict is
//! [`ErrorConfig::default`], so `parse_hgvs_with_config(input,
//! ErrorConfig::default())` is the strict entry. [`parse_hgvs`] applies **no**
//! mode at all — see its own documentation (#1632).

pub mod accession;
pub mod edit;
pub mod fast_path;
pub mod position;
pub mod variant;

use crate::error::FerroError;
use crate::error_handling::{
    CorrectionWarning, ErrorConfig, ErrorType, InputPreprocessor, ParseResultWithWarnings,
    ResolvedAction,
};
use crate::hgvs::alignment_symbols::alignment_only_symbol;
use crate::hgvs::noncoding_zones::{noncoding_zone_marker, NonCodingZone};
use crate::hgvs::HgvsVariant;

/// Parse an HGVS string into a variant
///
/// # Error handling: this entry applies the grammar only
///
/// `parse_hgvs` applies **no [`ErrorConfig`]**. It trims surrounding
/// whitespace, runs the fast path, and falls back to the grammar; the
/// input-hygiene ladder — both the repairs
/// lenient mode performs and the parse-stage refusals that make strict strict —
/// runs only in [`parse_hgvs_with_config`]. So this is not "strict mode": it is
/// a third behaviour that is neither mode, and it accepts inputs strict rejects.
/// For example `NM_024312.4:r.-6_-3g[6]` states a base label the range already
/// determines; `parse_hgvs` returns it unchanged, while
/// `parse_hgvs_with_config(_, ErrorConfig::strict())` refuses it by name.
///
/// Note the scope of that claim: it is the ladder that runs at the *parse*
/// stage. Strict mode also refuses at **normalize**, so switching to
/// [`parse_hgvs_with_config`] buys the parse-stage half only.
///
/// A bare-transcript intronic position (`W4007`) used to be the example here,
/// on the ground that "no parse-stage entry rejects it in any mode". That
/// stopped being true with #1630, which put `checklist.md:20`'s strict-mode
/// refusal at parse per
/// `rulings[absolute-prohibition-enforcement-stage]` — so it is now an example
/// of the *opposite*: a shape this entry accepts and the strict mode-aware
/// entry refuses. It still has a normalize-stage rung too, for callers that
/// arrive through this door.
///
/// Call [`parse_hgvs_with_config`] when the mode matters — validating untrusted
/// input, or matching what the CLI does with the same string. `ferro normalize`,
/// `ferro parse` and `ferro project` each take `--error-mode`, and all three
/// default to `strict`.
///
/// Do not generalise that default to the whole crate: it runs the other way one
/// layer down. [`ErrorConfig::default`] is strict, but `NormalizeConfig::default`
/// sets `ErrorConfig::lenient` for backwards compatibility, so a
/// `Normalizer::new(provider)` normalizes in **lenient** mode. Two `default()`s
/// pointing opposite ways is the reason to name the mode you want rather than
/// infer it.
///
/// This doc previously claimed strict handling "by default" (#1632). Whether
/// the function should be routed through `ErrorConfig::strict()` to make that
/// claim true is open; it would newly refuse inputs the crate's headline entry
/// point accepts today. `tests/it/issue_1632_parse_entry_applies_no_mode.rs`
/// pins the current behaviour so the two cannot drift apart again.
///
/// # Performance
///
/// Common substitutions and plain deletions/duplications (RefSeq / Ensembl /
/// LRG / Assembly `g.`/`c.` — ~72% of real ClinVar) take a specialized fast
/// path that is ~1.7x
/// faster end-to-end over the full ClinVar corpus; everything else falls back
/// to the full [`variant::parse_variant`] parser. The fast path is
/// observationally identical to the generic parser — both accept the same
/// inputs and produce the same [`HgvsVariant`] (guarded by the differential
/// test in `tests/fast_path_differential.rs`) — so it is transparent to
/// callers. To force the generic parser, call [`variant::parse_variant`]
/// directly.
///
/// # Example
///
/// ```
/// use ferro_hgvs::parse_hgvs;
///
/// let variant = parse_hgvs("NM_000088.3:c.459del").unwrap();
/// println!("Parsed: {}", variant);
/// ```
pub fn parse_hgvs(input: &str) -> Result<HgvsVariant, FerroError> {
    let trimmed = trim_hgvs(input);
    match fast_path::try_fast_path(trimmed) {
        fast_path::FastPathResult::Success(variant) => Ok(variant),
        // Fall back on the original (untrimmed) input so the generic parser
        // sees exactly what it did before the fast path existed.
        fast_path::FastPathResult::Fallback => variant::parse_variant(input),
    }
}

/// Trim surrounding whitespace, skipping the work when there is none.
///
/// `str::trim` is Unicode-aware, so even a string with no surrounding whitespace
/// pays for boundary `char` decoding + the Unicode-whitespace test on every call.
/// Virtually all HGVS strings begin and end with an ASCII non-whitespace byte
/// (`NM…`, `…A>G`, `…del`, `…=`, `…)`), so when both boundary bytes are ASCII and
/// non-whitespace the input cannot have surrounding whitespace (an ASCII byte is
/// never part of a multi-byte UTF-8 whitespace char) and is returned untouched.
/// Any ASCII-whitespace or non-ASCII (`>= 0x80`) boundary defers to `str::trim`,
/// so Unicode-whitespace handling is byte-for-byte unchanged.
#[inline]
fn trim_hgvs(input: &str) -> &str {
    let bytes = input.as_bytes();
    match (bytes.first(), bytes.last()) {
        (Some(&first), Some(&last))
            if first < 0x80
                && last < 0x80
                && !first.is_ascii_whitespace()
                && !last.is_ascii_whitespace() =>
        {
            input
        }
        _ => input.trim(),
    }
}

/// Parse an HGVS string using the fast path for common patterns.
///
/// **Equivalent to [`parse_hgvs`].** The fast path is now the default, so this
/// is a thin synonym retained for backward compatibility; new code should call
/// [`parse_hgvs`]. See [`parse_hgvs`] for the performance characteristics and
/// the identity guarantee versus the generic parser.
///
/// # Example
///
/// ```
/// use ferro_hgvs::parse_hgvs_fast;
///
/// let variant = parse_hgvs_fast("NC_000001.11:g.12345A>G").unwrap();
/// // Complex patterns transparently fall back to the full parser.
/// let variant = parse_hgvs_fast("NM_000088.3:c.100+5A>G").unwrap();
/// ```
#[inline]
pub fn parse_hgvs_fast(input: &str) -> Result<HgvsVariant, FerroError> {
    parse_hgvs(input)
}

/// Parse an HGVS string with configurable error handling.
///
/// This function applies preprocessing based on the error configuration,
/// then parses the (potentially corrected) input.
///
/// # Example
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
/// use ferro_hgvs::error_handling::{ErrorConfig, ErrorMode};
///
/// // Lenient mode: auto-correct common errors with warnings
/// let config = ErrorConfig::lenient();
/// let result = parse_hgvs_with_config("  NM_000088.3:c.459del  ", config);
/// assert!(result.is_ok());
///
/// let parsed = result.unwrap();
/// assert!(parsed.had_corrections()); // Whitespace was trimmed
/// ```
pub fn parse_hgvs_with_config(
    input: &str,
    config: ErrorConfig,
) -> Result<ParseResultWithWarnings<HgvsVariant>, FerroError> {
    // Resolve the bracket-cardinality action before `config` is consumed by the
    // preprocessor; the rule itself is applied post-parse on the AST below.
    let cardinality_action = config.action_for(ErrorType::NonConformantBracketCardinality);
    // Same, for the alignment-only-symbol rule (#1627).
    let alignment_symbol_action = config.action_for(ErrorType::AlignmentOnlySymbolInDescription);
    // Same, for the `n.-N` numbering-zone rule (#1748). Its `n.*N` sibling is
    // NOT here: that one is refused unconditionally, from `parse_variant`.
    let noncoding_zone_action = config.action_for(ErrorType::NonCodingPositionOutsideTranscript);
    // Same, for the genomic-family offset rule (#1628). `g.*10` and a leading
    // `g.-100` are NOT here: neither matches any production, so the grammar
    // refuses both in every mode before this point.
    let genomic_offset_action = config.action_for(ErrorType::GenomicPositionOffset);
    // Same, for the bare-transcript intronic rule (#1630). `checklist.md:20` is
    // a CONDITIONAL clause, so unlike the four above this one has no
    // unconditional arm at all: strict refuses, lenient and silent accept.
    let bare_intron_action = config.action_for(ErrorType::IntronicOnBareTranscript);

    // Create preprocessor and preprocess input
    let preprocessor = InputPreprocessor::new(config);
    let mut preprocess_result = preprocessor.preprocess(input);

    // Check if preprocessing failed
    if let Some(rejection) = preprocess_result.take_rejection_error() {
        return Err(rejection);
    }

    // Parse the preprocessed input
    let mut variant = variant::parse_variant(&preprocess_result.preprocessed)?;

    // Put back the reference run the preprocessor's textual rewrite dropped
    // (#1092). Lenient/silent turn `c.79_80GC>TT` into the canonical text
    // `c.79_80delinsTT` before the grammar runs, so the stated `GC` never
    // reaches the AST and a *false* claim about the reference would disappear
    // with no diagnostic. Restoring it into the provenance field — which
    // `Display` ignores, so the repaired output is byte-identical — routes the
    // claim through the same `validate_reference` check the neighbouring
    // `delGCinsTT` spelling already gets (#486).
    restore_stated_substitution_reference(&mut variant, &preprocess_result.original);

    // Apply the bracket/allele cardinality conformance rule on the parsed AST
    // (#493). This is a structural rule, identical across all coordinate
    // systems, so it is enforced once here rather than per-axis in the parser.
    let mut warnings = preprocess_result.warnings;
    let variant = apply_bracket_cardinality_rule(
        variant,
        cardinality_action,
        &preprocess_result.preprocessed,
        &mut warnings,
    )?;

    // Apply the alignment-only-symbol conformance rule on the parsed AST
    // (#1627). Like the cardinality rule this is structural and axis-independent,
    // so it is enforced once here rather than per-axis in the grammar — which is
    // also what keeps it out of the grammar, where it would refuse in every mode
    // and contradict `rulings[absolute-prohibition-enforcement-stage]`.
    apply_alignment_only_symbol_rule(
        &variant,
        alignment_symbol_action,
        &preprocess_result.preprocessed,
        &mut warnings,
    )?;

    // Apply the `n.-N` numbering-zone conformance rule on the parsed AST
    // (#1748). Structural and axis-keyed like the two above, and mode-gated for
    // the reason `rulings[absolute-prohibition-enforcement-stage]` gives: five
    // real ClinVar rows are `n.-N`, so lenient and silent must keep accepting
    // them. `n.*N` is the other half and is already gone by this point —
    // `parse_variant` refused it in every mode, on zero measured corpus cost.
    apply_noncoding_zone_rule(
        &variant,
        noncoding_zone_action,
        &preprocess_result.preprocessed,
        &mut warnings,
    )?;

    // Apply the genomic-family offset conformance rule on the parsed AST
    // (#1628). Same seam and same shape as the three above: `background/
    // numbering.md:6`/`:8`/`:11` say a `g.`/`o.`/`m.` nucleotide number does not
    // include a `+` or `-`, and there is no exon boundary on such an accession
    // for one to be measured from. Mode-gated here per
    // `rulings[absolute-prohibition-enforcement-stage]`; the OUTPUT half is
    // unconditional and lives in `Normalizer::normalize_core_checked`.
    apply_genomic_offset_rule(
        &variant,
        genomic_offset_action,
        &preprocess_result.preprocessed,
        &mut warnings,
    )?;

    // Apply the bare-transcript intronic conformance rule on the parsed AST
    // (#1630). Same seam as the four above, and the last of the three stage
    // defects `rulings[absolute-prohibition-enforcement-stage]`'s census names:
    // this clause was enforced with the right mode split at the wrong stage, so
    // a caller that parsed strictly and never normalized was told a bare
    // `NM_…:c.20+2del` conforms. The record's ground is general — "whether the
    // INPUT conforms is answered before the input is accepted, not part-way
    // through normalizing it."
    apply_bare_transcript_intron_rule(
        &variant,
        bare_intron_action,
        &preprocess_result.preprocessed,
        &mut warnings,
    )?;

    // Return result with warnings
    Ok(ParseResultWithWarnings::new(
        variant,
        warnings,
        preprocess_result.original,
        preprocess_result.preprocessed,
    ))
}

/// Enforce `docs/recommendations/checklist.md:20` on a parsed variant (#1630).
///
/// An `NM_`/`NR_` (or LRG) transcript reference "can only be used to describe
/// variants in introns … when a genomic reference sequence is given on which
/// the coding DNA reference sequence is annotated" — because a transcript
/// reference *is* the spliced sequence and contains no introns for an offset to
/// name. The full clause reading, its scope, and why there is no repair arm are
/// on [`crate::hgvs::bare_transcript_introns`].
///
/// # This is the CONDITIONAL clause, so there is no unconditional arm
///
/// The four rules above this one each pair a mode-gated *input* check with an
/// unconditional *output* refusal, because their shapes denote no sequence and
/// so rule 1 of the README ruleset bites in every mode. This one does not, and
/// must not: `rulings[bare-transcript-intronic-position]` decided that lenient
/// **accepts** the bare form and that "an input that already names one is still
/// left as authored". A bare intronic description denotes a perfectly good
/// sequence; what it lacks is the reference the clause conditions on. So the
/// only arm is the mode-gated one, and lenient's output is conformant.
///
/// # The normalize rung is kept
///
/// `Normalizer::normalize_core_checked` still refuses this in strict mode, with
/// the `EINTRONIC` tag the Mutalyzer conformance map keys off. That rung
/// answers for a caller who reaches the normalizer by another door — the
/// config-less `parse_hgvs` (#1632), or a lenient parse followed by a strict
/// `Normalizer`. Neither rung subsumes the other; both derive their scope from
/// [`crate::hgvs::bare_transcript_introns::bare_transcript_intronic_axis`], so
/// they cannot come to read the clause differently.
///
/// There is no repair arm, so this returns `Result<(), _>` like
/// [`apply_genomic_offset_rule`] rather than a rewritten variant: naming a
/// genomic parent needs an accession and an exon table the parser does not
/// hold, and the ruling says to leave an authored intronic offset as authored.
fn apply_bare_transcript_intron_rule(
    variant: &HgvsVariant,
    action: ResolvedAction,
    source: &str,
    warnings: &mut Vec<CorrectionWarning>,
) -> Result<(), FerroError> {
    let Some(found) = crate::hgvs::bare_transcript_introns::bare_transcript_intron(variant) else {
        return Ok(());
    };

    match action {
        ResolvedAction::Accept | ResolvedAction::SilentCorrect => Ok(()),
        ResolvedAction::WarnCorrect => {
            warnings.push(CorrectionWarning::new(
                ErrorType::IntronicOnBareTranscript,
                found.refusal(),
                None,
                source.to_string(),
                // No correction is derivable and none is wanted — the ruling
                // leaves an authored intronic offset as authored — so the
                // "corrected" text is the input itself.
                source.to_string(),
            ));
            Ok(())
        }
        ResolvedAction::Reject => Err(FerroError::Parse {
            pos: 0,
            msg: format!(
                "[{}] {}",
                ErrorType::IntronicOnBareTranscript.code(),
                found.refusal()
            ),
            diagnostic: None,
        }),
    }
}

/// Enforce `background/numbering.md:52` on a parsed `n.-N` description (#1748).
///
/// `numbering.md:50`–`:54` enumerates the non-coding DNA axis in full: `:52`
/// numbers it from the first nucleotide to the last, `:53` grants intronic
/// offsets *explicitly*, and `:54` forbids describing "variants in nucleotides
/// beyond the boundaries of a transcript reference sequence". `:53` is what
/// makes `:52`'s silence about `*` and `-` an exclusion rather than terseness —
/// the spec knew how to add a zone here and added exactly one — and
/// `numbering.md:45` records that the proposal to mark non-transcribed
/// nucleotides was rejected.
///
/// **Scoped to `n.` and nothing else.** `c.-1`/`c.*1` are anchored to the CDS
/// and are still inside the transcript, so `:54` does not reach them; and
/// `numbering.md:58` makes an `r.` description's zone set a property of the
/// underlying coding *or* non-coding reference, which this entry holds no
/// provider to resolve. See [`crate::hgvs::noncoding_zones`] for the reasoning.
///
/// # `-N` only. The `*N` half is refused before this runs.
///
/// This is the **mode-gated** arm, per the decided
/// `rulings[absolute-prohibition-enforcement-stage]`: strict validates input
/// conformance and so fails here, at parse; lenient and silent do not, and
/// accept — lenient with a `W4008` warning, silent with none. It stays that way
/// because the shape has measured real-world users: five ClinVar rows across
/// ferro's committed corpora are `n.-N`, three of them RMRP promoter variants.
///
/// `n.*N` has **none** (0 of 103,762 `n.`-axis rows), and is refused
/// unconditionally by `validate_noncoding_downstream_zone` in
/// [`variant::parse_variant`] — which runs inside the `?` on the line that
/// produced this function's `variant`, so a `Downstream` marker can never reach
/// here. The two are deliberately not symmetric; the measurement that split them
/// is on [`crate::hgvs::noncoding_zones`].
///
/// There is no repair arm, so this returns `Result<(), _>` rather than a
/// rewritten variant: re-expressing `n.-5` as an in-transcript coordinate needs
/// the transcript's length, which the parser does not hold, and `:52` denies the
/// zone rather than the spelling — so there is no conformant coordinate on this
/// reference to correct to.
///
/// Both messages come from [`crate::hgvs::noncoding_zones::NonCodingZoneMarker`]
/// so that `:54`'s scoping — it is conditioned on the transcript *being* the
/// reference sequence, which the selector form `NG_…(NR_…)` is not — cannot be
/// right in one and wrong in the other.
fn apply_noncoding_zone_rule(
    variant: &HgvsVariant,
    action: ResolvedAction,
    source: &str,
    warnings: &mut Vec<CorrectionWarning>,
) -> Result<(), FerroError> {
    let Some(found) = noncoding_zone_marker(variant, NonCodingZone::Upstream) else {
        return Ok(());
    };

    match action {
        ResolvedAction::Accept | ResolvedAction::SilentCorrect => Ok(()),
        ResolvedAction::WarnCorrect => {
            warnings.push(CorrectionWarning::new(
                ErrorType::NonCodingPositionOutsideTranscript,
                found.warning(),
                None,
                source.to_string(),
                // No correction is derivable without the transcript's length,
                // and `:52` denies the zone rather than the spelling.
                source.to_string(),
            ));
            Ok(())
        }
        ResolvedAction::Reject => Err(FerroError::Parse {
            pos: 0,
            msg: format!(
                "[{}] {}",
                ErrorType::NonCodingPositionOutsideTranscript.code(),
                found.refusal()
            ),
            diagnostic: None,
        }),
    }
}

/// Enforce `background/numbering.md:6`/`:8`/`:11` on a parsed variant (#1628).
///
/// The three genomic-family axes are numbered in three consecutive bullets and
/// each bullet ends the same way: nucleotide numbers based on a genomic (`:6`),
/// circular (`:8`) or mitochondrial (`:11`) reference sequence "do not include
/// `+`, `-`, `*`, or other prefixes". `checklist.md:16` says it of `g.` in the
/// checklist's register and `checklist.md:45` supplies the shape that actually
/// occurs — a range written with `-` where the separator is `_`.
///
/// **The stage is mode-dependent**, per the decided
/// `rulings[absolute-prohibition-enforcement-stage]`: strict validates input
/// conformance and so fails here, at parse; lenient and silent do not, and
/// accept — lenient with a `W4009` warning, silent without — and then fail at
/// *normalize*, because a genomic accession carries no exon boundary for the
/// offset to be measured from (see `Normalizer::normalize_core_checked`).
///
/// There is no repair arm, so this returns `Result<(), _>` like
/// [`apply_alignment_only_symbol_rule`] rather than a rewritten variant. Both
/// available repairs invent a variant the caller did not describe: dropping the
/// offset answers for `g.266del`, which is the silent flattening #1641 and
/// #1734 were filed to stop, and reading `g.266-268del` as a range guesses a
/// three-base deletion out of a one-base one.
fn apply_genomic_offset_rule(
    variant: &HgvsVariant,
    action: ResolvedAction,
    source: &str,
    warnings: &mut Vec<CorrectionWarning>,
) -> Result<(), FerroError> {
    let Some(found) = crate::hgvs::genomic_offsets::genomic_axis_offset(variant) else {
        return Ok(());
    };

    match action {
        ResolvedAction::Accept | ResolvedAction::SilentCorrect => Ok(()),
        ResolvedAction::WarnCorrect => {
            warnings.push(CorrectionWarning::new(
                ErrorType::GenomicPositionOffset,
                format!("{}; normalization will refuse it", found.refusal()),
                None,
                source.to_string(),
                // No correction is derivable, so the "corrected" text is the
                // input itself; the refusal comes at normalize.
                source.to_string(),
            ));
            Ok(())
        }
        ResolvedAction::Reject => Err(FerroError::Parse {
            pos: 0,
            msg: format!(
                "[{}] {}",
                ErrorType::GenomicPositionOffset.code(),
                found.refusal()
            ),
            diagnostic: None,
        }),
    }
}

/// Enforce `background/standards.md:39` on a parsed variant (#1627).
///
/// The table's two daggered symbols — `X` ("masked nucleotide", `:36`) and `-`
/// ("gap of indeterminate length", `:37`) — are footnoted "used in alignment
/// only", and `general.md:48` admits only IUPAC-IUBMB nucleotide symbols, of
/// which neither is one. The decided
/// `rulings[alignment-only-symbol-in-a-description]` rules that ferro must
/// refuse both.
///
/// **The stage is mode-dependent**, per the decided
/// `rulings[absolute-prohibition-enforcement-stage]`: strict validates input
/// conformance and so fails here, at parse; lenient and silent do not validate
/// input conformance and accept, and fail later at *normalize*, because a
/// masked base cannot be resolved to a sequence (see
/// `Normalizer::normalize_core_checked`). Lenient carries a `W3028` warning
/// saying so; silent is the same acceptance with no message.
///
/// There is no repair arm: `WarnCorrect` and `SilentCorrect` accept rather than
/// correct, since `X` names a base the alignment did not resolve and ferro will
/// not invent one. That is why this returns `Result<(), _>` rather than a
/// rewritten variant, unlike [`apply_bracket_cardinality_rule`].
fn apply_alignment_only_symbol_rule(
    variant: &HgvsVariant,
    action: ResolvedAction,
    source: &str,
    warnings: &mut Vec<CorrectionWarning>,
) -> Result<(), FerroError> {
    let Some(found) = alignment_only_symbol(variant) else {
        return Ok(());
    };

    match action {
        ResolvedAction::Accept | ResolvedAction::SilentCorrect => Ok(()),
        ResolvedAction::WarnCorrect => {
            warnings.push(CorrectionWarning::new(
                ErrorType::AlignmentOnlySymbolInDescription,
                format!(
                    "`{}` states `{}` ({}), which {} lists as used in alignment only; \
                     a description may not state it, and normalization will refuse it",
                    found.stated,
                    found.symbol,
                    found.meaning(),
                    found.clause(),
                ),
                None,
                source.to_string(),
                // No correction is derivable, so the "corrected" text is the
                // input itself; the refusal comes at normalize.
                source.to_string(),
            ));
            Ok(())
        }
        ResolvedAction::Reject => Err(FerroError::Parse {
            pos: 0,
            msg: format!(
                "[{}] `{}` states `{}` ({}), which {} lists as used in alignment only: \
                 it is not one of the IUPAC-IUBMB nucleotide symbols `{}` admits, \
                 so the description denotes no sequence. State the resolved bases instead \
                 (`{}` is the IUPAC symbol for an unknown base).",
                ErrorType::AlignmentOnlySymbolInDescription.code(),
                found.stated,
                found.symbol,
                found.meaning(),
                found.clause(),
                found.alphabet_clause(),
                found.unknown_base(),
            ),
            diagnostic: None,
        }),
    }
}

/// Restore the reference run stated by a deprecated `<ref>><alt>` description
/// onto the `Delins` the preprocessor's textual rewrite produced (#1092).
///
/// The rewrite happens on the input *string*, so the run is only recoverable
/// from the original text — this is the one place that is unavoidable, and it
/// is a repair path, not the conformance rule (which is AST-keyed, see
/// [`variant::validate_no_multibase_substitution`]).
///
/// Deliberately narrow: it fires only for a single non-allele nucleic-acid
/// variant whose edit is a bare `Delins` (no stated `deleted` / `deleted_length`
/// and no run already recorded), and only when the original names exactly one
/// such description. Anything more ambiguous is left alone rather than
/// attributing a run to the wrong edit.
fn restore_stated_substitution_reference(variant: &mut HgvsVariant, original: &str) {
    use crate::error_handling::corrections::stated_substitution_reference;
    use crate::hgvs::edit::{NaEdit, Sequence};
    use crate::hgvs::uncertainty::Mu;
    use std::str::FromStr;

    fn bare_delins_run(edit: &mut Mu<NaEdit>) -> Option<&mut Option<Sequence>> {
        match edit {
            Mu::Certain(NaEdit::Delins {
                deleted: None,
                deleted_length: None,
                substitution_reference: slot @ None,
                ..
            })
            | Mu::Uncertain(NaEdit::Delins {
                deleted: None,
                deleted_length: None,
                substitution_reference: slot @ None,
                ..
            }) => Some(slot),
            _ => None,
        }
    }

    let Some(run) = stated_substitution_reference(original) else {
        return;
    };
    let Ok(run) = Sequence::from_str(&run) else {
        return;
    };
    let edit = match variant {
        HgvsVariant::Genome(v) => &mut v.loc_edit.edit,
        HgvsVariant::Cds(v) => &mut v.loc_edit.edit,
        HgvsVariant::Tx(v) => &mut v.loc_edit.edit,
        HgvsVariant::Rna(v) => &mut v.loc_edit.edit,
        HgvsVariant::Mt(v) => &mut v.loc_edit.edit,
        HgvsVariant::Circular(v) => &mut v.loc_edit.edit,
        _ => return,
    };
    if let Some(slot) = bare_delins_run(edit) {
        *slot = Some(run);
    }
}

/// Enforce the HGVS bracket/allele cardinality rule on a parsed variant (#493).
///
/// `[ ]` is allele syntax. The spec admits only two conformant shapes, identically
/// across every coordinate system (`c/g/n/m/o/r/p`): one bracket group with two or
/// more cis members (`c.[76A>C;88G>T]`), or two or more trans groups
/// (`c.[76A>C];[88G>T]`). A standalone single-member bracket — `c.[76A>C]`,
/// `g.[1000G>A]`, `p.[=]`, `c.[?]`, `p.[(?)]` — is non-conformant
/// (`DNA/alleles.md`). It parses to an [`HgvsVariant::Allele`] with exactly one
/// member, which is the uniform signal this rule keys off.
///
/// The canonical repair drops the redundant brackets, unwrapping the singleton to
/// its sole member ([`AlleleVariant::Display`](crate::hgvs::variant::AlleleVariant)
/// already renders a singleton bracket-free, so the corrected output is identical).
/// `Reject` (strict) returns a `W3026` parse error; `WarnCorrect` (lenient) unwraps
/// and pushes a `W3026` warning; `SilentCorrect` (silent) unwraps quietly; `Accept`
/// leaves the wrapper intact. Conformant multi-member brackets pass through
/// untouched.
fn apply_bracket_cardinality_rule(
    variant: HgvsVariant,
    action: ResolvedAction,
    source: &str,
    warnings: &mut Vec<CorrectionWarning>,
) -> Result<HgvsVariant, FerroError> {
    // Only a standalone single-member allele bracket is non-conformant. Multi-member
    // cis groups and >=2 trans groups carry more than one member and pass through.
    let is_singleton = matches!(&variant, HgvsVariant::Allele(av) if av.variants.len() == 1);
    if !is_singleton {
        return Ok(variant);
    }

    match action {
        ResolvedAction::Accept => Ok(variant),
        ResolvedAction::Reject => {
            let canonical = variant.to_string();
            Err(FerroError::Parse {
                pos: 0,
                msg: format!(
                    "[{}] standalone single-member allele bracket is not HGVS-conformant: \
                     allele brackets require two or more cis members (`c.[a;b]`) or two or \
                     more trans groups (`c.[a];[b]`). Drop the brackets (canonical \
                     `{}`) or pair with a second allele (e.g. `[a];[=]`).",
                    ErrorType::NonConformantBracketCardinality.code(),
                    canonical,
                ),
                diagnostic: None,
            })
        }
        ResolvedAction::WarnCorrect | ResolvedAction::SilentCorrect => {
            // Canonical repair: unwrap the singleton wrapper to its sole member.
            let member = match variant {
                HgvsVariant::Allele(av) => av
                    .variants
                    .into_iter()
                    .next()
                    .expect("singleton allele has exactly one member"),
                // Unreachable: `is_singleton` already established this is an Allele.
                other => return Ok(other),
            };
            if action == ResolvedAction::WarnCorrect {
                let canonical = member.to_string();
                warnings.push(CorrectionWarning::new(
                    ErrorType::NonConformantBracketCardinality,
                    format!(
                        "standalone single-member allele bracket unwrapped to `{canonical}`; \
                         allele brackets require >=2 cis members or >=2 trans groups"
                    ),
                    None,
                    source.to_string(),
                    canonical,
                ));
            }
            Ok(member)
        }
    }
}

/// Parse an HGVS string with lenient error handling.
///
/// This is a convenience function that uses lenient mode, which auto-corrects
/// common errors and returns warnings.
///
/// # Examples
///
/// Auto-corrects an en-dash to a hyphen:
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// let result = parse_hgvs_lenient("NM_000088.3:c.100\u{2013}200del");
/// assert!(result.is_ok());
/// ```
///
/// Soft-validation warnings emitted in lenient mode include W1001 for
/// lowercase amino-acid codes:
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// let parsed = parse_hgvs_lenient("NP_000079.2:p.val600glu").unwrap();
/// assert_eq!(parsed.preprocessed_input, "NP_000079.2:p.Val600Glu");
/// assert!(parsed
///     .warnings
///     .iter()
///     .any(|w| w.error_type.code() == "W1001"));
/// ```
///
/// W1002 for one-letter amino-acid abbreviations:
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// let parsed = parse_hgvs_lenient("NP_000079.2:p.V600E").unwrap();
/// assert_eq!(parsed.preprocessed_input, "NP_000079.2:p.Val600Glu");
/// assert!(parsed
///     .warnings
///     .iter()
///     .any(|w| w.error_type.code() == "W1002"));
/// ```
///
/// W3001 for accessions missing a `.<version>` suffix (the warning fires
/// but the input is not auto-corrected — the version cannot be synthesised):
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// let parsed = parse_hgvs_lenient("NM_000088:c.100A>G").unwrap();
/// assert_eq!(parsed.preprocessed_input, "NM_000088:c.100A>G");
/// assert!(parsed
///     .warnings
///     .iter()
///     .any(|w| w.error_type.code() == "W3001"));
/// ```
///
/// ```
/// use ferro_hgvs::error_handling::ErrorType;
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// // SVA-008: single-position range is collapsed (W4003)
/// let result = parse_hgvs_lenient("NM_000088.3:c.123_123del").unwrap();
/// assert_eq!(result.preprocessed_input, "NM_000088.3:c.123del");
/// assert!(result
///     .warnings
///     .iter()
///     .any(|w| w.error_type == ErrorType::SinglePositionRange));
/// ```
///
/// ```
/// use ferro_hgvs::error_handling::ErrorType;
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// // SVA-010: empty delins is rewritten to del (W3012)
/// let result = parse_hgvs_lenient("NC_000001.11:g.100_102delins").unwrap();
/// assert_eq!(result.preprocessed_input, "NC_000001.11:g.100_102del");
/// assert!(result
///     .warnings
///     .iter()
///     .any(|w| w.error_type == ErrorType::EmptyDelinsInsert));
/// ```
///
/// ```
/// use ferro_hgvs::error_handling::ErrorType;
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// // SVA-007: deletion with size-count suffix warns but is not rewritten (W3011)
/// let result = parse_hgvs_lenient("NG_012232.1:g.123del6").unwrap();
/// assert!(result
///     .warnings
///     .iter()
///     .any(|w| w.error_type == ErrorType::DelSizeSuffix));
/// ```
pub fn parse_hgvs_lenient(input: &str) -> Result<ParseResultWithWarnings<HgvsVariant>, FerroError> {
    parse_hgvs_with_config(input, ErrorConfig::lenient())
}

/// Parse an HGVS string with silent error handling.
///
/// This is a convenience function that uses silent mode, which auto-corrects
/// common errors without generating warnings.
///
/// # Example
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_silent;
///
/// // This will silently auto-correct the en-dash to hyphen
/// let result = parse_hgvs_silent("NM_000088.3:c.100\u{2013}200del");
/// assert!(result.is_ok());
/// assert!(!result.unwrap().has_warnings()); // No warnings in silent mode
/// ```
pub fn parse_hgvs_silent(input: &str) -> Result<ParseResultWithWarnings<HgvsVariant>, FerroError> {
    parse_hgvs_with_config(input, ErrorConfig::silent())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error_handling::{ErrorOverride, ErrorType};

    #[test]
    fn test_parse_simple_substitution() {
        let result = parse_hgvs("NC_000001.11:g.12345A>G");
        assert!(result.is_ok());
    }

    #[test]
    fn test_trim_hgvs_matches_str_trim() {
        // `trim_hgvs` must be byte-for-byte identical to `str::trim` on every
        // input, including ASCII and multi-byte Unicode whitespace at either end.
        let cases = [
            "NM_000088.3:c.459A>G",                 // no whitespace (the fast path)
            "  NM_000088.3:c.459A>G  ",             // ASCII spaces
            "\tNM_000088.3:c.459A>G\n",             // ASCII tab / newline
            "\u{00A0}NM_000088.3:c.459A>G",         // leading non-breaking space (U+00A0)
            "NM_000088.3:c.459A>G\u{3000}",         // trailing ideographic space (U+3000)
            "\u{2028}NM_000088.3:c.459A>G\u{2029}", // line/paragraph separators
            "  ",                                   // all whitespace
            "",                                     // empty
            "x",                                    // single non-ws byte
        ];
        for input in cases {
            assert_eq!(
                trim_hgvs(input),
                input.trim(),
                "trim mismatch for {input:?}"
            );
        }
    }

    #[test]
    fn test_parse_deletion() {
        let result = parse_hgvs("NM_000088.3:c.459del");
        assert!(result.is_ok());
    }

    // Error handling mode tests
    #[test]
    fn test_parse_with_config_strict_valid() {
        let config = ErrorConfig::strict();
        let result = parse_hgvs_with_config("NM_000088.3:c.459del", config);
        assert!(result.is_ok());
        let parsed = result.unwrap();
        assert!(!parsed.had_corrections());
        assert!(!parsed.has_warnings());
    }

    #[test]
    fn test_parse_with_config_strict_rejects_en_dash() {
        let config = ErrorConfig::strict();
        let result = parse_hgvs_with_config("NM_000088.3:c.100\u{2013}200del", config);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_with_config_strict_rejects_whitespace() {
        let config = ErrorConfig::strict();
        let result = parse_hgvs_with_config("  NM_000088.3:c.459del  ", config);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_with_config_lenient_corrects_whitespace() {
        let config = ErrorConfig::lenient();
        let result = parse_hgvs_with_config("  NM_000088.3:c.459del  ", config);
        assert!(result.is_ok());
        let parsed = result.unwrap();
        assert!(parsed.had_corrections());
        assert!(parsed.has_warnings());
    }

    #[test]
    fn test_parse_with_config_silent_no_warnings() {
        let config = ErrorConfig::silent();
        let result = parse_hgvs_with_config("  NM_000088.3:c.459del  ", config);
        assert!(result.is_ok());
        let parsed = result.unwrap();
        assert!(parsed.had_corrections());
        assert!(!parsed.has_warnings()); // Silent mode = no warnings
    }

    #[test]
    fn test_parse_with_config_override() {
        // Lenient mode but override whitespace to reject
        let config =
            ErrorConfig::lenient().with_override(ErrorType::ExtraWhitespace, ErrorOverride::Reject);
        let result = parse_hgvs_with_config("  NM_000088.3:c.459del  ", config);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_lenient() {
        let result = parse_hgvs_lenient("  NM_000088.3:c.459del  ");
        assert!(result.is_ok());
        assert!(result.unwrap().had_corrections());
    }

    #[test]
    fn test_parse_silent() {
        let result = parse_hgvs_silent("  NM_000088.3:c.459del  ");
        assert!(result.is_ok());
        let parsed = result.unwrap();
        assert!(parsed.had_corrections());
        assert!(!parsed.has_warnings());
    }

    #[test]
    fn test_parse_lowercase_accession_lenient() {
        let result = parse_hgvs_lenient("nm_000088.3:c.459del");
        assert!(result.is_ok());
        let parsed = result.unwrap();
        assert!(parsed.had_corrections());
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.459del");
    }
}
