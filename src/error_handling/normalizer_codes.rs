//! The normalizer's SCREAMING_SNAKE diagnostic namespace, and how it resolves
//! against the `E`/`W` code registry (#2092).
//!
//! ferro emits diagnostics under two code namespaces, and `ferro normalize`
//! prints both to stderr in the same `warning[<CODE>]: <message>` shape:
//!
//! * the **preprocessor**'s `W`/`E`-numbered codes, one per
//!   [`ErrorType`](crate::error_handling::ErrorType), documented in
//!   [`registry`](super::registry) and configurable through `--error-mode` /
//!   `--ignore` / `--reject`; and
//! * the **normalizer**'s SCREAMING_SNAKE codes, one per
//!   [`NormalizationWarning`](crate::normalize::NormalizationWarning) /
//!   [`NormalizationInfo`](crate::normalize::NormalizationInfo) variant.
//!
//! Before #2092 `ferro explain` covered only the first, so a user who read
//! `warning[MEMBERS_COALESCED_FROM_REPORTED_FORM]:` out of `ferro normalize`
//! and asked about it got `Unknown code` — indistinguishable from typing a
//! code that does not exist.
//!
//! # Why most of these are aliases rather than new prose
//!
//! Nine of the fifteen normalizer codes name a condition the registry already
//! documents under a `W`-code: the two namespaces are two names for one thing,
//! not two things. `W5003`'s own explanation already says so in as many words
//! ("also published as `CanonicalSplitSkipped` in the NormalizationWarning
//! enum"). [`NORMALIZER_CODES`] makes that correspondence explicit and
//! machine-readable, so `explain` resolves through it instead of duplicating —
//! and drifting from — text that already exists.
//!
//! The remaining entries are normalizer-only: they have no `ErrorType`, so
//! there is no `W`-code to alias and nothing for `--ignore` / `--reject` to
//! act on. Their `CodeInfo` is written here, derived from the doc comment and
//! `Display` impl of the variant that emits it (each entry names its source).

use super::codes::{CodeCategory, CodeInfo, ModeBehavior};
use std::collections::HashMap;
use std::sync::OnceLock;

/// One code in the normalizer's SCREAMING_SNAKE namespace, and the registry
/// entry that documents it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NormalizerCode {
    /// The code as emitted — what
    /// [`NormalizationWarning::code`](crate::normalize::NormalizationWarning::code)
    /// or
    /// [`NormalizationInfo::code`](crate::normalize::NormalizationInfo::code)
    /// returns, and what the user reads out of `ferro normalize` stderr.
    pub code: &'static str,
    /// The enum variant that emits it. Carried so a single table can serve
    /// both `explain` (which keys on the code) and the error-code audit
    /// (which keys on the variant).
    pub variant: &'static str,
    /// The registry key documenting this condition: a `W`-code when the
    /// preprocessor namespace already names the same thing, otherwise
    /// [`Self::code`] itself, which resolves to the entry defined below.
    pub registry_code: &'static str,
}

impl NormalizerCode {
    /// Whether this code is documented by its own entry in this module rather
    /// than by an existing `W`-code entry in [`registry`](super::registry).
    pub fn is_normalizer_only(&self) -> bool {
        self.registry_code == self.code
    }
}

/// Every code the normalizer can emit, with the registry entry it resolves to.
///
/// Kept in the order the variants are declared in `src/normalize/mod.rs`
/// (warnings first, then infos) so the two read side by side.
///
/// **This table must cover every variant of both diagnostic enums.**
/// `tests/it/issue_2092_explain_normalizer_codes.rs` scans those enums out of
/// the source and fails if any code here is missing or unresolvable, so a new
/// variant cannot quietly ship without an explanation.
pub const NORMALIZER_CODES: &[NormalizerCode] = &[
    NormalizerCode {
        code: "REFSEQ_MISMATCH",
        variant: "RefSeqMismatch",
        registry_code: "W5001",
    },
    NormalizerCode {
        code: "OVERLAP_CONFLICTING_EDITS",
        variant: "OverlapConflict",
        registry_code: "W5002",
    },
    NormalizerCode {
        code: "MEMBERS_COALESCED_FROM_REPORTED_FORM",
        variant: "MembersCoalesced",
        registry_code: "W5005",
    },
    // The warning is generic over telomere/centromere markers; the registry
    // entry is specific to the only unresolvable one (`cen`).
    NormalizerCode {
        code: "UNRESOLVABLE_SPECIAL_POSITION",
        variant: "UnresolvableSpecialPosition",
        registry_code: "W4005",
    },
    NormalizerCode {
        code: "TRANSCRIPT_FLANK_NOT_DESCRIBABLE",
        variant: "TranscriptFlankNotDescribable",
        registry_code: "W4006",
    },
    // The warning names the procedural effect (the canonical split was
    // skipped); the registry entry names the underlying spec violation (the
    // variant is not encompassed by the reference). One condition, two names —
    // W5003's own explanation already records the pairing.
    NormalizerCode {
        code: "CANONICAL_SPLIT_SKIPPED",
        variant: "CanonicalSplitSkipped",
        registry_code: "W5003",
    },
    NormalizerCode {
        code: "CROSS_AXIS_VARIANT_NOT_SHUFFLED",
        variant: "CrossAxisVariantNotShuffled",
        registry_code: "CROSS_AXIS_VARIANT_NOT_SHUFFLED",
    },
    NormalizerCode {
        code: "AXIS_CLAMP_APPLIED",
        variant: "AxisClampApplied",
        registry_code: "AXIS_CLAMP_APPLIED",
    },
    NormalizerCode {
        code: "INITIATOR_MET_CANONICALIZATION",
        variant: "InitiatorMetCanonicalization",
        registry_code: "W3022",
    },
    NormalizerCode {
        code: "INSERTED_SEQUENCE_EXPANDED",
        variant: "InsertedSequenceExpanded",
        registry_code: "INSERTED_SEQUENCE_EXPANDED",
    },
    NormalizerCode {
        code: "POSITION_PAST_END",
        variant: "PositionPastEnd",
        registry_code: "W4004",
    },
    NormalizerCode {
        code: "INTRONIC_ON_BARE_TRANSCRIPT",
        variant: "IntronicOnBareTranscript",
        registry_code: "W4007",
    },
    NormalizerCode {
        code: "INCOMPLETE_CDS_START_REFERENCE",
        variant: "IncompleteCdsStartReference",
        registry_code: "W5004",
    },
    NormalizerCode {
        code: "REDUCED_CAPABILITY_NO_GENOME",
        variant: "ReducedCapabilityNoGenome",
        registry_code: "REDUCED_CAPABILITY_NO_GENOME",
    },
    // `NormalizationInfo`, not `NormalizationWarning`: info-grade, but it is
    // still a code a caller reads off a `NormalizeResult`.
    NormalizerCode {
        code: "SHUFFLE_APPLIED",
        variant: "ShuffleApplied",
        registry_code: "SHUFFLE_APPLIED",
    },
];

/// The registry key documenting `code`, or `None` when `code` is not a
/// normalizer code. `code` must already be uppercase.
pub fn registry_key_for(code: &str) -> Option<&'static str> {
    NORMALIZER_CODES
        .iter()
        .find(|entry| entry.code == code)
        .map(|entry| entry.registry_code)
}

/// The [`CodeInfo`] for a normalizer-only code — one with no `W`-code
/// counterpart. Returns `None` for aliased codes (look those up in
/// [`registry`](super::registry) under [`registry_key_for`]) and for anything
/// that is not a normalizer code at all.
pub fn normalizer_only_entry(code: &str) -> Option<&'static CodeInfo> {
    normalizer_only_registry().get(code)
}

/// Whether `code` is shaped like a normalizer code — SCREAMING_SNAKE with at
/// least one underscore — regardless of whether it is a *known* one.
///
/// Used to tell a user who typed an unrecognized SCREAMING_SNAKE code which
/// namespace they appear to be reaching for, rather than answering every miss
/// with the same undifferentiated "Unknown code".
pub fn has_normalizer_code_shape(code: &str) -> bool {
    code.contains('_')
        && !code.is_empty()
        && code
            .chars()
            .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || c == '_')
        && code.starts_with(|c: char| c.is_ascii_uppercase())
}

static NORMALIZER_ONLY_REGISTRY: OnceLock<HashMap<&'static str, CodeInfo>> = OnceLock::new();

fn normalizer_only_registry() -> &'static HashMap<&'static str, CodeInfo> {
    NORMALIZER_ONLY_REGISTRY.get_or_init(build_normalizer_only_registry)
}

/// Entries for the normalizer codes with no `ErrorType` and therefore no
/// `W`-code.
///
/// Every summary and explanation below is derived from the doc comment and
/// `Display` impl of the emitting variant in `src/normalize/mod.rs`; each entry
/// names its source. Nothing here states a normalization rule that is not
/// already stated there — CLI help text that invents semantics is worse than
/// no help text, because it drifts from the behaviour it claims to describe.
///
/// None of these carry a `mode_behavior` derived from `ErrorConfig`: without an
/// `ErrorType` there is nothing for `--error-mode` / `--ignore` / `--reject` to
/// resolve. The one exception is `REDUCED_CAPABILITY_NO_GENOME`, whose strict-
/// mode promotion is keyed on the mode directly
/// (`NormalizeConfig::should_reject_reduced_capability` asks
/// `ErrorMode::is_strict`, deliberately not `action_for`).
fn build_normalizer_only_registry() -> HashMap<&'static str, CodeInfo> {
    let mut map = HashMap::new();

    // Source: `NormalizationWarning::CrossAxisVariantNotShuffled` doc comment
    // and `Display` arm (`src/normalize/mod.rs`).
    map.insert(
        "CROSS_AXIS_VARIANT_NOT_SHUFFLED",
        CodeInfo {
            code: "CROSS_AXIS_VARIANT_NOT_SHUFFLED",
            name: "CrossAxisVariantNotShuffled",
            summary: "A c. variant spans two coordinate sub-axes, so the 3\u{2032}-rule shuffle \
                was skipped.",
            explanation: "The variant's start and end positions sit in different `c.` coordinate \
                sub-axes (5\u{2032}UTR / CDS / 3\u{2032}UTR). The 3\u{2032}-rule shuffle has no \
                well-defined semantics across an axis boundary \u{2014} the two ends would move \
                through differently-numbered coordinate spaces \u{2014} so ferro preserves the \
                canonical input position and reports the skip rather than guessing. The output is \
                otherwise canonical; only the shuffle step was declined. This code is emitted by \
                the normalizer and has no `ErrorType`, so `--error-mode`, `--ignore` and \
                `--reject` do not apply to it.",
            category: CodeCategory::NormalizerDiagnostic,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/background/numbering/"),
            related_codes: &["AXIS_CLAMP_APPLIED", "SHUFFLE_APPLIED"],
        },
    );

    // Source: `NormalizationWarning::AxisClampApplied` doc comment and
    // `Display` arm (`src/normalize/mod.rs`).
    map.insert(
        "AXIS_CLAMP_APPLIED",
        CodeInfo {
            code: "AXIS_CLAMP_APPLIED",
            name: "AxisClampApplied",
            summary: "A 3\u{2032}-rule shuffle was stopped at a CDS\u{2194}UTR sub-axis boundary.",
            explanation: "The 3\u{2032}-rule shuffle would have carried the variant across a \
                CDS\u{2194}UTR coordinate sub-axis boundary; the axis clamp constrained the \
                result to the boundary instead. The reported position is therefore the \
                boundary-most position on the input's own sub-axis, not the position an \
                unclamped shuffle would have produced. The warning names the direction that was \
                clamped (5\u{2032} or 3\u{2032}) and which bound did the clamping. This code is \
                emitted by the normalizer and has no `ErrorType`, so `--error-mode`, `--ignore` \
                and `--reject` do not apply to it.",
            category: CodeCategory::NormalizerDiagnostic,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/background/numbering/"),
            related_codes: &["CROSS_AXIS_VARIANT_NOT_SHUFFLED", "SHUFFLE_APPLIED"],
        },
    );

    // Source: `NormalizationWarning::InsertedSequenceExpanded` doc comment,
    // field docs, and `Display` arm (`src/normalize/mod.rs`).
    map.insert(
        "INSERTED_SEQUENCE_EXPANDED",
        CodeInfo {
            code: "INSERTED_SEQUENCE_EXPANDED",
            name: "InsertedSequenceExpanded",
            summary: "A bracketed or reference-range `ins[\u{2026}]` payload was expanded to a \
                flat literal sequence.",
            explanation: "The inserted-sequence payload was written in bracketed form \u{2014} a \
                literal (`ins[ATC]`), a reference range (`ins[100_120inv]`), or a compound of the \
                two (`ins[A;100_110]`) \u{2014} and canonicalization resolved it to the flat \
                literal sequence it denotes. The warning is emitted alongside the rewrite purely \
                for observability: it records the original payload and the literal it expanded \
                to, so a caller can audit which inputs were canonicalized and which were \
                preserved verbatim. It does not report a defect in the input. This code is \
                emitted by the normalizer and has no `ErrorType`, so `--error-mode`, `--ignore` \
                and `--reject` do not apply to it.",
            category: CodeCategory::NormalizerDiagnostic,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/insertion/",
            ),
            related_codes: &[],
        },
    );

    // Source: `NormalizationWarning::ReducedCapabilityNoGenome` doc comment and
    // `Display` arm, plus `NormalizeConfig::should_reject_reduced_capability`
    // for the mode behaviour (`src/normalize/config.rs`).
    map.insert(
        "REDUCED_CAPABILITY_NO_GENOME",
        CodeInfo {
            code: "REDUCED_CAPABILITY_NO_GENOME",
            name: "ReducedCapabilityNoGenome",
            summary: "A genome-requiring normalization step could not run, so the result is \
                best-effort rather than fully normalized.",
            explanation: "Normalization reached a step that needs genomic sequence, but the \
                configured reference provider carries none (a transcripts-only provider, for \
                example). The result is returned UNCHANGED from the point that step would have \
                refined it, and this warning marks the output as degraded so a reduced-capability \
                result is never mistaken for a fully-normalized one. This is an environmental \
                limitation, not a defect in the variant: the same input against a genome-backed \
                reference may normalize further \u{2014} the intronic and boundary-spanning paths \
                genuinely would, while an exon-junction 3\u{2032}-shuffle only might, since \
                whether the pattern extends into the intron is exactly what could not be checked. \
                Strict mode refuses rather than returning a knowingly-degraded result; lenient \
                and silent modes return the best-effort variant with this advisory. The strict \
                refusal is keyed on the mode itself, so it is deliberately not overridable \
                per-code.",
            category: CodeCategory::NormalizerDiagnostic,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: Some(ModeBehavior::always_warn_if_not_rejected()),
            hgvs_spec_url: None,
            related_codes: &["CANONICAL_SPLIT_SKIPPED"],
        },
    );

    // Source: `NormalizationInfo::ShuffleApplied` doc comment and `Display`
    // arm (`src/normalize/mod.rs`), plus `error_handling::info_map` for the
    // mutalyzer correspondence.
    map.insert(
        "SHUFFLE_APPLIED",
        CodeInfo {
            code: "SHUFFLE_APPLIED",
            name: "ShuffleApplied",
            summary: "The shuffle layer relocated the variant to its canonical position.",
            explanation: "An insertion, deletion or duplication in a repeated stretch can be \
                described at several positions; the HGVS arbitrary-position rule picks one. The \
                shuffle layer moved the variant to that position and recorded the move here, \
                carrying the input position and the normalized position so a caller can see \
                exactly what changed. Every shipped path shuffles 3\u{2032}, the rightmost form \
                the HGVS recommendations mandate; the direction is still reported explicitly \
                because ferro's own tests drive a 5\u{2032} arm as a differential oracle. This is \
                an info-grade signal rather than a warning \u{2014} nothing is wrong with the \
                input. The mutalyzer equivalent is `ICORRECTEDPOINT`, which is direction-agnostic \
                and therefore a lossy mapping. This code has no `ErrorType`, so `--error-mode`, \
                `--ignore` and `--reject` do not apply to it.",
            category: CodeCategory::NormalizerDiagnostic,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/general/"),
            related_codes: &["AXIS_CLAMP_APPLIED", "CROSS_AXIS_VARIANT_NOT_SHUFFLED"],
        },
    );

    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error_handling::ErrorType;

    #[test]
    fn every_table_entry_resolves_to_an_entry_or_an_alias() {
        for entry in NORMALIZER_CODES {
            if entry.registry_code == entry.code {
                assert!(
                    normalizer_only_entry(entry.code).is_some(),
                    "{} is marked normalizer-only but has no entry",
                    entry.code
                );
            } else {
                assert!(
                    entry.registry_code.starts_with('W'),
                    "{} aliases {}, which is not a W-code",
                    entry.code,
                    entry.registry_code
                );
            }
        }
    }

    #[test]
    fn normalizer_only_entries_are_all_reachable_from_the_table() {
        let tabled: Vec<&str> = NORMALIZER_CODES
            .iter()
            .filter(|e| e.registry_code == e.code)
            .map(|e| e.code)
            .collect();
        for code in normalizer_only_registry().keys() {
            assert!(
                tabled.contains(code),
                "{code} has an entry but is not listed in NORMALIZER_CODES, so `explain` cannot \
                 reach it"
            );
        }
    }

    #[test]
    fn code_and_key_agree_for_every_normalizer_only_entry() {
        for (key, info) in normalizer_only_registry() {
            assert_eq!(*key, info.code, "map key and CodeInfo::code disagree");
        }
    }

    /// The entries claim `--error-mode` / `--ignore` / `--reject` do not reach
    /// these codes. That claim is only true while they have no `ErrorType`, so
    /// it is asserted rather than trusted.
    #[test]
    fn normalizer_only_codes_have_no_error_type() {
        for entry in NORMALIZER_CODES {
            if entry.registry_code != entry.code {
                continue;
            }
            assert!(
                ErrorType::from_code(entry.code).is_none(),
                "{} now has an ErrorType; its explanation says the mode flags do not reach it",
                entry.code
            );
        }
    }

    #[test]
    fn table_codes_are_unique() {
        let mut seen: Vec<&str> = Vec::new();
        for entry in NORMALIZER_CODES {
            assert!(
                !seen.contains(&entry.code),
                "duplicate code {} in NORMALIZER_CODES",
                entry.code
            );
            seen.push(entry.code);
        }
    }

    #[test]
    fn shape_test_accepts_normalizer_codes_and_rejects_w_codes() {
        assert!(has_normalizer_code_shape("REFSEQ_MISMATCH"));
        assert!(has_normalizer_code_shape("NOT_A_REAL_CODE"));
        assert!(!has_normalizer_code_shape("W3003"));
        assert!(!has_normalizer_code_shape("E1001"));
        assert!(!has_normalizer_code_shape("W-LOAD-001"));
        assert!(!has_normalizer_code_shape(""));
        assert!(!has_normalizer_code_shape("refseq_mismatch"));
    }
}
