//! Reference sequence validation for normalization.
//!
//! Validates that stated reference bases in HGVS expressions match
//! the actual reference sequence.
//!
//! # Coordinate System
//!
//! | Parameter | Basis | Notes |
//! |-----------|-------|-------|
//! | `start`, `end` | 1-based | HGVS positions (inclusive) |
//! | Array indexing | 0-based | Converted via `hgvs_pos_to_index()` |
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::coords::hgvs_pos_to_index;
use crate::error::FerroError;
use crate::error_handling::ResolvedAction;
use crate::hgvs::edit::{Base, NaEdit, RepeatUnit, Sequence};
use crate::normalize::config::NormalizeConfig;

/// Result of reference validation
#[derive(Debug, Clone)]
pub struct ValidationResult {
    /// Whether validation passed
    pub valid: bool,
    /// Warning message if any (for lenient mode)
    pub warning: Option<String>,
    /// The actual reference sequence (for correction)
    pub actual_ref: Option<String>,
    /// The stated reference sequence
    pub stated_ref: Option<String>,
}

impl ValidationResult {
    /// Create a passing validation result
    pub fn ok() -> Self {
        Self {
            valid: true,
            warning: None,
            actual_ref: None,
            stated_ref: None,
        }
    }

    /// Create a mismatch result
    pub fn mismatch(stated: String, actual: String) -> Self {
        Self {
            valid: false,
            warning: Some(format!(
                "Reference mismatch: stated '{}' but actual is '{}'",
                stated, actual
            )),
            actual_ref: Some(actual),
            stated_ref: Some(stated),
        }
    }
}

/// Validate the reference bases in an edit against the actual sequence.
///
/// Returns a ValidationResult indicating whether the stated reference matches
/// the actual sequence, and provides correction information if not.
///
/// Stated bases (when present) take precedence. An explicit numeric length
/// (`del<N>`/`dup<N>`/`del<N>ins…`) is otherwise validated against the span
/// (#486); uncertain-extent dups (`dupN?`) are skipped.
///
/// # Arguments
/// * `edit` - The edit to validate
/// * `ref_seq` - The reference sequence (1-indexed, so position 1 is at index 0)
/// * `start` - Start position (1-indexed)
/// * `end` - End position (1-indexed, inclusive)
pub fn validate_reference(edit: &NaEdit, ref_seq: &[u8], start: u64, end: u64) -> ValidationResult {
    match edit {
        NaEdit::Substitution { reference, .. } => validate_single_base(reference, ref_seq, start),
        NaEdit::Deletion { sequence, length } => {
            if let Some(seq) = sequence {
                validate_sequence(seq.bases(), ref_seq, start, end)
            } else if let Some(n) = length {
                // Numeric `del<N>`: N must equal the position span (#486).
                validate_stated_length(*n, start, end)
            } else {
                // No sequence stated, nothing to validate
                ValidationResult::ok()
            }
        }
        NaEdit::Delins {
            deleted,
            deleted_length,
            ..
        } => {
            if let Some(seq) = deleted {
                // Non-recommended `delACinsGT` form: ferro parses the stated
                // deleted bases into `deleted`. Validate them against the
                // reference span (#486) — mirrors the `Deletion` arm.
                validate_sequence(seq.bases(), ref_seq, start, end)
            } else if let Some(n) = deleted_length {
                // Numeric `del<N>ins…`: N must equal the position span (#486).
                validate_stated_length(*n, start, end)
            } else {
                // Recommended short form `delinsXXX`: no stated deleted bases.
                ValidationResult::ok()
            }
        }
        NaEdit::Duplication {
            sequence,
            length,
            uncertain_extent,
        } => {
            if let Some(seq) = sequence {
                validate_sequence(seq.bases(), ref_seq, start, end)
            } else if let Some(n) = length {
                // Numeric `dup<N>`: N must equal the span (#486). Skipped when
                // `uncertain_extent` is set (`dupN?`/`dupN(a_b)`) — there the
                // span is intentionally not the duplicated length.
                if uncertain_extent.is_some() {
                    ValidationResult::ok()
                } else {
                    validate_stated_length(*n, start, end)
                }
            } else {
                // No sequence stated, nothing to validate
                ValidationResult::ok()
            }
        }
        NaEdit::Inversion { sequence, .. } => {
            if let Some(seq) = sequence {
                validate_sequence(seq.bases(), ref_seq, start, end)
            } else {
                // `inv<N>` numeric length is intentionally NOT span-checked here
                // (no ELENGTHMISMATCH corpus row; the preprocessor W3016 covers
                // it at the config-parse layer). See #486.
                ValidationResult::ok()
            }
        }
        NaEdit::Repeat {
            sequence: Some(unit),
            additional_counts,
            trailing,
            ..
        } => {
            // Mirror the skip-list used by `normalize_na_edit`'s repeat
            // arm: genotype notation (`A[6][1]`, `additional_counts`
            // non-empty) and VEP-style trailing sequences
            // (`c.212-18CTG[3]T`, `trailing.is_some()`) extend the
            // reference span semantics beyond `unit_len × k`. The
            // normalizer declines to rewrite these shapes — so running
            // the strict span × unit_len check on them would surface a
            // `ReferenceMismatch` without any compensating normalization
            // benefit. Defer to the existing pass-through behavior.
            if trailing.is_some() || !additional_counts.is_empty() {
                ValidationResult::ok()
            } else {
                validate_repeat_tract(unit, ref_seq, start, end)
            }
        }
        NaEdit::MultiRepeat { units } => validate_multirepeat_tract(units, ref_seq, start, end),
        // Other edit types don't have stated reference bases
        _ => ValidationResult::ok(),
    }
}

/// Validate that the reference span `[start, end]` is a clean tandem
/// repeat of `unit` (HGVS `repeated.md` invariant).
///
/// Per HGVS v21.0 the span `[start, end]` of a `unit[N]` description
/// covers the reference repeat tract. The two invariants are:
///   1. `span_len % unit_len == 0` (the span divides cleanly into whole
///      units), and
///   2. The reference bases at `[start, end]` equal `unit` repeated
///      `span_len / unit_len` times.
///
/// The variant count `N` is the alt-allele count and is independent of
/// the reference repeat count — `normalize_repeat` re-derives the
/// reference count by scanning. The consistency check is between the
/// *span* and the *reference bases* only.
fn validate_repeat_tract(
    unit: &Sequence,
    ref_seq: &[u8],
    start: u64,
    end: u64,
) -> ValidationResult {
    let unit_bases = unit.bases();
    if unit_bases.is_empty() {
        return ValidationResult::ok();
    }
    let start_idx = hgvs_pos_to_index(start);
    let end_idx = end as usize; // 1-based inclusive end = 0-based exclusive end

    if end_idx > ref_seq.len() || start_idx > end_idx {
        // Out-of-range or inverted span: leave to other validation gates
        // (range/position checks elsewhere catch this with a clearer
        // message than a "ref mismatch" would).
        return ValidationResult::ok();
    }
    let span = end_idx - start_idx;
    let unit_len = unit_bases.len();
    let unit_str: String = unit_bases.iter().map(|b| b.to_char()).collect();
    let actual_bytes = &ref_seq[start_idx..end_idx];
    let actual_str: String = actual_bytes
        .iter()
        .map(|&b| (b as char).to_ascii_uppercase())
        .collect();

    if !span.is_multiple_of(unit_len) {
        // Divisibility gate: span length must be a whole multiple of
        // unit_len. Report the unit and span explicitly so the
        // ReferenceMismatch error is actionable.
        return ValidationResult::mismatch(
            format!("{}[k] (k whole copies)", unit_str),
            format!(
                "{} ({} bp; unit_len {} does not divide span)",
                actual_str, span, unit_len
            ),
        );
    }
    let k = span / unit_len;
    // Build the expected reference tract: unit repeated k times.
    let mut expected = Vec::with_capacity(span);
    let unit_bytes: Vec<u8> = unit_bases.iter().map(|b| b.to_u8()).collect();
    for _ in 0..k {
        expected.extend_from_slice(&unit_bytes);
    }
    // Case-insensitive compare (matches `validate_sequence`).
    let matches = expected
        .iter()
        .zip(actual_bytes.iter())
        .all(|(a, b)| a.eq_ignore_ascii_case(b));
    if matches {
        return ValidationResult::ok();
    }
    let expected_str: String = expected
        .iter()
        .map(|&b| (b as char).to_ascii_uppercase())
        .collect();
    ValidationResult::mismatch(
        format!("{}[{}] ({})", unit_str, k, expected_str),
        actual_str,
    )
}

/// Validate that the reference span is the concatenation of the
/// declared mixed-repeat units (e.g. `CTG[2]TTG[1]CTG[11]`).
///
/// Multi-repeat notation describes the reference structure exactly: the
/// per-unit counts must reproduce the reference bases at `[start, end]`
/// when expanded and concatenated.
///
/// Validation mode is chosen by inspecting the per-unit counts:
///
/// - **All units `Exact`** — full validation. The expanded
///   concatenation has a deterministic length; we check both the
///   span-vs-sum length invariant and the base-by-base content match.
/// - **Mixed Exact / non-Exact** — *anchored partial validation*
///   (issue #279 + #395 item 1). The leading-`Exact` prefix and the
///   trailing-`Exact` suffix are deterministic (left-anchored and
///   right-anchored respectively); we validate both against the
///   reference and skip the ambiguous middle. If `prefix_len +
///   suffix_len > span_len`, that's a content/length mismatch — the
///   suffix would overlap the prefix or run off the start of the
///   tract, which is incompatible with any combination of middle-unit
///   counts.
/// - **First AND last unit non-`Exact`** — no validation. Neither end
///   is anchored.
///
/// Mismatches surface as `RefSeqMismatch`.
///
/// Missing reference data (out-of-range span, inverted range) is
/// handled defensively: those cases short-circuit to `ok()` so other
/// validation gates can surface a clearer error.
fn validate_multirepeat_tract(
    units: &[RepeatUnit],
    ref_seq: &[u8],
    start: u64,
    end: u64,
) -> ValidationResult {
    use crate::hgvs::edit::RepeatCount;

    if units.is_empty() {
        return ValidationResult::ok();
    }

    // Walk units left-to-right, building the deterministic prefix from
    // the leading `Exact` runs. Stop at the first non-`Exact` unit;
    // everything from there to the right anchor is the ambiguous middle.
    let mut prefix_bytes: Vec<u8> = Vec::new();
    let mut prefix_str = String::new();
    let mut prefix_units = 0usize;
    for u in units {
        let n = match u.count {
            RepeatCount::Exact(n) => n,
            _ => break,
        };
        let unit_bytes: Vec<u8> = u.sequence.bases().iter().map(|b| b.to_u8()).collect();
        let unit_str: String = u.sequence.bases().iter().map(|b| b.to_char()).collect();
        prefix_str.push_str(&format!("{}[{}]", unit_str, n));
        for _ in 0..n {
            prefix_bytes.extend_from_slice(&unit_bytes);
        }
        prefix_units += 1;
    }

    let all_exact = prefix_units == units.len();

    // Walk units right-to-left for the trailing-`Exact` suffix
    // (closes #395 item 1). Only meaningful when the left walk did not
    // consume every unit; otherwise the whole tract is the prefix and
    // there is no separate suffix to validate.
    let mut suffix_bytes: Vec<u8> = Vec::new();
    let mut suffix_str = String::new();
    let mut suffix_units = 0usize;
    if !all_exact {
        // Bound the right walk so prefix and suffix never overlap on
        // the same `Exact` units (would happen if the units in between
        // were all `Exact` too, but we've already established at least
        // one non-`Exact` unit exists in `[prefix_units, units.len())`).
        for u in units[prefix_units..].iter().rev() {
            let n = match u.count {
                RepeatCount::Exact(n) => n,
                _ => break,
            };
            let unit_bytes: Vec<u8> = u.sequence.bases().iter().map(|b| b.to_u8()).collect();
            let unit_str: String = u.sequence.bases().iter().map(|b| b.to_char()).collect();
            // Build suffix_str in left-to-right order so the diagnostic
            // reads naturally despite the right-to-left walk.
            suffix_str.insert_str(0, &format!("{}[{}]", unit_str, n));
            // Prepend bytes for the same reason — earlier units in the
            // tract appear earlier in the byte vector.
            let mut new_bytes =
                Vec::with_capacity(suffix_bytes.len() + unit_bytes.len() * n as usize);
            for _ in 0..n {
                new_bytes.extend_from_slice(&unit_bytes);
            }
            new_bytes.extend_from_slice(&suffix_bytes);
            suffix_bytes = new_bytes;
            suffix_units += 1;
        }
    }

    // No anchored units on either side → no validation. Preserves the
    // pre-#279 behavior for these descriptions.
    if prefix_units == 0 && suffix_units == 0 {
        return ValidationResult::ok();
    }

    let start_idx = hgvs_pos_to_index(start);
    let end_idx = end as usize;
    if end_idx > ref_seq.len() || start_idx > end_idx {
        return ValidationResult::ok();
    }
    let actual_bytes = &ref_seq[start_idx..end_idx];

    if all_exact {
        // Full validation: span length must equal expanded prefix
        // length, and bases must match.
        let actual_str: String = actual_bytes
            .iter()
            .map(|&b| (b as char).to_ascii_uppercase())
            .collect();
        if prefix_bytes.len() != actual_bytes.len() {
            return ValidationResult::mismatch(
                format!(
                    "{} ({} bp from declared multi-repeat units)",
                    prefix_str,
                    prefix_bytes.len()
                ),
                format!("{} ({} bp)", actual_str, actual_bytes.len()),
            );
        }
        let matches = prefix_bytes
            .iter()
            .zip(actual_bytes.iter())
            .all(|(a, b)| a.eq_ignore_ascii_case(b));
        if matches {
            return ValidationResult::ok();
        }
        let expected_bases_str: String = prefix_bytes
            .iter()
            .map(|&b| (b as char).to_ascii_uppercase())
            .collect();
        return ValidationResult::mismatch(
            format!("{} ({})", prefix_str, expected_bases_str),
            actual_str,
        );
    }

    // Anchored partial validation: check prefix and/or suffix; skip
    // the ambiguous middle.

    // Compute the minimum bp the middle units MUST occupy, given each
    // non-`Exact` unit's count semantics. Adding this to the prefix
    // and suffix lengths gives the tightest lower-bound length the
    // declared tract demands; any span shorter than that is incompatible
    // with EVERY count assignment and must reject. #428 (follow-up to
    // #395 item 1's `>` slack).
    //
    // Caveat: the HGVS spec
    // (`assets/hgvs-nomenclature/docs/recommendations/DNA/repeated.md`,
    // `syntax.yaml` `dna.rpt`) only formally defines `[n]` and
    // `[(lo_hi)]` for DNA repeat counts. ferro's parser additionally
    // accepts `[lo_hi]` (`Range`), `[lo_?]` (`MinUncertain`), `[?_hi]`
    // (`MaxUncertain`), and `[?]` (`Unknown`) as lenient extensions.
    // The per-variant minimums below are ferro's internal interpretation
    // of these out-of-spec uncertainty forms, not spec-mandated
    // semantics. Spec-defined forms (`Exact(n)`, `Range(lo, hi)`-style
    // explicit ranges visible in `[(lo_hi)]`) drive the tightening; the
    // lenient extensions are handled conservatively (Unknown /
    // MaxUncertain admit a 0-count interpretation; MinUncertain treats
    // `lo` as a hard floor).
    //
    // Per-variant minimum copy counts:
    //   - `Exact(n)`               → `n` (interior Exact units, sandwiched
    //                                 between non-Exact ones, contribute
    //                                 their exact bp).
    //   - `Range(lo, _)`           → `lo` (explicit fixed range; spec
    //                                 form `[lo_hi]`).
    //   - `UncertainRange(lo, _)`  → `lo` (uncertain-paren range
    //                                 `[(lo_hi)]`; `lo == 0` admits a
    //                                 0-count interpretation cleanly).
    //   - `MinUncertain(lo)`       → `lo` (`[lo_?]` — ferro convention:
    //                                 lower bound is a hard minimum).
    //   - `MaxUncertain(_)`        → `0` (`[?_hi]` — ferro convention:
    //                                 admits 0 copies).
    //   - `Unknown`                → `0` (`[?]` — ferro convention:
    //                                 admits 0 copies).
    debug_assert!(
        prefix_units + suffix_units <= units.len(),
        "prefix+suffix overlap invariant violated: prefix_units={} \
         suffix_units={} units.len()={}",
        prefix_units,
        suffix_units,
        units.len(),
    );
    let middle_units = &units[prefix_units..units.len() - suffix_units];
    let mut middle_min_bp: usize = 0;
    for u in middle_units {
        let min_copies = match u.count {
            RepeatCount::Exact(n) => n as usize,
            RepeatCount::Range(lo, _) => lo as usize,
            RepeatCount::UncertainRange(lo, _) => lo as usize,
            RepeatCount::MinUncertain(lo) => lo as usize,
            RepeatCount::MaxUncertain(_) | RepeatCount::Unknown => 0,
        };
        middle_min_bp = middle_min_bp.saturating_add(min_copies * u.sequence.bases().len());
    }

    // If prefix + middle_min + suffix together exceed the actual span,
    // the declared tract cannot fit at any middle-count assignment.
    // The previous `>` check (#395 item 1) accepted the `prefix +
    // suffix == span` boundary unconditionally; #428 tightens it by
    // requiring the middle's minimum bp to fit as well.
    //
    // The pin `interleaved_prefix_and_trailing_exact_consistent_accepted_strict`
    // (15-bp span = prefix(6) + middle_min(3, from `Range(1,3)`) + suffix(6))
    // continues to pass under the tightened check, because middle_min
    // exactly fills the gap.
    let declared_min_bp = prefix_bytes.len() + middle_min_bp + suffix_bytes.len();
    if declared_min_bp > actual_bytes.len() {
        let combined_str = if prefix_units == 0 {
            suffix_str.clone()
        } else if suffix_units == 0 {
            prefix_str.clone()
        } else {
            format!("{}…{}", prefix_str, suffix_str)
        };
        let actual_str: String = actual_bytes
            .iter()
            .map(|&b| (b as char).to_ascii_uppercase())
            .collect();
        // When middle_min_bp == 0 the overflow is purely prefix+suffix;
        // the "middle requires" note is vacuously true and confusing.
        // Mirror the conditional on both sides.
        let (declared_note, actual_note) = if middle_min_bp > 0 {
            (
                format!(
                    " (>= {} bp from declared anchored units; middle requires >= {} bp)",
                    declared_min_bp, middle_min_bp,
                ),
                format!(
                    " ({} bp; reference span too short for prefix+middle_min+suffix)",
                    actual_bytes.len()
                ),
            )
        } else {
            (
                format!(" ({} bp from declared anchored units)", declared_min_bp),
                format!(
                    " ({} bp; reference span too short for prefix+suffix)",
                    actual_bytes.len()
                ),
            )
        };
        return ValidationResult::mismatch(
            format!("{}{}", combined_str, declared_note),
            format!("{}{}", actual_str, actual_note),
        );
    }

    // Left-anchored prefix check (no-op when prefix_units == 0).
    if !prefix_bytes.is_empty() {
        let prefix_actual = &actual_bytes[..prefix_bytes.len()];
        let matches = prefix_bytes
            .iter()
            .zip(prefix_actual.iter())
            .all(|(a, b)| a.eq_ignore_ascii_case(b));
        if !matches {
            let expected_bases_str: String = prefix_bytes
                .iter()
                .map(|&b| (b as char).to_ascii_uppercase())
                .collect();
            let prefix_actual_str: String = prefix_actual
                .iter()
                .map(|&b| (b as char).to_ascii_uppercase())
                .collect();
            return ValidationResult::mismatch(
                format!("{} ({})", prefix_str, expected_bases_str),
                prefix_actual_str,
            );
        }
    }

    // Right-anchored suffix check (no-op when suffix_units == 0).
    if !suffix_bytes.is_empty() {
        let suffix_start = actual_bytes.len() - suffix_bytes.len();
        let suffix_actual = &actual_bytes[suffix_start..];
        let matches = suffix_bytes
            .iter()
            .zip(suffix_actual.iter())
            .all(|(a, b)| a.eq_ignore_ascii_case(b));
        if !matches {
            let expected_bases_str: String = suffix_bytes
                .iter()
                .map(|&b| (b as char).to_ascii_uppercase())
                .collect();
            let suffix_actual_str: String = suffix_actual
                .iter()
                .map(|&b| (b as char).to_ascii_uppercase())
                .collect();
            return ValidationResult::mismatch(
                format!("{} ({})", suffix_str, expected_bases_str),
                suffix_actual_str,
            );
        }
    }

    ValidationResult::ok()
}

/// Validate a single base against the reference
fn validate_single_base(stated: &Base, ref_seq: &[u8], pos: u64) -> ValidationResult {
    // Convert 1-based HGVS position to 0-based array index
    let idx = hgvs_pos_to_index(pos);

    if idx >= ref_seq.len() {
        return ValidationResult::mismatch(
            stated.to_char().to_string(),
            format!("(position {} out of range)", pos),
        );
    }

    let actual_byte = ref_seq[idx];
    let stated_byte = stated.to_u8();

    // Handle case-insensitive comparison
    if actual_byte.eq_ignore_ascii_case(&stated_byte) {
        ValidationResult::ok()
    } else {
        let actual_char = (actual_byte as char).to_ascii_uppercase();
        ValidationResult::mismatch(stated.to_char().to_string(), actual_char.to_string())
    }
}

/// Validate a sequence against the reference
fn validate_sequence(stated: &[Base], ref_seq: &[u8], start: u64, end: u64) -> ValidationResult {
    // Convert 1-based HGVS positions to 0-based array indices
    // start/end are inclusive in HGVS, so end_idx = end for half-open interval
    let start_idx = hgvs_pos_to_index(start);
    let end_idx = end as usize; // 1-based inclusive end = 0-based exclusive end

    if end_idx > ref_seq.len() {
        let stated_str: String = stated.iter().map(|b| b.to_char()).collect();
        return ValidationResult::mismatch(stated_str, format!("(position {} out of range)", end));
    }
    if start_idx > end_idx {
        let stated_str: String = stated.iter().map(|b| b.to_char()).collect();
        return ValidationResult::mismatch(
            stated_str,
            format!("(inverted range: start {} > end {})", start, end),
        );
    }

    let actual_bytes = &ref_seq[start_idx..end_idx];

    // Compare lengths first
    if stated.len() != actual_bytes.len() {
        let stated_str: String = stated.iter().map(|b| b.to_char()).collect();
        let actual_str: String = actual_bytes
            .iter()
            .map(|&b| (b as char).to_ascii_uppercase())
            .collect();
        return ValidationResult::mismatch(stated_str, actual_str);
    }

    // Compare each base
    for (stated_base, &actual_byte) in stated.iter().zip(actual_bytes.iter()) {
        let stated_byte = stated_base.to_u8();
        if !actual_byte.eq_ignore_ascii_case(&stated_byte) {
            let stated_str: String = stated.iter().map(|b| b.to_char()).collect();
            let actual_str: String = actual_bytes
                .iter()
                .map(|&b| (b as char).to_ascii_uppercase())
                .collect();
            return ValidationResult::mismatch(stated_str, actual_str);
        }
    }

    ValidationResult::ok()
}

/// Validate that an explicit numeric edit length (`del<N>`/`dup<N>`/`del<N>ins…`)
/// equals the position-interval span `[start, end]` (1-based inclusive). HGVS
/// requires `N == end - start + 1`; a disagreement is a `ReferenceMismatch`
/// (mutalyzer ELENGTHMISMATCH), independent of the reference bases.
fn validate_stated_length(declared: u64, start: u64, end: u64) -> ValidationResult {
    // Inverted ranges (start > end) are tolerated without panic, matching
    // `validate_repeat_tract` (the range gate elsewhere owns that diagnostic).
    // Guard before the u64 subtraction.
    if end < start {
        return ValidationResult::ok();
    }
    let span = end - start + 1;
    if declared == span {
        ValidationResult::ok()
    } else {
        ValidationResult::mismatch(
            format!("length {declared}"),
            format!("span {span} ({start}_{end})"),
        )
    }
}

/// Apply reference validation policy based on configuration.
///
/// Returns Ok(()) if validation passes or correction is allowed,
/// Returns Err if validation fails and strict mode rejects it.
///
/// Warnings are printed to stderr in lenient mode.
pub fn apply_validation_policy(
    result: &ValidationResult,
    config: &NormalizeConfig,
    variant_str: &str,
) -> Result<(), FerroError> {
    if result.valid {
        return Ok(());
    }

    let action = config.ref_mismatch_action();

    match action {
        ResolvedAction::Reject => Err(FerroError::ReferenceMismatch {
            location: variant_str.to_string(),
            expected: result.stated_ref.clone().unwrap_or_else(|| "?".to_string()),
            found: result.actual_ref.clone().unwrap_or_else(|| "?".to_string()),
        }),
        ResolvedAction::WarnCorrect => {
            // Print warning to stderr
            if let Some(ref warning) = result.warning {
                eprintln!("Warning: {} in '{}'", warning, variant_str);
            }
            Ok(())
        }
        ResolvedAction::SilentCorrect | ResolvedAction::Accept => {
            // Silently proceed
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::edit::{InsertedSequence, Sequence, UncertainDupExtent};
    use std::str::FromStr;

    #[test]
    fn test_validate_substitution_match() {
        let edit = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        };
        let ref_seq = b"ATGC";
        let result = validate_reference(&edit, ref_seq, 1, 1);
        assert!(result.valid);
    }

    #[test]
    fn test_validate_substitution_mismatch() {
        let edit = NaEdit::Substitution {
            reference: Base::G, // stated G
            alternative: Base::A,
        };
        let ref_seq = b"ATGC"; // actual is A at position 1
        let result = validate_reference(&edit, ref_seq, 1, 1);
        assert!(!result.valid);
        assert_eq!(result.stated_ref, Some("G".to_string()));
        assert_eq!(result.actual_ref, Some("A".to_string()));
    }

    #[test]
    fn test_validate_deletion_match() {
        let edit = NaEdit::Deletion {
            sequence: Some(Sequence::from_str("ATG").unwrap()),
            length: None,
        };
        let ref_seq = b"ATGC";
        let result = validate_reference(&edit, ref_seq, 1, 3);
        assert!(result.valid);
    }

    #[test]
    fn test_validate_deletion_mismatch() {
        let edit = NaEdit::Deletion {
            sequence: Some(Sequence::from_str("GGG").unwrap()),
            length: None,
        };
        let ref_seq = b"ATGC";
        let result = validate_reference(&edit, ref_seq, 1, 3);
        assert!(!result.valid);
        assert_eq!(result.stated_ref, Some("GGG".to_string()));
        assert_eq!(result.actual_ref, Some("ATG".to_string()));
    }

    #[test]
    fn test_validate_no_sequence() {
        // No stated sequence AND no stated length → nothing to validate.
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        let ref_seq = b"ATGC";
        let result = validate_reference(&edit, ref_seq, 1, 3);
        assert!(result.valid);
    }

    #[test]
    fn test_validate_delins_deleted_mismatch() {
        // `c.1_3delGGGinsAT`: stated deleted "GGG" but ref is "ATG".
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::from_str("AT").unwrap()),
            deleted: Some(Sequence::from_str("GGG").unwrap()),
            deleted_length: None,
        };
        let ref_seq = b"ATGC";
        let result = validate_reference(&edit, ref_seq, 1, 3);
        assert!(!result.valid);
        assert_eq!(result.stated_ref, Some("GGG".to_string()));
        assert_eq!(result.actual_ref, Some("ATG".to_string()));
    }

    #[test]
    fn test_validate_delins_deleted_match() {
        // `c.1_3delATGinsAT`: stated deleted "ATG" matches ref "ATG".
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::from_str("AT").unwrap()),
            deleted: Some(Sequence::from_str("ATG").unwrap()),
            deleted_length: None,
        };
        let ref_seq = b"ATGC";
        let result = validate_reference(&edit, ref_seq, 1, 3);
        assert!(result.valid);
    }

    #[test]
    fn test_validate_delins_no_deleted() {
        // Recommended short form `c.1_3delinsAT`: no stated deleted bases,
        // nothing to validate even when the span content differs.
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::from_str("AT").unwrap()),
            deleted: None,
            deleted_length: None,
        };
        let ref_seq = b"ATGC";
        let result = validate_reference(&edit, ref_seq, 1, 3);
        assert!(result.valid);
    }

    #[test]
    fn test_validate_delins_deleted_length_match() {
        // `c.1_3del3insAT`: deleted_length 3 == span 3, so it validates clean
        // (#486 checks deleted_length against the span; it is not compared as
        // bases). The mismatch case is covered separately below.
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::from_str("AT").unwrap()),
            deleted: None,
            deleted_length: Some(3),
        };
        let ref_seq = b"ATGC";
        let result = validate_reference(&edit, ref_seq, 1, 3);
        assert!(result.valid);
    }

    #[test]
    fn test_validate_case_insensitive() {
        let edit = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        };
        let ref_seq = b"atgc"; // lowercase
        let result = validate_reference(&edit, ref_seq, 1, 1);
        assert!(result.valid);
    }

    #[test]
    fn test_apply_policy_strict() {
        let result = ValidationResult::mismatch("G".to_string(), "A".to_string());
        let config = NormalizeConfig::strict();
        let err = apply_validation_policy(&result, &config, "c.1G>T");
        assert!(err.is_err());
    }

    #[test]
    fn test_apply_policy_lenient() {
        let result = ValidationResult::mismatch("G".to_string(), "A".to_string());
        let config = NormalizeConfig::lenient();
        let ok = apply_validation_policy(&result, &config, "c.1G>T");
        assert!(ok.is_ok()); // Should pass but emit warning
    }

    #[test]
    fn test_apply_policy_silent() {
        let result = ValidationResult::mismatch("G".to_string(), "A".to_string());
        let config = NormalizeConfig::silent();
        let ok = apply_validation_policy(&result, &config, "c.1G>T");
        assert!(ok.is_ok());
    }

    #[test]
    fn test_validate_sequence_inverted_range() {
        // Regression: inverted-range variants (start > end) like
        // NC_000011.10:g.5238138_5153222insTATTT must not panic.
        let edit = NaEdit::Deletion {
            sequence: Some(Sequence::from_str("ATG").unwrap()),
            length: None,
        };
        let ref_seq = b"ATGC";
        // start=3 end=1 → start_idx(2) > end_idx(1), should return mismatch not panic
        let result = validate_reference(&edit, ref_seq, 3, 1);
        assert!(!result.valid);
    }

    #[test]
    fn test_validate_deletion_length_mismatch() {
        // `c.45del4`: single position (span 1) but stated length 4.
        let edit = NaEdit::Deletion {
            sequence: None,
            length: Some(4),
        };
        let result = validate_reference(&edit, b"ATGC", 45, 45);
        assert!(!result.valid);
        assert_eq!(result.stated_ref, Some("length 4".to_string()));
    }

    #[test]
    fn test_validate_deletion_length_match() {
        // `c.1_4del4`: span 4 == stated 4.
        let edit = NaEdit::Deletion {
            sequence: None,
            length: Some(4),
        };
        let result = validate_reference(&edit, b"ATGC", 1, 4);
        assert!(result.valid);
    }

    #[test]
    fn test_validate_deletion_length_inverted_range_no_panic() {
        // start > end must be tolerated (guard), not panic on u64 underflow.
        let edit = NaEdit::Deletion {
            sequence: None,
            length: Some(4),
        };
        let result = validate_reference(&edit, b"ATGC", 3, 1);
        assert!(result.valid);
    }

    #[test]
    fn test_validate_duplication_length_mismatch() {
        let edit = NaEdit::Duplication {
            sequence: None,
            length: Some(4),
            uncertain_extent: None,
        };
        let result = validate_reference(&edit, b"ATGC", 45, 45);
        assert!(!result.valid);
        assert_eq!(result.stated_ref, Some("length 4".to_string()));
    }

    #[test]
    fn test_validate_duplication_uncertain_extent_skipped() {
        let edit = NaEdit::Duplication {
            sequence: None,
            length: Some(4),
            uncertain_extent: Some(UncertainDupExtent::Unknown),
        };
        let result = validate_reference(&edit, b"ATGC", 45, 45);
        assert!(result.valid);
    }

    #[test]
    fn test_validate_delins_deleted_length_mismatch() {
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::from_str("ATC").unwrap()),
            deleted: None,
            deleted_length: Some(4),
        };
        let result = validate_reference(&edit, b"ATGC", 45, 45);
        assert!(!result.valid);
        assert_eq!(result.stated_ref, Some("length 4".to_string()));
    }
}
