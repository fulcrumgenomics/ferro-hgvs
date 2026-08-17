//! Edit classification for HGVS generation.
//!
//! Classifies raw edits into HGVS-compatible types.

use super::align::RawEdit;
use crate::error::FerroError;

/// Classification of an edit for HGVS notation.
#[derive(Debug, Clone, PartialEq)]
pub enum EditClassification {
    /// Single nucleotide substitution (e.g., A>G).
    Substitution {
        /// Reference base.
        ref_base: char,
        /// Alternate base.
        alt_base: char,
    },
    /// Deletion of one or more bases.
    Deletion {
        /// Deleted sequence.
        deleted: String,
    },
    /// Insertion of one or more bases.
    Insertion {
        /// Inserted sequence.
        inserted: String,
    },
    /// Deletion-insertion (replacement).
    Delins {
        /// Deleted sequence.
        deleted: String,
        /// Inserted sequence.
        inserted: String,
    },
    /// Duplication of adjacent sequence.
    Duplication {
        /// Duplicated sequence.
        sequence: String,
    },
    /// Inversion of sequence.
    Inversion {
        /// Inverted sequence.
        sequence: String,
    },
    /// Repeated sequence with count.
    Repeat {
        /// Repeat unit.
        unit: String,
        /// Number of repeats.
        count: u32,
    },
}

/// Classify a raw edit into an HGVS-compatible type.
///
/// # Arguments
///
/// * `edit` - The raw edit to classify
/// * `reference` - The full reference sequence (for duplication detection)
/// * `detect_dup` - Whether to detect duplications
/// * `detect_inv` - Whether to detect inversions
pub fn classify_edit(
    edit: &RawEdit,
    reference: &str,
    detect_dup: bool,
    detect_inv: bool,
) -> Result<EditClassification, FerroError> {
    let ref_len = edit.ref_seq.len();
    let obs_len = edit.obs_seq.len();

    let classification = match (ref_len, obs_len) {
        // Single substitution
        (1, 1) => EditClassification::Substitution {
            ref_base: edit.ref_seq.chars().next().unwrap(),
            alt_base: edit.obs_seq.chars().next().unwrap(),
        },

        // Deletion
        (r, 0) if r > 0 => EditClassification::Deletion {
            deleted: edit.ref_seq.clone(),
        },

        // Insertion - check for duplication
        (0, a) if a > 0 => {
            if detect_dup && is_duplication(reference, edit.ref_start, &edit.obs_seq)? {
                EditClassification::Duplication {
                    sequence: edit.obs_seq.clone(),
                }
            } else {
                EditClassification::Insertion {
                    inserted: edit.obs_seq.clone(),
                }
            }
        }

        // Delins - check for inversion
        (r, a) if r > 0 && a > 0 => {
            if detect_inv && r == a && is_inversion(&edit.ref_seq, &edit.obs_seq) {
                EditClassification::Inversion {
                    sequence: edit.ref_seq.clone(),
                }
            } else if r == 1 && a == 1 {
                // Single base change
                EditClassification::Substitution {
                    ref_base: edit.ref_seq.chars().next().unwrap(),
                    alt_base: edit.obs_seq.chars().next().unwrap(),
                }
            } else {
                EditClassification::Delins {
                    deleted: edit.ref_seq.clone(),
                    inserted: edit.obs_seq.clone(),
                }
            }
        }

        _ => {
            // Fallback - shouldn't happen
            EditClassification::Delins {
                deleted: edit.ref_seq.clone(),
                inserted: edit.obs_seq.clone(),
            }
        }
    };

    Ok(classification)
}

/// Check if an insertion is a duplication of the preceding sequence.
///
/// According to HGVS, an insertion is a duplication if the inserted sequence
/// is identical to the sequence immediately 5' of the insertion point.
///
/// `position` is the alignment's 1-based `ref_start` for the insertion, so a
/// trailing insertion reports `position == reference.len() + 1`. Reading the
/// preceding bases at `reference[pos - ins_len..pos]` would then run past the
/// end of the reference and panic (#2128), which reached a user-facing CLI
/// entry point (`ferro describe`) as process exit 101 on well-formed input.
/// The out-of-range case now returns a typed [`FerroError::InvalidCoordinates`]
/// so the caller can decline cleanly instead.
pub fn is_duplication(reference: &str, position: u64, inserted: &str) -> Result<bool, FerroError> {
    if inserted.is_empty() {
        return Ok(false);
    }

    let pos = position as usize;
    let ins_len = inserted.len();

    // The insertion point must lie within the reference for the preceding-bases
    // read below to be in bounds. A point past the end is an unresolvable
    // coordinate, not a duplication — decline rather than panic.
    if pos > reference.len() {
        return Err(FerroError::InvalidCoordinates {
            msg: format!(
                "insertion point {pos} is past the end of the {}-base reference",
                reference.len()
            ),
        });
    }

    // Check if there's enough sequence before the insertion point
    if pos < ins_len {
        return Ok(false);
    }

    // Get the sequence immediately before the insertion point. Both `pos <=
    // reference.len()` and `pos >= ins_len` hold here, so the slice is in bounds.
    let preceding = &reference[pos - ins_len..pos];

    Ok(preceding == inserted)
}

/// Whether an equal-length delins is really an inversion — its replacement is
/// the reverse complement of what it deletes.
///
/// Three things this has to get right, and a local complement got all three
/// wrong. `classify_edit` calls it for *every* equal-length delins, so each was
/// a live misclassification rather than a latent one.
///
/// * **More than one nucleotide.** `inversion.md:15-16`: "the region inverted
///   contains **more than one nucleotide**. The description `g.234inv` is
///   therefore not allowed; a one-nucleotide inversion should be described as a
///   substitution." Without the guard every `A>T`, `T>A`, `G>C` and `C>G` came
///   back as an inversion instead of reaching the substitution arm below.
/// * **A code with no complement is not self-complementary.** The local helper's
///   fallback returned the character unchanged, so `complement(R) == R` held for
///   every code it did not model and `RR` to `RR` — no change at all — read as
///   an inversion.
/// * **Case.** Reference FASTAs are routinely soft-masked, so one inversion can
///   arrive as `ATG`/`cat` or `atg`/`CAT` and must classify the same way.
///
/// Built on [`crate::sequence::complement_base`], the crate's *comparison*
/// complement: it folds case and returns `None` for a byte it does not model, so
/// `X` cannot pass as its own complement. Deliberately **not**
/// [`crate::sequence::reverse_complement`], which is the *display* helper — it
/// preserves case and passes unmodelled characters through, because its job is
/// to render a sequence a caller reads back. Using it here would answer the `X`
/// case wrongly, which is the same conflation the local helper made.
///
/// Compared byte-wise against the reversed input so neither the fold nor the
/// reversal allocates.
pub fn is_inversion(deleted: &str, inserted: &str) -> bool {
    let (deleted, inserted) = (deleted.as_bytes(), inserted.as_bytes());
    if deleted.len() != inserted.len() || deleted.len() < 2 {
        return false;
    }
    deleted.iter().rev().zip(inserted).all(|(&base, &other)| {
        crate::sequence::complement_base(base).is_some_and(|c| c == other.to_ascii_uppercase())
    })
}

/// Find a repeat unit in a sequence.
///
/// Returns the unit and count if the sequence is composed of repeating units.
#[allow(dead_code)]
pub fn find_repeat_unit(sequence: &str) -> Option<(String, u32)> {
    if sequence.is_empty() {
        return None;
    }

    let len = sequence.len();

    // Try unit lengths from 1 to half the sequence length
    for unit_len in 1..=len / 2 {
        if !len.is_multiple_of(unit_len) {
            continue;
        }

        let unit = &sequence[0..unit_len];
        let count = len / unit_len;

        // Check if entire sequence is this unit repeated
        let is_repeat = (0..count).all(|i| {
            let start = i * unit_len;
            &sequence[start..start + unit_len] == unit
        });

        if is_repeat {
            return Some((unit.to_string(), count as u32));
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_edit(ref_start: u64, ref_seq: &str, obs_seq: &str) -> RawEdit {
        RawEdit {
            ref_start,
            ref_end: ref_start + ref_seq.len() as u64,
            obs_start: ref_start,
            obs_end: ref_start + obs_seq.len() as u64,
            ref_seq: ref_seq.to_string(),
            obs_seq: obs_seq.to_string(),
        }
    }

    #[test]
    fn test_classify_substitution() {
        let edit = make_edit(3, "G", "A");
        let class = classify_edit(&edit, "ATGC", true, true).unwrap();

        match class {
            EditClassification::Substitution { ref_base, alt_base } => {
                assert_eq!(ref_base, 'G');
                assert_eq!(alt_base, 'A');
            }
            _ => panic!("Expected Substitution"),
        }
    }

    #[test]
    fn test_classify_deletion() {
        let edit = make_edit(3, "GC", "");
        let class = classify_edit(&edit, "ATGCAT", true, true).unwrap();

        match class {
            EditClassification::Deletion { deleted } => {
                assert_eq!(deleted, "GC");
            }
            _ => panic!("Expected Deletion"),
        }
    }

    #[test]
    fn test_classify_insertion() {
        let edit = make_edit(3, "", "T");
        let class = classify_edit(&edit, "ATGC", true, true).unwrap();

        match class {
            EditClassification::Insertion { inserted } => {
                assert_eq!(inserted, "T");
            }
            _ => panic!("Expected Insertion"),
        }
    }

    #[test]
    fn test_classify_duplication() {
        // Insert "G" at position 3 (1-based) in "ATGC"
        // Position 3 is 'G', so inserting "G" after it is a duplication
        // The preceding character (at index 2, position 3) is 'G'
        let edit = make_edit(3, "", "G");
        let class = classify_edit(&edit, "ATGC", true, true).unwrap();

        match class {
            EditClassification::Duplication { sequence } => {
                assert_eq!(sequence, "G");
            }
            _ => panic!("Expected Duplication, got {:?}", class),
        }
    }

    #[test]
    fn test_classify_delins() {
        let edit = make_edit(3, "GC", "TT");
        let class = classify_edit(&edit, "ATGCAT", true, true).unwrap();

        match class {
            EditClassification::Delins { deleted, inserted } => {
                assert_eq!(deleted, "GC");
                assert_eq!(inserted, "TT");
            }
            _ => panic!("Expected Delins"),
        }
    }

    #[test]
    fn test_is_duplication() {
        // "ATGCAT" - inserting "CAT" at position 7 (after the full sequence)
        // Preceding "CAT" (positions 4-6) matches
        assert!(is_duplication("ATGCAT", 6, "CAT").unwrap());

        // Not a duplication - different sequence
        assert!(!is_duplication("ATGCAT", 6, "TTT").unwrap());

        // Not enough preceding sequence
        assert!(!is_duplication("AT", 2, "ATGC").unwrap());
    }

    /// An insertion point past the end of the reference is an unresolvable
    /// coordinate and must decline with a typed error rather than panic on the
    /// preceding-bases slice (#2128).
    #[test]
    fn is_duplication_past_the_reference_end_is_a_typed_error() {
        // pos = reference.len() + 1, as the alignment reports for a trailing
        // insertion.
        match is_duplication("ATG", 4, "G") {
            Err(FerroError::InvalidCoordinates { .. }) => {}
            other => panic!("expected InvalidCoordinates, got {other:?}"),
        }

        // The boundary `pos == reference.len()` stays in bounds and resolves.
        assert!(is_duplication("ATG", 3, "G").unwrap());
    }

    #[test]
    fn test_is_inversion() {
        // "ATG" inverted is "CAT" (reverse complement)
        assert!(is_inversion("ATG", "CAT"));

        // Not an inversion
        assert!(!is_inversion("ATG", "TTT"));

        // Different lengths
        assert!(!is_inversion("ATG", "CA"));
    }

    /// A one-nucleotide "inversion" is a substitution (`inversion.md:15-16`:
    /// "the region inverted contains **more than one nucleotide**. The
    /// description `g.234inv` is therefore not allowed").
    ///
    /// `classify_edit` reaches [`is_inversion`] for *any* equal-length delins,
    /// so without a length guard every `A>T`, `T>A`, `G>C` and `C>G` classified
    /// as an inversion instead of falling through to its own substitution arm.
    #[test]
    fn a_single_nucleotide_is_never_an_inversion() {
        for (deleted, inserted) in [("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")] {
            assert!(
                !is_inversion(deleted, inserted),
                "{deleted}>{inserted} is a substitution, not a 1 nt inversion"
            );
        }
    }

    /// A code with no modelled complement must not read as its own complement.
    ///
    /// This is #1249's trap in another module: a complement whose fallback
    /// returns the character unchanged makes `complement(R) == R` for every code
    /// it does not model, so `RR` to `RR` — no change at all — classified as an
    /// inversion. `R` complements to `Y`.
    #[test]
    fn a_code_without_a_complement_is_not_its_own_complement() {
        assert!(!is_inversion("RR", "RR"), "revcomp(RR) is YY, not RR");
        assert!(!is_inversion("XX", "XX"), "X has no complement to report");
        assert!(is_inversion("RY", "RY"), "revcomp(RY) really is RY");
    }

    /// The third class, and the one the two cases above cannot show: a code that
    /// genuinely *is* its own complement.
    ///
    /// `S` (G/C), `W` (A/T) and `N` all complement to themselves in
    /// [`crate::sequence::complement_base`], so an unchanged run of them really
    /// is an inversion and must be reported as one. Without this, a "fix" that
    /// simply refused every ambiguity code would satisfy
    /// `a_code_without_a_complement_is_not_its_own_complement` and be wrong —
    /// the three cases together pin unmodelled (`X`), modelled-but-not-self-
    /// complementary (`R`), and modelled-and-self-complementary as distinct.
    #[test]
    fn a_self_complementary_code_is_an_inversion_unchanged() {
        // Homogeneous runs only. Self-complementary is not the same as
        // palindromic: each base complementing to itself makes `revcomp` equal to
        // the plain *reverse*, so an unchanged run is an inversion exactly when
        // the run reads the same backwards.
        for code in ["SS", "WW", "NN"] {
            assert!(
                is_inversion(code, code),
                "`{code}` is a self-complementary palindrome, so it is its own \
                 reverse complement"
            );
        }
        // `SW` shows the distinction: both bases are self-complementary, but the
        // run is not a palindrome, so `revcomp(SW)` is `WS`.
        assert!(is_inversion("SW", "WS"), "revcomp(SW) is WS");
        assert!(
            !is_inversion("SW", "SW"),
            "revcomp(SW) is WS, so SW unchanged is not an inversion"
        );
    }

    /// Reference FASTAs are routinely soft-masked, so one inversion can arrive
    /// with either case on either side and must classify the same way.
    #[test]
    fn soft_masking_does_not_hide_an_inversion() {
        assert!(is_inversion("ATG", "cat"));
        assert!(is_inversion("atg", "CAT"));
        assert!(is_inversion("aTg", "cAt"));
    }

    #[test]
    fn test_find_repeat_unit() {
        assert_eq!(find_repeat_unit("ATAT"), Some(("AT".to_string(), 2)));
        assert_eq!(find_repeat_unit("AAA"), Some(("A".to_string(), 3)));
        assert_eq!(find_repeat_unit("ATCATCATC"), Some(("ATC".to_string(), 3)));
        assert_eq!(find_repeat_unit("ATGC"), None); // No repeating unit
        assert_eq!(find_repeat_unit(""), None);
    }
}
