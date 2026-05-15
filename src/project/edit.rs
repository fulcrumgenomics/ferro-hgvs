//! Strand-aware transformation of `NaEdit` from g. to c./n.

use crate::hgvs::edit::{Base, InsertedSequence, NaEdit, Sequence};
use crate::reference::Strand;

/// Transform a g.-coordinate edit into the equivalent edit on the transcript strand.
///
/// On the plus strand, the edit is returned unchanged. On the minus strand, the
/// ref/alt bases (and any embedded sequences) are reverse-complemented so the
/// resulting edit reads correctly on the transcript's sense strand.
#[allow(dead_code)] // TODO(issue-200): called from projector in Task 7
pub(crate) fn transform_edit_for_strand(edit: &NaEdit, strand: Strand) -> NaEdit {
    if strand == Strand::Plus {
        return edit.clone();
    }
    match edit {
        NaEdit::Substitution {
            reference,
            alternative,
        } => NaEdit::Substitution {
            reference: complement_base(*reference),
            alternative: complement_base(*alternative),
        },
        NaEdit::SubstitutionNoRef { alternative } => NaEdit::SubstitutionNoRef {
            alternative: complement_base(*alternative),
        },
        NaEdit::Deletion { sequence, length } => NaEdit::Deletion {
            sequence: sequence.as_ref().map(revcomp_sequence),
            length: *length,
        },
        NaEdit::Insertion { sequence } => NaEdit::Insertion {
            sequence: revcomp_inserted(sequence),
        },
        NaEdit::Delins {
            sequence,
            deleted,
            deleted_length,
        } => NaEdit::Delins {
            sequence: revcomp_inserted(sequence),
            deleted: deleted.as_ref().map(revcomp_sequence),
            deleted_length: *deleted_length,
        },
        NaEdit::Duplication {
            sequence,
            length,
            uncertain_extent,
        } => NaEdit::Duplication {
            sequence: sequence.as_ref().map(revcomp_sequence),
            length: *length,
            uncertain_extent: uncertain_extent.clone(),
        },
        NaEdit::DupIns { sequence } => NaEdit::DupIns {
            sequence: revcomp_inserted(sequence),
        },
        NaEdit::Inversion { sequence, length } => NaEdit::Inversion {
            sequence: sequence.as_ref().map(revcomp_sequence),
            length: *length,
        },
        // Identity, Unknown, Conversion, Repeat, MultiRepeat, Methylation, CopyNumber,
        // Splice, NoProduct, PositionOnly: no revcomp needed at the edit level.
        other => other.clone(),
    }
}

/// Complement a single IUPAC base (no reversal — reversal is handled by revcomp_sequence).
fn complement_base(b: Base) -> Base {
    match b {
        Base::A => Base::T,
        Base::T => Base::A,
        Base::G => Base::C,
        Base::C => Base::G,
        Base::U => Base::A,
        Base::R => Base::Y,
        Base::Y => Base::R,
        Base::S => Base::S,
        Base::W => Base::W,
        Base::K => Base::M,
        Base::M => Base::K,
        Base::B => Base::V,
        Base::V => Base::B,
        Base::D => Base::H,
        Base::H => Base::D,
        Base::N => Base::N,
    }
}

/// Reverse-complement a `Sequence` value.
fn revcomp_sequence(s: &Sequence) -> Sequence {
    use crate::sequence::reverse_complement;
    reverse_complement(&s.to_string())
        .parse()
        .expect("reverse_complement produces only valid IUPAC bases")
}

/// Reverse-complement an `InsertedSequence`, handling only the `Literal` variant.
///
/// For non-literal variants (counts, ranges, repeats, named elements, etc.)
/// there is no embedded sequence to revcomp, so they are cloned unchanged.
fn revcomp_inserted(ins: &InsertedSequence) -> InsertedSequence {
    match ins {
        InsertedSequence::Literal(seq) => InsertedSequence::Literal(revcomp_sequence(seq)),
        other => other.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::edit::{InsertedSequence, NaEdit, Sequence};
    use crate::reference::Strand;

    fn seq(s: &str) -> Sequence {
        s.parse().unwrap()
    }

    #[test]
    fn substitution_plus_strand_unchanged() {
        let edit = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        };
        assert_eq!(transform_edit_for_strand(&edit, Strand::Plus), edit);
    }

    #[test]
    fn substitution_minus_strand_revcomps() {
        let edit = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        };
        assert_eq!(
            transform_edit_for_strand(&edit, Strand::Minus),
            NaEdit::Substitution {
                reference: Base::T,
                alternative: Base::C
            }
        );
    }

    #[test]
    fn deletion_minus_strand_revcomps_sequence() {
        let edit = NaEdit::Deletion {
            sequence: Some(seq("ATG")),
            length: None,
        };
        assert_eq!(
            transform_edit_for_strand(&edit, Strand::Minus),
            NaEdit::Deletion {
                sequence: Some(seq("CAT")),
                length: None
            }
        );
    }

    #[test]
    fn insertion_minus_strand_revcomps_sequence() {
        // InsertedSequence::Literal is the correct variant for literal sequences.
        let edit = NaEdit::Insertion {
            sequence: InsertedSequence::Literal(seq("ATG")),
        };
        let out = transform_edit_for_strand(&edit, Strand::Minus);
        assert!(
            out.to_string().contains("CAT"),
            "expected revcomp 'CAT' in display, got {}",
            out
        );
    }

    #[test]
    fn delins_minus_strand_revcomps_both() {
        // Plus-strand delATGinsCC → minus-strand del(revcomp ATG = CAT) ins(revcomp CC = GG).
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(seq("CC")),
            deleted: Some(seq("ATG")),
            deleted_length: None,
        };
        let out = transform_edit_for_strand(&edit, Strand::Minus);
        let s = out.to_string();
        assert!(
            s.contains("CAT"),
            "expected deleted 'CAT' in display, got {}",
            s
        );
        assert!(
            s.contains("GG"),
            "expected inserted 'GG' in display, got {}",
            s
        );
    }

    #[test]
    fn duplication_minus_strand_revcomps_sequence() {
        let edit = NaEdit::Duplication {
            sequence: Some(seq("ATG")),
            length: None,
            uncertain_extent: None,
        };
        let out = transform_edit_for_strand(&edit, Strand::Minus);
        let s = out.to_string();
        assert!(
            s.contains("CAT"),
            "expected revcomp 'CAT' in dup display, got {}",
            s
        );
    }

    #[test]
    fn inversion_minus_strand_length_unchanged() {
        // Length-only inversions don't carry sequence — should be returned as-is.
        let edit = NaEdit::Inversion {
            sequence: None,
            length: Some(5),
        };
        assert_eq!(transform_edit_for_strand(&edit, Strand::Minus), edit);
    }

    #[test]
    fn inversion_minus_strand_revcomps_sequence() {
        let edit = NaEdit::Inversion {
            sequence: Some(seq("ATG")),
            length: None,
        };
        let out = transform_edit_for_strand(&edit, Strand::Minus);
        let s = out.to_string();
        assert!(
            s.contains("CAT"),
            "expected revcomp 'CAT' in inv display, got {}",
            s
        );
    }
}
