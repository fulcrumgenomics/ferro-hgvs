//! Strand-aware transformation of `NaEdit` from g. to c./n.

use crate::hgvs::edit::{Base, NaEdit, Sequence};
use crate::reference::Strand;

/// Transform a g.-coordinate edit into the equivalent edit on the transcript strand.
///
/// On the plus strand, the edit is returned unchanged. On the minus strand, the
/// ref/alt bases (and any embedded sequences) are reverse-complemented so the
/// resulting edit reads correctly on the transcript's sense strand.
pub(crate) fn transform_edit_for_strand(edit: &NaEdit, strand: Strand) -> NaEdit {
    // Implementation comes next (Task 4).
    let _ = strand;
    edit.clone()
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
