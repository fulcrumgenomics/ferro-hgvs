//! Strand-aware transformation of `NaEdit` between g. and c./n./r.

use crate::hgvs::edit::{Base, InsertedSequence, NaEdit, Sequence};
use crate::reference::Strand;

/// Transform an edit between the genomic axis and the transcript axis.
///
/// On the **plus strand**, no reverse-complement is needed — the
/// transcript reads the same as the genome. However, r. inputs carry
/// `Base::U` which is not a valid DNA letter; when projecting
/// r.→g. we still need to translate every `U` to `T` so the output
/// is valid DNA. Closes #395 item 4: previously the plus-strand path
/// returned the edit unchanged, which left RNA `u` substitutions
/// emitting invalid `g.<N>U>A` shapes. Non-U bases are passed
/// through unchanged.
///
/// On the **minus strand**, ref/alt bases (and any embedded sequences)
/// are reverse-complemented. The `Base::U → Base::A` complement in
/// `complement_base` already handles RNA `u` correctly on this path.
pub fn transform_edit_for_strand(edit: &NaEdit, strand: Strand) -> NaEdit {
    if strand == Strand::Plus {
        return u_to_t_edit(edit);
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

/// U→T translation for an edit, used by `transform_edit_for_strand` on
/// `Strand::Plus`. Translates RNA `U` to DNA `T` in single bases and
/// embedded sequences without reversing or complementing other bases.
/// Closes #395 item 4.
fn u_to_t_edit(edit: &NaEdit) -> NaEdit {
    match edit {
        NaEdit::Substitution {
            reference,
            alternative,
        } => NaEdit::Substitution {
            reference: u_to_t_base(*reference),
            alternative: u_to_t_base(*alternative),
        },
        NaEdit::SubstitutionNoRef { alternative } => NaEdit::SubstitutionNoRef {
            alternative: u_to_t_base(*alternative),
        },
        NaEdit::Deletion { sequence, length } => NaEdit::Deletion {
            sequence: sequence.as_ref().map(u_to_t_sequence),
            length: *length,
        },
        NaEdit::Insertion { sequence } => NaEdit::Insertion {
            sequence: u_to_t_inserted(sequence),
        },
        NaEdit::Delins {
            sequence,
            deleted,
            deleted_length,
        } => NaEdit::Delins {
            sequence: u_to_t_inserted(sequence),
            deleted: deleted.as_ref().map(u_to_t_sequence),
            deleted_length: *deleted_length,
        },
        NaEdit::Duplication {
            sequence,
            length,
            uncertain_extent,
        } => NaEdit::Duplication {
            sequence: sequence.as_ref().map(u_to_t_sequence),
            length: *length,
            uncertain_extent: uncertain_extent.clone(),
        },
        NaEdit::DupIns { sequence } => NaEdit::DupIns {
            sequence: u_to_t_inserted(sequence),
        },
        NaEdit::Inversion { sequence, length } => NaEdit::Inversion {
            sequence: sequence.as_ref().map(u_to_t_sequence),
            length: *length,
        },
        // Other edit shapes have no embedded sequence to translate.
        other => other.clone(),
    }
}

fn u_to_t_base(b: Base) -> Base {
    match b {
        Base::U => Base::T,
        other => other,
    }
}

fn u_to_t_sequence(s: &Sequence) -> Sequence {
    let translated: String = s
        .to_string()
        .chars()
        .map(|c| match c {
            'U' => 'T',
            'u' => 't',
            other => other,
        })
        .collect();
    translated
        .parse()
        .expect("U→T translation produces only valid IUPAC bases")
}

fn u_to_t_inserted(ins: &InsertedSequence) -> InsertedSequence {
    match ins {
        InsertedSequence::Literal(seq) => InsertedSequence::Literal(u_to_t_sequence(seq)),
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
