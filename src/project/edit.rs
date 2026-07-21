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
            substitution_reference,
        } => NaEdit::Delins {
            sequence: revcomp_inserted(sequence),
            deleted: deleted.as_ref().map(revcomp_sequence),
            deleted_length: *deleted_length,
            substitution_reference: substitution_reference.as_ref().map(revcomp_sequence),
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
        // A tandem-repeat unit reads on the transcript's sense strand, so on a
        // minus-strand transcript it must be reverse-complemented to land on the
        // genomic reference (e.g. `CA[5]` → `TG[5]`, `T[8]` → `A[8]`) — #852.
        // The repeat counts are orientation-invariant and carry through unchanged.
        //
        // CAVEAT: `trailing` AND `additional_counts` are reverse-complemented in
        // place only. On the opposite strand a *trailing* partial unit becomes a
        // *leading* one, and a mixed repeat's unit *order* reverses — neither is
        // re-ordered here. Latent-only: `normalize` bails when `trailing.is_some()`
        // or `additional_counts` is non-empty (normalize/mod.rs), and no projected
        // repeat in scope has them. The first mixed-repeat minus-strand projection
        // must revisit this (#852 notes). Mirrors the in-place handling in
        // `u_to_t_edit`.
        NaEdit::Repeat {
            sequence,
            count,
            additional_counts,
            trailing,
        } => NaEdit::Repeat {
            sequence: sequence.as_ref().map(revcomp_sequence),
            count: count.clone(),
            additional_counts: additional_counts.clone(),
            trailing: trailing.as_ref().map(revcomp_sequence),
        },
        // Identity, Unknown, Conversion, MultiRepeat, Methylation, CopyNumber,
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
///
/// KNOWN GAP (#856): a multi-part `Complex` payload whose parts include literals
/// (e.g. `ins[ATG;CCC]`) is cloned unchanged here, so its literal members are
/// left in transcript orientation on minus-strand projection. A *single*
/// bracketed literal no longer reaches this path — it is collapsed to `Literal`
/// at parse (see `parser::edit::parse_bracketed_inserted_sequence`) — so the
/// common case is correct; the residual multi-part case is a separate, pre-
/// existing limitation. Reverse-complementing a `Complex` correctly would also
/// require reversing the part order and projecting any position-range members,
/// which is out of scope here.
fn revcomp_inserted(ins: &InsertedSequence) -> InsertedSequence {
    match ins {
        InsertedSequence::Literal(seq) => InsertedSequence::Literal(revcomp_sequence(seq)),
        other => other.clone(),
    }
}

/// U→T translation for an edit. Translates RNA `U` to DNA `T` in single bases
/// and embedded sequences without reversing or complementing other bases.
///
/// Used by `transform_edit_for_strand` on `Strand::Plus` (closes #395 item 4)
/// and by the RNA→VCF converter, which lowers `r.` edits onto DNA-coordinate
/// records and must not leak uracil into the VCF REF/ALT.
pub(crate) fn u_to_t_edit(edit: &NaEdit) -> NaEdit {
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
            substitution_reference,
        } => NaEdit::Delins {
            sequence: u_to_t_inserted(sequence),
            deleted: deleted.as_ref().map(u_to_t_sequence),
            deleted_length: *deleted_length,
            substitution_reference: substitution_reference.as_ref().map(u_to_t_sequence),
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
        NaEdit::Repeat {
            sequence,
            count,
            additional_counts,
            trailing,
        } => NaEdit::Repeat {
            sequence: sequence.as_ref().map(u_to_t_sequence),
            count: count.clone(),
            additional_counts: additional_counts.clone(),
            trailing: trailing.as_ref().map(u_to_t_sequence),
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
    // `Base::to_char()` emits uppercase, so `s.to_string()` is always
    // uppercase. The lowercase `'u'` arm is defensive only — kept in
    // case a future `Sequence::from_str` accepts lowercase input
    // without normalizing.
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

// Same KNOWN GAP as `revcomp_inserted` (#856): literal members of a multi-part
// `Complex` payload are not U→T translated; a single bracketed literal is
// collapsed to `Literal` at parse, so the common case is handled.
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
    fn repeat_minus_strand_revcomps_unit() {
        use crate::hgvs::edit::RepeatCount;
        // A tandem-repeat unit must be reverse-complemented on the minus strand:
        // revcomp(CA) = TG (#852). Counts are orientation-invariant.
        let edit = NaEdit::Repeat {
            sequence: Some(seq("CA")),
            count: RepeatCount::Exact(5),
            additional_counts: Vec::new(),
            trailing: None,
        };
        match transform_edit_for_strand(&edit, Strand::Minus) {
            NaEdit::Repeat {
                sequence, count, ..
            } => {
                assert_eq!(sequence.unwrap().to_string(), "TG");
                assert_eq!(count, RepeatCount::Exact(5));
            }
            other => panic!("expected Repeat, got {other:?}"),
        }
    }

    #[test]
    fn substitution_plus_strand_non_u_unchanged() {
        // A>G has no U, so plus-strand pass-through is bit-identical.
        // (With U the plus-strand path translates U→T — covered by
        // `tests/issue_395_rna_u_to_t_plus_strand.rs`.)
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
            substitution_reference: None,
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
