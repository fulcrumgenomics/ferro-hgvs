//! Issue #395 item 4 — RNA `u` to DNA `t` translation on plus-strand
//! projection.
//!
//! `src/project/edit.rs::transform_edit_for_strand` previously
//! returned `edit.clone()` unchanged on `Strand::Plus`. For r. inputs
//! (RNA-coordinate variants) that carry `Base::U`, projecting to the
//! genomic axis (DNA) requires translating `u → t` independently of
//! strand: the genomic reference is DNA, not RNA. The `Base::U →
//! Base::A` complement map at `src/project/edit.rs:71` only ran via
//! `complement_base` on the `Strand::Minus` path, leaving plus-strand
//! r.→g. projections emitting invalid DNA like `g.<N>U>A`.
//!
//! # Spec basis
//!
//! `assets/hgvs-nomenclature/docs/recommendations/general.md` defines
//! r./n./c. variants on RNA sequences with lowercase nucleotides
//! including `u`. The DNA reference sequences (`g.`, `m.`, `o.`) use
//! `T`. A round-trip r.→g. projection must therefore translate `u`
//! to `t` regardless of strand orientation.

use ferro_hgvs::hgvs::edit::{Base, InsertedSequence, NaEdit, Sequence};
use ferro_hgvs::project::edit::transform_edit_for_strand;
use ferro_hgvs::reference::Strand;

#[test]
fn plus_strand_substitution_u_to_a_translates_to_t_to_a_dna() {
    // r.100u>a — RNA `u` substitution to `a`. On plus strand projecting
    // to g. (DNA), the `u` becomes `t` (no reverse complement on Plus).
    let r_edit = NaEdit::Substitution {
        reference: Base::U,
        alternative: Base::A,
    };
    let g_edit = transform_edit_for_strand(&r_edit, Strand::Plus);
    assert_eq!(
        g_edit,
        NaEdit::Substitution {
            reference: Base::T,
            alternative: Base::A,
        },
        "plus-strand r.u>a must project to g.T>A (DNA), not g.U>A (invalid DNA)",
    );
}

#[test]
fn plus_strand_substitution_a_to_u_translates_alt_u_to_t() {
    // r.100a>u — RNA `a` substitution to `u`. Plus-strand projection:
    // alt `u` → `t`.
    let r_edit = NaEdit::Substitution {
        reference: Base::A,
        alternative: Base::U,
    };
    let g_edit = transform_edit_for_strand(&r_edit, Strand::Plus);
    assert_eq!(
        g_edit,
        NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::T,
        },
        "plus-strand r.a>u must project to g.A>T (DNA)",
    );
}

#[test]
fn plus_strand_substitution_non_u_unchanged() {
    // Negative control: a substitution with no `u` is returned
    // bit-for-bit identical on plus strand.
    let r_edit = NaEdit::Substitution {
        reference: Base::A,
        alternative: Base::G,
    };
    let g_edit = transform_edit_for_strand(&r_edit, Strand::Plus);
    assert_eq!(
        g_edit, r_edit,
        "plus-strand non-u edit must pass through unchanged"
    );
}

#[test]
fn plus_strand_deletion_u_in_sequence_translates_to_t() {
    // r.100_102delauu — deleted sequence carries `u`. Plus-strand
    // projection: the stated deleted DNA sequence is `ATT`.
    let r_edit = NaEdit::Deletion {
        sequence: Some("AUU".parse::<Sequence>().unwrap()),
        length: None,
    };
    let g_edit = transform_edit_for_strand(&r_edit, Strand::Plus);
    let expected_seq: Sequence = "ATT".parse().unwrap();
    match g_edit {
        NaEdit::Deletion {
            sequence: Some(s), ..
        } => assert_eq!(
            s, expected_seq,
            "plus-strand del sequence with U must translate U→T"
        ),
        other => panic!("expected Deletion with sequence, got {other:?}"),
    }
}

#[test]
fn plus_strand_insertion_u_in_inserted_translates_to_t() {
    // r.100_101insauu — inserted sequence carries `u`. Plus-strand
    // projection: `ATT`.
    let r_edit = NaEdit::Insertion {
        sequence: InsertedSequence::Literal("AUU".parse::<Sequence>().unwrap()),
    };
    let g_edit = transform_edit_for_strand(&r_edit, Strand::Plus);
    let expected_seq: Sequence = "ATT".parse().unwrap();
    match g_edit {
        NaEdit::Insertion {
            sequence: InsertedSequence::Literal(s),
        } => assert_eq!(
            s, expected_seq,
            "plus-strand ins sequence with U must translate U→T"
        ),
        other => panic!("expected Insertion with Literal sequence, got {other:?}"),
    }
}

#[test]
fn plus_strand_delins_u_in_both_translates_both() {
    // r.100_102delauuinscu — both deleted and inserted sequences have
    // `u`. Plus-strand: deleted `ATT`, inserted `CT`.
    let r_edit = NaEdit::Delins {
        sequence: InsertedSequence::Literal("CU".parse::<Sequence>().unwrap()),
        deleted: Some("AUU".parse::<Sequence>().unwrap()),
        deleted_length: None,
    };
    let g_edit = transform_edit_for_strand(&r_edit, Strand::Plus);
    let expected_del: Sequence = "ATT".parse().unwrap();
    let expected_ins: Sequence = "CT".parse().unwrap();
    match g_edit {
        NaEdit::Delins {
            sequence: InsertedSequence::Literal(ins),
            deleted: Some(del),
            ..
        } => {
            assert_eq!(del, expected_del, "delins deleted U→T");
            assert_eq!(ins, expected_ins, "delins inserted U→T");
        }
        other => panic!("expected Delins, got {other:?}"),
    }
}

#[test]
fn minus_strand_u_in_substitution_still_revcomps() {
    // Regression for the minus-strand path: r.U>A on minus strand
    // becomes g.A>T (U complement = A, A complement = T). The
    // existing complement_base map already handles U → A; this test
    // just pins that the minus-strand path is unaffected by the
    // item 4 fix.
    let r_edit = NaEdit::Substitution {
        reference: Base::U,
        alternative: Base::A,
    };
    let g_edit = transform_edit_for_strand(&r_edit, Strand::Minus);
    assert_eq!(
        g_edit,
        NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::T,
        },
        "minus-strand r.u>a must project to g.A>T (existing complement_base behavior)",
    );
}
