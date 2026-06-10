//! Positionless breakpoint insertion within a composite cis allele (#546).
//!
//! HGVS `DNA/complex.md` describes a de novo pericentric inversion with a
//! deletion and a 1 bp insertion at the break point, where the insertion is
//! written WITHOUT a coordinate because the junction position is implied:
//! `NC_000002.12:g.[32310435_32310710del;32310711_171827243inv;insG]`.

use ferro_hgvs::hgvs::edit::NaEdit;
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::python_helpers::get_indel_length;

/// The canonical composite SV from `DNA/complex.md` round-trips, including the
/// positionless `insG` breakpoint member.
#[test]
fn composite_allele_with_positionless_insertion_round_trips() {
    let s = "NC_000002.12:g.[32310435_32310710del;32310711_171827243inv;insG]";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
}

/// A standalone positionless insertion is NOT spec-valid (the form only exists
/// inside a composite breakpoint context) and must stay rejected.
#[test]
fn standalone_positionless_insertion_is_rejected() {
    assert!(
        parse_hgvs("NC_000002.12:g.insG").is_err(),
        "standalone positionless `g.insG` must be rejected"
    );
}

/// A positionless `ins` is valid ONLY as a member of a multi-member cis
/// composite allele (an implied breakpoint between other edits). It must stay
/// rejected in every other bracket/join context — a lone single-member bracket,
/// a ring `::` segment, a trans `];[` member, and an unknown-phase `(;)` member.
#[test]
fn positionless_insertion_rejected_outside_multi_member_cis_group() {
    for s in [
        "NC_000022.11:g.[insG]",              // lone single-member bracket
        "NC_000022.11:g.100_200del::insG",    // ring `::` segment
        "NC_000002.12:g.[insG];[100_200del]", // trans `];[` member
        "NC_000002.12:g.[100_200del(;)insG]", // unknown-phase `(;)` member
    ] {
        assert!(
            parse_hgvs(s).is_err(),
            "positionless `ins` outside a multi-member cis group must be rejected: `{s}`"
        );
    }
}

/// A genuine multi-member cis breakpoint group with a positionless `ins`
/// round-trips even with only two members (the spec form has three, but the
/// "implied breakpoint between edits" rule holds for >= 2).
#[test]
fn two_member_cis_breakpoint_insertion_round_trips() {
    let s = "NC_000002.12:g.[100_200del;insG]";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
}

/// Helper: extract the `insG` sub-variant from the composite allele.
fn get_breakpoint_member() -> HgvsVariant {
    let s = "NC_000002.12:g.[32310435_32310710del;32310711_171827243inv;insG]";
    let allele = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    let HgvsVariant::Allele(a) = allele else {
        panic!("expected Allele, got {allele:?}");
    };
    a.variants
        .into_iter()
        .find(|v| {
            matches!(v, HgvsVariant::Genome(g)
                if matches!(g.loc_edit.edit.inner(), Some(NaEdit::BreakpointInsertion { .. })))
        })
        .expect("composite allele must contain a BreakpointInsertion member")
}

/// `NaEdit::is_positionless()` returns `true` for the breakpoint insertion
/// member and `false` for ordinary positioned edits.
#[test]
fn breakpoint_insertion_is_positionless() {
    let ins_member = get_breakpoint_member();
    let HgvsVariant::Genome(g) = &ins_member else {
        panic!("expected Genome variant");
    };
    let edit = g
        .loc_edit
        .edit
        .inner()
        .expect("breakpoint member must have a NaEdit");
    assert!(
        edit.is_positionless(),
        "BreakpointInsertion must report is_positionless() == true; got false for {edit:?}"
    );
}

/// `get_indel_length` returns `Some(1)` for the single-base `insG` breakpoint
/// insertion — i.e. the accessor no longer falls through to `_ => None`.
#[test]
fn breakpoint_insertion_indel_length_is_sequence_length() {
    let ins_member = get_breakpoint_member();
    let len = get_indel_length(&ins_member);
    assert_eq!(
        len,
        Some(1),
        "indel_length for `insG` breakpoint insertion must be Some(1); got {len:?}"
    );
}
