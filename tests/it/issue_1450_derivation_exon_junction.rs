//! #1450 — the sequence-first derivation must not move a member across an
//! exon/exon junction that the per-member pipeline would clamp.
//!
//! `general.md:44` exempts deletions and duplications around exon/exon junctions
//! from the 3' rule on the `c.`, `r.` and `n.` axes. The per-member pipeline
//! honours it. The sequence-first cis path did not, so the *same* member was
//! clamped when described alone and shifted across two junctions when described
//! beside a sibling — a confluence defect in which whether a member is clamped
//! depended on whether it had a sibling.
//!
//! `NM_001234.1` from `MockProvider::with_test_data` has exons at transcript
//! 1-15 / 16-30 / 31-44 and `cds_start = 5`, with a G-run spanning c.9-c.33 —
//! one homopolymer crossing both junctions, which is what makes the shift
//! reachable at all.

use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

fn normalized(input: &str) -> String {
    let normalizer = Normalizer::new(MockProvider::with_test_data());
    let variant = parse_hgvs(input).expect("input parses");
    normalizer
        .normalize(&variant)
        .expect("normalization succeeds")
        .to_string()
}

/// The lone spellings, which never left the per-member pipeline and are the
/// reference answer the pair has to agree with. Pinned here so a regression in
/// the clamp itself cannot be mistaken for a change in the pair below.
#[test]
fn a_lone_member_still_clamps_at_its_exon_junction() {
    assert_eq!(normalized("NM_001234.1:c.10del"), "NM_001234.1:c.11del");
    assert_eq!(normalized("NM_001234.1:c.13del"), "NM_001234.1:c.26del");
}

/// The defect. Each member must land where it lands alone, rather than both
/// being carried past their junctions and merged into one 2-base deletion in the
/// third exon.
#[test]
fn a_cis_pair_does_not_shift_its_members_across_exon_junctions() {
    assert_eq!(
        normalized("NM_001234.1:c.[10del;13del]"),
        "NM_001234.1:c.[11del;26del]"
    );
}

/// The same pair spelled from a different starting position. `c.9del` and
/// `c.10del` are the same variant (both delete one G from the run), so this must
/// reach the same answer — the confluence half of the defect.
#[test]
fn the_same_pair_spelled_from_another_position_agrees() {
    assert_eq!(
        normalized("NM_001234.1:c.[9del;13del]"),
        "NM_001234.1:c.[11del;26del]"
    );
}

/// `general.md:44` names the `r.` axis too, and it went through the same path.
#[test]
fn the_rna_axis_carries_the_same_exception() {
    assert_eq!(
        normalized("NM_001234.1:r.[10del;13del]"),
        "NM_001234.1:r.[11del;26del]"
    );
}

/// The `n.` axis, which `general.md:44` names alongside `c.` and `r.` and which
/// the guard covers through the same `CisKind` arm.
///
/// Spelled in transcript coordinates, so these are the `c.` cases above shifted
/// by the delta of 4: `n.14`/`n.17` are `c.10`/`c.13`, and the junctions sit at
/// `n.15|n.16` and `n.30|n.31`. Pinned rather than assumed equal to the `c.`
/// arm — the two reach the guard through different frames (`axis_frame` gives
/// `delta: 0` here against `4` there), which is exactly the kind of difference
/// that has hidden an axis-specific defect before.
#[test]
fn the_noncoding_axis_carries_the_same_exception() {
    // The lone spellings first, as the reference answer.
    assert_eq!(normalized("NM_001234.1:n.14del"), "NM_001234.1:n.15del");
    assert_eq!(normalized("NM_001234.1:n.17del"), "NM_001234.1:n.30del");

    assert_eq!(
        normalized("NM_001234.1:n.[14del;17del]"),
        "NM_001234.1:n.[15del;30del]",
        "each member must land where it lands alone, rather than being carried \
         past its junction and merged with its sibling (#1450)"
    );
    // And the confluence half: a different spelling of the same variant.
    assert_eq!(
        normalized("NM_001234.1:n.[13del;17del]"),
        "NM_001234.1:n.[15del;30del]"
    );
}
