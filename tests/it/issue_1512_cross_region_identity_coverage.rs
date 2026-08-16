//! #1512 — `drop_identity_members_covered_by_siblings` was blind to a covering
//! sibling in a different region, so a cancelled identity member that a sibling
//! demonstrably claims survived, leaving two members on one position.
//!
//! The pass drops an identity member (`c.X=`) when a `blocks_sibling_shift`
//! sibling's span overlaps it — an identity beside a sibling that claims the same
//! bases is a contradiction, two members on one position, which the apply oracle
//! declines. After #1506 gave `MemberSpan` sequence coordinates so that any two
//! members of one molecule are comparable across a region boundary, the coverage
//! predicate still gated on `s.region == span.region`. That is a **reader
//! declining a cross-region interval**, the #1512 family: the pass only compares
//! two spans (its one mutation is `members.retain(...)`, it never writes an axis
//! coordinate), so it must use the sequence coordinates it now has, not refuse
//! the comparison whenever the two members sit in different regions.
//!
//! Measured on `CORE` under a `1..=15` CDS. `c.15` is the last CDS base and
//! `c.*1` the first 3'UTR base, so a sibling starting at `c.15` and reaching
//! `c.*1` — or one authored at `c.*1` — crosses the boundary and its start
//! region differs from a `c.*1=` identity's:
//!
//! ```text
//! c.[15=;15_*1dup]    ->  c.15_*1dup     both members in the CDS; already dropped
//! c.[*1=;15_*1dup]    ->  c.[15_*1dup;*1=]   PRE-FIX: identity in the 3'UTR survives
//! ```
//!
//! The `15_*1dup` copies `c.15` and `c.*1`, so under `blocks_sibling_shift` it
//! names `c.*1`; the `*1=` names it too. Same accession, overlapping sequence
//! coordinates, different start regions — so the identity was kept and the
//! rendered allele claimed `c.*1` twice.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Under a `1..=15` CDS, `c.15` is `C` and `c.*1` is `T` (positions 15 and 16 of
/// this core).
const CORE: &str = "GCGCTAGTCTCGCCCTGTTA";

fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let length = CORE.len() as u64;
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        CORE.to_string(),
        Some(1),
        Some(15),
        vec![Exon::new(1, 1, length)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

fn normalize(input: &str, direction: ShuffleDirection) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_direction(direction),
    )
    .normalize(&variant)
    .expect("normalization must succeed")
    .to_string()
}

/// The bug. A `c.*1=` identity in the 3'UTR is covered by a `c.15_*1dup` whose
/// span starts one base earlier, in the CDS. The two members claim `c.*1`, so the
/// identity is a cancellation residue to drop — exactly as the in-region control
/// below already does — and the sequence-preserving output is the duplication
/// alone.
#[test]
fn a_cross_region_identity_covered_by_a_duplication_is_dropped() {
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize("NM_TEST.1:c.[*1=;15_*1dup]", direction),
            "NM_TEST.1:c.15_*1dup",
            "the 3'UTR identity is covered by the boundary-crossing dup and must \
             be dropped ({direction:?})"
        );
    }
}

/// The in-region control. Both members sit in the CDS, so this dropped correctly
/// before the fix and must keep doing so — it is what makes the cross-region case
/// above a boundary problem rather than a coverage-logic one.
#[test]
fn an_in_region_identity_covered_by_a_duplication_is_dropped() {
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize("NM_TEST.1:c.[15=;15_*1dup]", direction),
            "NM_TEST.1:c.15_*1dup",
            "the in-region identity coverage must be unaffected ({direction:?})"
        );
    }
}

/// The discriminating control against over-generalising the fix. Removing the
/// region gate must not drop an identity that no sibling *overlaps*: the `c.*3=`
/// identity sits in the 3'UTR and the `c.1_2dup` in the CDS — a cross-region
/// pair, now visible to each other, whose sequence coordinates do not intersect
/// (`*3` is well 3' of the duplicated CDS bases). The overlap test must decline,
/// so the identity survives. A fix that dropped it would be deleting the overlap
/// check, not the region check. The identity's survival is the whole property;
/// the dup's own 3'/5' landing spot is orthogonal shuffle behaviour, so this
/// asserts the member is present rather than pinning the full allele.
#[test]
fn a_cross_region_identity_nothing_overlaps_survives() {
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        let out = normalize("NM_TEST.1:c.[*3=;1_2dup]", direction);
        assert!(
            out.contains("*3="),
            "an uncovered identity must survive however its cross-region sibling \
             shuffles; got {out} ({direction:?})"
        );
    }
}
