//! CLI surface + axis-logic tests for `ferro project` (#626).
//!
//! Two kinds of test: spawned-binary surface tests (clap validation, exit
//! codes — no reference data needed) and in-process axis/decline tests that
//! call the library `project_axis` with a `MockProvider` fixture (so they need
//! no on-disk manifest).

use std::process::Command;

use ferro_hgvs::cli::project::{project_axis, Axis, AxisOutcome};
use ferro_hgvs::data::{CdotMapper, CdotTranscript, Projector};
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::project::VariantProjector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;

fn ferro() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ferro"))
}

// ===== spawned-binary surface tests =====

#[test]
fn project_rejects_invalid_axis() {
    let out = ferro()
        .args([
            "project",
            "NM_000088.3:c.589G>T",
            "--axis",
            "q",
            "--reference",
            "/tmp",
        ])
        .output()
        .unwrap();
    assert!(!out.status.success(), "invalid --axis must fail");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("axis"), "stderr: {stderr}");
}

#[test]
fn project_requires_axis_and_reference() {
    // Missing --axis (required) -> clap error.
    let out = ferro()
        .args(["project", "NM_000088.3:c.589G>T", "--reference", "/tmp"])
        .output()
        .unwrap();
    assert!(!out.status.success(), "missing --axis must fail");
    // Missing --reference (required) -> clap error.
    let out = ferro()
        .args(["project", "NM_000088.3:c.589G>T", "--axis", "p"])
        .output()
        .unwrap();
    assert!(!out.status.success(), "missing --reference must fail");
}

#[test]
fn project_bad_manifest_exits_nonzero() {
    // A nonexistent manifest dir surfaces a hard error (nonzero exit), not a panic.
    let out = ferro()
        .args([
            "project",
            "NM_000088.3:c.589G>T",
            "--axis",
            "g",
            "--reference",
            "/nonexistent-xyz-ferro",
        ])
        .output()
        .unwrap();
    assert!(!out.status.success());
}

// ===== in-process axis-logic tests (MockProvider fixture) =====

/// A single plus-strand coding transcript NM_TEST.1 on NC_000001.11
/// [1000..1008], CDS = the whole 9-base exon "ATGCGCTAA". The contig is the
/// chromosome accession (NOT "chr1") because a bare `NC_000001.11:g.` input
/// resolves its contig from the accession.
fn fixture() -> VariantProjector<MockProvider> {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("TESTGENE".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            exons: vec![[1000, 1009, 0, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_TEST.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        TxStrand::Plus,
        "ATGCGCTAA".to_string(),
        Some(1),
        Some(9),
        vec![Exon::new(1, 1, 9)],
        Some("NC_000001.11".to_string()),
        Some(1000),
        Some(1008),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    provider.add_genomic_sequence(
        "NC_000001.11",
        format!("{}ATGCGCTAA{}", "N".repeat(1000), "N".repeat(100)),
    );

    VariantProjector::new(projector, provider)
}

#[test]
fn coding_input_n_axis_renders() {
    let vp = fixture();
    let variant = parse_hgvs("NM_TEST.1:c.4C>A").unwrap();
    let outcome = project_axis(&vp, &variant, Axis::Noncoding, None).unwrap();
    assert!(
        matches!(outcome, AxisOutcome::Rendered { .. }),
        "got {outcome:?}"
    );
}

#[test]
fn bare_genomic_with_transcript_renders_coding() {
    let vp = fixture();
    let variant = parse_hgvs("NC_000001.11:g.1003C>A").unwrap();
    let outcome = project_axis(&vp, &variant, Axis::Coding, Some("NM_TEST.1")).unwrap();
    match outcome {
        AxisOutcome::Rendered {
            transcript_id,
            output,
            ..
        } => {
            assert_eq!(transcript_id, "NM_TEST.1");
            assert!(output.contains(":c."), "got {output}");
        }
        other => panic!("expected Rendered, got {other:?}"),
    }
}

/// #1198: an unavailable axis must report the engine's *own* explanation, not a
/// string synthesized from the axis code alone.
///
/// `NM_TEST.1` begins its CDS at transcript base 1, so `c.-1` lies 5' of the
/// transcript's first base and `noncoding_from_coding` refuses it (#1193) with a
/// precise message. That message was computed and then dropped, leaving the user
/// with the generic "no n. representation for this variant".
///
/// A deletion, not a substitution: `c.-1` is not a base of this transcript at
/// all, so a substitution would have to assert a reference base that cannot be
/// checked. The issue's own reproducer (`c.-238del`) has the same shape.
#[test]
fn unavailable_n_axis_reports_the_engines_explanation() {
    let vp = fixture();
    let variant = parse_hgvs("NM_TEST.1:c.-1del").unwrap();
    let outcome = project_axis(&vp, &variant, Axis::Noncoding, None).unwrap();
    match outcome {
        AxisOutcome::Unavailable {
            reason, warnings, ..
        } => {
            // Asserted by substance, not by exact string: the load-bearing claim
            // is that the reason names the refused position and the transcript,
            // which is what the synthesized string could never do. Pinning the
            // whole sentence would couple this test to `FerroError`'s Display
            // prefix and to #1193's prose, neither of which it has an opinion on.
            assert_ne!(
                reason, "no n. representation for this variant",
                "the synthesized string must not survive when the engine explained itself"
            );
            // The `FerroError` class label must not come with it: `c.-1` is a
            // valid `c.` coordinate, so "Invalid coordinates:" beside
            // `status: unavailable` would assert two contradictory things.
            assert!(
                !reason.starts_with("Invalid coordinates:"),
                "the error's class label must be stripped; got {reason:?}"
            );
            for expected in ["n.0", "NM_TEST.1", "5' of the transcript's first base"] {
                assert!(
                    reason.contains(expected),
                    "the reason must carry the engine's explanation ({expected:?} missing); \
                     got {reason:?}"
                );
            }
            // #1182's contract still holds on this arm: an axis can be
            // unavailable *because* of what a warning describes, so the reason
            // must not have displaced the warnings.
            assert_eq!(
                warnings.iter().map(|w| w.code.as_str()).collect::<Vec<_>>(),
                ["POSITION_PAST_END"],
                "the reason and the warnings must coexist"
            );
        }
        other => panic!("expected Unavailable, got {other:?}"),
    }
}

/// The converse of the test above: a derivable axis must not acquire a decline
/// reason. A reason that is always populated would be as uninformative as the
/// generic string it replaces.
#[test]
fn a_derivable_n_axis_records_no_decline_reason() {
    let vp = fixture();
    let variant = parse_hgvs("NM_TEST.1:c.4C>A").unwrap();
    let projection = vp.project_variant(&variant, "NM_TEST.1").unwrap();
    assert!(
        projection.noncoding.is_some(),
        "c.4 is inside the CDS, so the n. axis must derive"
    );
    assert_eq!(projection.axis_decline_reasons.noncoding, None);
}

/// #1198: an allele's `n.` axis is all-or-nothing, so when it is absent at least
/// one member's was — and the member's recorded reason must be carried up rather
/// than the aggregate falling back to the generic string. Dropping it there would
/// reintroduce the same defect one level out.
#[test]
fn an_allele_carries_up_a_members_decline_reason() {
    let vp = fixture();
    // The second member lies 5' of the transcript's first base, so its own `n.`
    // derivation is refused and the whole allele's `n.` axis goes absent.
    let variant = parse_hgvs("NM_TEST.1:c.[4C>A;-1del]").unwrap();
    let outcome = project_axis(&vp, &variant, Axis::Noncoding, None).unwrap();
    match outcome {
        AxisOutcome::Unavailable { reason, .. } => {
            assert!(
                reason.contains("allele member") && reason.contains("n.0"),
                "the member's explanation must reach the aggregate, labelled as a \
                 member's; got {reason:?}"
            );
        }
        other => panic!("expected Unavailable, got {other:?}"),
    }
}
