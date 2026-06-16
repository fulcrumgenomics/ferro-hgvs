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
        } => {
            assert_eq!(transcript_id, "NM_TEST.1");
            assert!(output.contains(":c."), "got {output}");
        }
        other => panic!("expected Rendered, got {other:?}"),
    }
}
