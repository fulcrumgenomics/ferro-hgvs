//! CLI surface test for `--strict-reference` on `ferro normalize` (#1001).
//!
//! Builds a minimal stamped reference (a manifest plus one small stamped
//! artifact, `ng_hosted_transcripts.json`), then drifts the artifact in place
//! after stamping so the recomputed content identity no longer matches the
//! recorded one. `MultiFastaProvider::from_manifest_inner` runs the identity
//! check right after manifest validation, before any cdot/genome loading, so
//! this minimal (cdot-less) reference still reaches the verify path — a
//! genomic-only normalize input needs no cdot.
//!
//! Asserts the two load-time behaviors: `--strict-reference` hard-fails with
//! the identity-mismatch error text; the default path only warns (and does
//! not hard-fail with that same error text — the run may still fail
//! downstream for lack of cdot/genome data, which is unrelated to this
//! check).

use std::path::Path;
use std::process::Command;

use ferro_hgvs::prepare::manifest::ReferenceManifest;

fn ferro() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ferro"))
}

/// Write a minimal stamped reference dir (manifest + one small stamped
/// artifact), then mutate the artifact in place so the recorded identity no
/// longer matches the recomputed one. Mirrors
/// `reference::multi_fasta::tests::stamped_reference` (Task 4).
fn drifted_stamped_reference(dir: &Path) {
    std::fs::write(
        dir.join("ng_hosted_transcripts.json"),
        b"{\"schema_version\":1}",
    )
    .unwrap();
    let mut manifest = ReferenceManifest {
        reference_dir: dir.to_path_buf(),
        transcript_count: 1,
        available_prefixes: vec!["NM_".to_string()],
        ng_hosted_transcripts: Some(std::path::PathBuf::from("ng_hosted_transcripts.json")),
        ..Default::default()
    };
    manifest.save().unwrap();

    // Drift the stamped artifact after it was stamped, so the identity
    // recomputed at load time no longer matches the recorded one.
    std::fs::write(
        dir.join("ng_hosted_transcripts.json"),
        b"{\"schema_version\":1,\"drifted\":true}",
    )
    .unwrap();
}

#[test]
fn normalize_strict_reference_fails_on_drift_but_default_warns() {
    let dir = tempfile::tempdir().unwrap();
    drifted_stamped_reference(dir.path());
    let reference = dir.path().to_str().unwrap();

    // `--strict-reference`: the load-time identity mismatch is a hard error.
    let strict_out = ferro()
        .args([
            "normalize",
            "--strict-reference",
            "--reference",
            reference,
            "NC_000001.11:g.100A>G",
        ])
        .output()
        .unwrap();
    assert!(
        !strict_out.status.success(),
        "--strict-reference must exit non-zero on a drifted reference"
    );
    // `main` returns `Result<(), Box<dyn Error>>`, so a propagated error is
    // printed via the process's default `{:?}` (Debug) termination handler,
    // not the `FerroError` `Display` message -- assert on the variant name
    // rather than the Display wording.
    let strict_stderr = String::from_utf8_lossy(&strict_out.stderr);
    assert!(
        strict_stderr.contains("ReferenceIdentityMismatch"),
        "strict stderr: {strict_stderr}"
    );

    // Default (no flag): the mismatch only warns. Don't assert the overall
    // exit code -- this minimal reference has no cdot/genome data, so
    // normalization itself may still fail downstream for unrelated reasons.
    // The distinguishing signal is the "warning:" wording rather than the
    // hard "does not match its recorded identity" error used above.
    let default_out = ferro()
        .args([
            "normalize",
            "--reference",
            reference,
            "NC_000001.11:g.100A>G",
        ])
        .output()
        .unwrap();
    let default_stderr = String::from_utf8_lossy(&default_out.stderr);
    assert!(
        default_stderr.contains("warning: reference content does not match its recorded identity"),
        "default stderr: {default_stderr}"
    );
}
