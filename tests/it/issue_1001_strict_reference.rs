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

/// Write a minimal reference dir with `manifest.json` present but the
/// `reference_identity` field absent — the "legacy"/unstamped case.
///
/// `ReferenceManifest::save` always stamps, so building this fixture goes
/// around `save()`: stamp normally, then strip `reference_identity` back out
/// of the on-disk JSON.
fn unstamped_reference(dir: &Path) {
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

    let manifest_path = dir.join("manifest.json");
    let mut value: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&manifest_path).unwrap()).unwrap();
    value.as_object_mut().unwrap().remove("reference_identity");
    std::fs::write(&manifest_path, serde_json::to_vec_pretty(&value).unwrap()).unwrap();
}

/// Write a minimal *stamped* (non-drifted) reference dir — same shape as
/// [`drifted_stamped_reference`] but without the post-stamp mutation.
fn stamped_reference(dir: &Path) {
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
}

#[test]
fn check_reports_unstamped_then_write_identity_stamps_it() {
    let dir = tempfile::tempdir().unwrap();
    unstamped_reference(dir.path());
    let reference = dir.path().to_str().unwrap();

    // Plain `check` on an unstamped reference: reports "unstamped", exits 0.
    let out = ferro()
        .args(["check", "--reference", reference])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "check on an unstamped reference must exit 0; stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("Reference identity: unstamped"),
        "stdout: {stdout}"
    );

    // `--write-identity` stamps it.
    let write_out = ferro()
        .args(["check", "--reference", reference, "--write-identity"])
        .output()
        .unwrap();
    assert!(
        write_out.status.success(),
        "write-identity on a fresh unstamped reference must exit 0; stderr: {}",
        String::from_utf8_lossy(&write_out.stderr)
    );
    let write_stdout = String::from_utf8_lossy(&write_out.stdout);
    assert!(
        write_stdout.contains("Wrote reference identity:"),
        "stdout: {write_stdout}"
    );

    let manifest_value: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.path().join("manifest.json")).unwrap()).unwrap();
    assert!(
        manifest_value
            .get("reference_identity")
            .and_then(|v| v.as_str())
            .is_some(),
        "manifest should now have a reference_identity: {manifest_value}"
    );

    // A second plain check now reports "verified".
    let verify_out = ferro()
        .args(["check", "--reference", reference])
        .output()
        .unwrap();
    assert!(
        verify_out.status.success(),
        "check after stamping must exit 0; stderr: {}",
        String::from_utf8_lossy(&verify_out.stderr)
    );
    let verify_stdout = String::from_utf8_lossy(&verify_out.stdout);
    assert!(
        verify_stdout.contains("Reference identity: verified"),
        "stdout: {verify_stdout}"
    );
}

#[test]
fn write_identity_refuses_to_overwrite_a_differing_stamp_without_force() {
    let dir = tempfile::tempdir().unwrap();
    stamped_reference(dir.path());
    let reference = dir.path().to_str().unwrap();

    // Drift a stamped artifact after stamping.
    std::fs::write(
        dir.path().join("ng_hosted_transcripts.json"),
        b"{\"schema_version\":1,\"drifted\":true}",
    )
    .unwrap();

    // Without --force: refuses, exits non-zero.
    let refuse_out = ferro()
        .args(["check", "--reference", reference, "--write-identity"])
        .output()
        .unwrap();
    assert!(
        !refuse_out.status.success(),
        "write-identity over a differing stamp without --force must exit non-zero"
    );

    // With --force: re-stamps, exits 0.
    let force_out = ferro()
        .args([
            "check",
            "--reference",
            reference,
            "--write-identity",
            "--force",
        ])
        .output()
        .unwrap();
    assert!(
        force_out.status.success(),
        "write-identity --force must exit 0; stderr: {}",
        String::from_utf8_lossy(&force_out.stderr)
    );
    let force_stdout = String::from_utf8_lossy(&force_out.stdout);
    assert!(
        force_stdout.contains("Wrote reference identity:"),
        "stdout: {force_stdout}"
    );

    // A subsequent plain check reports "verified" against the drifted content.
    let verify_out = ferro()
        .args(["check", "--reference", reference])
        .output()
        .unwrap();
    assert!(verify_out.status.success());
    let verify_stdout = String::from_utf8_lossy(&verify_out.stdout);
    assert!(
        verify_stdout.contains("Reference identity: verified"),
        "stdout: {verify_stdout}"
    );
}
