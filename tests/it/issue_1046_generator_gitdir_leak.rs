//! Issue #1046 — `generate_spec_fixture` must resolve the spec submodule's
//! `commit_sha` from the submodule itself, even when git's hook environment
//! (`GIT_DIR` / `GIT_WORK_TREE`) points at the OUTER repo.
//!
//! Root cause: the generator resolved the SHA with `git -C assets/hgvs-nomenclature
//! rev-parse HEAD`. An ambient `GIT_DIR` overrides the `-C` directory's repo
//! discovery, so under the pre-push `generate_spec_fixture --check` hook — which
//! git runs with `GIT_DIR`/`GIT_WORK_TREE` exported for the outer repo — the
//! embedded `commit_sha` resolved to the outer branch HEAD instead of the
//! submodule pin. The one differing field then failed `--check` on every push.

use std::path::{Path, PathBuf};
use std::process::Command;

fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// `git -C <dir> rev-parse <rev>` with the ambient GIT_* env cleared so the
/// `-C` directory's own repo is the one discovered.
fn rev_parse(dir: &Path, rev: &str) -> Option<String> {
    let out = Command::new("git")
        .arg("-C")
        .arg(dir)
        .args(["rev-parse", rev])
        .env_remove("GIT_DIR")
        .env_remove("GIT_WORK_TREE")
        .env_remove("GIT_INDEX_FILE")
        .output()
        .ok()?;
    if !out.status.success() {
        return None;
    }
    Some(String::from_utf8_lossy(&out.stdout).trim().to_string())
}

/// Path to the built `generate_spec_fixture` example, derived from the running
/// test binary's location (`<target>/<profile>/deps/<test-exe>` →
/// `<target>/<profile>/examples/generate_spec_fixture`). This is robust to
/// `CARGO_TARGET_DIR` and the active profile, and avoids invoking the example
/// through `cargo run` (which would leak our hostile GIT env into cargo too).
fn example_binary() -> PathBuf {
    let mut p = std::env::current_exe().expect("current_exe");
    p.pop(); // drop the test-exe filename → .../deps/
    p.pop(); // drop deps → .../<profile>/
    p.push("examples");
    p.push("generate_spec_fixture");
    p
}

#[test]
fn commit_sha_is_submodule_pin_under_hook_git_env() {
    let root = manifest_dir();
    let submodule = root.join("assets/hgvs-nomenclature");

    // Requires the spec submodule (CI checks it out; local dev inits it). The
    // generator can't run without it, so skip rather than hard-fail.
    if !submodule.join("docs").is_dir() {
        eprintln!("skipping issue_1046: assets/hgvs-nomenclature submodule not initialized");
        return;
    }

    let expected = rev_parse(&submodule, "HEAD").expect("resolve submodule HEAD");
    let outer_head = rev_parse(&root, "HEAD").expect("resolve outer repo HEAD");

    // The regression is only meaningful when the outer HEAD differs from the
    // submodule pin (otherwise the bug and the fix are indistinguishable).
    assert_ne!(
        outer_head, expected,
        "test precondition: outer HEAD must differ from the submodule pin",
    );

    // Discover the outer repo's real git dir with a clean env, then feed it back
    // as GIT_DIR/GIT_WORK_TREE to reproduce exactly what git exports for a hook.
    let outer_git_dir = {
        let out = Command::new("git")
            .arg("-C")
            .arg(&root)
            .args(["rev-parse", "--absolute-git-dir"])
            .env_remove("GIT_DIR")
            .env_remove("GIT_WORK_TREE")
            .output()
            .expect("git rev-parse --absolute-git-dir");
        assert!(out.status.success(), "resolve outer git dir");
        String::from_utf8_lossy(&out.stdout).trim().to_string()
    };

    // Build the example with a clean env (a fast no-op when already built).
    let build = Command::new(env!("CARGO"))
        .current_dir(&root)
        .args([
            "build",
            "--quiet",
            "--features",
            "dev",
            "--example",
            "generate_spec_fixture",
        ])
        .status()
        .expect("build generate_spec_fixture example");
    assert!(build.success(), "example build failed");

    let bin = example_binary();
    assert!(
        bin.exists(),
        "example binary not found at {}",
        bin.display()
    );

    let tmp = std::env::temp_dir().join(format!(
        "ferro_issue_1046_fixture_{}.json",
        std::process::id()
    ));

    // Run the example directly (no cargo) with the hostile GIT env set.
    let status = Command::new(&bin)
        .current_dir(&root)
        .arg("--output")
        .arg(&tmp)
        .env("GIT_DIR", &outer_git_dir)
        .env("GIT_WORK_TREE", &root)
        .status()
        .expect("run generate_spec_fixture example");
    assert!(status.success(), "generator exited with failure");

    let text = std::fs::read_to_string(&tmp).expect("read generated fixture");
    let _ = std::fs::remove_file(&tmp);
    let json: serde_json::Value = serde_json::from_str(&text).expect("parse fixture json");
    let got = json["spec"]["commit_sha"]
        .as_str()
        .expect("spec.commit_sha field present");

    assert_eq!(
        got, expected,
        "commit_sha must be the submodule pin ({expected}), not the outer HEAD ({outer_head})",
    );
}
