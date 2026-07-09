//! Verify (and stamp) the committed `src/arbitrate/spec_citations.json`
//! against the vendored HGVS spec submodule. Run with `--check` in CI.
//!
//!   cargo run --features dev --example generate_spec_citations -- --check
//!
//! Curated, not auto-extracted: a human chooses the excerpt per operation;
//! this tool proves the excerpt still appears verbatim in its spec file and
//! stamps the current spec version. See design spec §5.3.

use std::fs;
use std::path::Path;
use std::process::Command;

fn main() -> anyhow::Result<()> {
    let check = std::env::args().any(|a| a == "--check");
    let spec_dir = Path::new("assets/hgvs-nomenclature");
    let map_path = Path::new("src/arbitrate/spec_citations.json");

    let raw = fs::read_to_string(map_path)
        .map_err(|e| anyhow::anyhow!("read {}: {e}", map_path.display()))?;
    let mut map: serde_json::Value = serde_json::from_str(&raw)?;

    // Stamp the spec version from `git describe --tags` in the submodule.
    // Never downgrade a known version to "unknown" — if git fails for any
    // reason, preserve whatever version is already committed.
    if let Some(version) = detect_spec_version(spec_dir) {
        map["spec_version"] = serde_json::Value::String(version);
    }

    // Verify every excerpt still exists verbatim in its named file.
    let entries = map["entries"].as_object().cloned().unwrap_or_default();
    let mut problems = Vec::new();
    for (op, cites) in &entries {
        for c in cites.as_array().cloned().unwrap_or_default() {
            let file = c["file"].as_str().unwrap_or_default();
            let excerpt = c["excerpt"].as_str().unwrap_or_default();
            let content = fs::read_to_string(spec_dir.join(file)).unwrap_or_default();
            // Normalize away whitespace differences and markdown emphasis
            // markers (`*`/`_`) so a curated plain-text excerpt still matches
            // spec source that renders it with **bold**/_italic_ emphasis.
            let norm = |s: &str| {
                s.chars()
                    .filter(|c| *c != '*' && *c != '_')
                    .collect::<String>()
                    .split_whitespace()
                    .collect::<Vec<_>>()
                    .join(" ")
            };
            if !norm(&content).contains(&norm(excerpt)) {
                problems.push(format!("[{op}] excerpt not found verbatim in {file}"));
            }
        }
    }
    if !problems.is_empty() {
        anyhow::bail!("spec-citation drift:\n  {}", problems.join("\n  "));
    }

    let rendered = serde_json::to_string_pretty(&map)? + "\n";
    if check {
        if rendered != raw {
            anyhow::bail!(
                "src/arbitrate/spec_citations.json is out of date; \
                 re-run without --check to regenerate"
            );
        }
        println!("spec_citations.json up to date");
    } else {
        fs::write(map_path, rendered)?;
        println!("wrote {}", map_path.display());
    }
    Ok(())
}

/// Detect the vendored spec's version via `git describe --tags` inside the
/// submodule checkout. Returns `None` (rather than "unknown") if the git
/// invocation fails for any reason, so callers can preserve the previously
/// committed version instead of downgrading it.
fn detect_spec_version(spec_dir: &Path) -> Option<String> {
    let output = Command::new("git")
        .arg("-C")
        .arg(spec_dir)
        .arg("describe")
        .arg("--tags")
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let version = String::from_utf8(output.stdout).ok()?.trim().to_string();
    if version.is_empty() {
        None
    } else {
        Some(version)
    }
}
