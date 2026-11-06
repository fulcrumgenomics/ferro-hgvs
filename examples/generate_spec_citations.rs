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

use ferro_hgvs::conformance::completeness::CaptureLedger;

fn main() -> anyhow::Result<()> {
    let check = std::env::args().any(|a| a == "--check");
    if let Some(message) = ferro_hgvs::normalize::partition_switch_startup_error() {
        anyhow::bail!("{message}");
    }
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

    // Verify every excerpt still exists verbatim in its named file. Each
    // citation record is routed through the ledger, so a malformed or drifted
    // record is counted as a drop and `finish` refuses before anything is
    // written or approved.
    let entries = map["entries"]
        .as_object()
        .ok_or_else(|| anyhow::anyhow!("`entries` must be an object"))?;
    let mut ledger = CaptureLedger::new("spec citations");
    for (op, cites) in entries {
        let Some(cites) = cites.as_array() else {
            ledger.record_drop(op.as_str(), "citations must be an array");
            continue;
        };
        for (i, c) in cites.iter().enumerate() {
            ledger.record(format!("{op}[{i}]"), validate_citation(spec_dir, c));
        }
    }
    ledger
        .finish()
        .map_err(|e| anyhow::anyhow!("spec-citation drift: {e}"))?;

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

/// Check one citation record: `file` and `excerpt` must be non-empty strings,
/// and the excerpt must still appear in that spec file. Whitespace and
/// markdown emphasis markers (`*`/`_`) are normalized away on both sides, so a
/// curated plain-text excerpt still matches source that renders it with
/// **bold**/_italic_ emphasis.
fn validate_citation(spec_dir: &Path, c: &serde_json::Value) -> Result<(), String> {
    let field = |name: &str| -> Result<&str, String> {
        match c[name].as_str() {
            Some(v) if !v.trim().is_empty() => Ok(v),
            _ => Err(format!("`{name}` must be a non-empty string")),
        }
    };
    let file = field("file")?;
    let excerpt = field("excerpt")?;
    let content =
        fs::read_to_string(spec_dir.join(file)).map_err(|e| format!("read {file}: {e}"))?;
    let norm = |s: &str| {
        s.chars()
            .filter(|c| *c != '*' && *c != '_')
            .collect::<String>()
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ")
    };
    if norm(&content).contains(&norm(excerpt)) {
        Ok(())
    } else {
        Err(format!("excerpt not found verbatim in {file}"))
    }
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

#[cfg(test)]
mod tests {
    use super::validate_citation;
    use std::path::PathBuf;

    fn spec_dir_with(file: &str, body: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "ferro-spec-citations-test-{}-{file}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join(file), body).unwrap();
        dir
    }

    #[test]
    fn a_present_excerpt_validates_through_emphasis() {
        let dir = spec_dir_with("a.md", "it **must** be described as a duplication");
        let c = serde_json::json!({"file": "a.md", "excerpt": "it must be described"});
        assert_eq!(validate_citation(&dir, &c), Ok(()));
    }

    #[test]
    fn an_empty_or_missing_excerpt_is_refused() {
        // `contains("")` is always true, so an empty excerpt would otherwise
        // approve a citation that quotes nothing.
        let dir = spec_dir_with("b.md", "anything");
        for c in [
            serde_json::json!({"file": "b.md", "excerpt": ""}),
            serde_json::json!({"file": "b.md"}),
            serde_json::json!({"file": "b.md", "excerpt": 7}),
        ] {
            assert!(validate_citation(&dir, &c).is_err(), "{c}");
        }
    }

    #[test]
    fn a_missing_file_field_is_refused() {
        let dir = spec_dir_with("c.md", "text");
        let c = serde_json::json!({"excerpt": "text"});
        assert!(validate_citation(&dir, &c).is_err());
    }
}
