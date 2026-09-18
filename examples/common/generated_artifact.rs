//! Shared `--check` semantics for the gitignored generated spec artifacts.
//!
//! `hgvs_spec_normalization.json` and `hgvs_spec_enumeration.json` are build
//! artifacts, not committed files (see `docs/TESTING.md`), so `--check` answers
//! exactly one question: **is my local artifact current?**
//!
//! It deliberately does *not* validate the committed inputs — the plain
//! (non-`--check`) generation run does that, by harvesting the spec checkout and
//! resolving the curated overrides against it. That is why CI and the pre-push
//! hooks run plain generation rather than `--check`: they care about the inputs,
//! which are committed, not about a local artifact's freshness.
//!
//! Keeping those two questions apart is the point. Conflating them is what made
//! `--check` fail on a fresh worktree: it required the artifact to pre-exist in
//! order to answer a question about the inputs, so a never-generated file — the
//! normal state of a gitignored artifact on a new checkout — aborted the run
//! with a bare `No such file or directory`.

use std::path::Path;

/// Compare `rendered` against the artifact at `path`.
///
/// * identical            → `Ok(())`
/// * present and differs  → error naming the regeneration command; `path` is
///   left untouched, so the stale content remains available for inspection
/// * absent               → **not** an error. A gitignored artifact that has
///   never been generated is a cold cache, not drift: there is no committed
///   baseline it could have drifted from. Materialise it and say so.
///
/// `label` names the artifact in messages; `generator` is the generator to
/// rerun. Both generators are `[[bin]]` targets, so the printed command must say
/// `--bin`; `--example` would name a target that no longer exists.
pub fn check_up_to_date(
    path: &Path,
    rendered: &str,
    label: &str,
    generator: &str,
) -> anyhow::Result<()> {
    match std::fs::read_to_string(path) {
        Ok(on_disk) if on_disk == rendered => Ok(()),
        Ok(_) => {
            eprintln!(
                "{label} {} is out of date; rerun: cargo run --features dev --bin {generator}",
                path.display()
            );
            Err(anyhow::anyhow!("{label} out of date"))
        }
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            std::fs::write(path, rendered)
                .map_err(|e| anyhow::anyhow!("write {}: {e}", path.display()))?;
            eprintln!(
                "{label} {} was absent (gitignored build artifact); generated it",
                path.display()
            );
            Ok(())
        }
        Err(e) => Err(anyhow::anyhow!("read {}: {e}", path.display())),
    }
}
