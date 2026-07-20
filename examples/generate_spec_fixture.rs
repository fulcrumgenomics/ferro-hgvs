//! Generator for the HGVS v21.0 spec normalization fixture.
//!
//! Walks the vendored `assets/hgvs-nomenclature/` markdown + syntax.yaml,
//! extracts every HGVS string, runs each through ferro's `Normalizer`, merges
//! hand overrides, and writes `tests/fixtures/grammar/hgvs_spec_normalization.json`.
//!
//! Two modes:
//!   - default:  regenerate the fixture in-place
//!   - --check:  regenerate in-memory and exit non-zero if the on-disk fixture
//!     differs (used in CI to catch drift)
//!
//! Run: `cargo run --features dev --example generate_spec_fixture`
//! Check: `cargo run --features dev --example generate_spec_fixture -- --check`
//!
//! See: docs/superpowers/specs/2026-05-03-issue-84-hgvs-v21-normalization-fixture-design.md

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use serde::{Deserialize, Deserializer, Serialize};

#[derive(Parser, Debug)]
#[command(about = "Generate the HGVS v21.0 spec normalization fixture")]
struct Cli {
    /// Path to the vendored spec checkout (must be at tag 21.0.0).
    #[arg(long, default_value = "assets/hgvs-nomenclature")]
    spec_dir: PathBuf,

    /// Path to the hand-curated overrides file.
    #[arg(
        long,
        default_value = "tests/fixtures/grammar/hgvs_spec_normalization_overrides.json"
    )]
    overrides: PathBuf,

    /// Path to the output fixture file.
    #[arg(
        long,
        default_value = "tests/fixtures/grammar/hgvs_spec_normalization.json"
    )]
    output: PathBuf,

    /// Check that the on-disk fixture is byte-identical to a fresh regeneration.
    /// Exits 1 on mismatch. Does not modify the file.
    #[arg(long)]
    check: bool,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: &Cli) -> anyhow::Result<()> {
    let candidates = sources::discover(&cli.spec_dir)?;
    let overrides = overrides::load(&cli.overrides)?;
    let rows = runner::build_rows(&candidates, &overrides)?;
    let rendered = render::render(&rows, &cli.spec_dir)?;

    if cli.check {
        let on_disk = std::fs::read_to_string(&cli.output)
            .map_err(|e| anyhow::anyhow!("read {}: {e}", cli.output.display()))?;
        if on_disk != rendered {
            eprintln!(
                "fixture {} is out of date; rerun: cargo run --features dev --example generate_spec_fixture",
                cli.output.display()
            );
            return Err(anyhow::anyhow!("fixture out of date"));
        }
    } else {
        std::fs::write(&cli.output, rendered)?;
    }
    Ok(())
}

#[path = "common/spec_harvest.rs"]
mod spec_harvest;

use spec_harvest::{classify, prefix, sources};

// ---------- overrides ----------

mod overrides {
    use super::*;

    #[derive(Debug, Default, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct Overrides {
        #[serde(default)]
        pub by_input: BTreeMap<String, OverrideEntry>,
    }

    /// `spec_expected` uses the doubly-optional pattern so the override file
    /// can distinguish three states:
    ///   * field absent     → `None`             (use the auto-default)
    ///   * field set to null → `Some(None)`       (force "spec rejects this")
    ///   * field set to "X" → `Some(Some("X"))`  (force a canonical output)
    #[derive(Debug, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct OverrideEntry {
        #[serde(default)]
        pub status: Option<String>,
        #[serde(default, deserialize_with = "deserialize_some")]
        pub spec_expected: Option<Option<String>>,
        /// Free-form note attached to the row in the emitted fixture.
        /// Used to capture nuances where the auto-classified `status`
        /// doesn't fully convey ferro's behavior — e.g. rows whose
        /// `current` is the spec's recommended canonical form even
        /// though the original `input` is `<code class="invalid">`.
        /// Closes-after: #353.
        #[serde(default)]
        pub note: Option<String>,
        /// Force a specific prefixed form for bare fragments. Wins over the
        /// auto-default in `prefix::default_prefixed`.
        #[serde(default)]
        pub input_prefixed: Option<String>,
        /// Mark this row as requiring a real reference sequence to evaluate
        /// (e.g. §2.1 3'-rule shifting examples). The test consumer skips
        /// `Some(true)` rows until #82's `from_manifest` loader lands. There
        /// is no auto-detection — auditors flip this field by hand during
        /// the #83 sweep when they see a row whose correctness depends on
        /// reference data.
        #[serde(default)]
        pub requires_reference: Option<bool>,
    }

    fn deserialize_some<'de, T, D>(d: D) -> Result<Option<T>, D::Error>
    where
        T: Deserialize<'de>,
        D: Deserializer<'de>,
    {
        T::deserialize(d).map(Some)
    }

    pub fn load(path: &Path) -> anyhow::Result<Overrides> {
        if !path.exists() {
            return Ok(Overrides::default());
        }
        let text = std::fs::read_to_string(path)?;
        let parsed: Overrides = serde_json::from_str(&text)
            .map_err(|e| anyhow::anyhow!("parse {}: {e}", path.display()))?;
        Ok(parsed)
    }
}

// ---------- runner ----------

mod runner {
    use super::*;
    use ferro_hgvs::reference::mock::MockProvider;
    use ferro_hgvs::{parse_hgvs, Normalizer};

    #[derive(Debug, Serialize)]
    pub struct Row {
        pub input: String,
        /// Default-prefixed form for bare illustrative fragments
        /// (e.g. `c.1083A>C` → `NM_004006.2:c.1083A>C`). When `Some`, this is
        /// the string ferro was actually run against, and the test consumer
        /// asserts against this form. Verbatim spec text stays in `input`.
        #[serde(skip_serializing_if = "Option::is_none")]
        pub input_prefixed: Option<String>,
        pub current: String,
        /// `None` (serialized as JSON `null`) means the spec rejects this input
        /// — either via an explicit `<code class="invalid">…</code>` marker in
        /// the spec text or via an override entry with `spec_expected: null`.
        pub spec_expected: Option<String>,
        pub status: String,
        pub coordinate_system: String,
        pub source_kind: String,
        pub source_paths: Vec<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        pub working_group: Option<String>,
        /// Free-form clarifying note. Surfaced from the override file
        /// for rows where the `status` bucket doesn't fully capture
        /// ferro's behavior (e.g. canonical form emitted from a
        /// spec-rejected input). Closes-after: #353.
        #[serde(skip_serializing_if = "Option::is_none")]
        pub note: Option<String>,
        /// Mirrors the override flag — `Some(true)` means the row needs a
        /// real reference sequence to evaluate and is skipped by the test
        /// consumer. Surfaced on the row so consumers don't have to read
        /// the override file.
        #[serde(skip_serializing_if = "Option::is_none")]
        pub requires_reference: Option<bool>,
        /// Warning codes ferro currently emits for this row, sorted alphabetically.
        /// Empty for rows with no warnings.
        #[serde(skip_serializing_if = "Vec::is_empty")]
        pub expected_warnings: Vec<String>,
    }

    fn source_kind_priority(k: sources::SourceKind) -> u8 {
        use sources::SourceKind::*;
        match k {
            Recommendation => 5,
            SyntaxYaml => 4,
            Background => 3,
            Consultation => 2,
            PortedLegacyProbe => 1,
        }
    }

    fn source_kind_str(k: sources::SourceKind) -> &'static str {
        use sources::SourceKind::*;
        match k {
            Recommendation => "recommendation",
            SyntaxYaml => "syntax_yaml",
            Background => "background",
            Consultation => "consultation",
            PortedLegacyProbe => "ported_legacy_probe",
        }
    }

    fn coord_system(input: &str) -> String {
        // Find the first occurrence of `:<sys>.` and return the single-letter
        // system. Strings without a colon may be bare-prefix forms (e.g.
        // `c.123del`) used in some spec examples — handle that too.
        for prefix in ["c.", "g.", "m.", "n.", "o.", "p.", "r."] {
            if input.contains(&format!(":{prefix}")) || input.starts_with(prefix) {
                return prefix[..1].to_string();
            }
        }
        "?".to_string()
    }

    /// Chromosome landmarks the spec uses in place of a numeric position.
    /// They carry no edit meaning, so a coordinate built only from them (e.g.
    /// `g.pter`) is still edit-less.
    const POSITION_LANDMARKS: [&str; 3] = ["pter", "qter", "cen"];

    /// True when `input` carries no edit operator: a coordinate-system prefix
    /// (`c.`/`g.`/`m.`/`n.`/`o.`/`p.`/`r.`) followed only by position tokens —
    /// position characters (digits, `+ - * ?`) or a chromosome landmark
    /// (`pter`/`qter`/`cen`) — joined by `_`.
    ///
    /// The spec formats position references in prose (e.g. "the A-stretch
    /// running from position `c.5690` to `c.5697`") with the same inline-code
    /// markup as variant examples, so the harvester scrapes them; a bare
    /// coordinate is not a describable variant. Edit-less strings are not all
    /// prose references, though: the spec's whole-variant markers (`p.?`, `p.0`,
    /// `r.?`, …) and its edit-less negatives (`c.123-65_-50`, an incomplete
    /// range the spec requires be *rejected*) share this shape. So the caller
    /// pairs this test with evidence that the string is only a position — see
    /// [`is_position_reference`].
    fn is_edit_less_coordinate(input: &str) -> bool {
        // Strip a leading accession: the coordinate follows the *first* colon.
        // Splitting on the last colon instead would look past the junctions of a
        // chimeric description (`NC_1:g.pter_100::NC_2:g.200_cen_qter`) and judge
        // only its tail, wrongly calling the whole string edit-less.
        let core = input.split_once(':').map_or(input, |(_, rest)| rest);
        let rest = match core.get(0..2) {
            Some("c." | "g." | "m." | "n." | "o." | "p." | "r.") => &core[2..],
            _ => return false,
        };
        // A bare prefix (`c.`) is edit-less; but past that, an empty token from
        // `_`-splitting means a malformed range (`c.1_`, `c.1__2`) — not a clean
        // position, so it is not treated as an edit-less coordinate.
        rest.is_empty()
            || rest.split('_').all(|token| {
                !token.is_empty()
                    && (POSITION_LANDMARKS.contains(&token)
                        || token
                            .chars()
                            .all(|ch| matches!(ch, '0'..='9' | '+' | '-' | '*' | '?')))
            })
    }

    /// True when a harvested string is a position reference lifted out of spec
    /// prose rather than a variant, and so should not become a fixture row.
    ///
    /// Two signals, both requiring [`is_edit_less_coordinate`] (no edit
    /// operator) to hold:
    ///
    /// * ferro cannot parse it at all — a position cited in a sentence, e.g.
    ///   `c.` or a truncated fragment; or
    /// * ferro normalizes it by merely appending `=` (e.g. `g.12345678` →
    ///   `NC_000023.11:g.12345678=`), i.e. it invented the "identity" edit the
    ///   string never had. Scoring that as a divergence is noise.
    ///
    /// Edit-less strings that ferro parses *and* leaves alone are kept: the
    /// whole-variant markers (`p.?`, `p.0`) and the spec's edit-less negatives
    /// such as `c.123-65_-50`, which the checklist requires be rejected.
    ///
    /// `parse_ok`/`current` come from `target` (the accession-prefixed form)
    /// while the shape test reads the raw `input`. These are equivalent here:
    /// prefixing only prepends an accession, and `is_edit_less_coordinate`
    /// strips any accession before testing, so both forms share the shape.
    fn is_position_reference(input: &str, target: &str, current: &str, parse_ok: bool) -> bool {
        is_edit_less_coordinate(input) && (!parse_ok || current == format!("{target}="))
    }

    pub fn build_rows(
        candidates: &[sources::Candidate],
        overrides: &overrides::Overrides,
    ) -> anyhow::Result<Vec<Row>> {
        #[derive(Default)]
        struct Agg {
            kinds: Vec<sources::SourceKind>,
            paths: Vec<String>,
            wg: Option<String>,
            spec_rejects: bool,
        }
        let mut by_input: BTreeMap<String, Agg> = BTreeMap::new();
        for c in candidates {
            let a = by_input.entry(c.input.clone()).or_default();
            a.kinds.push(c.source_kind);
            a.paths.push(c.source_path.clone());
            if c.working_group.is_some() {
                a.wg = c.working_group.clone();
            }
            if c.intent == sources::SpanIntent::SpecRejects {
                a.spec_rejects = true;
            }
        }

        let normalizer = Normalizer::new(MockProvider::new());

        let mut rows = Vec::with_capacity(by_input.len());
        for (input, agg) in by_input {
            let kind = agg
                .kinds
                .iter()
                .max_by_key(|k| source_kind_priority(**k))
                .copied()
                .unwrap();
            let mut paths: Vec<String> = agg.paths.clone();
            paths.sort();
            paths.dedup();

            let ov = overrides.by_input.get(&input);

            // input_prefixed resolution. The spec routinely writes bare
            // illustrative fragments like `c.1083A>C` whose accession is
            // implied by the surrounding paragraph. Pin a default accession
            // per coord system so ferro can actually parse them. Override
            // wins for rows that need a non-default prefix.
            let input_prefixed: Option<String> = ov
                .and_then(|o| o.input_prefixed.clone())
                .or_else(|| prefix::default_prefixed(&input));

            // What ferro is run against: the prefixed form when present,
            // else the raw input. This is also the anchor for the
            // spec_expected default (the spec wrote the bare form as
            // shorthand for the prefixed form).
            let target = input_prefixed.as_deref().unwrap_or(&input);
            let (current, parse_ok, normalize_ok, expected_warnings) =
                run_ferro(&normalizer, target);

            // Drop position references harvested from spec prose (e.g.
            // "position `c.5690`"): they are not describable variants, just
            // positions cited in a sentence. An explicit override keeps the row
            // regardless (the auditor deliberately classified it).
            if ov.is_none() && is_position_reference(&input, target, &current, parse_ok) {
                continue;
            }

            // spec_expected resolution. Source order: override > structural
            // extraction > default-to-target.
            //   * override.spec_expected = Some(value) → use value verbatim
            //     (Some(None) forces "spec rejects"; Some(Some(s)) forces a
            //     canonical output)
            //   * spec marks the input via `<code class="invalid">…</code>`
            //     → None (spec rejects)
            //   * otherwise → Some(target) (spec offered the target as
            //     a canonical form; ferro is expected to round-trip it)
            let spec_expected: Option<String> = match ov.and_then(|o| o.spec_expected.clone()) {
                Some(value) => value,
                None => {
                    if agg.spec_rejects {
                        None
                    } else {
                        Some(target.to_string())
                    }
                }
            };

            // Status is derived from the taxonomy table in `classify::classify`.
            // Overrides may pin a non-default status (e.g. for rows the
            // auditor has hand-classified into a sub-bucket).
            let auto_status =
                classify::classify(parse_ok, normalize_ok, &current, spec_expected.as_deref());
            let status = ov
                .and_then(|o| o.status.clone())
                .unwrap_or_else(|| auto_status.to_string());

            rows.push(Row {
                coordinate_system: coord_system(&input),
                input: input.clone(),
                input_prefixed,
                current,
                spec_expected,
                status,
                source_kind: source_kind_str(kind).to_string(),
                source_paths: paths,
                working_group: agg.wg,
                note: ov.and_then(|o| o.note.clone()),
                requires_reference: ov.and_then(|o| o.requires_reference),
                expected_warnings,
            });
        }

        let row_inputs: std::collections::BTreeSet<&str> =
            rows.iter().map(|r| r.input.as_str()).collect();
        let unknown: Vec<&str> = overrides
            .by_input
            .keys()
            .filter(|k| !row_inputs.contains(k.as_str()))
            .map(String::as_str)
            .collect();
        if !unknown.is_empty() {
            anyhow::bail!(
                "overrides reference unknown inputs (typo or spec drift): {:?}",
                unknown
            );
        }

        rows.sort_by(|a, b| {
            a.coordinate_system
                .cmp(&b.coordinate_system)
                .then(a.status.cmp(&b.status))
                .then(a.input.cmp(&b.input))
        });
        Ok(rows)
    }

    fn run_ferro(
        normalizer: &Normalizer<MockProvider>,
        input: &str,
    ) -> (String, bool, bool, Vec<String>) {
        match parse_hgvs(input) {
            Err(e) => (format!("parse error: {e}"), false, false, Vec::new()),
            Ok(v) => match normalizer.normalize_with_diagnostics(&v) {
                Err(e) => (format!("normalize error: {e}"), true, false, Vec::new()),
                Ok(n) => {
                    let mut codes: Vec<String> =
                        n.warnings.iter().map(|w| w.code().to_string()).collect();
                    codes.sort();
                    codes.dedup();
                    (format!("{}", n.result), true, true, codes)
                }
            },
        }
    }

    #[cfg(test)]
    mod tests {
        use super::{is_edit_less_coordinate, is_position_reference};

        #[test]
        fn flags_bare_position_references() {
            for bare in [
                "c.5690",
                "c.147",
                "c.1210-33",
                "NM_004006.1:c.5690",
                "AB053210.2:r.1289-365_1289-73",
                "c.", // bare prefix
                "c.*",
                "g.12345678",
                "g.pter", // chromosome landmark, still no edit operator
                "p.?",    // edit-less marker — shape matches; caller's gate keeps it
            ] {
                assert!(
                    is_edit_less_coordinate(bare),
                    "{bare:?} should be an edit-less coordinate"
                );
            }
        }

        #[test]
        fn does_not_flag_real_variants() {
            for variant in [
                "c.5697del",
                "c.76A>T",
                "NM_004006.2:c.1210-33del",
                "g.123_456dup",
                "LRG_199t1:c.11T>G",
                "p.Met1ext-5", // has an edit token (`ext`), not bare
                "g.1234=",
                // Malformed ranges: an empty `_`-token is not a clean position,
                // so these are not classified as edit-less (would otherwise be
                // wrongly dropped as prose if ferro rejects them).
                "c.1_",
                "c.1__2",
                // Chimeric junctions: edit-less-looking tail, but the whole
                // string is a (spec-rejected) description, not a position.
                "NC_000002.12:g.pter_8247756::NC_000011.10:g.15825273_cen_qter",
                "g.[chr11:pter_15825272::chr2:8247757_cen_qter]",
            ] {
                assert!(
                    !is_edit_less_coordinate(variant),
                    "{variant:?} is a real variant, not a bare coordinate"
                );
            }
        }

        #[test]
        fn drops_coordinates_ferro_merely_stamps_with_equals() {
            assert!(is_position_reference(
                "g.12345678",
                "NC_000023.11:g.12345678",
                "NC_000023.11:g.12345678=",
                true,
            ));
            assert!(is_position_reference(
                "g.pter",
                "NC_000023.11:g.pter",
                "NC_000023.11:g.pter=",
                true,
            ));
        }

        #[test]
        fn drops_unparsable_prose_fragments() {
            assert!(is_position_reference(
                "c.",
                "NM_004006.2:c.",
                "parse error",
                false
            ));
        }

        #[test]
        fn keeps_edit_less_spec_negatives_and_markers() {
            // `c.123-65_-50` is an incomplete range the spec requires be
            // rejected: edit-less, but ferro leaves it untouched, so it stays.
            assert!(!is_position_reference(
                "c.123-65_-50",
                "NM_004006.2:c.123-65_-50",
                "NM_004006.2:c.123-65_-50",
                true,
            ));
            // Whole-variant markers parse and round-trip unchanged.
            assert!(!is_position_reference(
                "p.?",
                "NP_000079.2:p.?",
                "NP_000079.2:p.?",
                true
            ));
            // Real variants never match, whatever ferro does with them.
            assert!(!is_position_reference(
                "g.123del",
                "NC_000023.11:g.123del",
                "NC_000023.11:g.123del",
                true,
            ));
        }
    }
}

// ---------- render ----------

mod render {
    use super::*;

    #[derive(Serialize)]
    struct Document<'a> {
        description: &'a str,
        spec: SpecBlock<'a>,
        generated_utc: String,
        summary: Summary,
        rows: &'a [runner::Row],
    }

    #[derive(Serialize)]
    struct SpecBlock<'a> {
        source: &'a str,
        tag: &'a str,
        commit_sha: String,
    }

    #[derive(Serialize, Default)]
    struct Summary {
        total: usize,
        by_status: BTreeMap<String, usize>,
        by_coordinate_system: BTreeMap<String, usize>,
        by_source_kind: BTreeMap<String, usize>,
    }

    fn read_submodule_commit(spec_dir: &Path) -> anyhow::Result<String> {
        // Clear the ambient GIT_* env before resolving the submodule pin. Git
        // exports GIT_DIR/GIT_WORK_TREE (pointing at the OUTER repo) when it
        // runs a hook, and those override the `-C spec_dir` repo discovery — so
        // under the pre-push `--check` hook this would otherwise resolve the
        // outer branch HEAD instead of the submodule's commit (#1046).
        let out = std::process::Command::new("git")
            .arg("-C")
            .arg(spec_dir)
            .args(["rev-parse", "HEAD"])
            .env_remove("GIT_DIR")
            .env_remove("GIT_WORK_TREE")
            .env_remove("GIT_INDEX_FILE")
            .output()
            .map_err(|e| anyhow::anyhow!("git rev-parse: {e}"))?;
        if !out.status.success() {
            anyhow::bail!(
                "git rev-parse failed: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        }
        Ok(String::from_utf8_lossy(&out.stdout).trim().to_string())
    }

    pub fn render(rows: &[runner::Row], spec_dir: &Path) -> anyhow::Result<String> {
        let commit_sha = read_submodule_commit(spec_dir)?;
        let mut summary = Summary {
            total: rows.len(),
            ..Default::default()
        };
        for r in rows {
            *summary.by_status.entry(r.status.clone()).or_default() += 1;
            *summary
                .by_coordinate_system
                .entry(r.coordinate_system.clone())
                .or_default() += 1;
            *summary
                .by_source_kind
                .entry(r.source_kind.clone())
                .or_default() += 1;
        }
        let doc = Document {
            description:
                "HGVS v21.0 spec normalization fixture - pins ferro's current normalize() output. \
                          Companion to issue #84.",
            spec: SpecBlock {
                source: "https://github.com/HGVSnomenclature/hgvs-nomenclature",
                tag: "21.0.0",
                commit_sha,
            },
            generated_utc: "fixture-byte-stable".to_string(),
            summary,
            rows,
        };
        let mut out = serde_json::to_string_pretty(&doc)?;
        out.push('\n');
        Ok(out)
    }
}
