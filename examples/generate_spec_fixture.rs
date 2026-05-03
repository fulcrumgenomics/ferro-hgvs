//! Generator for the HGVS v21.0 spec normalization fixture.
//!
//! Walks the vendored `assets/hgvs-nomenclature/` markdown + syntax.yaml,
//! extracts every HGVS string, runs each through ferro's `Normalizer`, merges
//! hand overrides, and writes `tests/fixtures/grammar/hgvs_spec_normalization.json`.
//!
//! Two modes:
//!   - default:  regenerate the fixture in-place
//!   - --check:  regenerate in-memory and exit non-zero if the on-disk fixture
//!               differs (used in CI to catch drift)
//!
//! Run: `cargo run --features dev --example generate_spec_fixture`
//! Check: `cargo run --features dev --example generate_spec_fixture -- --check`
//!
//! See: docs/superpowers/specs/2026-05-03-issue-84-hgvs-v21-normalization-fixture-design.md

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use serde::{Deserialize, Serialize};

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

// ---------- sources ----------

mod sources {
    use super::*;
    use std::ffi::OsStr;

    #[derive(Debug)]
    pub struct Candidate {
        pub input: String,
        pub source_kind: SourceKind,
        pub source_path: String,
        pub working_group: Option<String>,
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum SourceKind {
        Recommendation,
        Background,
        Consultation,
        SyntaxYaml,
        PortedLegacyProbe,
    }

    pub fn discover(spec_dir: &Path) -> anyhow::Result<Vec<Candidate>> {
        let mut all = Vec::new();
        for (path, kind, wg) in whitelisted_files(spec_dir) {
            let text = std::fs::read_to_string(&path)?;
            let rel = path
                .strip_prefix(spec_dir)
                .unwrap_or(&path)
                .to_string_lossy()
                .replace('\\', "/");
            for cand_text in extract_codespans(&text) {
                if let Some(input) = canonicalize(&cand_text) {
                    all.push(Candidate {
                        input,
                        source_kind: kind,
                        source_path: rel.clone(),
                        working_group: wg.clone(),
                    });
                }
            }
        }
        all.extend(extract_from_syntax_yaml(spec_dir)?);
        all.extend(ported_legacy_probes());
        Ok(all)
    }

    fn whitelisted_files(spec_dir: &Path) -> Vec<(PathBuf, SourceKind, Option<String>)> {
        let mut out = Vec::new();
        let docs = spec_dir.join("docs");

        let recs = docs.join("recommendations");
        for entry in walkdir_md(&recs) {
            out.push((entry, SourceKind::Recommendation, None));
        }

        for name in ["numbering.md", "refseq.md", "simple.md"] {
            let p = docs.join("background").join(name);
            if p.exists() {
                out.push((p, SourceKind::Background, None));
            }
        }

        let cons = docs.join("consultation");
        for n in 1..=10 {
            if n == 8 {
                continue; // SVD-WG008 is meta about reference sequences; no new HGVS forms.
            }
            let name = format!("SVD-WG{:03}.md", n);
            let p = cons.join(&name);
            if p.exists() {
                let wg = format!("SVD-WG{:03}", n);
                out.push((p, SourceKind::Consultation, Some(wg)));
            }
        }
        out
    }

    fn walkdir_md(dir: &Path) -> Vec<PathBuf> {
        let mut out = Vec::new();
        if !dir.exists() {
            return out;
        }
        let mut stack = vec![dir.to_path_buf()];
        while let Some(d) = stack.pop() {
            let entries = match std::fs::read_dir(&d) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("warning: cannot read dir {}: {e}", d.display());
                    continue;
                }
            };
            for e in entries.flatten() {
                let p = e.path();
                if p.is_dir() {
                    stack.push(p);
                } else if p.extension() == Some(OsStr::new("md")) {
                    out.push(p);
                }
            }
        }
        out.sort();
        out
    }

    /// Return every code span (inline `…` and fenced ``` …``` block) text in the markdown.
    fn extract_codespans(md: &str) -> Vec<String> {
        use pulldown_cmark::{Event, Parser, Tag, TagEnd};
        let mut out = Vec::new();
        let mut in_code_block = false;
        let mut buf = String::new();
        for ev in Parser::new(md) {
            match ev {
                Event::Start(Tag::CodeBlock(_)) => {
                    in_code_block = true;
                    buf.clear();
                }
                Event::End(TagEnd::CodeBlock) => {
                    if !buf.is_empty() {
                        for tok in buf.split_whitespace() {
                            out.push(tok.to_string());
                        }
                    }
                    in_code_block = false;
                    buf.clear();
                }
                Event::Text(t) if in_code_block => buf.push_str(&t),
                Event::Code(c) => out.push(c.into_string()),
                _ => {}
            }
        }
        out
    }

    /// Decide whether a candidate code-span text is a plausible HGVS string.
    /// Filtering is liberal — anything that looks shaped like HGVS is kept and
    /// fed to ferro for the real parse decision. We do reject obvious template
    /// placeholders and tokens with unbalanced brackets/parens so the fixture
    /// isn't cluttered with prose-template fragments.
    ///
    /// We deliberately do NOT strip outer parentheses: predicted variants like
    /// `c.(75_77)` and `p.(Trp24Cys)` carry semantically meaningful parens, so
    /// rewriting raw spec tokens would collapse distinct inputs and lose
    /// coverage. Bracket/paren balance below catches genuine prose splits
    /// (e.g. `g.[variant1` from a sentence that wraps across spans).
    fn canonicalize(s: &str) -> Option<String> {
        let s = s.trim();
        if s.is_empty() {
            return None;
        }
        if s.contains(char::is_whitespace) {
            return None;
        }
        const PREFIXES: &[&str] = &["c.", "g.", "m.", "n.", "o.", "p.", "r."];
        let looks_like = PREFIXES.iter().any(|p| s.starts_with(p))
            || PREFIXES.iter().any(|p| s.contains(&format!(":{p}")));
        if !looks_like {
            return None;
        }
        if s.contains('\\') || s.contains('`') {
            return None;
        }
        // Drop placeholder template fragments the spec uses in prose
        // (e.g. `g.[variant1;variant2]`, `p.[`, etc.).
        if s.contains("variant1") || s.contains("variant2") || s.contains("variantN") {
            return None;
        }
        // Drop tokens with unbalanced brackets or parens — these are always
        // prose fragments captured because the spec splits an example across
        // multiple backtick spans. (Brackets/parens balance only across the
        // single span.)
        if s.matches('[').count() != s.matches(']').count()
            || s.matches('(').count() != s.matches(')').count()
        {
            return None;
        }
        Some(s.to_string())
    }

    fn extract_from_syntax_yaml(spec_dir: &Path) -> anyhow::Result<Vec<Candidate>> {
        let path = spec_dir.join("docs/syntax.yaml");
        if !path.exists() {
            return Ok(Vec::new());
        }
        let text = std::fs::read_to_string(&path)?;
        let v: serde_yaml::Value = serde_yaml::from_str(&text)?;
        let mut harvested = Vec::new();
        walk_yaml_examples(&v, &mut harvested);
        let rel = "docs/syntax.yaml".to_string();
        Ok(harvested
            .into_iter()
            .filter_map(|s| {
                canonicalize(&s).map(|input| Candidate {
                    input,
                    source_kind: SourceKind::SyntaxYaml,
                    source_path: rel.clone(),
                    working_group: None,
                })
            })
            .collect())
    }

    fn walk_yaml_examples(v: &serde_yaml::Value, out: &mut Vec<String>) {
        match v {
            serde_yaml::Value::Mapping(m) => {
                for (k, val) in m {
                    if k.as_str() == Some("examples") {
                        if let Some(seq) = val.as_sequence() {
                            for item in seq {
                                if let Some(s) = item.as_str() {
                                    out.push(s.to_string());
                                }
                            }
                        }
                    } else {
                        walk_yaml_examples(val, out);
                    }
                }
            }
            serde_yaml::Value::Sequence(s) => {
                for x in s {
                    walk_yaml_examples(x, out);
                }
            }
            _ => {}
        }
    }

    /// Canonical-shape probes ported from the legacy hgvs_spec_compliance_tests.rs.
    /// These cover spec-conformant shapes that aren't literal v21 spec text.
    const PORTED_LEGACY_PROBES: &[&str] = &[
        // canonical DNA
        "NC_000001.11:g.1234del",
        "NC_000001.11:g.1234_2345del",
        "NC_000001.11:g.123delinsAC",
        "NC_000001.11:g.123_129delinsAC",
        "NC_000001.11:g.1234dup",
        "NC_000001.11:g.1234_2345dup",
        "NC_000001.11:g.1234_1235insACGT",
        "NC_000001.11:g.1234_2345inv",
        "NC_000001.11:g.1234=",
        "NC_000023.10:g.33038255C>A",
        // canonical RNA
        "NM_004006.3:r.123_127del",
        "NM_004006.3:r.123_127delinsag",
        "NM_004006.3:r.123_345dup",
        "NM_004006.3:r.123_124insauc",
        "NM_004006.3:r.123_345inv",
        "NM_004006.3:r.123c>g",
        // canonical protein
        "NP_003997.2:p.Val7del",
        "NP_003997.2:p.Lys23_Val25del",
        "NP_003997.2:p.Val7dup",
        "NP_003997.2:p.Lys23_Val25dup",
        "NP_003997.1:p.Trp24Cys",
        "NP_003997.1:p.Trp24Ter",
        "NP_003997.1:p.W24*",
        // canonical coding DNA
        "NM_004006.2:c.145_147delinsTGG",
        "NM_004006.2:c.20dup",
        "NM_004006.2:c.20_23dup",
        "NM_004006.2:c.5657_5660inv",
        "NM_004006.2:c.4145_4160inv",
        "NM_004006.2:c.4375C>T",
    ];

    fn ported_legacy_probes() -> Vec<Candidate> {
        PORTED_LEGACY_PROBES
            .iter()
            .map(|s| Candidate {
                input: (*s).to_string(),
                source_kind: SourceKind::PortedLegacyProbe,
                source_path: "tests/hgvs_spec_compliance_tests.rs (legacy)".to_string(),
                working_group: None,
            })
            .collect()
    }
}

// ---------- classify ----------

mod classify {
    pub fn classify_from_run(parse_ok: bool, normalize_ok: bool) -> &'static str {
        match (parse_ok, normalize_ok) {
            (false, _) => "reject",
            (true, false) => "needs-reference",
            (true, true) => "compliant",
        }
    }
}

// ---------- overrides ----------

mod overrides {
    use super::*;

    #[derive(Debug, Default, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct Overrides {
        #[serde(default)]
        pub by_input: BTreeMap<String, OverrideEntry>,
    }

    #[derive(Debug, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct OverrideEntry {
        #[serde(default)]
        pub status: Option<String>,
        #[serde(default)]
        pub spec_expected: Option<String>,
        #[serde(default)]
        pub todo: Option<String>,
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
        pub current: String,
        pub spec_expected: String,
        pub status: String,
        pub coordinate_system: String,
        pub source_kind: String,
        pub source_paths: Vec<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        pub working_group: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        pub todo: Option<String>,
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

    pub fn build_rows(
        candidates: &[sources::Candidate],
        overrides: &overrides::Overrides,
    ) -> anyhow::Result<Vec<Row>> {
        #[derive(Default)]
        struct Agg {
            kinds: Vec<sources::SourceKind>,
            paths: Vec<String>,
            wg: Option<String>,
        }
        let mut by_input: BTreeMap<String, Agg> = BTreeMap::new();
        for c in candidates {
            let a = by_input.entry(c.input.clone()).or_default();
            a.kinds.push(c.source_kind);
            a.paths.push(c.source_path.clone());
            if c.working_group.is_some() {
                a.wg = c.working_group.clone();
            }
        }

        let normalizer = Normalizer::new(MockProvider::new());

        let mut rows = Vec::with_capacity(by_input.len());
        for (input, agg) in by_input {
            let (current, parse_ok, normalize_ok) = run_ferro(&normalizer, &input);
            let auto_status = classify::classify_from_run(parse_ok, normalize_ok);

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
            let status = ov
                .and_then(|o| o.status.clone())
                .unwrap_or_else(|| auto_status.to_string());
            // spec_expected default depends on status:
            //   - reject:    default = current. Reject rows are dominated by
            //                spec-conformant rejections (malformed prose
            //                fragments, unsupported syntax that the spec also
            //                rejects). Defaulting spec_expected to ferro's
            //                current rejection means "by default, assume ferro
            //                is doing the right thing per spec." Reviewers
            //                override (e.g. set spec_expected to input) for
            //                the rejects ferro is genuinely wrong about (e.g.
            //                #83 A.3 chr1: aliases).
            //   - other:     default = input. The spec offered the input as
            //                a canonical form; if ferro rewrites it (current
            //                != input), that's a candidate pair row.
            let spec_expected = ov.and_then(|o| o.spec_expected.clone()).unwrap_or_else(|| {
                if status == "reject" {
                    current.clone()
                } else {
                    input.clone()
                }
            });
            // Auto-attach a #83 todo whenever current diverges from
            // spec_expected. With the per-status default above, this
            // attaches todos only to genuine audit items.
            const ISSUE_83: &str = "https://github.com/fulcrumgenomics/ferro-hgvs/issues/83";
            let todo = ov.and_then(|o| o.todo.clone()).or_else(|| {
                if current != spec_expected {
                    Some(ISSUE_83.to_string())
                } else {
                    None
                }
            });

            rows.push(Row {
                coordinate_system: coord_system(&input),
                input: input.clone(),
                current,
                spec_expected,
                status,
                source_kind: source_kind_str(kind).to_string(),
                source_paths: paths,
                working_group: agg.wg,
                todo,
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

    fn run_ferro(normalizer: &Normalizer<MockProvider>, input: &str) -> (String, bool, bool) {
        match parse_hgvs(input) {
            Err(e) => (format!("parse error: {e}"), false, false),
            Ok(v) => match normalizer.normalize(&v) {
                Err(e) => (format!("normalize error: {e}"), true, false),
                Ok(n) => (format!("{}", n), true, true),
            },
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
        ferro_version: &'a str,
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
        let out = std::process::Command::new("git")
            .arg("-C")
            .arg(spec_dir)
            .args(["rev-parse", "HEAD"])
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
            ferro_version: env!("CARGO_PKG_VERSION"),
            generated_utc: "fixture-byte-stable".to_string(),
            summary,
            rows,
        };
        let mut out = serde_json::to_string_pretty(&doc)?;
        out.push('\n');
        Ok(out)
    }
}
