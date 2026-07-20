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
        pub intent: SpanIntent,
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum SourceKind {
        Recommendation,
        Background,
        Consultation,
        SyntaxYaml,
        PortedLegacyProbe,
    }

    /// Whether the spec marks this codespan as a rejection example.
    /// `<code class="invalid">…</code>` → `SpecRejects`; backticked codespans,
    /// fenced blocks, syntax.yaml entries, and ported probes → `Plain`.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum SpanIntent {
        Plain,
        SpecRejects,
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
                        intent: SpanIntent::Plain,
                    });
                }
            }
            for cand_text in extract_class_invalid_codespans(&text) {
                if let Some(input) = canonicalize(&cand_text) {
                    all.push(Candidate {
                        input,
                        source_kind: kind,
                        source_path: rel.clone(),
                        working_group: wg.clone(),
                        intent: SpanIntent::SpecRejects,
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

    /// Return every variant string the spec explicitly marks as invalid via
    /// `<code class="invalid">…</code>`. The marker is the spec's own
    /// typographic convention for "this string is wrong" — used 35 times
    /// across `recommendations/` in v21.0. Captures only flat content (no
    /// nested HTML); nested cases would need bespoke parsing and don't occur
    /// in v21.0.
    fn extract_class_invalid_codespans(md: &str) -> Vec<String> {
        use regex::Regex;
        let re = Regex::new(r#"<code\s+class="invalid">([^<]*)</code>"#).unwrap();
        re.captures_iter(md).map(|c| c[1].to_string()).collect()
    }

    /// Return every code span (inline `…` and fenced ``` …``` block) text in the markdown.
    ///
    /// Directly-adjacent inline-code spans are glued back together, including the
    /// inner text of any `<code class="spotN">…</code>` highlight wedged between
    /// them. The spec typesets a single variant across multiple code spans when it
    /// wants to highlight a fragment — an edit (`` `g.123_456`<code…>dup</code> ``),
    /// an LRG transcript/protein suffix
    /// (`` `LRG_199`<code…>t1</code>`:c.11T>G` ``) — so without reassembly the
    /// harvester would see severed pieces (`LRG_199`, `:c.11T>G`) instead of the
    /// intact variant. A plain-text run, line break, or block boundary ends a run.
    fn extract_codespans(md: &str) -> Vec<String> {
        use pulldown_cmark::{Event, Parser, Tag, TagEnd};
        let mut out = Vec::new();
        let mut in_code_block = false;
        let mut block_buf = String::new();
        // Current run of adjacent inline-code + spot-highlight text.
        let mut run = String::new();
        let mut in_spot = false;

        macro_rules! flush_run {
            () => {{
                if !run.is_empty() {
                    out.push(std::mem::take(&mut run));
                }
                in_spot = false;
            }};
        }

        for ev in Parser::new(md) {
            match ev {
                Event::Start(Tag::CodeBlock(_)) => {
                    flush_run!();
                    in_code_block = true;
                    block_buf.clear();
                }
                Event::End(TagEnd::CodeBlock) => {
                    if !block_buf.is_empty() {
                        for tok in block_buf.split_whitespace() {
                            out.push(tok.to_string());
                        }
                    }
                    in_code_block = false;
                    block_buf.clear();
                }
                Event::Text(t) if in_code_block => block_buf.push_str(&t),
                // Highlighted fragment inside a `<code class="spotN">…</code>` span
                // continues the current variant run.
                Event::Text(t) if in_spot => run.push_str(&t),
                // Adjacent inline-code spans accumulate into one run.
                Event::Code(c) => run.push_str(&c),
                Event::InlineHtml(h) => {
                    let h = h.as_ref();
                    if h.starts_with("<code") && h.contains("class=\"spot") {
                        in_spot = true; // continue the run through the highlight
                    } else if in_spot && h.trim() == "</code>" {
                        in_spot = false; // highlight closed; run continues
                    } else {
                        flush_run!(); // any other inline HTML ends the run
                    }
                }
                // Plain prose text, line breaks, and structural events end a run.
                _ => flush_run!(),
            }
        }
        if !run.is_empty() {
            out.push(run);
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
        // HTML remnants leaked from the spec's `<code class="invalid">…</code>`
        // negative examples inside indented (code-block) prose: no real HGVS
        // string contains `<` or `"`. (`>` is legitimate — substitutions.)
        if s.contains('<') || s.contains('"') {
            return None;
        }
        // Drop placeholder template fragments the spec uses in prose
        // (e.g. `g.[variant1;variant2]`, `p.[`, etc.).
        if s.contains("variant1") || s.contains("variant2") || s.contains("variantN") {
            return None;
        }
        // Drop the spec's syntax-template diagrams: a `#` count placeholder
        // (`delN[#]`) or an English placeholder word joined by a hyphen
        // (`fragment-start`, `last-normal`). Real HGVS never contains `#`, and
        // its hyphens are digit-flanked intronic offsets (`c.1210-33`), never a
        // lowercase-word-hyphen-word run.
        if s.contains('#') || contains_word_hyphen_word(s) {
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

    /// True if `s` contains a `<lowercase-word>-<lowercase-word>` run — the
    /// spec's placeholder-word convention (`fragment-start`, `last-normal`).
    /// HGVS hyphens are always digit-flanked intronic offsets, so this only
    /// fires on syntax-template diagrams, never on real variants.
    fn contains_word_hyphen_word(s: &str) -> bool {
        let b = s.as_bytes();
        for i in 0..b.len() {
            if b[i] == b'-' {
                let left = i >= 2 && b[i - 1].is_ascii_lowercase() && b[i - 2].is_ascii_lowercase();
                let right = i + 2 < b.len()
                    && b[i + 1].is_ascii_lowercase()
                    && b[i + 2].is_ascii_lowercase();
                if left && right {
                    return true;
                }
            }
        }
        false
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
                    intent: SpanIntent::Plain,
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
                intent: SpanIntent::Plain,
            })
            .collect()
    }

    #[cfg(test)]
    mod tests {
        use super::{canonicalize, extract_codespans};

        /// The spec highlights part of a variant with `<code class="spotN">…</code>`
        /// wedged between adjacent inline-code spans, e.g.
        /// `` `LRG_199`<code class="spot1">t1</code>`:c.11T>G` ``. The harvester must
        /// reassemble the intact variant, not emit the severed `:c.11T>G` tail.
        #[test]
        fn reassembles_lrg_spot_split_variant() {
            let md = "- e.g., `LRG_199`<code class=\"spot1\">t1</code>`:c.11T>G`.\n";
            let spans = extract_codespans(md);
            assert!(
                spans.iter().any(|s| s == "LRG_199t1:c.11T>G"),
                "expected reassembled `LRG_199t1:c.11T>G`, got {spans:?}"
            );
            assert!(
                !spans.iter().any(|s| s == ":c.11T>G"),
                "severed orphan tail `:c.11T>G` must not be emitted, got {spans:?}"
            );
        }

        /// A spot-highlighted *edit* (e.g. `dup`) trailing an inline-code position
        /// span must fold back into the position to form the whole variant.
        #[test]
        fn reassembles_spot_highlighted_edit() {
            let md = "text `g.123_456`<code class=\"spot1\">dup</code> more\n";
            let spans = extract_codespans(md);
            assert!(
                spans.iter().any(|s| s == "g.123_456dup"),
                "expected reassembled `g.123_456dup`, got {spans:?}"
            );
        }

        /// Plain prose text between two inline-code spans is a boundary: the spans
        /// must NOT be glued across it.
        #[test]
        fn does_not_join_across_prose() {
            let md = "the `g.123del` and `g.456dup` variants\n";
            let spans = extract_codespans(md);
            assert!(spans.iter().any(|s| s == "g.123del"), "got {spans:?}");
            assert!(spans.iter().any(|s| s == "g.456dup"), "got {spans:?}");
            assert!(
                !spans.iter().any(|s| s.contains("g.123delg.456dup")),
                "adjacent-but-prose-separated spans must not glue, got {spans:?}"
            );
        }

        /// HTML remnants leaked from indented (code-block) `<code class="invalid">`
        /// examples must be rejected by `canonicalize` — no real HGVS string
        /// contains `<` or `"`.
        #[test]
        fn canonicalize_rejects_html_remnants() {
            assert_eq!(
                canonicalize("class=\"invalid\">NM_000109.3:c.-401C>T</code>"),
                None
            );
            assert_eq!(canonicalize("NM_004006.1:c.5690</code>"), None);
        }

        /// The spec's syntax-template diagrams use hyphenated English placeholder
        /// words (`fragment-start`, `last-normal`) and a `#` count placeholder —
        /// these are typeset examples of *syntax*, not variants, and must be
        /// dropped. Real HGVS never contains a lowercase-word-hyphen-word run
        /// (intronic offsets are digit-hyphenated) nor `#`.
        #[test]
        fn canonicalize_rejects_syntax_templates() {
            assert_eq!(canonicalize("g.(fragment-start_fragment-end)delN[#]"), None);
            assert_eq!(
                canonicalize("g.(last-normal_first-duplicated)_(last-duplicated_first-normal)dup"),
                None
            );
        }

        /// The template guards must NOT reject legitimate variants: lowercase RNA
        /// bases and digit-hyphenated intronic offsets are fine.
        #[test]
        fn canonicalize_keeps_legit_rna_and_intronic() {
            assert_eq!(
                canonicalize("r.-125_-123cug[4]"),
                Some("r.-125_-123cug[4]".to_string())
            );
            assert_eq!(
                canonicalize("NM_004006.2:c.1210-33del"),
                Some("NM_004006.2:c.1210-33del".to_string())
            );
        }
    }
}

// ---------- prefix ----------

mod prefix {
    /// Default reference accessions for bare illustrative fragments
    /// (e.g. `c.1083A>C` with no accession). The spec routinely omits the
    /// accession when the surrounding paragraph implies one; ferro requires a
    /// full accession to parse, so we prepend a sensible default before
    /// feeding the input to the normalizer.
    ///
    /// The choices match the most-cited accessions in the v21.0 spec corpus.
    /// They are load-bearing constants — re-validate when bumping the spec
    /// submodule (see CONTRIBUTING.md).
    pub const DEFAULTS: &[(char, &str)] = &[
        ('c', "NM_004006.2"),
        ('n', "NR_002196.1"),
        ('r', "NM_004006.3"),
        ('g', "NC_000023.11"),
        ('p', "NP_003997.1"),
        ('m', "NC_012920.1"),
        ('o', "NC_000023.11"),
    ];

    /// Compute the default-prefixed form for an input. Returns `None` when
    /// the input already carries an accession (`<acc>:<sys>.…`), uses an
    /// allele-list shape we can't unambiguously rewrite, or has an unknown
    /// coord system.
    pub fn default_prefixed(input: &str) -> Option<String> {
        // Allele lists like `[c.X;c.Y]` are wrapped; we don't auto-prefix
        // because the right form is `<acc>:c.[X;Y]`, not `[<acc>:c.X;c.Y]`.
        // Auditors can force a value via overrides.input_prefixed if needed.
        if input.starts_with('[') {
            return None;
        }
        // Only the top-level head (before any `[...]` inserted-sequence
        // reference) determines the coordinate system and whether an accession
        // is already present. An inner reference like `ins[L37425.1:r.23_361]`
        // carries its own `:r.` that must not be mistaken for the variant's own
        // accession — otherwise a bare `r.123_124ins[…:r.…]` is wrongly judged
        // already-prefixed and never gets a default accession, so it fails to
        // parse and lands in `parse-error`. (#955)
        let head = &input[..input.find('[').unwrap_or(input.len())];
        let coord = coord_system_char(head)?;
        // Already prefixed: any colon before the coord-system letter in the
        // head implies a `<accession>:<sys>.` form already in place.
        let needle = format!(":{coord}.");
        if head.contains(&needle) {
            return None;
        }
        // Bare form: must start with `<sys>.`.
        let bare = format!("{coord}.");
        if !input.starts_with(&bare) {
            return None;
        }
        let prefix = DEFAULTS
            .iter()
            .find(|(c, _)| *c == coord)
            .map(|(_, p)| *p)?;
        Some(format!("{prefix}:{input}"))
    }

    fn coord_system_char(input: &str) -> Option<char> {
        for (c, _) in DEFAULTS {
            let starts = format!("{c}.");
            let mid = format!(":{c}.");
            if input.starts_with(&starts) || input.contains(&mid) {
                return Some(*c);
            }
        }
        None
    }

    #[cfg(test)]
    mod tests {
        use super::default_prefixed;

        #[test]
        fn prefixes_bare_fragment_with_inner_reference() {
            // The inner `[…:r.…]` reference must NOT be read as the variant's own
            // accession; the bare `r.` fragment still gets a default accession.
            assert_eq!(
                default_prefixed("r.123_124ins[L37425.1:r.23_361]").as_deref(),
                Some("NM_004006.3:r.123_124ins[L37425.1:r.23_361]")
            );
            assert_eq!(
                default_prefixed("g.?_?ins[NC_000023.10:g.(12345_23456)_(34567_45678)]").as_deref(),
                Some("NC_000023.11:g.?_?ins[NC_000023.10:g.(12345_23456)_(34567_45678)]")
            );
        }

        #[test]
        fn leaves_already_prefixed_input_untouched() {
            // A genuine top-level accession still short-circuits (no double prefix).
            assert_eq!(default_prefixed("NM_004006.3:r.123_124insauc"), None);
        }
    }
}

// ---------- classify ----------

mod classify {
    /// Status taxonomy is a function of (parse_ok, normalize_ok, spec_expected,
    /// current == spec_expected). `spec_expected: None` is the spec-rejects-this
    /// sentinel (set by override or by `<code class="invalid">` extraction).
    ///
    /// | spec_expected | parse | normalize | current == spec_expected | status              |
    /// |---------------|-------|-----------|--------------------------|---------------------|
    /// | `None`        | err   | —         | —                        | `correctly-rejected`|
    /// | `None`        | ok    | —         | —                        | `false-acceptance`  |
    /// | `Some`        | err   | —         | —                        | `parse-error`       |
    /// | `Some`        | ok    | err       | —                        | `needs-reference`   |
    /// | `Some`        | ok    | ok        | true                     | `preserved`         |
    /// | `Some`        | ok    | ok        | false                    | `diverges`          |
    pub fn classify(
        parse_ok: bool,
        normalize_ok: bool,
        current: &str,
        spec_expected: Option<&str>,
    ) -> &'static str {
        match (parse_ok, normalize_ok, spec_expected) {
            (false, _, None) => "correctly-rejected",
            (false, _, Some(_)) => "parse-error",
            // Spec rejects the input: any successful parse is a false-acceptance,
            // regardless of whether ferro's normalize() then errored.
            (true, _, None) => "false-acceptance",
            (true, false, Some(_)) => "needs-reference",
            (true, true, Some(expected)) if current == expected => "preserved",
            (true, true, Some(_)) => "diverges",
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
