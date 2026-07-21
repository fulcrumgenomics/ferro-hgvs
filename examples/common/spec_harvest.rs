//! Shared HGVS-spec harvesting machinery for the `examples/` generators.
//!
//! Extracted verbatim from `generate_spec_fixture.rs` so a second generator
//! (`generate_spec_enumeration.rs`) can reuse the same discovery, prefixing and
//! classification logic instead of duplicating it. Behaviour is unchanged:
//! `generate_spec_fixture --check` must stay byte-identical across this move.
//!
//! This file lives under `examples/common/` (a subdirectory without a
//! `main.rs`), so Cargo does not build it as an example target of its own;
//! each generator pulls it in with `#[path = "common/spec_harvest.rs"] mod`.
#![allow(dead_code)]

use std::path::{Path, PathBuf};

pub mod sources {
    use super::*;
    use anyhow::Context as _;
    use std::ffi::OsStr;

    #[derive(Debug)]
    pub struct Candidate {
        pub input: String,
        pub source_kind: SourceKind,
        pub source_path: String,
        /// 1-based line in `source_path` the span was harvested from. `None`
        /// for candidates with no markdown origin (syntax.yaml, ported probes).
        /// Provenance only — never serialized into the normalization fixture,
        /// which must stay byte-stable across this field's introduction.
        pub source_line: Option<usize>,
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
        for (path, kind, wg) in whitelisted_files(spec_dir)? {
            let text = std::fs::read_to_string(&path)?;
            let rel = path
                .strip_prefix(spec_dir)
                .unwrap_or(&path)
                .to_string_lossy()
                .replace('\\', "/");
            for (cand_text, offset) in extract_codespans(&text) {
                if let Some(input) = canonicalize(&cand_text) {
                    all.push(Candidate {
                        input,
                        source_kind: kind,
                        source_path: rel.clone(),
                        source_line: Some(line_of(&text, offset)),
                        working_group: wg.clone(),
                        intent: SpanIntent::Plain,
                    });
                }
            }
            for (cand_text, offset) in extract_class_invalid_codespans(&text) {
                if let Some(input) = canonicalize(&cand_text) {
                    all.push(Candidate {
                        input,
                        source_kind: kind,
                        source_path: rel.clone(),
                        source_line: Some(line_of(&text, offset)),
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

    fn whitelisted_files(
        spec_dir: &Path,
    ) -> anyhow::Result<Vec<(PathBuf, SourceKind, Option<String>)>> {
        let mut out = Vec::new();
        let docs = spec_dir.join("docs");

        let recs = docs.join("recommendations");
        for entry in walkdir_md(&recs)? {
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
        Ok(out)
    }

    /// Public alias of [`walkdir_md`] for the `negatives` harvester.
    pub fn walkdir_md_pub(dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
        walkdir_md(dir)
    }

    /// Public alias of [`canonicalize`] for the `negatives` harvester.
    pub fn canonicalize_pub(s: &str) -> Option<String> {
        canonicalize(s)
    }

    fn walkdir_md(dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
        let mut out = Vec::new();
        if !dir.exists() {
            return Ok(out);
        }
        let mut stack = vec![dir.to_path_buf()];
        while let Some(d) = stack.pop() {
            // Propagate read failures rather than warn-and-continue: a directory
            // the harvester cannot read yields a silently incomplete candidate
            // set, and this feeds a conformance fixture, so an incomplete scan
            // must fail generation (matching `discover`'s error propagation).
            let entries = std::fs::read_dir(&d)
                .with_context(|| format!("read spec documentation directory {}", d.display()))?;
            for e in entries {
                let e =
                    e.with_context(|| format!("read a directory entry under {}", d.display()))?;
                let p = e.path();
                if p.is_dir() {
                    stack.push(p);
                } else if p.extension() == Some(OsStr::new("md")) {
                    out.push(p);
                }
            }
        }
        out.sort();
        Ok(out)
    }

    /// Return every variant string the spec explicitly marks as invalid via
    /// `<code class="invalid">…</code>`. The marker is the spec's own
    /// typographic convention for "this string is wrong" — used 35 times
    /// across `recommendations/` in v21.0. Captures only flat content (no
    /// nested HTML); nested cases would need bespoke parsing and don't occur
    /// in v21.0.
    pub fn extract_class_invalid_codespans(md: &str) -> Vec<(String, usize)> {
        use regex::Regex;
        let re = Regex::new(r#"<code\s+class="invalid">([^<]*)</code>"#).unwrap();
        re.captures_iter(md)
            .map(|c| {
                let whole = c.get(0).expect("group 0 always present");
                (c[1].to_string(), whole.start())
            })
            .collect()
    }

    /// 1-based line number of a byte offset within `md`.
    pub fn line_of(md: &str, offset: usize) -> usize {
        md[..offset.min(md.len())]
            .bytes()
            .filter(|b| *b == b'\n')
            .count()
            + 1
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
    pub fn extract_codespans(md: &str) -> Vec<(String, usize)> {
        use pulldown_cmark::{Event, Parser, Tag, TagEnd};
        let mut out: Vec<(String, usize)> = Vec::new();
        let mut in_code_block = false;
        let mut block_buf = String::new();
        let mut block_start = 0usize;
        // Current run of adjacent inline-code + spot-highlight text, and the
        // byte offset in `md` where the run began (for file:line provenance).
        let mut run = String::new();
        let mut run_start = 0usize;
        let mut in_spot = false;

        macro_rules! flush_run {
            () => {{
                if !run.is_empty() {
                    out.push((std::mem::take(&mut run), run_start));
                }
                in_spot = false;
            }};
        }

        for (ev, range) in Parser::new(md).into_offset_iter() {
            match ev {
                Event::Start(Tag::CodeBlock(_)) => {
                    flush_run!();
                    in_code_block = true;
                    block_buf.clear();
                    block_start = range.start;
                }
                Event::End(TagEnd::CodeBlock) => {
                    if !block_buf.is_empty() {
                        for tok in block_buf.split_whitespace() {
                            out.push((tok.to_string(), block_start));
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
                Event::Code(c) => {
                    if run.is_empty() {
                        run_start = range.start;
                    }
                    run.push_str(&c);
                }
                Event::InlineHtml(h) => {
                    let h = h.as_ref();
                    if h.starts_with("<code") && h.contains("class=\"spot") {
                        in_spot = true; // continue the run through the highlight
                        if run.is_empty() {
                            run_start = range.start;
                        }
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
            out.push((run, run_start));
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
                    source_line: None,
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
                source_line: None,
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
                spans.iter().any(|(s, _)| s == "LRG_199t1:c.11T>G"),
                "expected reassembled `LRG_199t1:c.11T>G`, got {spans:?}"
            );
            assert!(
                !spans.iter().any(|(s, _)| s == ":c.11T>G"),
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
                spans.iter().any(|(s, _)| s == "g.123_456dup"),
                "expected reassembled `g.123_456dup`, got {spans:?}"
            );
        }

        /// Plain prose text between two inline-code spans is a boundary: the spans
        /// must NOT be glued across it.
        #[test]
        fn does_not_join_across_prose() {
            let md = "the `g.123del` and `g.456dup` variants\n";
            let spans = extract_codespans(md);
            assert!(spans.iter().any(|(s, _)| s == "g.123del"), "got {spans:?}");
            assert!(spans.iter().any(|(s, _)| s == "g.456dup"), "got {spans:?}");
            assert!(
                !spans.iter().any(|(s, _)| s.contains("g.123delg.456dup")),
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

pub mod prefix {
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

pub mod classify {
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

// ---------- negatives ----------

/// Harvester for the spec's **prose** negatives — sentences that forbid a form
/// in words rather than by wrapping it in `<code class="invalid">`.
///
/// The spec follows RFC 2119 (`docs/recommendations/general/style.md`), and is
/// deliberate about its keywords, so the matched phrase determines the
/// normative level: only MUST-level phrasing may ever hard-fail a test.
pub mod negatives {
    use super::sources::{extract_class_invalid_codespans, extract_codespans, line_of};
    use super::*;

    /// RFC 2119 strength of a prose rule.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum NormativeLevel {
        /// "is not correct", "not allowed", "can not be …" — may hard-fail.
        Must,
        /// "should not be", "the recommendation is not to", "preferably" —
        /// report-only, never a hard failure.
        Should,
    }

    impl NormativeLevel {
        pub fn as_str(self) -> &'static str {
            match self {
                NormativeLevel::Must => "must",
                NormativeLevel::Should => "should",
            }
        }
    }

    /// MUST-level prose markers, longest-first so the more specific phrase wins.
    const MUST_MARKERS: &[&str] = &[
        "is not correct",
        "are not correct",
        "not allowed",
        "can not be",
        "cannot be",
        "not described as",
        "is not a valid",
        "is incorrect",
    ];

    /// SHOULD-level prose markers (RFC 2119 SHOULD NOT / recommendation).
    const SHOULD_MARKERS: &[&str] = &[
        "should not be",
        "the recommendation is not",
        "recommendation is not to",
        "should preferably",
        "is discouraged",
    ];

    /// One markdown line that states a negative rule, with every HGVS-shaped
    /// code span found on that line.
    #[derive(Debug)]
    pub struct ProseNegative {
        pub source_path: String,
        pub source_line: usize,
        /// The prose marker that matched (verbatim).
        pub marker: &'static str,
        pub normative_level: NormativeLevel,
        /// The full markdown line, trimmed.
        pub sentence: String,
        /// HGVS-shaped strings on the line that the spec marked invalid.
        pub invalid_spans: Vec<String>,
        /// HGVS-shaped plain code spans on the line — repair-target candidates.
        /// These are *candidates only*: the spec frequently cites a bare
        /// coordinate ("`c.-244`") on the same line, so a repair target is only
        /// ever asserted when a hand-curated override names it.
        pub plain_spans: Vec<String>,
    }

    /// Classify a line by its strongest negative marker, if any.
    fn marker_for(line: &str) -> Option<(&'static str, NormativeLevel)> {
        // SHOULD is checked first: "should not be" contains no MUST marker, but
        // keeping the order explicit documents that SHOULD never upgrades.
        for m in SHOULD_MARKERS {
            if line.contains(m) {
                return Some((m, NormativeLevel::Should));
            }
        }
        for m in MUST_MARKERS {
            if line.contains(m) {
                return Some((m, NormativeLevel::Must));
            }
        }
        None
    }

    /// Scan the whitelisted recommendation files for prose negatives.
    pub fn discover(spec_dir: &Path) -> anyhow::Result<Vec<ProseNegative>> {
        let mut out = Vec::new();
        let recs = spec_dir.join("docs/recommendations");
        for path in super::sources::walkdir_md_pub(&recs)? {
            let text = std::fs::read_to_string(&path)?;
            let rel = path
                .strip_prefix(spec_dir)
                .unwrap_or(&path)
                .to_string_lossy()
                .replace('\\', "/");
            // Line-scoped harvesting: the spec writes one rule per list item /
            // paragraph line, so a line is the natural rule unit. Spans are
            // re-extracted per line (rather than sliced out of the whole-file
            // pass) so the offsets stay local and unambiguous.
            for (idx, line) in text.split('\n').enumerate() {
                let Some((marker, level)) = marker_for(line) else {
                    continue;
                };
                let invalid_spans: Vec<String> = extract_class_invalid_codespans(line)
                    .into_iter()
                    .filter_map(|(s, _)| super::sources::canonicalize_pub(&s))
                    .collect();
                let plain_spans: Vec<String> = extract_codespans(line)
                    .into_iter()
                    .filter_map(|(s, _)| super::sources::canonicalize_pub(&s))
                    .collect();
                let _ = line_of; // offsets are line-local here
                out.push(ProseNegative {
                    source_path: rel.clone(),
                    source_line: idx + 1,
                    marker,
                    normative_level: level,
                    sentence: line.trim().to_string(),
                    invalid_spans,
                    plain_spans,
                });
            }
        }
        out.sort_by(|a, b| {
            a.source_path
                .cmp(&b.source_path)
                .then(a.source_line.cmp(&b.source_line))
        });
        Ok(out)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn must_and_should_markers_are_distinguished() {
            assert_eq!(
                marker_for("a one-nucleotide inversion is not allowed").map(|m| m.1),
                Some(NormativeLevel::Must)
            );
            assert_eq!(
                marker_for("the recommendation is not to use this form").map(|m| m.1),
                Some(NormativeLevel::Should)
            );
            assert_eq!(marker_for("a plain sentence about numbering"), None);
        }

        /// "should not be" is RFC 2119 SHOULD NOT — it must never be promoted to
        /// MUST just because a MUST marker also appears later in the sentence.
        #[test]
        fn should_not_be_never_upgrades_to_must() {
            assert_eq!(
                marker_for("this should not be used, it is not correct either").map(|m| m.1),
                Some(NormativeLevel::Should)
            );
        }
    }
}
