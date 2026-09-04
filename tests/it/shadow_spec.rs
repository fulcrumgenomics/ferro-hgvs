//! Shadow-spec checker.
//!
//! Verifies that every authored shadow page under `docs/src/shadow-spec/` stays true to the
//! pinned HGVS spec and to ferro's real behavior. The page is the source of truth and the
//! published (mdBook) deliverable; the parts the machine checks are the same bytes the reader
//! sees ("visible is checked"), so there is no second copy to drift.
//!
//! Per anchor (a heading beginning with a `path:line` inline-code token):
//!   1. quote     — the block-quote is a (whitespace-collapsed) substring of the cited spec span
//!   2. examples  — every `| Input | Verdict | Normalizes to | Notes |` row is executed
//!   3. links     — every ruling-id link resolves in the ledger; every `file.md:line` cross-ref resolves
//!   4. coverage  — reported (not gated): which spec lines the page anchors
//!
//! Verdicts describe ferro's OUTPUT against the HGVS recommendations:
//!   recommended -> ferro's output is the recommended form (input parses; `Normalizes to` = self,
//!                  or the recommended form ferro rewrites the input into)
//!   conformant  -> ferro's output is valid HGVS but not yet the recommended form — a limitation
//!                  (input parses; `Normalizes to` shows the current output)
//!   refused     -> the input is not valid HGVS; ferro rejects it in strict mode (correct)
//!   bug         -> ferro's output is not valid HGVS (input parses; `Normalizes to` pins the
//!                  invalid output, which must not re-parse). None expected on a healthy page.
//!
//! Reference tiers for `Normalizes to`: the committed GRCh38 slice gates every PR (hermetic);
//! when `FERRO_MANIFEST` is set the same assertions cross-check the full live reference.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, MultiFastaProvider, NormalizeConfig, Normalizer, ReferenceProvider};

const SHADOW_ROOT: &str = "docs/src/shadow-spec";
const SPEC_ROOT: &str = "assets/hgvs-nomenclature/docs";
const LEDGER: &str = "tests/fixtures/grammar/hgvs_spec_normalization_overrides.json";
/// Committed real-GRCh38 transcript slice — the gating reference for `Normalizes to`, so
/// normalization is checked on every PR without a multi-GB reference. Regenerate with
/// `tests/fixtures/shadow_spec/build_fixture.py`.
const SLICE: &str = "tests/fixtures/shadow_spec/transcripts.json";
/// The single shadow-spec example corpus — committed + blessed, generated from the pages.
/// Kept out of `docs/src/` so mdBook does not try to render it.
const CORPUS: &str = "tests/fixtures/shadow_spec/corpus.jsonl";
/// Blessed per-input alignment sidecar, served to the page JS by mdBook (kept in `docs/src`).
const ALIGNMENTS: &str = "docs/src/shadow-spec/alignments.json";

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn read(path: &Path) -> String {
    std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

/// Collapse all whitespace runs to a single space and trim — the same normalization the
/// ledger generator uses so a re-wrapped spec clause still matches.
fn collapse(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}

/// GitHub URL of the rendered ledger, linked from each transcluded "Why" block.
const CONTRACT_URL: &str =
    "https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md";

/// Ledger ruling ids -> `summary` (None where the record has no summary yet).
fn ledger_summaries() -> BTreeMap<String, Option<String>> {
    let text = read(&repo_root().join(LEDGER));
    let v: serde_json::Value = serde_json::from_str(&text).expect("parse ledger json");
    let mut out = BTreeMap::new();
    for r in v["rulings"].as_array().expect("rulings array") {
        let id = r["id"].as_str().expect("ruling id").to_string();
        let summary = r["summary"].as_str().map(|s| s.to_string());
        out.insert(id, summary);
    }
    out
}

/// The rendered "Why" body for a comma-separated list of ruling ids: one blockquote line
/// per ruling, each carrying the ledger's own `summary` and a link to the full record.
/// Returns `Err(id)` for the first cited id that is missing or has no summary.
fn render_why_body(
    ids: &[String],
    summaries: &BTreeMap<String, Option<String>>,
) -> Result<String, String> {
    let mut lines = Vec::new();
    for id in ids {
        match summaries.get(id) {
            Some(Some(summary)) => {
                lines.push(format!("> **[{id}]({CONTRACT_URL})** — {summary}"));
            }
            _ => return Err(id.clone()),
        }
    }
    Ok(lines.join("\n>\n"))
}

/// Transclude every `<!-- why:START -->…<!-- why:END:id,id -->` block in `text` from the
/// ledger summaries. In bless mode the returned text has each block's body rewritten; otherwise
/// the text is unchanged and a stale/missing block is reported in `errors`.
fn transclude_why_blocks(
    text: &str,
    rel: &str,
    summaries: &BTreeMap<String, Option<String>>,
    bless: bool,
    errors: &mut Vec<String>,
) -> String {
    let start = "<!-- why:START -->";
    let end_prefix = "<!-- why:END:";
    let end_suffix = "-->";
    let mut out = String::new();
    let mut rest = text;
    while let Some(s) = rest.find(start) {
        out.push_str(&rest[..s]);
        out.push_str(start);
        rest = &rest[s + start.len()..];

        let Some(e) = rest.find(end_prefix) else {
            errors.push(format!("{rel}: `why:START` without a matching `why:END`"));
            out.push_str(rest);
            return out;
        };
        let body = &rest[..e];
        let after_prefix = &rest[e + end_prefix.len()..];
        let Some(sfx) = after_prefix.find(end_suffix) else {
            errors.push(format!("{rel}: `why:END` marker is not closed with `-->`"));
            out.push_str(rest);
            return out;
        };
        let ids: Vec<String> = after_prefix[..sfx]
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();
        let end_marker = format!("{end_prefix}{}{end_suffix}", &after_prefix[..sfx]);

        match render_why_body(&ids, summaries) {
            Err(id) => {
                errors.push(format!(
                    "{rel}: `why` block cites '{id}', which is missing from the ledger or has no `summary`"
                ));
                out.push_str(body);
            }
            Ok(rendered) => {
                let want = format!("\n{rendered}\n");
                if bless {
                    out.push_str(&want);
                } else {
                    if body != want {
                        errors.push(format!(
                            "{rel}: `why` block for [{}] is stale — regenerate with BLESS_SHADOW_CORPUS=1",
                            ids.join(", ")
                        ));
                    }
                    out.push_str(body);
                }
            }
        }
        out.push_str(&end_marker);
        rest = &after_prefix[sfx + end_suffix.len()..];
    }
    out.push_str(rest);
    out
}

/// Ledger ruling ids -> status.
fn ledger_ids() -> BTreeMap<String, String> {
    let text = read(&repo_root().join(LEDGER));
    let v: serde_json::Value = serde_json::from_str(&text).expect("parse ledger json");
    let mut out = BTreeMap::new();
    for r in v["rulings"].as_array().expect("rulings array") {
        let id = r["id"].as_str().expect("ruling id").to_string();
        let status = r["status"].as_str().unwrap_or("").to_string();
        out.insert(id, status);
    }
    out
}

/// What an example row asserts about `normalize(input)`.
enum Norm {
    /// No `Normalizes to` cell / `—` — normalization is not asserted for this row.
    Unset,
    /// `self` — the input is a fixed point: `normalize(input) == input`.
    Fixed,
    /// A spelling — `normalize(input)` must equal it (a rewrite, or a pinned gap).
    To(String),
}

struct Example {
    spelling: String,
    verdict: String,
    normalizes_to: Norm,
    row: usize,
}

/// Parse a `Normalizes to` cell into a [`Norm`].
fn parse_norm(cell: &str) -> Norm {
    let t = cell.trim();
    if t.is_empty() || t == "—" || t == "-" {
        return Norm::Unset;
    }
    if t.to_lowercase().contains("self") {
        return Norm::Fixed;
    }
    // The expected output is the first backtick token, else the arrow-stripped text.
    if let Some(tok) = backtick_tokens(t).into_iter().next() {
        return Norm::To(tok);
    }
    Norm::To(t.trim_start_matches('→').trim().to_string())
}

/// Flank width (bases each side) shown around an edit in a shadow-spec alignment.
const ALIGN_FLANK: u64 = 12;
/// Edited spans wider than this are windowed to head + ellipsis + tail.
const ALIGN_WINDOW_MAX: usize = 40;
/// Bases kept at each end of a windowed edit.
const ALIGN_WINDOW_KEEP: usize = 15;

/// Render a base string in the axis's own casing: RNA lowercase with `u`, DNA uppercase.
fn align_cast(s: &str, is_rna: bool) -> String {
    if is_rna {
        s.chars()
            .map(|c| match c {
                'T' | 't' => 'u',
                o => o.to_ascii_lowercase(),
            })
            .collect()
    } else {
        s.to_ascii_uppercase()
    }
}

/// Window a padded allele to `head … tail` when it exceeds the cap; else return it whole.
fn align_window(chars: &[char]) -> String {
    if chars.len() > ALIGN_WINDOW_MAX {
        let head: String = chars[..ALIGN_WINDOW_KEEP].iter().collect();
        let tail: String = chars[chars.len() - ALIGN_WINDOW_KEEP..].iter().collect();
        format!("{head} … {tail}")
    } else {
        chars.iter().collect()
    }
}

/// Per-column identity marker for two equal-length spans: `|` where the base is unchanged
/// (same reference position before and after the edit), a space where it changed.
fn align_connector(a: &[char], b: &[char]) -> String {
    a.iter()
        .zip(b)
        .map(|(x, y)| if x == y { '|' } else { ' ' })
        .collect()
}

/// The pure formatting half of [`build_alignment`]: given the already-cast reference and
/// alternate alleles and their flanks, produce the 3-line alignment JSON
/// (`ref` / `connector` / `alt`, plus the edited-column span `editStart`/`editEnd` and a
/// size `note`). Split out from the provider I/O so it can be unit-tested against
/// hand-computed expectations without a reference. Honest by construction: the connector
/// marks per-column identity only for equal-length edits; a length-changing edit has no
/// per-base correspondence across the changed span, so its connector there is blank.
fn render_alignment(del: &[char], ins: &[char], left: &str, right: &str) -> serde_json::Value {
    let equal = del.len() == ins.len();
    let width = del.len().max(ins.len());
    // Pad both alleles to a common width with `-` so the flanks line up on both rows.
    let pad = |bases: &[char]| -> Vec<char> {
        bases
            .iter()
            .copied()
            .chain(std::iter::repeat_n('-', width - bases.len()))
            .collect()
    };
    let ref_full = pad(del);
    let alt_full = pad(ins);
    let windowed = width > ALIGN_WINDOW_MAX;

    let ref_edit = align_window(&ref_full);
    let alt_edit = align_window(&alt_full);
    let connector_edit: String = if !equal {
        " ".repeat(ref_edit.chars().count())
    } else if !windowed {
        align_connector(&ref_full, &alt_full)
    } else {
        // Head and tail columns correspond position-for-position; the middle is elided.
        let n = ALIGN_WINDOW_KEEP;
        let head = align_connector(&ref_full[..n], &alt_full[..n]);
        let tail = align_connector(
            &ref_full[ref_full.len() - n..],
            &alt_full[alt_full.len() - n..],
        );
        format!("{head}   {tail}") // three columns under " … "
    };

    let note = if equal {
        let changed = del.iter().zip(ins).filter(|(r, a)| r != a).count();
        format!("{width} nt · {changed} changed")
    } else if del.is_empty() {
        format!("insertion · {} nt", ins.len())
    } else if ins.is_empty() {
        format!("deletion · {} nt", del.len())
    } else {
        format!("delins · {} → {} nt", del.len(), ins.len())
    };

    let edit_start = left.chars().count();
    let edit_end = edit_start + ref_edit.chars().count();
    serde_json::json!({
        "ref": format!("{left}{ref_edit}{right}"),
        "connector": format!(
            "{}{connector_edit}{}",
            "|".repeat(edit_start),
            "|".repeat(right.chars().count())
        ),
        "alt": format!("{left}{alt_edit}{right}"),
        "editStart": edit_start,
        "editEnd": edit_end,
        "note": note,
    })
}

/// Build a biopython-style 3-line alignment for one example spelling, computed from its SPDI
/// against the committed slice. Returns `None` for inputs that do not resolve on the slice
/// (foreign accessions, parse-only rows) or that denote no change.
fn build_alignment<P: ReferenceProvider + ?Sized>(
    input: &str,
    provider: &P,
) -> Option<serde_json::Value> {
    let variant = parse_hgvs(input).ok()?;
    let spdi = hgvs_to_spdi(&variant, provider).ok()?;
    let is_rna = input.contains(":r.");
    let del: Vec<char> = align_cast(&spdi.deletion, is_rna).chars().collect();
    let ins: Vec<char> = align_cast(&spdi.insertion, is_rna).chars().collect();
    if del == ins {
        return None; // identity — nothing to show
    }

    let seq = &spdi.sequence;
    let seq_len = provider.get_sequence_length(seq).ok()?;
    let pos = spdi.position;
    let del_len = del.len() as u64;
    let left = align_cast(
        &provider
            .get_sequence(seq, pos.saturating_sub(ALIGN_FLANK), pos)
            .ok()?,
        is_rna,
    );
    // Right flank, clamped to the sequence end: an edit ending exactly at the 3' end has an
    // empty flank rather than a dropped alignment (the 5' end is symmetric via saturating_sub).
    let rstart = pos + del_len;
    let rend = (rstart + ALIGN_FLANK).min(seq_len);
    let right = if rend > rstart {
        align_cast(&provider.get_sequence(seq, rstart, rend).ok()?, is_rna)
    } else {
        String::new()
    };

    Some(render_alignment(&del, &ins, &left, &right))
}

/// Bless (`BLESS_SHADOW_CORPUS=1`) or check a committed, generated artifact: on bless, write
/// `rendered` to `path`; otherwise push a stale-artifact error when the on-disk bytes differ.
fn bless_or_check(path: &Path, label: &str, rendered: &str, bless: bool, errors: &mut Vec<String>) {
    if bless {
        std::fs::write(path, rendered).unwrap_or_else(|e| panic!("write {label}: {e}"));
        eprintln!(
            "shadow-spec: blessed {label} ({} lines)",
            rendered.lines().count()
        );
    } else {
        let on_disk = std::fs::read_to_string(path).unwrap_or_default();
        if on_disk != rendered {
            errors.push(format!(
                "{label} is stale ({} committed lines vs {} rendered) — regenerate with \
                 BLESS_SHADOW_CORPUS=1",
                on_disk.lines().count(),
                rendered.lines().count()
            ));
        }
    }
}

struct Anchor {
    /// spec path relative to `SPEC_ROOT`, e.g. `recommendations/DNA/substitution.md`
    spec_rel: String,
    lstart: usize,
    lend: usize,
    heading_lineno: usize,
    quote: String,
    examples: Vec<Example>,
    ruling_ids: Vec<String>,
    crossrefs: Vec<(String, usize, usize)>, // (spec_rel, lstart, lend)
}

/// Parse `name.md:line` or `name.md:a-b`. Returns (basename, lstart, lend).
fn parse_ref_token(tok: &str) -> Option<(String, usize, usize)> {
    let (name, lines) = tok.split_once(':')?;
    if !name.ends_with(".md") {
        return None;
    }
    let (a, b) = match lines.split_once('-') {
        Some((a, b)) => (a.parse().ok()?, b.parse().ok()?),
        None => {
            let n = lines.parse().ok()?;
            (n, n)
        }
    };
    Some((name.to_string(), a, b))
}

/// All text runs that sit between backticks (inline code), across the whole slice.
fn backtick_tokens(s: &str) -> Vec<String> {
    s.split('`')
        .enumerate()
        .filter(|(i, _)| i % 2 == 1)
        .map(|(_, part)| part.to_string())
        .collect()
}

/// Labels of every `[label](url)` markdown link in the slice.
fn md_link_labels(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut i = 0;
    while let Some(rel) = s[i..].find('[') {
        let open = i + rel;
        if let Some(crel) = s[open + 1..].find(']') {
            let close = open + 1 + crel;
            if s[close + 1..].starts_with('(') {
                out.push(s[open + 1..close].to_string());
                i = close + 1;
            } else {
                i = open + 1;
            }
        } else {
            break;
        }
    }
    out
}

fn looks_like_ruling_id(label: &str) -> bool {
    label.len() >= 8
        && label.contains('-')
        && label
            .chars()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '-')
}

/// Parse one shadow page into anchors. `spec_dir_rel` is the page's mirrored spec directory,
/// e.g. `recommendations/DNA`, used to resolve same-directory basenames.
fn parse_page(text: &str, spec_dir_rel: &str) -> Vec<Anchor> {
    let lines: Vec<&str> = text.lines().collect();

    // Find anchor heading line indices: a heading whose text starts with a backtick `name.md:lines`.
    let mut heads: Vec<(usize, String, usize, usize)> = Vec::new(); // (line_idx, spec_rel, lstart, lend)
    for (idx, raw) in lines.iter().enumerate() {
        let t = raw.trim_start();
        if !t.starts_with('#') {
            continue;
        }
        let after = t.trim_start_matches('#').trim_start();
        if !after.starts_with('`') {
            continue;
        }
        let toks = backtick_tokens(after);
        if let Some(tok) = toks.first() {
            if let Some((name, a, b)) = parse_ref_token(tok) {
                heads.push((idx, format!("{spec_dir_rel}/{name}"), a, b));
            }
        }
    }

    let mut anchors = Vec::new();
    for (h, (idx, spec_rel, a, b)) in heads.iter().enumerate() {
        let body_start = idx + 1;
        let body_end = heads.get(h + 1).map(|n| n.0).unwrap_or(lines.len());
        let body_lines = &lines[body_start..body_end];
        let body = body_lines.join("\n");

        // quote: first contiguous run of `>` lines
        let mut quote = String::new();
        let mut in_q = false;
        for l in body_lines {
            let t = l.trim_start();
            if t.starts_with('>') {
                in_q = true;
                let frag = t.trim_start_matches('>').trim();
                if !quote.is_empty() {
                    quote.push(' ');
                }
                quote.push_str(frag);
            } else if in_q && !t.is_empty() {
                break;
            } else if in_q {
                // blank line inside/after quote -> stop after first run
                break;
            }
        }

        // examples table: `| Input | Verdict | Normalizes to | Notes |`
        // (`Normalizes to` is optional; older 3-column tables still parse.)
        let mut examples = Vec::new();
        let mut cols: Option<(usize, usize, Option<usize>)> = None;
        for (row_off, l) in body_lines.iter().enumerate() {
            let t = l.trim();
            if !t.starts_with('|') {
                cols = None;
                continue;
            }
            let cells: Vec<String> = t
                .trim_matches('|')
                .split('|')
                .map(|c| c.trim().to_string())
                .collect();
            let lowered: Vec<String> = cells.iter().map(|c| c.to_lowercase()).collect();
            if cols.is_none() {
                let sp = lowered.iter().position(|c| c == "input" || c == "spelling");
                let vd = lowered.iter().position(|c| c == "verdict");
                let nm = lowered.iter().position(|c| c == "normalizes to");
                if let (Some(sp), Some(vd)) = (sp, vd) {
                    cols = Some((sp, vd, nm));
                }
                continue;
            }
            // separator row like |---|---|
            if cells
                .iter()
                .all(|c| c.chars().all(|ch| ch == '-' || ch == ':') && !c.is_empty())
            {
                continue;
            }
            let (sp, vd, nm) = cols.unwrap();
            let spelling = cells
                .get(sp)
                .map(|c| c.trim_matches('`').to_string())
                .unwrap_or_default();
            let verdict_cell = lowered.get(vd).cloned().unwrap_or_default();
            let verdict = ["recommended", "conformant", "refused", "bug"]
                .into_iter()
                .find(|k| verdict_cell.split_whitespace().any(|w| w == *k))
                .unwrap_or("")
                .to_string();
            let normalizes_to = nm
                .and_then(|i| cells.get(i))
                .map(|c| parse_norm(c))
                .unwrap_or(Norm::Unset);
            examples.push(Example {
                spelling,
                verdict,
                normalizes_to,
                row: body_start + row_off + 1,
            });
        }

        // ruling ids: link labels that look like ledger ids
        let ruling_ids: Vec<String> = md_link_labels(&body)
            .into_iter()
            .filter(|l| looks_like_ruling_id(l))
            .collect();

        // cross-refs: backtick tokens shaped `name.md:lines`; the basename is resolved at
        // check time by walking up the spec tree (a page cites same-dir siblings and the
        // shared parent-level pages like `general.md`).
        let crossrefs: Vec<(String, usize, usize)> = backtick_tokens(&body)
            .iter()
            .filter_map(|tok| parse_ref_token(tok))
            .collect();

        anchors.push(Anchor {
            spec_rel: spec_rel.clone(),
            lstart: *a,
            lend: *b,
            heading_lineno: idx + 1,
            quote,
            examples,
            ruling_ids,
            crossrefs,
        });
    }
    anchors
}

/// `normalize(input)` as a rendered string, generic over the provider so the committed slice
/// (`MockProvider`) and the full reference (`MultiFastaProvider`) share one call site.
fn normalize_str<P: ferro_hgvs::ReferenceProvider>(
    n: &Normalizer<P>,
    input: &str,
) -> Option<String> {
    parse_hgvs(input)
        .ok()
        .and_then(|v| n.normalize(&v).ok())
        .map(|o| o.to_string())
}

/// Resolve a cross-ref basename to a spec-relative path by walking up from the citing page's
/// spec directory to the spec root; the nearest ancestor holding the file wins. Same-dir
/// siblings resolve immediately; shared pages like `general.md` resolve at `recommendations/`.
fn resolve_crossref(spec_root: &Path, spec_dir_rel: &str, name: &str) -> Option<String> {
    let mut dir = PathBuf::from(spec_dir_rel);
    loop {
        let cand = if dir.as_os_str().is_empty() {
            PathBuf::from(name)
        } else {
            dir.join(name)
        };
        let rel = cand.to_string_lossy().replace('\\', "/");
        if spec_root.join(&rel).is_file() {
            return Some(rel);
        }
        if dir.as_os_str().is_empty() {
            return None;
        }
        dir = dir.parent().map(|p| p.to_path_buf()).unwrap_or_default();
    }
}

/// Lines `a..=b` (1-based) of a spec file, joined and whitespace-collapsed.
fn spec_span(spec_abs: &Path, a: usize, b: usize) -> Option<String> {
    let text = std::fs::read_to_string(spec_abs).ok()?;
    let lines: Vec<&str> = text.lines().collect();
    if a == 0 || b > lines.len() || a > b {
        return None;
    }
    Some(collapse(&lines[a - 1..b].join(" ")))
}

fn all_shadow_pages(root: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    let mut stack = vec![root.to_path_buf()];
    while let Some(dir) = stack.pop() {
        let Ok(rd) = std::fs::read_dir(&dir) else {
            continue;
        };
        for e in rd.flatten() {
            let p = e.path();
            if p.is_dir() {
                stack.push(p);
            } else if p.extension().and_then(|x| x.to_str()) == Some("md") {
                out.push(p);
            }
        }
    }
    out.sort();
    out
}

#[test]
fn shadow_spec_pages_are_current() {
    let root = repo_root();
    let shadow_root = root.join(SHADOW_ROOT);
    let spec_root = root.join(SPEC_ROOT);
    let ledger = ledger_ids();

    let pages = all_shadow_pages(&shadow_root);
    assert!(
        !pages.is_empty(),
        "no shadow pages found under {} — a silent zero is not a pass",
        shadow_root.display()
    );

    // Gating reference tier: the committed real-GRCh38 slice. `Normalizes to` is asserted
    // against it on every PR — no manifest, no network, hermetic.
    let slice_normalizer = Normalizer::with_config(
        MockProvider::from_json(&root.join(SLICE))
            .unwrap_or_else(|e| panic!("load committed slice {SLICE}: {e}")),
        NormalizeConfig::default(),
    );
    // Cross-check tier (nightly): the full live reference. When FERRO_MANIFEST is present the
    // same assertions run against it too, catching drift between the slice and current hg38.
    // Under FERRO_REQUIRE_MANIFEST its absence is a failure, not a silent skip.
    let manifest_normalizer = std::env::var("FERRO_MANIFEST").ok().map(|p| {
        Normalizer::with_config(
            MultiFastaProvider::from_manifest(&p)
                .unwrap_or_else(|e| panic!("load FERRO_MANIFEST {p}: {e}")),
            NormalizeConfig::default(),
        )
    });
    let require_manifest = std::env::var("FERRO_REQUIRE_MANIFEST").is_ok();
    let bless = std::env::var("BLESS_SHADOW_CORPUS").is_ok();
    let summaries = ledger_summaries();

    let mut errors: Vec<String> = Vec::new();
    let mut examples_run = 0usize;
    let mut anchors_total = 0usize;
    let mut norm_checked = 0usize;
    let mut norm_crosschecked = 0usize;
    // The single shadow-spec example corpus: one JSONL row per example, in page order. It is
    // the ONLY thing built from these pages and is committed + blessed (BLESS_SHADOW_CORPUS=1).
    let mut corpus: Vec<String> = Vec::new();

    // Biopython-style alignments for each executable example, computed off the committed slice
    // and served to the page's JS via a blessed sidecar under docs/. Prototype: currently only
    // the RNA/inversion page is generated, until the visual format is approved.
    let align_provider = MockProvider::from_json(&root.join(SLICE)).ok();
    let mut alignments: BTreeMap<String, serde_json::Value> = BTreeMap::new();

    for page in &pages {
        let rel = page.strip_prefix(&shadow_root).unwrap();
        let spec_dir_rel = rel
            .parent()
            .map(|p| p.to_string_lossy().replace('\\', "/"))
            .unwrap_or_default();
        let raw = read(page);
        // Transclude the "Why" blocks from the ledger summaries. In bless mode the file is
        // rewritten in place; otherwise a stale or missing block is reported.
        let rel_str = rel.display().to_string();
        let text = transclude_why_blocks(&raw, &rel_str, &summaries, bless, &mut errors);
        if bless && text != raw {
            std::fs::write(page, &text).expect("write transcluded page");
        }
        let anchors = parse_page(&text, &spec_dir_rel);
        // A page with no `path:line` anchor headings is a prose landing page (e.g. the
        // section `index.md`), not a shadow of a spec page — nothing to check.
        if anchors.is_empty() {
            continue;
        }

        for a in &anchors {
            anchors_total += 1;
            let where_ = format!("{}:{}", rel.display(), a.heading_lineno);
            let spec_abs = spec_root.join(&a.spec_rel);

            // 1. quote
            if a.quote.trim().is_empty() {
                errors.push(format!("{where_}: anchor has no spec block-quote"));
            } else {
                match spec_span(&spec_abs, a.lstart, a.lend) {
                    None => errors.push(format!(
                        "{where_}: cannot read spec span {}:{}-{}",
                        a.spec_rel, a.lstart, a.lend
                    )),
                    Some(span) => {
                        if !span.contains(&collapse(&a.quote)) {
                            errors.push(format!(
                                "{where_}: quote not found at {}:{}-{}\n    quote: {}\n    spec : {}",
                                a.spec_rel,
                                a.lstart,
                                a.lend,
                                collapse(&a.quote),
                                span
                            ));
                        }
                    }
                }
            }

            // 3a. ruling links resolve in the ledger
            for id in &a.ruling_ids {
                if !ledger.contains_key(id) {
                    errors.push(format!("{where_}: ruling link '{id}' not in ledger"));
                }
            }
            // 3b. cross-refs resolve in the spec (basename walked up the spec tree)
            let anchor_dir = Path::new(&a.spec_rel)
                .parent()
                .map(|p| p.to_string_lossy().replace('\\', "/"))
                .unwrap_or_default();
            for (name, xa, xb) in &a.crossrefs {
                let resolved = resolve_crossref(&spec_root, &anchor_dir, name)
                    .filter(|rel| spec_span(&spec_root.join(rel), *xa, *xb).is_some());
                if resolved.is_none() {
                    errors.push(format!(
                        "{where_}: cross-ref {name}:{xa}-{xb} does not resolve"
                    ));
                }
            }

            // 2. execute examples
            for ex in &a.examples {
                let loc = format!("{}:{}", rel.display(), ex.row);
                let parses = parse_hgvs(&ex.spelling).is_ok();
                let strict_refused =
                    parse_hgvs_with_config(&ex.spelling, ErrorConfig::strict()).is_err();
                examples_run += 1;
                // The verdict describes ferro's OUTPUT against the recommendations. `recommended`
                // and `conformant` both require the input to parse; their difference (is the output
                // the recommended form, or a valid-but-not-yet-recommended one?) is the author's
                // reading, grounded by the spec quote and the executed `Normalizes to` assertion
                // below. `refused` and `bug` are mechanically checked here.
                match ex.verdict.as_str() {
                    "recommended" | "conformant" => {
                        if !parses {
                            errors.push(format!(
                                "{loc}: '{}' marked {} but parse_hgvs refused it",
                                ex.spelling, ex.verdict
                            ));
                        }
                    }
                    // The input is not valid HGVS; ferro correctly rejects it in strict mode.
                    "refused" => {
                        if !strict_refused {
                            errors.push(format!(
                                "{loc}: '{}' marked refused but strict mode accepted it",
                                ex.spelling
                            ));
                        }
                    }
                    // ferro accepts the input (parses) but its output is not valid HGVS — a defect.
                    // The `Normalizes to` cell pins the invalid output, and the assertion below
                    // verifies it does not re-parse; fixing the bug flips this row.
                    "bug" => {
                        if !parses {
                            errors.push(format!(
                                "{loc}: '{}' marked bug but parse_hgvs refused the INPUT (a refused input is not a bug)",
                                ex.spelling
                            ));
                        }
                        match &ex.normalizes_to {
                            Norm::To(out) if parse_hgvs(out).is_ok() => errors.push(format!(
                                "{loc}: '{}' marked bug but its pinned output '{out}' is valid HGVS — not a bug",
                                ex.spelling
                            )),
                            Norm::Unset | Norm::Fixed => errors.push(format!(
                                "{loc}: '{}' marked bug but names no invalid output in `Normalizes to`",
                                ex.spelling
                            )),
                            _ => {}
                        }
                    }
                    other => errors.push(format!(
                        "{loc}: '{}' has unknown verdict '{other}' (expected recommended/conformant/refused/bug)",
                        ex.spelling
                    )),
                }

                // reference tier: assert normalize(input) == the page's `Normalizes to`
                let expected = match &ex.normalizes_to {
                    Norm::Unset => None,
                    Norm::Fixed => Some(ex.spelling.clone()),
                    Norm::To(s) => Some(s.clone()),
                };
                if let Some(expected) = &expected {
                    // gating: the committed slice
                    norm_checked += 1;
                    match normalize_str(&slice_normalizer, &ex.spelling) {
                        Some(g) if &g == expected => {}
                        Some(g) => errors.push(format!(
                            "{loc}: '{}' normalizes to '{g}', page says '{expected}' (committed slice)",
                            ex.spelling
                        )),
                        None => errors.push(format!(
                            "{loc}: '{}' could not be normalized against the committed slice",
                            ex.spelling
                        )),
                    }
                    // cross-check: the full live reference, when present
                    if let Some(n) = &manifest_normalizer {
                        norm_crosschecked += 1;
                        match normalize_str(n, &ex.spelling) {
                            Some(g) if &g == expected => {}
                            Some(g) => errors.push(format!(
                                "{loc}: '{}' normalizes to '{g}' on the full reference, page says '{expected}' — slice drift from current hg38",
                                ex.spelling
                            )),
                            None => errors.push(format!(
                                "{loc}: '{}' could not be normalized against the full reference",
                                ex.spelling
                            )),
                        }
                    }
                }

                // Prototype: a visual alignment for each executable row on the inversion page.
                if expected.is_some() && rel.ends_with("RNA/inversion.md") {
                    if let Some(p) = &align_provider {
                        if let Some(al) = build_alignment(&ex.spelling, p) {
                            alignments.entry(ex.spelling.clone()).or_insert(al);
                        }
                    }
                }

                // one corpus row per example (JSON, stable field order)
                let norm_json = match &expected {
                    Some(e) => format!("{:?}", e),
                    None => "null".to_string(),
                };
                corpus.push(format!(
                    "{{\"page\":{:?},\"anchor\":{:?},\"input\":{:?},\"verdict\":{:?},\"normalizes_to\":{}}}",
                    rel.display().to_string(),
                    format!("{}:{}-{}", a.spec_rel, a.lstart, a.lend),
                    ex.spelling,
                    ex.verdict,
                    norm_json,
                ));
            }

            // 4. coverage (reported, not gated)
        }

        // coverage line per page
        let spec_abs = spec_root.join(anchors[0].spec_rel.clone());
        let spec_lines = std::fs::read_to_string(&spec_abs)
            .map(|t| t.lines().count())
            .unwrap_or(0);
        let covered: usize = anchors.iter().map(|a| a.lend - a.lstart + 1).sum();
        eprintln!(
            "shadow-spec coverage: {} — {} anchors, ~{}/{} spec lines annotated",
            rel.display(),
            anchors.len(),
            covered,
            spec_lines
        );
    }

    assert!(
        examples_run > 0,
        "no examples executed — a silent zero is not a pass"
    );
    eprintln!(
        "shadow-spec: {} pages, {anchors_total} anchors, {examples_run} examples executed",
        pages.len()
    );
    eprintln!(
        "shadow-spec: normalization checked={norm_checked} (committed slice), cross-checked={norm_crosschecked} (full reference)"
    );
    if require_manifest && manifest_normalizer.is_none() {
        errors.push(
            "FERRO_REQUIRE_MANIFEST is set but FERRO_MANIFEST is absent — the cross-check tier could not run".into(),
        );
    }

    // Committed + blessed corpus: BLESS_SHADOW_CORPUS=1 rewrites it, otherwise it must match.
    let corpus_path = root.join(CORPUS);
    let rendered = if corpus.is_empty() {
        String::new()
    } else {
        format!("{}\n", corpus.join("\n"))
    };
    bless_or_check(&corpus_path, CORPUS, &rendered, bless, &mut errors);

    // Blessed alignment sidecar, served to the page JS by mdBook. Gated by the same flag as the
    // corpus so one bless command refreshes both. Kept in docs/src so mdBook copies it into the
    // book output (it carries no `path:line` heading, so parse_page never treats it as a page).
    // Only run when the slice provider loaded: otherwise `alignments` is empty for a reason
    // unrelated to the sidecar's contents, and comparing empty-vs-committed would falsely report
    // it stale. (In practice the slice always loads — `slice_normalizer` panics above if it does
    // not — but do not couple the sidecar's verdict to that.)
    if align_provider.is_some() {
        let align_path = root.join(ALIGNMENTS);
        let align_rendered =
            serde_json::to_string_pretty(&alignments).expect("serialize alignments") + "\n";
        bless_or_check(&align_path, ALIGNMENTS, &align_rendered, bless, &mut errors);
    }

    assert!(
        errors.is_empty(),
        "shadow-spec check failed ({} problem(s)):\n{}",
        errors.len(),
        errors.join("\n")
    );
}

/// Grounds [`render_alignment`]'s presentation logic against hand-computed expectations —
/// independent of the bless-vs-disk self-consistency the page gate relies on, and exercising
/// the substitution / deletion / insertion / delins / windowed branches (the single-page
/// prototype only ever feeds it inversions and one insertion).
mod alignment_render_tests {
    use super::render_alignment;

    fn cv(s: &str) -> Vec<char> {
        s.chars().collect()
    }
    fn s<'a>(v: &'a serde_json::Value, k: &str) -> &'a str {
        v[k].as_str().unwrap()
    }

    #[test]
    fn substitution_marks_the_single_changed_column() {
        let v = render_alignment(&cv("A"), &cv("G"), "ACC", "GGT");
        assert_eq!(s(&v, "ref"), "ACCAGGT");
        assert_eq!(s(&v, "alt"), "ACCGGGT");
        assert_eq!(s(&v, "connector"), "||| |||");
        assert_eq!(v["editStart"], 3);
        assert_eq!(v["editEnd"], 4);
        assert_eq!(s(&v, "note"), "1 nt · 1 changed");
    }

    #[test]
    fn deletion_gaps_the_alternate_row_with_a_blank_connector() {
        let v = render_alignment(&cv("AT"), &cv(""), "CCC", "GGG");
        assert_eq!(s(&v, "ref"), "CCCATGGG");
        assert_eq!(s(&v, "alt"), "CCC--GGG");
        assert_eq!(s(&v, "connector"), "|||  |||");
        assert_eq!(s(&v, "note"), "deletion · 2 nt");
    }

    #[test]
    fn insertion_gaps_the_reference_row() {
        let v = render_alignment(&cv(""), &cv("GC"), "CCC", "GGG");
        assert_eq!(s(&v, "ref"), "CCC--GGG");
        assert_eq!(s(&v, "alt"), "CCCGCGGG");
        assert_eq!(s(&v, "note"), "insertion · 2 nt");
    }

    #[test]
    fn unequal_delins_pads_the_shorter_allele_and_blanks_the_connector() {
        let v = render_alignment(&cv("A"), &cv("GCT"), "CC", "GG");
        assert_eq!(s(&v, "ref"), "CCA--GG");
        assert_eq!(s(&v, "alt"), "CCGCTGG");
        assert_eq!(s(&v, "connector"), "||   ||");
        assert_eq!(s(&v, "note"), "delins · 1 → 3 nt");
    }

    #[test]
    fn a_wide_span_is_windowed_and_the_note_reports_the_full_size() {
        let del: String = std::iter::repeat_n('A', 50).collect();
        let mut ins = del.clone();
        ins.replace_range(0..1, "C"); // 3 changes, one of them
        ins.replace_range(25..26, "C"); // elided in the windowed middle
        ins.replace_range(49..50, "C");
        let v = render_alignment(&cv(&del), &cv(&ins), "TT", "GG");
        assert!(
            s(&v, "ref").contains(" … "),
            "wide span should be windowed: {}",
            s(&v, "ref")
        );
        assert_eq!(s(&v, "note"), "50 nt · 3 changed"); // counts the elided change too
        assert_eq!(v["editStart"], 2); // 2-base left flank
        assert_eq!(v["editEnd"], 2 + 15 + 3 + 15); // head + " … " + tail
    }
}
