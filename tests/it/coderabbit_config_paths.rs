//! Liveness check for `.coderabbit.yaml`.
//!
//! The review config targets instructions at modules and names specific test
//! files. Both rot silently: a renamed directory turns a `path:` glob into a
//! pattern that matches nothing, and a moved test turns a "confirm the test
//! still covers it" pointer into a dangling reference. Neither failure is
//! visible — CodeRabbit does not report an instruction that never fired, so the
//! config degrades into misdirection while still looking maintained. (This
//! happened: the instructions pointed at `tests/<name>.rs` for a year after the
//! integration suite moved to `tests/it/`.)
//!
//! This test makes that failure loud:
//!
//!   1. every `path:` glob under `reviews.path_instructions` matches at least
//!      one file in the repository,
//!   2. every repo-relative path named in backticks anywhere in the config
//!      exists on disk, and
//!   3. no backticked path is split across lines — the `>-` folded scalars fold
//!      a line break into a space, so a wrapped path reaches CodeRabbit as
//!      `tests/it/ foo.rs`, which names nothing, and
//!   4. every `path:` glob uses only syntax the matcher here implements, so a
//!      verdict is never issued about a pattern that was silently read as a
//!      literal.
//!
//! The file universe is the git index, not the working tree — see
//! [`repository_files`] — so the verdict is the same on a developer host with
//! build artifacts lying around as on a clean checkout.
//!
//! `reviews.path_filters` are deliberately not checked: they are exclusions, so
//! one that matches nothing is inert rather than misleading, and several of them
//! point into the `assets/` submodule, which is empty in a checkout that has not
//! initialized it.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

/// Directory names skipped wherever they appear — build output and tooling
/// caches, never a legitimate target for a review instruction.
const SKIPPED_DIR_NAMES: &[&str] = &[".git", ".venv", "node_modules", "target"];

/// Directories skipped only at the repository root: the vendored spec submodule
/// (empty unless initialized) and the top-level reference-data and vendor trees.
/// Matched by root-relative path, not by name, so `src/data/` is still walked.
const SKIPPED_ROOT_DIRS: &[&str] = &["assets", "data", "vendor"];

/// Paths that may be named in the config while legitimately absent from a
/// checkout: both are gitignored build artifacts that the config exists to keep
/// *out* of diffs, so requiring them to exist would invert the intent.
const ALLOWED_ABSENT: &[&str] = &[
    "tests/fixtures/grammar/hgvs_spec_normalization.json",
    "tests/fixtures/grammar/hgvs_spec_enumeration.json",
];

/// File extensions worth treating as repo-relative path references when they
/// appear in backticks. Anything else (prose, command fragments) is ignored.
const PATH_EXTENSIONS: &[&str] = &[
    ".rs", ".json", ".pyi", ".py", ".txt", ".md", ".yml", ".yaml", ".toml",
];

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn config_text() -> String {
    let path = repo_root().join(".coderabbit.yaml");
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()))
}

/// Every repository file, as a path relative to the repo root with `/`
/// separators.
///
/// Tracked files are the authority: they are the only files CodeRabbit can ever
/// see in a diff, and they are the same on every checkout. The directory walk
/// below also counts untracked and gitignored files — the two generated spec
/// fixtures under `tests/fixtures/grammar/` among them — so a glob whose only
/// match is a local build artifact would look live on a developer host and be
/// dead on a clean one, which is the determinism this test exists to provide.
/// The walk is kept as the fallback for a checkout with no usable `git`.
fn repository_files() -> Vec<String> {
    let root = repo_root();
    if let Some(tracked) = tracked_files(&root) {
        return tracked;
    }
    let mut files = Vec::new();
    collect_files(&root, &root, &mut files);
    files.sort();
    files
}

/// Files known to the git index, `/`-separated and sorted. `None` when `git` is
/// unavailable or the directory is not a work tree.
fn tracked_files(root: &Path) -> Option<Vec<String>> {
    let output = std::process::Command::new("git")
        .arg("-C")
        .arg(root)
        .args(["ls-files", "-z"])
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    // `-z` emits raw, `/`-separated pathnames with no quoting or escaping, so
    // they need no separator rewrite — unlike the `collect_files` walk.
    let mut files: Vec<String> = String::from_utf8(output.stdout)
        .ok()?
        .split('\0')
        .filter(|entry| !entry.is_empty())
        .map(str::to_string)
        .collect();
    if files.is_empty() {
        return None;
    }
    files.sort();
    Some(files)
}

fn collect_files(root: &Path, dir: &Path, out: &mut Vec<String>) {
    let entries = match std::fs::read_dir(dir) {
        Ok(entries) => entries,
        Err(_) => return,
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let name = entry.file_name().to_string_lossy().into_owned();
        let Ok(relative) = path.strip_prefix(root) else {
            continue;
        };
        let relative = relative.to_string_lossy().replace('\\', "/");
        if path.is_dir() {
            let skipped = SKIPPED_DIR_NAMES.contains(&name.as_str())
                || SKIPPED_ROOT_DIRS.contains(&relative.as_str());
            if !skipped {
                collect_files(root, &path, out);
            }
        } else {
            out.push(relative);
        }
    }
}

/// The `path:` value of every entry under `reviews.path_instructions`.
fn instruction_globs(text: &str) -> Vec<String> {
    let config: serde_yaml::Value =
        serde_yaml::from_str(text).expect(".coderabbit.yaml is not valid YAML");
    let instructions = config
        .get("reviews")
        .and_then(|reviews| reviews.get("path_instructions"))
        .and_then(|value| value.as_sequence())
        .expect("reviews.path_instructions is missing or not a list");
    instructions
        .iter()
        .map(|entry| {
            entry
                .get("path")
                .and_then(|path| path.as_str())
                .expect("a path_instructions entry has no `path`")
                .to_string()
        })
        .collect()
}

/// The `knowledge_base.code_guidelines.filePatterns` entries, which name the
/// files CodeRabbit loads as project guidelines.
///
/// Covered for the same reason as `path_instructions`: these are plain paths,
/// and renaming one silently stops the guideline being loaded at all. That is a
/// quieter failure than a dead instruction glob, since nothing downstream ever
/// mentions the file again.
fn code_guideline_patterns(text: &str) -> Vec<String> {
    let config: serde_yaml::Value =
        serde_yaml::from_str(text).expect(".coderabbit.yaml is not valid YAML");
    let Some(patterns) = config
        .get("knowledge_base")
        .and_then(|kb| kb.get("code_guidelines"))
        .and_then(|guidelines| guidelines.get("filePatterns"))
    else {
        return Vec::new();
    };
    patterns
        .as_sequence()
        .expect("knowledge_base.code_guidelines.filePatterns is not a list")
        .iter()
        .map(|entry| {
            entry
                .as_str()
                .expect("a filePatterns entry is not a string")
                .to_string()
        })
        .collect()
}

/// Rejoins a span that a `>-` folded scalar wrapped across lines, by dropping
/// each line break together with the indentation that follows it. Spaces within
/// a line are preserved, so prose is still recognisable as prose.
fn unfold(span: &str) -> String {
    let mut out = String::with_capacity(span.len());
    let mut folding = false;
    for ch in span.chars() {
        match ch {
            '\n' | '\r' => folding = true,
            ' ' | '\t' if folding => {}
            _ => {
                folding = false;
                out.push(ch);
            }
        }
    }
    out
}

/// Whether a backticked span names a repo-relative path. Bare file names
/// (`merge.rs`) are skipped — they are module references, not paths — as are
/// patterns containing wildcards, which are covered by the glob check, and
/// `<name>`-style placeholders, which stand for a family of files.
fn is_path_reference(candidate: &str) -> bool {
    candidate.contains('/')
        && !candidate.contains('*')
        && !candidate.contains('<')
        && !candidate.contains(' ')
        && PATH_EXTENSIONS.iter().any(|ext| candidate.ends_with(ext))
}

/// The backtick-delimited spans of `text`, in order.
///
/// `split('`')` alternates prose, code, prose, code — so the code spans are the
/// odd indices, but only while the backticks are balanced. An odd count flips
/// that parity and every consumer then reads prose as code and code as prose,
/// finding nothing and reporting success. That silent pass is the exact failure
/// mode this file exists to remove, so the parity is asserted rather than
/// assumed.
fn code_spans(text: &str) -> Vec<&str> {
    assert!(
        text.matches('`').count().is_multiple_of(2),
        "unbalanced backtick: the span scan would read every code span as prose and vice versa, \
         so both backtick checks would pass while checking nothing"
    );
    text.split('`').skip(1).step_by(2).collect()
}

/// Backtick-quoted repo-relative paths mentioned anywhere in the config.
fn named_paths(text: &str) -> BTreeSet<String> {
    let mut paths = BTreeSet::new();
    for raw in code_spans(text) {
        // Classify the rejoined span. A folded scalar can wrap a path across
        // lines, and the raw span then carries a newline plus indentation —
        // which the "no interior space" rule rejects, silently skipping exactly
        // the dangling pointer this test exists to catch.
        let candidate = unfold(raw);
        if is_path_reference(&candidate) && !ALLOWED_ABSENT.contains(&candidate.as_str()) {
            paths.insert(candidate);
        }
    }
    paths
}

/// Expands `{a,b}` alternations into one pattern per alternative, recursively so
/// nested groups are handled.
fn expand_braces(pattern: &str) -> Vec<String> {
    let Some(open) = pattern.find('{') else {
        return vec![pattern.to_string()];
    };
    let mut depth = 0usize;
    let mut close = None;
    let mut alternatives = Vec::new();
    let mut current = String::new();
    for (index, ch) in pattern[open..].char_indices() {
        let absolute = open + index;
        match ch {
            '{' => {
                depth += 1;
                if depth > 1 {
                    current.push(ch);
                }
            }
            '}' => {
                depth -= 1;
                if depth == 0 {
                    close = Some(absolute);
                    alternatives.push(std::mem::take(&mut current));
                    break;
                }
                current.push(ch);
            }
            ',' if depth == 1 => alternatives.push(std::mem::take(&mut current)),
            _ => current.push(ch),
        }
    }
    let close = close.unwrap_or_else(|| panic!("unbalanced `{{` in glob: {pattern}"));
    let (prefix, suffix) = (&pattern[..open], &pattern[close + 1..]);
    alternatives
        .iter()
        .flat_map(|alternative| expand_braces(&format!("{prefix}{alternative}{suffix}")))
        .collect()
}

/// Minimatch metacharacters this matcher does not implement. Treating one as a
/// literal would silently misclassify the glob — a `path:` instruction using it
/// could be reported live while matching nothing in the real review, or the
/// reverse. Reject them instead, so adding one to the config fails loudly here
/// and names what has to be implemented.
const UNSUPPORTED_GLOB_CHARS: &[char] = &['?', '[', ']', '(', ')', '!'];

/// The first unsupported metacharacter in `pattern`, if any. `*`, `**` and
/// `{a,b}` are the supported forms.
fn unsupported_glob_syntax(pattern: &str) -> Option<char> {
    pattern
        .chars()
        .find(|ch| UNSUPPORTED_GLOB_CHARS.contains(ch))
}

/// Matches a brace-free glob against a `/`-separated path. `**` matches any
/// number of whole segments (including none); `*` matches within one segment;
/// and neither crosses into a dot-leading component unless the pattern spells
/// the dot out, matching the minimatch semantics CodeRabbit applies.
///
/// Panics on syntax this matcher does not implement — see
/// [`unsupported_glob_syntax`].
fn glob_matches(pattern: &str, path: &str) -> bool {
    assert!(
        unsupported_glob_syntax(pattern).is_none(),
        "glob `{pattern}` uses syntax this matcher does not implement; it would be matched as a \
         literal and the liveness verdict would be wrong"
    );
    let pattern_segments: Vec<&str> = pattern.split('/').collect();
    let path_segments: Vec<&str> = path.split('/').collect();
    match_segments(&pattern_segments, &path_segments)
}

fn match_segments(pattern: &[&str], path: &[&str]) -> bool {
    match pattern.first() {
        None => path.is_empty(),
        // `**` spans whole segments, but like minimatch it does not descend
        // through a dot-leading directory — see `segment_matches`.
        Some(&"**") => (0..=path.len())
            .take_while(|&skip| skip == 0 || !path[skip - 1].starts_with('.'))
            .any(|skip| match_segments(&pattern[1..], &path[skip..])),
        Some(segment) => match path.first() {
            Some(candidate) if segment_matches(segment, candidate) => {
                match_segments(&pattern[1..], &path[1..])
            }
            _ => false,
        },
    }
}

/// Wildcard match within a single path segment (`*` only; `/` never matches).
///
/// CodeRabbit evaluates these globs with minimatch, which matches a dot-leading
/// component only when the pattern component starts with a literal dot — a
/// leading wildcard never does. Without that rule a glob whose only match is a
/// hidden file looks live here while firing on nothing in the actual review.
fn segment_matches(pattern: &str, segment: &str) -> bool {
    if segment.starts_with('.') && !pattern.starts_with('.') {
        return false;
    }
    let parts: Vec<&str> = pattern.split('*').collect();
    if parts.len() == 1 {
        return pattern == segment;
    }
    let mut rest = segment;
    if let Some(first) = parts.first() {
        match rest.strip_prefix(first) {
            Some(stripped) => rest = stripped,
            None => return false,
        }
    }
    if let Some(last) = parts.last() {
        match rest.strip_suffix(last) {
            // `rest` shrank from both ends; `last` may overlap `first` on a
            // short segment, which `strip_suffix` on the already-stripped rest
            // correctly rejects.
            Some(stripped) => rest = stripped,
            None => return false,
        }
    }
    for middle in parts.iter().take(parts.len() - 1).skip(1) {
        match rest.find(middle) {
            Some(index) => rest = &rest[index + middle.len()..],
            None => return false,
        }
    }
    true
}

/// Brace alternatives of `glob` that select no file. Checked per alternative
/// rather than per glob: `src/{batch,parallel.rs}` "matches something" as a
/// whole even when `src/batch` is dead, so a whole-glob check would miss
/// exactly the case this test exists to catch.
fn dead_alternatives(glob: &str, files: &[String]) -> Vec<String> {
    expand_braces(glob)
        .into_iter()
        .filter(|pattern| !files.iter().any(|file| glob_matches(pattern, file)))
        .collect()
}

/// Every `path:` glob — and every brace alternative within it — must select at
/// least one file. One that matches nothing is an instruction that can never
/// fire.
#[test]
fn every_path_instruction_matches_a_file() {
    let files = repository_files();
    let orphans: Vec<String> = instruction_globs(&config_text())
        .iter()
        .flat_map(|glob| dead_alternatives(glob, &files))
        .collect();
    assert!(
        orphans.is_empty(),
        ".coderabbit.yaml has path_instructions globs that match no file — the instruction \
         can never fire. Repoint or remove: {orphans:#?}"
    );
}

/// The same liveness question for `knowledge_base.code_guidelines.filePatterns`.
///
/// A dead entry here is quieter than a dead `path_instructions` glob: the
/// guideline simply never loads, and nothing in a review says so.
#[test]
fn every_code_guideline_pattern_matches_a_file() {
    let files = repository_files();
    let patterns = code_guideline_patterns(&config_text());
    assert!(
        !patterns.is_empty(),
        "knowledge_base.code_guidelines.filePatterns is missing or empty — if the key was \
         renamed or removed, this test is no longer checking anything"
    );
    let orphans: Vec<String> = patterns
        .iter()
        .flat_map(|pattern| dead_alternatives(pattern, &files))
        .collect();
    assert!(
        orphans.is_empty(),
        ".coderabbit.yaml names code-guideline files that do not exist, so those guidelines \
         are never loaded. Repoint or remove: {orphans:#?}"
    );
}

/// The matcher implements `*`, `**` and `{a,b}`. A `path:` using anything else
/// would be matched as a literal, so the liveness verdict above would be about
/// a pattern CodeRabbit never evaluates.
#[test]
fn every_path_instruction_uses_supported_glob_syntax() {
    let unsupported: Vec<String> = instruction_globs(&config_text())
        .iter()
        .flat_map(|glob| expand_braces(glob))
        .filter_map(|pattern| {
            unsupported_glob_syntax(&pattern).map(|ch| format!("{pattern} (uses `{ch}`)"))
        })
        .collect();
    assert!(
        unsupported.is_empty(),
        ".coderabbit.yaml uses glob syntax this test's matcher does not implement, so the \
         liveness check cannot speak for these patterns. Implement them in `segment_matches` or \
         rewrite the globs: {unsupported:#?}"
    );
}

#[test]
fn unsupported_glob_syntax_is_detected() {
    assert_eq!(unsupported_glob_syntax("src/**/*.rs"), None);
    assert_eq!(unsupported_glob_syntax("tests/it/mutalyzer_*.rs"), None);
    assert_eq!(unsupported_glob_syntax("src/lib.?s"), Some('?'));
    assert_eq!(unsupported_glob_syntax("src/[ab]/*.rs"), Some('['));
    assert_eq!(unsupported_glob_syntax("src/!(batch)/*.rs"), Some('!'));
}

#[test]
#[should_panic(expected = "does not implement")]
fn glob_matches_rejects_unsupported_syntax() {
    glob_matches("src/[ab].rs", "src/a.rs");
}

#[test]
#[should_panic(expected = "unbalanced backtick")]
fn code_spans_rejects_an_unbalanced_backtick() {
    code_spans("a `tests/it/foo.rs` and a stray `");
}

#[test]
fn dead_alternatives_flags_a_bare_directory_alternative() {
    let files = vec![
        "src/batch/mod.rs".to_string(),
        "src/parallel.rs".to_string(),
    ];
    assert_eq!(
        dead_alternatives("src/{batch,parallel.rs}", &files),
        vec!["src/batch".to_string()]
    );
    assert!(dead_alternatives("src/{batch/**/*.rs,parallel.rs}", &files).is_empty());
}

/// Every backticked repo-relative path must exist. These are the "confirm the
/// test still covers it" pointers; a dangling one misdirects the reviewer.
#[test]
fn every_named_path_exists() {
    let root = repo_root();
    let missing: Vec<String> = named_paths(&config_text())
        .into_iter()
        .filter(|path| !root.join(path).exists())
        .collect();
    assert!(
        missing.is_empty(),
        ".coderabbit.yaml names paths that do not exist — a dangling pointer silently \
         misdirects the reviewer. Repoint or remove: {missing:#?}"
    );
}

/// A backticked path must stay on one line. The `path_instructions` values are
/// `>-` folded scalars, and YAML folds a line break into a space — so a wrapped
/// path reaches CodeRabbit as `tests/it/ foo.rs`, a pointer to nothing, even
/// though the file it meant to name exists.
#[test]
fn no_backticked_path_is_split_across_lines() {
    let text = config_text();
    let split: Vec<String> = code_spans(&text)
        .into_iter()
        .filter(|span| span.contains('\n'))
        .map(unfold)
        .filter(|candidate| is_path_reference(candidate))
        .collect();
    assert!(
        split.is_empty(),
        ".coderabbit.yaml wraps a backticked path across lines. The folded scalar turns the \
         break into a space, so CodeRabbit reads a path that names nothing. Rewrap the line so \
         each of these stays intact: {split:#?}"
    );
}

#[test]
fn glob_matcher_handles_globstar_and_wildcards() {
    assert!(glob_matches("src/**/*.rs", "src/normalize/shuffle.rs"));
    assert!(glob_matches("src/**/*.rs", "src/lib.rs"));
    assert!(glob_matches(
        "**/*.proptest-regressions",
        "tests/x.proptest-regressions"
    ));
    assert!(glob_matches(
        "tests/it/mutalyzer_*.rs",
        "tests/it/mutalyzer_tests.rs"
    ));
    assert!(!glob_matches("src/**/*.rs", "src/lib.py"));
    assert!(!glob_matches("tests/it/*.rs", "tests/it/common/mod.rs"));
}

/// A directory named without a trailing `/**` selects only the directory entry
/// itself, never the files inside it — the failure mode that left `src/batch`
/// uncovered.
#[test]
fn glob_matcher_does_not_match_files_under_a_bare_directory() {
    assert!(!glob_matches("src/batch", "src/batch/mod.rs"));
    assert!(glob_matches("src/batch/**/*.rs", "src/batch/mod.rs"));
}

#[test]
fn brace_expansion_covers_every_alternative() {
    let expanded = expand_braces("src/{a,b/**/*.rs,c.rs}");
    assert_eq!(expanded, vec!["src/a", "src/b/**/*.rs", "src/c.rs"]);
}

#[test]
fn named_paths_ignores_bare_file_names_and_globs() {
    let text = "`merge.rs` `tests/it/projection.rs` `tests/it/*.rs` `cargo run --bin x`";
    let paths = named_paths(text);
    assert_eq!(
        paths.into_iter().collect::<Vec<_>>(),
        vec!["tests/it/projection.rs".to_string()]
    );
}

/// A span the folded scalar wrapped is still a path reference, not prose: the
/// break plus indentation must not be mistaken for an interior space.
#[test]
fn named_paths_rejoins_a_wrapped_backtick_span() {
    let text = "confirm `tests/it/\n        idempotency_tests.rs` covers it";
    assert_eq!(
        named_paths(text).into_iter().collect::<Vec<_>>(),
        vec!["tests/it/idempotency_tests.rs".to_string()]
    );
    // Rejoining must not swallow a genuine intra-line space.
    assert!(named_paths("`cargo run\n        --bin x`").is_empty());
}

/// Minimatch — CodeRabbit's matcher — only lets a literal dot match a
/// dot-leading component, so a glob that "matches" only hidden files is dead in
/// the review even though a naive matcher calls it live.
#[test]
fn glob_matcher_uses_minimatch_dotfile_semantics() {
    assert!(!glob_matches("src/*.rs", "src/.hidden.rs"));
    assert!(!glob_matches("**/*.yml", ".github/workflows/ci.yml"));
    assert!(!glob_matches("*/workflows/**", ".github/workflows/ci.yml"));
    // An explicit literal dot still matches, which is why the real
    // `.github/workflows/**` instruction stays live.
    assert!(glob_matches(
        ".github/workflows/**",
        ".github/workflows/ci.yml"
    ));
}
