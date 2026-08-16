//! Every `ferro-benchmark …` invocation printed to a user or written in a doc
//! comment must name a subcommand — and spell flags — that the real clap command
//! tree still accepts.
//!
//! # Why this guard exists
//!
//! The benchmark CLI was reorganised from flat, hyphenated subcommands
//! (`setup-uta`, `extract-clinvar`, `compare-normalize`, …) into grouped ones
//! (`setup uta`, `extract clinvar`, `benchmark normalize`, …). The old spellings
//! then survived as string literals — most damagingly in *runtime* "next step"
//! hints such as the line `setup uta` prints on success, which a contributor
//! copy-pastes and watches fail with `error: unrecognized subcommand`. This is
//! the second time these strings have gone stale (#1893, then #2065), so the
//! class is made self-policing here rather than re-swept by hand.
//!
//! # What is checked
//!
//! The set is **derived from source**, not maintained as a list: every
//! occurrence of the literal `ferro-benchmark` in the benchmark sources, docs
//! and `pixi.toml` is scanned, and two things are lifted out of the text that
//! follows it — the subcommand path, and the flags spelled after it on the same
//! line. Both are resolved against `Cli::command()`, the same clap tree the
//! binary parses with:
//!
//! - a **subcommand path** fails when, at some point, the current command has no
//!   matching subcommand *and* takes no positional argument that could consume
//!   the token (so the token could only have been a subcommand);
//! - a **flag** fails when the command the path resolved to declares no argument
//!   with that long name or short name.
//!
//! The flag half is not decoration. The #2065 sweep corrected nine subcommand
//! spellings and left `setup seqrepo --seqrepo-dir <path>` behind on a line it
//! had just rewritten: the subcommand was now right and the flag was dead, so
//! `check`'s own printed remedy still died with `error: unexpected argument`.
//! `--seqrepo-dir` is a real flag on `prepare`, which is exactly why reading the
//! line does not catch it — only resolving it does.
//!
//! # What it does not reach
//!
//! It is a **token scan**, a floor and not a proof. Specifically:
//!
//! - It reads one line at a time, so flags written on a shell **continuation
//!   line** (`… prepare mutalyzer \` + `    --ferro-reference …`) are invisible;
//!   most of `docs/BENCHMARK_GUIDE.md` is written that way.
//! - It recognises a flag only in its bare `--long-name` / `-s` spelling. A
//!   `--flag=value` spelling, a flag assembled at runtime, or one reached
//!   through a variable is not seen.
//! - It checks that a flag **exists**, never that its *value* makes sense; and
//!   it does not check that required flags are present.
//! - When a `{}` placeholder truncates the path (`prepare {} --seqrepo-dir …`),
//!   flags are resolved against the command reached so far. That is exact for
//!   this CLI, because `prepare`, `check`, `parse` and `normalize` are leaf
//!   commands carrying a tool *positional* rather than subcommand groups — the
//!   flags genuinely live on the command the placeholder follows.
//! - A subcommand or flag reached only through a variable is invisible to it.
//!
//! What it buys is that a dead form written the way all twelve #2065 sites were
//! written — a literal, on one line, in one of the scanned files — cannot be
//! reintroduced in silence, in either half.
#![cfg(feature = "benchmark")]

use std::path::{Path, PathBuf};

use clap::{Command, CommandFactory};
use ferro_hgvs::benchmark::cli::Cli;

/// An invocation lifted out of a source string, with where it was found.
#[derive(Debug)]
struct Form {
    /// The subcommand tokens following `ferro-benchmark` (e.g. `["setup", "uta"]`).
    path: Vec<String>,
    /// The flag tokens spelled after the path on the same line, as written
    /// (e.g. `["--output-dir", "-o"]`).
    flags: Vec<String>,
    /// `file:line`, for a legible failure message.
    origin: String,
}

/// Whether a token can begin (or continue) a subcommand path: a lowercase word,
/// optionally hyphenated, as clap subcommand and tool-positional names are
/// spelled (`setup`, `start-uta`, `hgvs-rs`).
fn is_command_word(tok: &str) -> bool {
    let mut chars = tok.chars();
    match chars.next() {
        Some(c) if c.is_ascii_lowercase() => {}
        _ => return false,
    }
    tok.chars()
        .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '-')
}

/// Whether a token ends the invocation rather than continuing it: a shell
/// operator, a redirection, a comment, or a line-continuation backslash.
///
/// This is what keeps flags attributed to the right binary. `pixi.toml` spells
/// `… --bin ferro-benchmark -- check biocommons && cargo run --release …` on one
/// line; without the boundary, cargo's own `--release` would be read as a flag of
/// `ferro-benchmark check` and reported as dead.
fn is_invocation_boundary(tok: &str) -> bool {
    matches!(tok, "&&" | "||" | ";" | "|" | "&" | ">" | ">>" | "<" | "<<")
        || tok.starts_with('#')
        || tok.contains(">&")
        // A trailing `\` (or `\\`, inside a Rust string literal) continues the
        // command onto a line this scan cannot see.
        || tok.chars().all(|c| c == '\\')
}

/// The flag a token spells, if it spells one: a bare `--long-name` or `-s`.
/// A placeholder (`--{}`), a bare `--`, an `=`-joined spelling and a negative
/// number are all refused rather than guessed at.
fn as_flag_token(tok: &str) -> Option<&str> {
    if let Some(long) = tok.strip_prefix("--") {
        let starts_lower = long.starts_with(|c: char| c.is_ascii_lowercase());
        let body_ok = long
            .chars()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '-');
        return (starts_lower && body_ok).then_some(tok);
    }
    let short = tok.strip_prefix('-')?;
    let mut chars = short.chars();
    match (chars.next(), chars.next()) {
        (Some(c), None) if c.is_ascii_alphabetic() => Some(tok),
        _ => None,
    }
}

/// Extract the subcommand path and the flags from the text immediately following
/// one `ferro-benchmark` occurrence. Returns `None` when the occurrence is a bare
/// mention rather than an invocation (an inline `` `ferro-benchmark` ``, a
/// possessive, `ferro-benchmark --help`, and so on).
fn extract_form(after: &str) -> Option<(Vec<String>, Vec<String>)> {
    // The name must be followed by whitespace to be an invocation; a backtick,
    // apostrophe or period means it is being named, not run.
    let first = after.chars().next()?;
    if !first.is_whitespace() {
        return None;
    }

    let mut tokens = after.split_whitespace().peekable();

    // A `cargo run … --bin ferro-benchmark -- <cmd>` separator: step over it.
    if tokens.peek() == Some(&"--") {
        tokens.next();
    }

    let mut path = Vec::new();
    let mut flags = Vec::new();
    // The path runs until the first token that cannot be a subcommand; flags are
    // collected from there to the end of the invocation.
    let mut path_open = true;

    for raw in tokens {
        if is_invocation_boundary(raw) {
            break;
        }
        // Strip surrounding markup/punctuation so a backtick- or quote-wrapped
        // command (`` `ferro-benchmark generate-biocommons-settings` ``) or a
        // trailing comma is still recognised.
        let tok = raw.trim_matches(|c: char| "`'\".,;:()[]".contains(c));

        if path_open && is_command_word(tok) {
            path.push(tok.to_string());
            // A subcommand path is at most a group plus a leaf plus a tool
            // positional; nothing legitimate runs deeper.
            if path.len() == 4 {
                path_open = false;
            }
            continue;
        }
        path_open = false;
        if let Some(flag) = as_flag_token(tok) {
            flags.push(flag.to_string());
        }
    }

    if path.is_empty() {
        None
    } else {
        Some((path, flags))
    }
}

/// Scan one file's text for every `ferro-benchmark` invocation form.
fn scan_text(text: &str, rel_path: &str, out: &mut Vec<Form>) {
    const NEEDLE: &str = "ferro-benchmark";
    for (line_no, line) in text.lines().enumerate() {
        let mut search_from = 0;
        while let Some(idx) = line[search_from..].find(NEEDLE) {
            let start = search_from + idx;
            let after = &line[start + NEEDLE.len()..];
            if let Some((path, flags)) = extract_form(after) {
                out.push(Form {
                    path,
                    flags,
                    origin: format!("{}:{}", rel_path, line_no + 1),
                });
            }
            search_from = start + NEEDLE.len();
        }
    }
}

/// What a subcommand path resolved to.
enum Resolved<'a> {
    /// The path names live commands. `Command` is the one that owns any flags
    /// spelled after it — either the leaf the path reached, or the command whose
    /// positional absorbed the remaining tokens.
    Command(&'a Command),
    /// A token could only have been a subcommand, and no such subcommand exists.
    Dead,
}

/// Resolve a subcommand path against the clap tree.
fn resolve_path<'a>(root: &'a Command, path: &[String]) -> Resolved<'a> {
    let mut cur = root;
    for tok in path {
        if let Some(sub) = cur.find_subcommand(tok) {
            cur = sub;
        } else if cur.get_positionals().next().is_some() {
            // `cur` takes a positional; this token (and the rest) are values, so
            // `cur` is still the command any following flags belong to.
            return Resolved::Command(cur);
        } else {
            return Resolved::Dead;
        }
    }
    Resolved::Command(cur)
}

/// Whether `cmd` declares an argument the flag token names. `flag` is as written
/// (`--output-dir` or `-o`).
fn command_accepts_flag(cmd: &Command, flag: &str) -> bool {
    if let Some(long) = flag.strip_prefix("--") {
        return cmd.get_arguments().any(|arg| {
            arg.get_long() == Some(long)
                || arg
                    .get_all_aliases()
                    .is_some_and(|aliases| aliases.contains(&long))
        });
    }
    let short = flag.trim_start_matches('-').chars().next();
    cmd.get_arguments()
        .any(|arg| arg.get_short().is_some() && arg.get_short() == short)
}

/// The clap tree, `build`-ed so that auto-generated (`--help`) and propagated
/// arguments are present on every subcommand — otherwise a legitimate
/// `--help` would read as a dead flag.
fn built_cli() -> Command {
    let mut root = Cli::command();
    root.build();
    root
}

/// Repository-relative source paths whose `ferro-benchmark` strings this guard
/// polices. Missing files are skipped (docs may be reorganised); the scan is a
/// floor, so a dropped path only narrows coverage, never breaks the build.
///
/// Not every file mentioning `ferro-benchmark` belongs here. The list is
/// deliberately confined to files whose mentions are *invocations* — the repo's
/// own `CLAUDE.md` and `tests/fixtures/README.md`, for instance, name the tool
/// mid-sentence in prose the token scan would read as a subcommand path.
const SCAN_TARGETS: &[&str] = &[
    "src/bin/benchmark.rs",
    "src/benchmark/mod.rs",
    "src/benchmark/cli.rs",
    "src/benchmark/biocommons.rs",
    "src/benchmark/compare.rs",
    "src/benchmark/hgvs_rs.rs",
    "pixi.toml",
    "README.md",
    "docs/BENCHMARK_GUIDE.md",
    "docs/BENCHMARK_RUNBOOK.md",
    "docs/BIOCOMMONS_LOCAL_SETUP.md",
];

fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn collect_forms() -> Vec<Form> {
    let root = manifest_dir();
    let mut forms = Vec::new();
    for rel in SCAN_TARGETS {
        let path: &Path = Path::new(rel);
        let full = root.join(path);
        if let Ok(text) = std::fs::read_to_string(&full) {
            scan_text(&text, rel, &mut forms);
        }
    }
    forms
}

/// The whole point: every scanned invocation form resolves against the CLI —
/// subcommand path and flags alike.
#[test]
fn every_documented_invocation_form_resolves() {
    let root = built_cli();
    let forms = collect_forms();

    // A scanner that finds nothing would pass vacuously; the sources carry many
    // forms, so guard the floor. The flag floor is separate because the flag half
    // could go blind (a scanner change that stops at the first flag again) while
    // the subcommand half stayed healthy — which is the state this guard shipped
    // in, and the state #2074's dead flag survived.
    assert!(
        forms.len() >= 15,
        "scanner found only {} ferro-benchmark forms — the scan is probably broken",
        forms.len()
    );
    // Measured on the tree this shipped against: 185 flag tokens across 146
    // forms. The floor is a vacuity guard, not a census — it is set well below
    // the measurement so that reorganising a doc does not redden it, and well
    // above zero so that a scanner which stops collecting flags does.
    let flags_seen: usize = forms.iter().map(|f| f.flags.len()).sum();
    assert!(
        flags_seen >= 50,
        "scanner found only {flags_seen} flag tokens across {} forms — the flag half is probably blind",
        forms.len()
    );

    let mut dead: Vec<String> = Vec::new();
    for form in &forms {
        let spelling = form.path.join(" ");
        match resolve_path(&root, &form.path) {
            Resolved::Dead => dead.push(format!(
                "  {}  ->  ferro-benchmark {spelling}  (no such subcommand)",
                form.origin,
            )),
            Resolved::Command(cmd) => {
                for flag in &form.flags {
                    if !command_accepts_flag(cmd, flag) {
                        dead.push(format!(
                            "  {}  ->  ferro-benchmark {spelling} {flag}  (`{}` declares no {flag})",
                            form.origin,
                            cmd.get_name(),
                        ));
                    }
                }
            }
        }
    }

    assert!(
        dead.is_empty(),
        "dead ferro-benchmark invocation form(s) found — each names a subcommand or flag the CLI no longer accepts:\n{}",
        dead.join("\n")
    );
}

/// The resolver must distinguish a live grouped form from a retired flat one —
/// otherwise the guard above could pass while asserting nothing. These are the
/// exact #2065 before/after spellings.
#[test]
fn resolver_separates_live_forms_from_dead_ones() {
    let root = built_cli();
    let resolves = |form: &[&str]| {
        let path: Vec<String> = form.iter().map(|s| s.to_string()).collect();
        matches!(resolve_path(&root, &path), Resolved::Command(_))
    };

    let live: &[&[&str]] = &[
        &["benchmark", "normalize"],
        &["extract", "clinvar"],
        &["extract", "sample"],
        &["setup", "uta"],
        &["setup", "seqrepo"],
        &["setup", "start-uta"],
        &["parse", "ferro"], // `ferro` is the tool positional, not a subcommand
        &["prepare", "biocommons"],
        &["prepare", "mutalyzer"],
        &["benchmark", "matrix"],
    ];
    for form in live {
        assert!(
            resolves(form),
            "expected live form to resolve: ferro-benchmark {}",
            form.join(" ")
        );
    }

    let dead: &[&[&str]] = &[
        &["compare", "normalize"],
        &["extract-clinvar"],
        &["setup-uta"],
        &["setup-seqrepo"],
        &["start-uta"],
        &["populate-protein-cache"],
        &["generate-biocommons-settings"],
        &["parse-ferro"],
        &["run"],
        &["sample"], // `sample` lives under `extract`, never at the top level
    ];
    for form in dead {
        assert!(
            !resolves(form),
            "expected dead form to be rejected: ferro-benchmark {}",
            form.join(" ")
        );
    }
}

/// The flag half of the resolver, pinned the same way. Every `dead` row below is
/// a flag that **exists somewhere in this CLI** and is simply not on the command
/// it is paired with — which is the only shape worth guarding, and the shape
/// #2074 shipped: `--seqrepo-dir` is real, on `prepare`, and dead on
/// `setup seqrepo`.
#[test]
fn resolver_separates_live_flags_from_dead_ones() {
    let root = built_cli();
    let accepts = |form: &[&str], flag: &str| match resolve_path(
        &root,
        &form.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
    ) {
        Resolved::Command(cmd) => command_accepts_flag(cmd, flag),
        Resolved::Dead => panic!("fixture path does not resolve: ferro-benchmark {form:?}"),
    };

    let live: &[(&[&str], &str)] = &[
        (&["setup", "seqrepo"], "--output-dir"),
        (&["setup", "seqrepo"], "-o"),
        (&["setup", "seqrepo"], "--instance"),
        (&["setup", "uta"], "--uta-dump"),
        (&["prepare", "biocommons"], "--seqrepo-dir"),
        (&["prepare", "mutalyzer"], "--ferro-reference"),
        (&["check", "biocommons"], "--seqrepo-path"),
        (&["benchmark", "normalize"], "--validator"),
        (&["benchmark", "normalize"], "--uta-db-url"),
        (&["normalize", "mutalyzer"], "-i"),
        (&["compare", "results", "normalize"], "-o"),
        // Auto-generated arguments must resolve too, or a documented
        // `--help` would read as dead.
        (&["setup", "seqrepo"], "--help"),
    ];
    for (form, flag) in live {
        assert!(
            accepts(form, flag),
            "expected live flag to resolve: ferro-benchmark {} {flag}",
            form.join(" ")
        );
    }

    let dead: &[(&[&str], &str)] = &[
        // The #2074 defect itself, in the exact spelling that shipped.
        (&["setup", "seqrepo"], "--seqrepo-dir"),
        (&["setup", "seqrepo"], "--seqrepo-path"),
        (&["check", "biocommons"], "--seqrepo-dir"),
        (&["prepare", "biocommons"], "--seqrepo-path"),
        (&["extract", "clinvar"], "--seed"),
        (&["prepare", "ferro"], "--validator"),
        (&["setup", "seqrepo"], "-x"),
    ];
    for (form, flag) in dead {
        assert!(
            !accepts(form, flag),
            "expected dead flag to be rejected: ferro-benchmark {} {flag}",
            form.join(" ")
        );
    }
}

/// A flag belongs to the binary it was written for. `pixi.toml` chains two
/// invocations with `&&` on one line, so without a shell boundary the scan would
/// attribute `cargo run`'s own `--release` / `--bin` to `ferro-benchmark check`
/// and report a defect that is not there.
#[test]
fn the_scanner_attributes_flags_to_the_right_binary() {
    let (path, flags) = extract_form(
        " -- check biocommons && cargo run --release --features benchmark \
         --bin ferro-benchmark -- check mutalyzer",
    )
    .expect("the chained form is an invocation");
    assert_eq!(path, ["check", "biocommons"]);
    assert!(
        flags.is_empty(),
        "cargo's flags must not be attributed to ferro-benchmark; got {flags:?}"
    );

    // A line-continuation ends the invocation for the same reason: the flags are
    // on a line this scan cannot see, so none may be attributed here.
    let (path, flags) = extract_form(" prepare mutalyzer \\").expect("an invocation");
    assert_eq!(path, ["prepare", "mutalyzer"]);
    assert!(flags.is_empty(), "got {flags:?}");

    // And the flags that ARE on the line are collected, placeholders and
    // values skipped.
    let (path, flags) =
        extract_form(" prepare {} --seqrepo-dir {} --uta-dump <path>").expect("an invocation");
    assert_eq!(path, ["prepare"]);
    assert_eq!(flags, ["--seqrepo-dir", "--uta-dump"]);
}
