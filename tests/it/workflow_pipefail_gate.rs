//! A workflow step whose failure is meant to fail the job must not be able to
//! exit 0 with its real work broken.
//!
//! # The failure this exists for
//!
//! `nightly-mutalyzer.yml`'s `Diff FAIL set against xfail baseline` step — the
//! gate written to close #1998, whose entire purpose is to turn the nightly red
//! on drift — ended in
//!
//! ```sh
//! python3 scripts/diff_nightly_xfail_baseline.py --junit "$JUNIT" --baseline "$BASELINE" \
//!   | tee -a "$GITHUB_STEP_SUMMARY"
//! ```
//!
//! A step that names no `shell:` runs under `bash -e {0}`. `pipefail` is **not**
//! in that, and a pipeline's exit status under `-e` is its **last** command's —
//! so `tee` reported 0 while the diff exited 1. Measured over a synthetic JUnit
//! carrying the 38 baselined failures plus one new one, the step printed
//! `### NEW FAILURE (1)` and exited **0**; the same script under
//! `bash --noprofile --norc -eo pipefail {0}` exited **1**. The job concluded
//! `success`, `report-failure` (`if: failure()`) never fired, and no issue
//! opened. #1998's own failure mode, reproduced inside its own fix (#2082).
//!
//! The convention was already there to be slipped from: `report-failure.yml`,
//! `nightly-perf.yml`, `nightly-soak.yml`, `ci.yml` and `release-binaries.yml`
//! all set `pipefail` explicitly in the blocks that need it. One step did not,
//! and **nothing could see that** — `actionlint` exits 0 on the file either way.
//! This test is the thing that can see it.
//!
//! # Where the line is drawn, and why it is drawn there
//!
//! **Every `run:` block that contains a pipeline must enable `pipefail`.** Not
//! "every gating block", not "every block whose pipe is the last command".
//!
//! Three narrower rules were considered and all three are rejected for the same
//! reason — each keys the guard on a fact that can change in a *later* PR
//! without anyone re-reading the pipeline it exempted:
//!
//! - **"only steps without `continue-on-error: true`."** GitHub's own semantics
//!   make that the exact definition of a gating step, so it is tempting. But
//!   `continue-on-error` is a one-line field on an unrelated axis: deleting it
//!   promotes the step to a gate, and nothing in that diff mentions the pipe.
//!   The nightly job in this very repository carries a `continue-on-error` step
//!   and a gating step ten lines apart.
//! - **"only where the pipe is the last command."** The status of *any*
//!   pipeline under `-e` is its last command's, so `curl … | sh` in the middle of
//!   a block swallows a failed download exactly as thoroughly. Two of the four
//!   pre-existing violations this guard found on its first run were that shape.
//! - **"only where the pipeline's status is unobserved"** — i.e. exempt a
//!   pipeline used as an `if`/`while` condition, or as the left operand of
//!   `&&`/`||`, where the shell reads the status rather than discarding it.
//!   That is a real distinction, and it is still the wrong line: `pipefail`
//!   changes *which* status those constructs read, so the exemption trades a
//!   silent exit-0 for a silently-taken branch. It also requires a shell parser
//!   good enough to tell a condition from a command, which is precisely the kind
//!   of judgement a guard should not be making at 2am on someone else's PR.
//!
//! So the rule is the loud one, and the escape hatch is explicit rather than
//! inferred: a block may carry a [`EXEMPT_MARKER`] comment naming a reason. That
//! is a decision recorded at the site, in the diff, by the person making it —
//! which is what the three rules above cannot be. Nothing is exempt today.
//!
//! The cost of the loud rule is small and was measured before choosing it: the
//! whole of `.github/workflows/` holds nine `run:` blocks containing a pipeline.
//!
//! # The set of workflows is DERIVED, following #2060
//!
//! `release_tag_gate.rs`'s `every_release_reachable_root_job_is_tag_gated` is
//! the model. A guard that iterates a hardcoded roster of the gates it already
//! knows about is blind to the gate someone adds next week — measured there,
//! where deleting a gate outright left four assertions passing at EXIT=0. So
//! every workflow, every job and every step here comes out of
//! `.github/workflows/`, and the floors below exist so that a scanner which has
//! gone blind fails loudly instead of passing over an empty set.
//!
//! # What it does NOT catch
//!
//! - **It is a lexer, not a shell.** [`ScriptScan`] tracks single quotes, double
//!   quotes, backslash escapes, `#` comments and here-documents, which is enough
//!   for every script in this repository and for the shapes that recur. It does
//!   not parse. A `case` pattern (`case $x in a|b)`) reads as a pipeline to it —
//!   a **false positive**, whose remedy is to add `pipefail` (harmless) or the
//!   exemption marker (deliberate). False positives are the direction chosen.
//! - **A pipeline inside a double-quoted command substitution is missed.**
//!   `line="$(a | b)"` is code the scanner reads as string content, because it
//!   does not track `$(…)` nesting. A **false negative**, and the one this file
//!   would fix first if the tree grew a case that mattered — today the two such
//!   lines (`nightly-perf.yml:147`, `prune-ci-caches.yml:88`) are in blocks that
//!   set `pipefail` anyway. Backticks are *not* treated as quotes, so a pipeline
//!   inside `` `…` `` is seen.
//! - **`shell:` templates are read as text.** A custom `shell:` string counts as
//!   safe if the word `pipefail` appears in it, which is what the two forms in
//!   use here (`bash`, and an explicit `bash -eo pipefail {0}`) reduce to.
//! - **`shell: sh` is judged, and the remedy the message names is a bashism.**
//!   `sh -e {0}` carries no `pipefail` and `set -o pipefail` is not POSIX, so
//!   the fix there is `shell: bash`, not a `set`. Nothing uses `shell: sh`
//!   today; if something does, read the message with that in mind.
//! - **A Windows job's un-shelled steps are out of scope.** GitHub's default
//!   shell there is `pwsh`, whose `|` is an object pipeline and not a
//!   status-dropping construct, so the class does not exist. A job is treated as
//!   possibly-Windows if the word appears anywhere in its `runs-on:` or its
//!   `strategy.matrix`, which is deliberately over-broad — the failure direction
//!   is a missed check on a Windows job, never a false demand on a Linux one.
//! - **It says nothing about what a step *does* with a non-zero status.** A
//!   block can still `|| true` its own gate away, or `exit 0` at the end. This
//!   guard closes one seam: the status the shell *computes* for a pipeline.
//! - **Ordering is read by line, not by control flow.** A `set -o pipefail`
//!   after a pipeline correctly fails to protect it, but one inside a branch
//!   that is never taken (`if false; then set -o pipefail; fi`) is credited, and
//!   a `set` later on the *same line* as a pipeline is too. Both are false
//!   negatives; neither is a shape anyone writes, and closing them needs an
//!   evaluator rather than a lexer.
//! - **`composite` actions are not scanned.** There are none in this repository
//!   (`.github/actions/` does not exist); if one is added, its `run:` steps are
//!   invisible here. [`WORKFLOW_DIR`] is the derivation's boundary and that
//!   boundary is a real limit, not an oversight.
//! - **Committed `*.sh` files are out of scope, and that boundary was
//!   measured.** A `run:` block is where GitHub reads a shell's exit status as a
//!   step's verdict, which is the seam this file is about; a standalone script
//!   has no such contract until a workflow invokes it. Of the eight committed
//!   shell scripts, seven already set `pipefail` — including every one a
//!   workflow calls — and the eighth (`tests/validation_test.sh`) is a manual
//!   developer script that no workflow, document or source file references and
//!   that sets no `-e` either, so no status it computes gates anything.

use std::path::PathBuf;

use serde_yaml::Value;

/// Where the workflows live. The set is derived from this directory rather than
/// named; see the module docs for why a roster is the wrong shape here.
const WORKFLOW_DIR: &str = ".github/workflows";

/// The marker that exempts one `run:` block, written as a shell comment inside
/// the block itself:
///
/// ```sh
/// # pipefail-exempt: <why this pipeline's left-hand status may be discarded>
/// ```
///
/// A reason is **required** — a bare marker fails — and an exemption on a block
/// with no pipeline in it also fails, so a marker cannot outlive the pipeline it
/// was written for. That makes the exemption set shrink-only in the same sense
/// as `generator_completeness.rs`'s `LEDGER_EXEMPT`, without any list to
/// maintain: the marker lives at the site it excuses.
const EXEMPT_MARKER: &str = "pipefail-exempt:";

/// A floor on the number of workflows scanned, so a directory read or an `on:`
/// parse that has silently stopped matching fails loudly rather than passing
/// vacuously over an empty set. The analogue of
/// `release_tag_gate.rs`'s `MINIMUM_RELEASE_WORKFLOWS`.
///
/// Comfortably below the current count, because the floor's job is to catch
/// *blindness*, not to pin the repository's size.
const MINIMUM_WORKFLOWS: usize = 10;

/// A floor on the number of `run:` blocks reached, for the same reason. A job
/// or step traversal that has gone blind reports zero violations.
///
/// 117 today; the floor is set well below so that deleting a workflow is not a
/// failure, because the floor's subject is the *scanner*, not the repository.
const MINIMUM_RUN_STEPS: usize = 80;

/// A floor on the number of `run:` blocks in which the scanner found a pipeline.
///
/// **This is the one that matters**, and it is the reason the other two are not
/// enough. Every assertion in this file is driven by [`ScriptScan::pipelines`];
/// a scanner that stops recognising `|` — a quote-state bug, a comment-stripping
/// bug, an over-eager heredoc — finds no pipelines anywhere and reports a clean
/// tree. That is indistinguishable from success, which is the failure mode this
/// whole file is about.
const MINIMUM_PIPELINE_STEPS: usize = 5;

/// A floor on the number of `run:` blocks resolved to a bash-family shell, so a
/// shell-resolution bug that classifies everything as "not bash" — which would
/// also report a clean tree — fails instead.
///
/// 116 of the 117 today: exactly one step runs on a Windows runner without an
/// explicit `shell:`, and is out of scope for the reason the module docs give.
const MINIMUM_BASH_STEPS: usize = 80;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// How a `run:` block's script relates to `pipefail`, and where its pipelines
/// are.
///
/// Built by a single left-to-right pass so that quote state, comment state and
/// here-document state are all resolved once, consistently. Reading the script
/// with three independent greps is how the defect that prompted this file
/// survived review: the step's own comment contained the literal text
/// `set -o pipefail` — asserting, falsely, that the shell supplied it — so a
/// grep for that string over the raw block answers **yes** for the one block in
/// the repository where the answer is no.
#[derive(Debug, Default)]
struct ScriptScan {
    /// 1-based line numbers, within the `run:` block, carrying a pipeline
    /// operator outside quotes, comments and here-documents.
    pipelines: Vec<usize>,
    /// The line of the first `set` at a command position that enables
    /// `pipefail`, if any.
    ///
    /// A line, not a flag, because **order matters**: `curl … | sh` followed by
    /// `set -o pipefail` is not protected by it, and a boolean cannot tell the
    /// two apart.
    sets_pipefail_at: Option<usize>,
    /// The reason given on an [`EXEMPT_MARKER`] comment, if one is present.
    /// `Some("")` distinguishes a marker with no reason from no marker at all.
    exemption: Option<String>,
}

/// Whether a `#` at `index` in `line` opens a comment, given that the scanner is
/// not inside quotes there. Shell only treats `#` as a comment at the start of a
/// word, so `foo#bar` and `${x#y}` are not comments.
fn opens_comment(line: &[char], index: usize) -> bool {
    index == 0 || line[index - 1].is_whitespace()
}

/// The here-document delimiter a line opens, if any.
///
/// Handles `<<WORD`, `<< WORD`, `<<-WORD`, `<<'WORD'` and `<<"WORD"`. `<<<` is a
/// here-*string* and opens nothing. Only the first on a line is tracked; two
/// here-documents opened by one command line is a shape this repository does not
/// contain and the scanner does not claim to handle.
fn heredoc_delimiter(code: &str) -> Option<String> {
    let bytes: Vec<char> = code.chars().collect();
    let mut i = 0;
    while i + 1 < bytes.len() {
        if bytes[i] != '<' || bytes[i + 1] != '<' {
            i += 1;
            continue;
        }
        let mut j = i + 2;
        if j < bytes.len() && bytes[j] == '<' {
            // A here-string (`<<<`), not a here-document.
            i = j + 1;
            continue;
        }
        if j < bytes.len() && bytes[j] == '-' {
            j += 1;
        }
        while j < bytes.len() && bytes[j] == ' ' {
            j += 1;
        }
        let quote = bytes.get(j).copied().filter(|c| *c == '\'' || *c == '"');
        if quote.is_some() {
            j += 1;
        }
        let start = j;
        while j < bytes.len()
            && (bytes[j].is_alphanumeric() || bytes[j] == '_' || bytes[j] == '-' || bytes[j] == '.')
        {
            j += 1;
        }
        if j > start {
            return Some(bytes[start..j].iter().collect());
        }
        i += 2;
    }
    None
}

/// Whether `word_start` in one line of code is a command position — the start of
/// the line, or immediately after a `;`, `&`, `|`, `(`, `{`, or a `then`/`else`/
/// `do` keyword.
///
/// The second line of defence behind quote-stripping: it is what keeps
/// `git log --grep set -o pipefail` from reading as a `set`. Quoted prose about
/// `pipefail` — which this repository's scripts genuinely contain — never
/// reaches here, because the caller passes the line with quoted spans removed.
fn is_command_position(code: &str, word_start: usize) -> bool {
    let trimmed = code[..word_start].trim_end();
    trimmed.is_empty()
        || trimmed.ends_with(';')
        || trimmed.ends_with('&')
        || trimmed.ends_with('|')
        || trimmed.ends_with('(')
        || trimmed.ends_with('{')
        || trimmed.ends_with("then")
        || trimmed.ends_with("else")
        || trimmed.ends_with("do")
}

/// Scan one `run:` block.
///
/// The pass carries quote state **across** lines, because a shell string may
/// span them and a per-line reset would read the second half of one as code. It
/// carries here-document state for the same reason.
fn scan_script(script: &str) -> ScriptScan {
    let mut scan = ScriptScan::default();
    // `None` outside quotes; `Some('\'')` or `Some('"')` inside one.
    let mut quote: Option<char> = None;
    let mut heredoc: Option<String> = None;

    for (lineno, raw) in script.lines().enumerate() {
        let lineno = lineno + 1;
        if let Some(delimiter) = &heredoc {
            if raw.trim() == delimiter.as_str() {
                heredoc = None;
            }
            continue;
        }

        let chars: Vec<char> = raw.chars().collect();
        // The line with comments stripped, quotes and their contents RETAINED:
        // what the here-document detector reads, since `<<'EOF'` quotes its own
        // delimiter.
        let mut line_code = String::new();
        // The same line with quoted spans and the quote characters themselves
        // REMOVED: what the pipeline search and the `set … pipefail` search
        // read. Stripping quotes is what makes `echo "set -o pipefail"` prose
        // rather than a `set`, and a `|` inside a string not a pipeline.
        //
        // Built during the same pass rather than re-derived afterwards, because
        // quote state carries ACROSS lines — a per-line re-derivation reads the
        // tail of a multi-line string as code, which is a shape
        // `nightly-soak.yml` actually contains.
        let mut unquoted: Vec<char> = Vec::new();
        let mut i = 0;
        while i < chars.len() {
            let c = chars[i];
            match quote {
                None => {
                    if c == '\\' {
                        // Escapes the next character, whatever it is.
                        i += 2;
                        continue;
                    }
                    if c == '\'' || c == '"' {
                        quote = Some(c);
                        line_code.push(c);
                        i += 1;
                        continue;
                    }
                    if c == '#' && opens_comment(&chars, i) {
                        let comment: String = chars[i..].iter().collect();
                        if let Some(at) = comment.find(EXEMPT_MARKER) {
                            let reason = comment[at + EXEMPT_MARKER.len()..].trim().to_string();
                            scan.exemption = Some(reason);
                        }
                        break;
                    }
                    line_code.push(c);
                    unquoted.push(c);
                    i += 1;
                }
                Some(q) => {
                    if q == '"' && c == '\\' {
                        i += 2;
                        continue;
                    }
                    if c == q {
                        quote = None;
                    }
                    line_code.push(c);
                    i += 1;
                }
            }
        }

        let mut k = 0;
        while k < unquoted.len() {
            if unquoted[k] != '|' {
                k += 1;
                continue;
            }
            if unquoted.get(k + 1) == Some(&'|') {
                // `||` — a logical OR, not a pipeline. Skip both characters so a
                // third `|` beside them cannot read as a pipeline either.
                k += 2;
                continue;
            }
            if k > 0 && (unquoted[k - 1] == '>' || unquoted[k - 1] == '|') {
                // `>|` is the noclobber override, and a `|` whose predecessor is
                // a `|` is the tail of an operator already read.
                k += 1;
                continue;
            }
            scan.pipelines.push(lineno);
            k += 1;
        }

        if scan.sets_pipefail_at.is_none() {
            let line: String = unquoted.iter().collect();
            if sets_pipefail(&line) {
                scan.sets_pipefail_at = Some(lineno);
            }
        }

        if quote.is_none() {
            if let Some(delimiter) = heredoc_delimiter(&line_code) {
                heredoc = Some(delimiter);
            }
        }
    }

    scan
}

/// Whether one line of comment-stripped, **unquoted** script enables `pipefail`
/// with a `set` at a command position.
///
/// Two refusals, both deliberate and both erring the same way — toward demanding
/// `pipefail` rather than crediting it:
///
/// - `set +o pipefail` **disables** it, and is read as such.
/// - The text must survive quote-stripping, so `echo "set -o pipefail"` is
///   prose. The corollary is that `set -o "pipefail"` is not recognised; nobody
///   writes it, and being wrong in that direction costs a spurious demand rather
///   than a missed one.
fn sets_pipefail(code: &str) -> bool {
    let mut from = 0;
    while let Some(at) = code[from..].find("set") {
        let start = from + at;
        from = start + 3;
        let after = code[start + 3..].chars().next();
        if after.is_some_and(|c| !c.is_whitespace()) {
            continue;
        }
        if !is_command_position(code, start) {
            continue;
        }
        // The rest of this simple command: up to a `;`, `&`, `|` or newline.
        let rest_end = code[start..]
            .find(['\n', ';', '&', '|'])
            .map_or(code.len(), |offset| start + offset);
        let words: Vec<&str> = code[start + 3..rest_end].split_whitespace().collect();
        if !words.contains(&"pipefail") {
            continue;
        }
        // `set -o pipefail` enables, `set +o pipefail` disables.
        let disabled = words
            .iter()
            .take_while(|word| **word != "pipefail")
            .any(|word| word.starts_with('+'));
        if !disabled {
            return true;
        }
    }
    false
}

/// Every `.github/workflows/*.{yml,yaml}`, as `(repo-relative path, parsed
/// document)` pairs sorted by path.
fn workflows() -> Vec<(String, Value)> {
    let dir = repo_root().join(WORKFLOW_DIR);
    let mut found = Vec::new();
    let entries = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("read {}: {e}", dir.display()))
        .filter_map(Result::ok);
    for entry in entries {
        let path = entry.path();
        if path
            .extension()
            .is_none_or(|extension| extension != "yml" && extension != "yaml")
        {
            continue;
        }
        let name = path
            .file_name()
            .expect("a file has a name")
            .to_string_lossy()
            .into_owned();
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        let document: Value = serde_yaml::from_str(&text)
            .unwrap_or_else(|e| panic!("{WORKFLOW_DIR}/{name} is not valid YAML: {e}"));
        found.push((format!("{WORKFLOW_DIR}/{name}"), document));
    }
    found.sort_by(|(left, _), (right, _)| left.cmp(right));

    assert!(
        found.len() >= MINIMUM_WORKFLOWS,
        "found {} workflow file(s) in {WORKFLOW_DIR}; at least {MINIMUM_WORKFLOWS} are there \
         today, so the directory read has gone blind rather than the repository having shrunk. \
         Every assertion in this file iterates over that set, so an empty one passes them all.",
        found.len(),
    );
    found
}

/// The `defaults.run.shell` declared at some level, if any.
fn declared_shell(node: &Value) -> Option<&str> {
    node.get("defaults")?.get("run")?.get("shell")?.as_str()
}

/// Whether a job might run on a Windows runner, read over-broadly from its
/// `runs-on:` and its whole `strategy.matrix`.
///
/// Over-broad on purpose: the guard's failure direction here must be a missed
/// check on a Windows job (where the class does not exist — `pwsh` has no
/// `pipefail` and its `|` is not a status-dropping construct), never a demand
/// for `set -o pipefail` in a PowerShell script.
fn may_run_on_windows(job: &Value) -> bool {
    fn mentions_windows(node: &Value) -> bool {
        match node {
            Value::String(text) => text.to_ascii_lowercase().contains("windows"),
            Value::Sequence(items) => items.iter().any(mentions_windows),
            Value::Mapping(entries) => entries.values().any(mentions_windows),
            _ => false,
        }
    }
    job.get("runs-on").is_some_and(mentions_windows)
        || job
            .get("strategy")
            .and_then(|strategy| strategy.get("matrix"))
            .is_some_and(mentions_windows)
}

/// One `run:` block, located and scanned.
#[derive(Debug)]
struct RunStep {
    file: String,
    job: String,
    /// The step's `name:`, or a positional label when it has none.
    step: String,
    /// The resolved `shell:`, or `None` for GitHub's default.
    shell: Option<String>,
    /// Whether the resolved shell is one where `pipefail` is the question — a
    /// bash-family shell that does not already supply it.
    bash_family: bool,
    /// Whether the resolved shell itself supplies `pipefail`.
    shell_supplies_pipefail: bool,
    scan: ScriptScan,
}

impl RunStep {
    /// The location, in the form a failure message should name it.
    fn locus(&self) -> String {
        format!(
            "{} :: job `{}` :: step `{}`",
            self.file, self.job, self.step
        )
    }

    /// The pipelines in this block that `pipefail` does not cover.
    ///
    /// Empty when the shell supplies `pipefail`, or when an in-block `set`
    /// reaches them. **Order is part of the question**: a `set -o pipefail`
    /// written *after* a pipeline does not protect it, and a block-level boolean
    /// cannot see that.
    fn unprotected_pipelines(&self) -> Vec<usize> {
        if !self.bash_family || self.shell_supplies_pipefail {
            return Vec::new();
        }
        self.scan
            .pipelines
            .iter()
            .copied()
            .filter(|line| self.scan.sets_pipefail_at.is_none_or(|set| *line < set))
            .collect()
    }

    /// Whether this block's pipelines can silently discard a left-hand failure.
    fn is_violation(&self) -> bool {
        self.scan.exemption.is_none() && !self.unprotected_pipelines().is_empty()
    }
}

/// Every `run:` block in every workflow, scanned.
fn run_steps() -> Vec<RunStep> {
    let mut steps = Vec::new();
    for (file, workflow) in &workflows() {
        let workflow_shell = declared_shell(workflow);
        let Some(jobs) = workflow.get("jobs").and_then(Value::as_mapping) else {
            // A reusable-workflow document with no `jobs:` is not a shape this
            // repository has, but an absent mapping must not read as "no work".
            panic!("{file} declares no jobs mapping");
        };
        for (job_name, job) in jobs {
            let Some(job_name) = job_name.as_str() else {
                continue;
            };
            let job_shell = declared_shell(job).or(workflow_shell);
            let windows = may_run_on_windows(job);
            let Some(job_steps) = job.get("steps").and_then(Value::as_sequence) else {
                // A `uses:` job calls a reusable workflow and has no steps of
                // its own; that workflow is scanned in its own right.
                continue;
            };
            for (index, step) in job_steps.iter().enumerate() {
                let Some(script) = step.get("run").and_then(Value::as_str) else {
                    continue;
                };
                let shell = step
                    .get("shell")
                    .and_then(Value::as_str)
                    .or(job_shell)
                    .map(str::to_string);
                let name = step
                    .get("name")
                    .and_then(Value::as_str)
                    .map(str::to_string)
                    .unwrap_or_else(|| format!("#{index} (unnamed)"));
                let (bash_family, supplies) = classify_shell(shell.as_deref(), windows);
                steps.push(RunStep {
                    file: file.clone(),
                    job: job_name.to_string(),
                    step: name,
                    shell,
                    bash_family,
                    shell_supplies_pipefail: supplies,
                    scan: scan_script(script),
                });
            }
        }
    }
    steps
}

/// `(is a bash-family shell, that shell already supplies pipefail)`.
///
/// GitHub expands `shell: bash` to `bash --noprofile --norc -eo pipefail {0}`
/// and the **default** (no `shell:`, non-Windows runner) to `bash -e {0}` —
/// which is the whole of this file's subject. `sh` expands to `sh -e {0}`, also
/// without `pipefail`; it is in scope, and `set -o pipefail` is a bashism that
/// `dash` rejects, so the remedy there is `shell: bash` rather than a `set`.
/// Nothing in this repository uses `shell: sh`.
fn classify_shell(shell: Option<&str>, windows_runner: bool) -> (bool, bool) {
    match shell {
        None => (!windows_runner, false),
        Some(text) => {
            let first = text.split_whitespace().next().unwrap_or("");
            let bash_family = matches!(first, "bash" | "sh");
            if !bash_family {
                return (false, false);
            }
            // `shell: bash` is expanded by GitHub with `-eo pipefail`; a custom
            // template says for itself.
            let supplies = first == "bash" && (text.trim() == "bash" || text.contains("pipefail"));
            (true, supplies)
        }
    }
}

/// **The assertion.** No `run:` block may contain a pipeline whose left-hand
/// failure the shell will discard.
///
/// This is the one that would have failed on the tree that shipped #2079's
/// first revision, and it names the file and the step so the failure is
/// actionable without opening a workflow.
#[test]
fn no_run_block_pipes_without_pipefail() {
    let steps = run_steps();
    let violations: Vec<&RunStep> = steps.iter().filter(|step| step.is_violation()).collect();

    let report: Vec<String> = violations
        .iter()
        .map(|step| {
            format!(
                "\n  {}\n    shell: {} (GitHub expands a step naming no `shell:` to `bash -e {{0}}` \
                 — no `pipefail`)\n    unprotected pipeline(s) at line(s) {:?} of the `run:` \
                 block; `set … pipefail` at {:?}",
                step.locus(),
                step.shell.as_deref().unwrap_or("<default>"),
                step.unprotected_pipelines(),
                step.scan.sets_pipefail_at,
            )
        })
        .collect();

    assert!(
        violations.is_empty(),
        "{} `run:` block(s) contain a pipeline and do not enable `pipefail`:{}\n\n\
         Under `bash -e`, a pipeline's exit status is its LAST command's, so a failure on the \
         left of the `|` is discarded and the step exits 0 — which is how the gate written to \
         close #1998 came to report success while printing `NEW FAILURE` (#2082).\n\n\
         Fix by adding `set -euo pipefail` as the first line of the `run:` block (the convention \
         this repository already follows in `report-failure.yml`, `nightly-perf.yml`, \
         `nightly-soak.yml`, `ci.yml` and `release-binaries.yml`), or `shell: bash` on the step, \
         which GitHub expands with `-eo pipefail`.\n\n\
         If the discarded status is genuinely fine here, say so at the site:\n    \
         # {EXEMPT_MARKER} <reason>\n\
         A reason is required, and the marker fails if the block stops containing a pipeline, so \
         it cannot outlive what it excuses.",
        violations.len(),
        report.join(""),
    );

    eprintln!(
        "{} `run:` block(s) scanned; {} contain a pipeline, all with `pipefail`",
        steps.len(),
        steps
            .iter()
            .filter(|step| !step.scan.pipelines.is_empty())
            .count(),
    );
}

/// The scan reached enough of the repository to be able to fail.
///
/// Split from the assertion above rather than folded into it because the two ask
/// opposite questions: that one fails when it finds something, this one fails
/// when it finds *nothing*, and a guard whose subject has silently emptied
/// passes the first while asserting nothing at all. `spec_fixture_setup_filter`
/// (#2033) and `release_tag_gate` (#2060) both carry this pair for the same
/// reason.
#[test]
fn the_scan_is_not_vacuous() {
    let steps = run_steps();
    let bash: Vec<&RunStep> = steps.iter().filter(|step| step.bash_family).collect();
    let piping: Vec<&RunStep> = steps
        .iter()
        .filter(|step| !step.scan.pipelines.is_empty())
        .collect();

    assert!(
        steps.len() >= MINIMUM_RUN_STEPS,
        "found {} `run:` block(s) across {WORKFLOW_DIR}; at least {MINIMUM_RUN_STEPS} are there \
         today, so the job/step traversal has gone blind. A traversal that reaches no step \
         reports no violation.",
        steps.len(),
    );
    assert!(
        bash.len() >= MINIMUM_BASH_STEPS,
        "only {} of {} `run:` block(s) resolved to a bash-family shell; at least \
         {MINIMUM_BASH_STEPS} do today. `classify_shell` has stopped recognising the default \
         shell, and a guard that believes nothing is bash demands `pipefail` of nothing.",
        bash.len(),
        steps.len(),
    );
    assert!(
        piping.len() >= MINIMUM_PIPELINE_STEPS,
        "the scanner found a pipeline in only {} of {} `run:` block(s); at least \
         {MINIMUM_PIPELINE_STEPS} contain one today.\n\
         This is the floor that matters: every assertion in this file is driven by the pipeline \
         scan, so a quote-state, comment or here-document bug in `scan_script` reports a clean \
         tree — which is exactly what a clean tree looks like.",
        piping.len(),
        steps.len(),
    );

    eprintln!(
        "{} run steps, {} bash-family, {} containing a pipeline",
        steps.len(),
        bash.len(),
        piping.len(),
    );
}

/// An exemption marker must carry a reason, and must still have a pipeline to
/// excuse.
///
/// The second half is what keeps the exemption set shrink-only without a list:
/// a marker left behind after its pipeline is deleted fails, so it is removed in
/// the change that removes the pipeline rather than accumulating as a licence
/// nobody re-reads. Nothing carries the marker today, which is why this test
/// also states its own denominator.
#[test]
fn every_pipefail_exemption_is_reasoned_and_live() {
    let steps = run_steps();
    let exempt: Vec<&RunStep> = steps
        .iter()
        .filter(|step| step.scan.exemption.is_some())
        .collect();

    for step in &exempt {
        let reason = step.scan.exemption.as_deref().unwrap_or_default();
        assert!(
            !reason.is_empty(),
            "{} carries a bare `# {EXEMPT_MARKER}` with no reason. The marker's whole value is \
             that it records a decision someone made; an unreasoned one records only that the \
             check was silenced.",
            step.locus(),
        );
        assert!(
            !step.scan.pipelines.is_empty(),
            "{} carries `# {EXEMPT_MARKER} {reason}` but its `run:` block contains no pipeline. \
             Delete the marker in the change that deleted the pipeline — a stale exemption is a \
             standing licence for whatever pipeline is written there next.",
            step.locus(),
        );
    }

    eprintln!("{} pipefail exemption(s) in force", exempt.len());
}

// ---------------------------------------------------------------------------
// The scanner's own tests.
//
// The floors above catch a scanner that has gone blind across the board. These
// catch the narrower, more likely failure: one construct read wrongly. Each case
// is a shape that occurs in this repository's workflows, and the first is the
// one that made the defect invisible to review.
// ---------------------------------------------------------------------------

/// The literal text `set -o pipefail` inside a **comment** does not enable it.
///
/// This is not a hypothetical corner: the block that prompted this file carried
/// exactly that comment, asserting that GitHub's default shell supplied
/// `pipefail`. A guard that greps the raw `run:` text answers "yes, it's there"
/// for the one block in the repository where it is not.
#[test]
fn a_comment_mentioning_pipefail_does_not_supply_it() {
    let scan = scan_script(
        "# `set -o pipefail` (GitHub's bash default) propagates the exit code.\n\
         python3 report.py | tee -a \"$SUMMARY\"\n",
    );
    assert_eq!(
        scan.sets_pipefail_at, None,
        "a comment must not count as a `set`"
    );
    assert_eq!(scan.pipelines, vec![2]);
}

/// `echo "set -o pipefail"` is prose, not a `set`.
#[test]
fn pipefail_inside_a_string_argument_does_not_supply_it() {
    let scan = scan_script("echo \"set -o pipefail is not on by default\" | tee log\n");
    assert_eq!(scan.sets_pipefail_at, None);
    assert_eq!(scan.pipelines, vec![1]);
}

/// The real thing, in each of the spellings this repository uses, located on
/// the line it is actually written on.
#[test]
fn a_real_set_enables_pipefail() {
    for (script, expected_line) in [
        ("set -euo pipefail\ncurl x | sh\n", 1),
        ("set -o pipefail\ncurl x | sh\n", 1),
        ("set -e\nset -o pipefail\ncurl x | sh\n", 2),
        ("if true; then set -o pipefail; fi\ncurl x | sh\n", 1),
        ("  set -o pipefail\ncurl x | sh\n", 1),
    ] {
        assert_eq!(
            scan_script(script).sets_pipefail_at,
            Some(expected_line),
            "not detected on the right line in {script:?}"
        );
    }
}

/// `set +o pipefail` turns it off and is not credited as turning it on.
#[test]
fn disabling_pipefail_is_not_enabling_it() {
    let scan = scan_script("set +o pipefail\ncurl x | sh\n");
    assert_eq!(scan.sets_pipefail_at, None);
    assert_eq!(scan.pipelines, vec![2]);
}

/// `||` is a logical OR, and `>|` a redirect. Neither is a pipeline.
#[test]
fn logical_or_and_noclobber_are_not_pipelines() {
    let scan = scan_script("a || b\nc >| out\nd ||| e\n");
    assert!(scan.pipelines.is_empty(), "{:?}", scan.pipelines);
}

/// A `|` inside quotes, inside a comment, or inside a here-document body is
/// data, not a pipeline.
#[test]
fn quoted_commented_and_heredoc_pipes_are_not_pipelines() {
    let scan = scan_script(
        "grep -E 'a|b' file\n\
         echo \"x|y\"\n\
         # a | b\n\
         cat <<'EOF'\n\
         literal | text\n\
         EOF\n\
         awk '{print $1 \"|\" $2}' file\n",
    );
    assert!(scan.pipelines.is_empty(), "{:?}", scan.pipelines);
}

/// An unquoted here-document delimiter still opens a here-document, and `<<<` is
/// a here-string that opens nothing.
#[test]
fn heredoc_state_is_tracked_in_both_spellings() {
    let unquoted = scan_script("cat <<EOF\nliteral | text\nEOF\nreal | pipe\n");
    assert_eq!(unquoted.pipelines, vec![4]);

    let herestring = scan_script("grep x <<<\"a|b\"\nreal | pipe\n");
    assert_eq!(herestring.pipelines, vec![2]);
}

/// A quoted string may span lines, and the scanner must not treat the second
/// line as code. `nightly-soak.yml` contains exactly this shape.
#[test]
fn quote_state_carries_across_lines() {
    let scan = scan_script("echo \"first | second\n third | fourth\"\nreal | pipe\n");
    assert_eq!(scan.pipelines, vec![3]);
}

/// A pipeline anywhere counts, not only as the last command — `curl … | sh` in
/// the middle of a block discards a failed download exactly as thoroughly.
#[test]
fn a_pipeline_that_is_not_the_last_command_still_counts() {
    let scan = scan_script("curl -sSf https://example/install | sh\necho done\n");
    assert_eq!(scan.pipelines, vec![1]);
}

/// A `set -o pipefail` written *after* a pipeline does not protect it, and the
/// guard reads the order rather than a block-level boolean.
///
/// Worth its own test because the boolean version of this scan — which is what
/// the file was first written with — passes on the "after" case, and passing
/// there is the same shape of wrong as the defect the file exists for: a check
/// that looks satisfied by text that does nothing.
#[test]
fn pipefail_must_be_established_before_the_pipeline() {
    let before = step_from("set -euo pipefail\ncurl x | sh\n");
    assert_eq!(before.unprotected_pipelines(), Vec::<usize>::new());
    assert!(!before.is_violation());

    let after = step_from("curl x | sh\nset -euo pipefail\necho done\n");
    assert_eq!(after.unprotected_pipelines(), vec![1]);
    assert!(after.is_violation());
}

/// A block whose resolved shell supplies `pipefail` needs no `set`, and a
/// non-shell step is not judged at all.
#[test]
fn the_shell_can_supply_what_the_block_does_not() {
    let mut supplied = step_from("curl x | sh\n");
    supplied.shell = Some("bash".to_string());
    supplied.shell_supplies_pipefail = true;
    assert_eq!(supplied.unprotected_pipelines(), Vec::<usize>::new());

    let mut not_a_shell = step_from("curl x | sh\n");
    not_a_shell.bash_family = false;
    assert_eq!(not_a_shell.unprotected_pipelines(), Vec::<usize>::new());
}

/// An unreasoned or stale exemption is caught by
/// [`every_pipefail_exemption_is_reasoned_and_live`] over the real workflows;
/// this pins the predicate the marker feeds, since nothing carries one today.
#[test]
fn an_exemption_marker_suppresses_the_violation() {
    let exempt = step_from("# pipefail-exempt: the left side cannot fail\ncurl x | sh\n");
    assert_eq!(exempt.unprotected_pipelines(), vec![2]);
    assert!(!exempt.is_violation());
}

/// A synthetic [`RunStep`] carrying `script`, resolved as GitHub's default
/// shell on a non-Windows runner — the case this whole file is about.
fn step_from(script: &str) -> RunStep {
    RunStep {
        file: "<synthetic>".to_string(),
        job: "<job>".to_string(),
        step: "<step>".to_string(),
        shell: None,
        bash_family: true,
        shell_supplies_pipefail: false,
        scan: scan_script(script),
    }
}

/// The exemption marker is read with its reason, and a bare one is
/// distinguishable from an absent one.
#[test]
fn the_exemption_marker_is_read_with_its_reason() {
    let reasoned = scan_script("# pipefail-exempt: the left side cannot fail\na | b\n");
    assert_eq!(
        reasoned.exemption.as_deref(),
        Some("the left side cannot fail")
    );

    let bare = scan_script("# pipefail-exempt:\na | b\n");
    assert_eq!(bare.exemption.as_deref(), Some(""));

    let absent = scan_script("a | b\n");
    assert_eq!(absent.exemption, None);
}

/// `shell:` resolution, in each form that decides whether `pipefail` is even a
/// question.
#[test]
fn shell_resolution_decides_where_pipefail_is_a_question() {
    // The default on a non-Windows runner: bash, WITHOUT pipefail. The subject.
    assert_eq!(classify_shell(None, false), (true, false));
    // The default on a Windows runner is `pwsh`, where the class does not exist.
    assert_eq!(classify_shell(None, true), (false, false));
    // `shell: bash` is expanded by GitHub with `-eo pipefail`.
    assert_eq!(classify_shell(Some("bash"), false), (true, true));
    // A custom template says for itself.
    assert_eq!(
        classify_shell(Some("bash --noprofile -eo pipefail {0}"), false),
        (true, true)
    );
    assert_eq!(classify_shell(Some("bash -e {0}"), false), (true, false));
    // `sh -e {0}` carries no pipefail and is in scope.
    assert_eq!(classify_shell(Some("sh"), false), (true, false));
    // Not a shell script at all.
    assert_eq!(classify_shell(Some("python"), false), (false, false));
    assert_eq!(classify_shell(Some("pwsh"), false), (false, false));
}

/// A `#` that does not start a word is not a comment — `${x#y}` and `a#b` are
/// ordinary text, and a `|` after one is still a pipeline.
#[test]
fn a_mid_word_hash_does_not_open_a_comment() {
    let scan = scan_script("echo \"${VAR#prefix}\" | tee log\n");
    assert_eq!(scan.pipelines, vec![1]);
}
