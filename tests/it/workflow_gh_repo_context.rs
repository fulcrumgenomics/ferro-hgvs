//! Every repo-scoped `gh` call in a workflow must be able to name its repo.
//!
//! `gh issue`, `gh pr`, `gh run` and friends operate on "the current
//! repository". With no `--repo` flag and no `GH_REPO` in the environment, `gh`
//! works out which one that is by shelling out to `git` in the working
//! directory — so a job that never runs `actions/checkout` has nothing to
//! resolve against and the command dies before it does anything:
//!
//! ```text
//! failed to run git: fatal: not a git repository (or any of the parent directories): .git
//! ```
//!
//! **This is the guard for a defect that shipped and stayed invisible.**
//! `report-failure.yml` is the reusable job six workflows call to open a
//! tracking issue when they break after merge — the alerting that exists
//! precisely because those workflows cannot fail a pull request. It had no
//! checkout and no `GH_REPO`, so it exited 1 on every invocation without filing
//! anything, and the only symptom was a second red X underneath the red X it was
//! supposed to be reporting. A reporter that fails exactly when it is needed
//! looks identical to one that is never needed.
//!
//! Note what a test of `report-failure.yml` alone would have been worth: nothing
//! much, since a guard written against one known-broken file passes the moment
//! that file is fixed and says nothing about the next one. The class is what is
//! checked here — six other workflows shell out to `gh`, and two of them
//! (`audit-pr-runs`, `prune-ci-caches`) are correct today only because they use
//! `gh api` with the repo already interpolated into the path.
//!
//! Parsed rather than grepped, for the reason `generator_completeness.rs` gives
//! about `Cargo.toml`: the question is a property of the workflow's structure —
//! which steps sit in which job, and which `env` blocks are in scope for them —
//! and a substring scan cannot answer it. It would also match `gh` commands
//! quoted inside a `with:` input, which is real text in `audit-pr-runs.yml`'s
//! reproduction hint and is documentation rather than a command anything runs.

use serde_yaml::Value;
use std::path::PathBuf;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// `gh` subcommands that act on "the current repository" and therefore need one
/// named.
///
/// Deliberately a list of the repo-scoped command *groups* rather than every
/// leaf: `gh issue list` and `gh issue create` fail identically, and enumerating
/// leaves would leave the next one added silently unguarded. `gh api` is not
/// here because an absolute path (`/repos/OWNER/REPO/...`) needs no resolution —
/// it is handled separately, since the placeholder form does.
///
/// The list is deliberately wider than what the workflows currently call —
/// measured, they use `issue`, `api`, `pr`, `run` and `release` — because a name
/// missing here fails in the silent direction: a future checkout-less step using
/// it passes the scan and dies at runtime, which is the defect this file exists to
/// catch. `secret`, `variable` and `ruleset` are the likeliest additions (a release
/// or publish job setting one), so they are listed before anyone needs them.
const REPO_SCOPED: &[&str] = &[
    "issue", "pr", "release", "run", "workflow", "cache", "label", "browse", "secret", "variable",
    "ruleset",
];

fn workflow_files() -> Vec<PathBuf> {
    let dir = repo_root().join(".github/workflows");
    let mut files: Vec<PathBuf> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("read {}: {e}", dir.display()))
        .map(|entry| entry.expect("readable dir entry").path())
        .filter(|path| {
            matches!(
                path.extension().and_then(|e| e.to_str()),
                Some("yml") | Some("yaml")
            )
        })
        .collect();
    files.sort();
    assert!(
        !files.is_empty(),
        ".github/workflows holds no workflow files; this guard would pass \
         vacuously, which is the failure mode it exists to prevent"
    );
    files
}

/// A `run:` script reduced to the shell lines worth scanning.
///
/// Two normalisations, both load-bearing. Comment lines are dropped, because
/// `report-failure.yml` explains its own `gh` usage in a comment directly above
/// the call and that prose must not be read as a second invocation. Backslash
/// continuations are joined, because `--repo` habitually lands on a later line
/// than the subcommand — treating those as separate lines would report a
/// correctly-flagged call as unflagged.
fn shell_lines(run: &str) -> Vec<String> {
    let mut lines: Vec<String> = Vec::new();
    let mut pending = String::new();
    for raw in run.lines() {
        let line = raw.trim();
        if pending.is_empty() && line.starts_with('#') {
            continue;
        }
        if let Some(head) = line.strip_suffix('\\') {
            pending.push_str(head.trim_end());
            pending.push(' ');
            continue;
        }
        pending.push_str(line);
        lines.push(std::mem::take(&mut pending));
    }
    if !pending.is_empty() {
        lines.push(pending);
    }
    lines
}

/// One shell line split at the operators that end a command.
///
/// `--repo` is a property of an invocation, not of a line, so a line carrying two
/// commands has to be judged in two pieces: `gh pr view --repo x/y && gh issue
/// create` supplies a repository to the first and nothing to the second, and a
/// whole-line `contains("--repo ")` clears both. Splitting on `&&`, `||`, `;` and
/// `|` is coarse — it does not understand quoting — but it errs toward *more*
/// scrutiny, which is the safe direction for a guard: a split inside a quoted
/// string can only produce a piece that fails the check, never one that passes it
/// wrongly.
fn shell_segments(line: &str) -> Vec<&str> {
    line.split("&&")
        .flat_map(|piece| piece.split("||"))
        .flat_map(|piece| piece.split(';'))
        .flat_map(|piece| piece.split('|'))
        .map(str::trim)
        .filter(|piece| !piece.is_empty())
        .collect()
}

/// The subcommand of each `gh` invocation on one shell line.
///
/// A line can carry more than one (`gh issue comment … && gh issue edit …`), so
/// every occurrence is returned rather than the first. `gh` must start a command
/// — preceded by nothing, whitespace, or a shell operator — so that a word
/// merely *ending* in "gh" does not match.
fn gh_subcommands(line: &str) -> Vec<&str> {
    let bytes = line.as_bytes();
    let mut found = Vec::new();
    for (at, _) in line.match_indices("gh ") {
        let starts_command = at == 0
            || matches!(
                bytes[at - 1],
                b' ' | b'\t' | b'|' | b'(' | b'&' | b';' | b'$' | b'`'
            );
        if !starts_command {
            continue;
        }
        if let Some(word) = line[at + 3..].split_whitespace().next() {
            found.push(word);
        }
    }
    found
}

/// Whether an invocation resolves its own repo, either by flag or by not needing
/// one.
fn names_its_repo(line: &str, subcommand: &str) -> bool {
    if line.contains("--repo ") || line.contains("--repo=") || line.contains(" -R ") {
        return true;
    }
    // `gh api` substitutes `{owner}`/`{repo}` from the resolved repository, so
    // the placeholder form needs context exactly as much as `gh issue` does;
    // a path with the repo already written into it needs none.
    if subcommand == "api" {
        return !line.contains("{owner}") && !line.contains("{repo}");
    }
    !REPO_SCOPED.contains(&subcommand)
}

/// Whether `GH_REPO` is set anywhere in scope for a step, to a value that could
/// actually resolve a repository.
///
/// Presence alone is not enough: `GH_REPO: ""` sets the key and leaves `gh` with
/// nothing to resolve, so a guard keyed on `contains_key` would call the very
/// state it exists to catch a pass. An unexpanded expression (`${{ … }}`) is
/// accepted — its value is not knowable from the YAML, and refusing it would
/// reject the fix this file was written to pin.
fn env_sets_gh_repo(scopes: &[Option<&Value>]) -> bool {
    scopes.iter().flatten().any(|env| {
        env.as_mapping()
            .and_then(|m| m.get(Value::String("GH_REPO".into())))
            .is_some_and(|value| match value {
                Value::String(s) => !s.trim().is_empty(),
                // A non-string scalar is odd here but is not empty.
                other => !other.is_null(),
            })
    })
}

#[test]
fn repo_scoped_gh_calls_can_resolve_their_repository() {
    let mut checked = 0usize;
    let mut offences: Vec<String> = Vec::new();

    for path in workflow_files() {
        let text = std::fs::read_to_string(&path).expect("workflow is readable");
        let doc: Value =
            serde_yaml::from_str(&text).unwrap_or_else(|e| panic!("parse {}: {e}", path.display()));
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("?");
        let workflow_env = doc.get("env");

        let Some(jobs) = doc.get("jobs").and_then(Value::as_mapping) else {
            continue;
        };
        for (job_id, job) in jobs {
            let job_id = job_id.as_str().unwrap_or("?");
            let job_env = job.get("env");
            // A job that delegates to a reusable workflow has no steps of its
            // own; the callee is scanned as its own file.
            let Some(steps) = job.get("steps").and_then(Value::as_sequence) else {
                continue;
            };
            // Tracked in step ORDER, not as `steps.iter().any(...)`. A checkout
            // later in the job has not run yet when an earlier step shells out to
            // `gh`, so an aggregate would clear exactly the arrangement that
            // fails: `gh` first, checkout after. Measured: no job is in that
            // state today, so this closes it before it happens rather than after.
            let mut checked_out = false;

            for step in steps {
                if step
                    .get("uses")
                    .and_then(Value::as_str)
                    .is_some_and(|u| u.starts_with("actions/checkout"))
                {
                    checked_out = true;
                    continue;
                }
                let Some(run) = step.get("run").and_then(Value::as_str) else {
                    continue;
                };
                let scopes = [workflow_env, job_env, step.get("env")];
                for line in shell_lines(run) {
                    // Per INVOCATION, not per line: `gh pr view --repo x/y && gh
                    // issue create` carries `--repo` on the line but not on the
                    // second command, and judging the line would clear it.
                    for segment in shell_segments(&line) {
                        for subcommand in gh_subcommands(segment) {
                            if names_its_repo(segment, subcommand) {
                                continue;
                            }
                            checked += 1;
                            if checked_out || env_sets_gh_repo(&scopes) {
                                continue;
                            }
                            offences.push(format!(
                                "{name} / job `{job_id}`: `gh {subcommand}` has no `--repo`, \
                                 no `GH_REPO` in scope, and the repository is not checked out \
                                 at that point in the job, so `gh` cannot resolve a target \
                                 and the step exits 1"
                            ));
                        }
                    }
                }
            }
        }
    }

    // A guard that scanned nothing would pass, and would keep passing after a
    // refactor moved every `gh` call out of reach of the scan.
    assert!(
        checked > 0,
        "no repo-scoped `gh` invocations were found in any workflow. Either they \
         are all flagged with `--repo` now — in which case relax this assertion \
         deliberately — or the scan stopped matching the shape they are written in"
    );
    assert!(
        offences.is_empty(),
        "workflow `gh` calls cannot resolve their repository:\n  {}",
        offences.join("\n  ")
    );
}

/// The reporter is the case this guard was written for, so it is pinned by name.
///
/// The scan above is the general claim and would go green if `report-failure.yml`
/// were deleted. It is the one workflow whose whole purpose is to be reachable
/// when something else has already failed, and it is called by six others, so
/// losing it silently is worth its own assertion.
#[test]
fn the_failure_reporter_names_its_repository() {
    let path = repo_root().join(".github/workflows/report-failure.yml");
    let text = std::fs::read_to_string(&path).expect("report-failure.yml is readable");
    let doc: Value = serde_yaml::from_str(&text).expect("report-failure.yml parses");

    let step = doc
        .get("jobs")
        .and_then(|j| j.get("report"))
        .and_then(|j| j.get("steps"))
        .and_then(Value::as_sequence)
        .and_then(|steps| steps.iter().find(|s| s.get("run").is_some()))
        .expect("the `report` job still has a `run` step");

    let env = step.get("env").expect("the `run` step sets env");
    for key in ["GH_TOKEN", "GH_REPO"] {
        assert!(
            env.get(key).is_some(),
            "`report-failure.yml` must set `{key}`: without a token it cannot \
             authenticate and without `GH_REPO` it cannot resolve the repository, \
             and either way it fails exactly when a caller has already failed"
        );
    }
}

/// `GH_REPO` present but empty must not satisfy the scan.
///
/// Keyed on `contains_key`, the scan called `GH_REPO: ""` a pass — the exact
/// state it exists to catch, since `gh` has nothing to resolve from it and the
/// step dies the same way the unset case did. Pinned on synthetic YAML because
/// the real workflow cannot hold the broken value and stay green.
#[test]
fn an_empty_gh_repo_does_not_satisfy_the_scan() {
    let cases = [
        ("GH_REPO: \"\"", false),
        ("GH_REPO: \"   \"", false),
        ("GH_REPO: fulcrumgenomics/ferro-hgvs", true),
        ("GH_REPO: ${{ github.repository }}", true),
    ];
    for (yaml, want) in cases {
        let doc: Value = serde_yaml::from_str(yaml).expect("case parses");
        assert_eq!(
            env_sets_gh_repo(&[Some(&doc)]),
            want,
            "`{yaml}` should {} a resolvable repository",
            if want { "supply" } else { "NOT supply" }
        );
    }

    // And the absent case, which is what the whole guard was written for.
    let doc: Value = serde_yaml::from_str("GH_TOKEN: x").expect("case parses");
    assert!(
        !env_sets_gh_repo(&[Some(&doc)]),
        "an env block with no GH_REPO must not satisfy the scan"
    );
}

/// A line carrying two commands is judged per command, not per line.
///
/// Both of these were latent rather than live — measured, no workflow has a line
/// with two `gh` calls in a `run:`, and no job runs `gh` before a later checkout
/// — so they are pinned on synthetic input. A guard whose whole purpose is to
/// catch the *next* occurrence must not be satisfied by there being none today.
#[test]
fn a_repo_flag_on_one_command_does_not_clear_its_neighbour() {
    let line = "gh pr view 1 --repo owner/name && gh issue create --title x";
    let segments = shell_segments(line);
    assert_eq!(segments.len(), 2, "the line must split into two commands");

    // First command: flagged, so it resolves.
    assert!(
        names_its_repo(segments[0], "pr"),
        "`--repo` on the first command must count for the first command"
    );
    // Second command: not flagged, and must NOT inherit the first one's flag.
    assert!(
        !names_its_repo(segments[1], "issue"),
        "`gh issue create` supplies no repository and must not be cleared by the \
         `--repo` on the command before it"
    );
    // And the whole-line reading is the bug this replaces: it clears both.
    assert!(
        names_its_repo(line, "issue"),
        "sanity: judged as one line, the unflagged command is wrongly cleared — \
         which is why the scan works on segments"
    );
}
