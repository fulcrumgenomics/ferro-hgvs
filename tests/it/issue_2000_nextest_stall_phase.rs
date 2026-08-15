//! #2000: a stalled archive run must name the phase it stopped in, and be bounded.
//!
//! # What was measured, and what it refutes
//!
//! `Test (2/6)` stalled twice in three hours and was stopped by hand both times;
//! the second instance held the merge queue. Both job logs were read against a
//! healthy run of the same shard, and they stop at the same place:
//!
//! ```text
//! healthy   19:57:40  Extracted 66 files ... in 1.09s
//!           19:57:40  Starting 1794 tests across 24 binaries      <- 0.101s later
//!           ...       PASS [...] ...
//!
//! hung #1   19:57:40  Extracted 66 files ... in 0.82s
//!           21:10:01  ##[error]The operation was canceled.        <- 72 minutes, nothing between
//!
//! hung #2   21:44:29  Extracted 66 files ... in 1.08s
//!           22:05:38  ##[error]The operation was canceled.        <- 21 minutes, nothing between
//! ```
//!
//! `Starting N tests across M binaries` is printed once nextest has enumerated
//! every test binary, before any test body runs. Neither hung run reached it, so
//! **no test ever started**: the stall is in nextest's binary-list phase, in
//! which each archived binary is executed with `--list` to enumerate its tests.
//! The runner's orphan reaping says the same thing and says which binary —
//! `cargo-nextest` plus exactly one test binary in each case, and a *different*
//! one each time (`it-f2431c9a8062b86c`, then `ferro_hgvs-686e43b9c82c0b52`).
//!
//! Three things the issue proposed follow from that, and all three are refuted:
//!
//! - **A per-test `slow-timeout` with `terminate-after` cannot fire.** It governs
//!   test execution, and execution never began. (`.config/nextest.toml`'s comment
//!   already warns that a `[profile.ci]` sized by a round number kills the
//!   censuses; this is the separate reason a per-test bound is not the answer
//!   *here*.)
//! - **`--status-level` reports on tests**, of which there were none.
//! - **"Which tests land in shard 2" is the wrong question.** Listing enumerates
//!   every binary identically on every shard; `--partition` and `-E` are applied
//!   to the resulting list. No shard does less listing than any other, and the
//!   two stalls were on different binaries, so it is not a property of one
//!   binary's contents either.
//!
//! # What this guard pins
//!
//! The root cause is not established — one child of `cargo-nextest` stops
//! returning, transiently, and the next attempt runs the same tree in 91 s. What
//! *is* established is the discriminator, so the next occurrence costs one log
//! line instead of two investigations. `scripts/run_nextest_archive.sh` bounds
//! each archive run and, on a timeout, classifies it by whether that line was
//! ever printed.
//!
//! The classifier is checked against the **recorded bytes of the real runs**,
//! not against a paraphrase, because the naive pattern does not match them:
//! `ci.yml` sets `CARGO_TERM_COLOR: always`, so the line arrives as
//! `\e[32;1m    Starting\e[0m \e[1m1794\e[0m tests across ...` and a
//! `grep 'Starting [0-9]* tests across'` over the raw log returns **0 on a
//! healthy run** — which would classify every execution-phase stall as a
//! list-phase one, the exact misdiagnosis this file exists to prevent.

use std::path::PathBuf;
use std::process::Command;

/// The workflow whose archive runs this guard is about.
const CI_YML: &str = ".github/workflows/ci.yml";

/// The wrapper every archive run goes through.
const WRAPPER: &str = "scripts/run_nextest_archive.sh";

/// One line, byte-for-byte, from the healthy `Test (2/6)` of run 31905255274
/// attempt 2 (job 95070402725) — colour codes included, because they are what
/// defeats the obvious pattern.
const HEALTHY_STARTING_LINE: &str = concat!(
    "\u{1b}[32;1m    Starting\u{1b}[0m \u{1b}[1m1794\u{1b}[0m tests across ",
    "\u{1b}[1m24\u{1b}[0m binaries (\u{1b}[1m9391\u{1b}[0m tests \u{1b}[33;1mskipped\u{1b}[0m)"
);

/// The last line either hung run printed, byte-for-byte, from `Test (2/6)` of
/// run 31905255274 attempt 1 (job 95062135173). Everything after it is the
/// cancellation, 72 minutes later.
const HUNG_LAST_LINE: &str = concat!(
    "\u{1b}[32;1m   Extracted\u{1b}[0m \u{1b}[1m66\u{1b}[0m files to ",
    "\u{1b}[1m/home/runner/work/ferro-hgvs/ferro-hgvs\u{1b}[0m in 0.82s"
);

/// A log that stopped before extraction finished.
///
/// **Constructed, not recorded** — and that is worth stating, because every
/// other fixture in this file is the recorded bytes of a real run and the
/// distinction is the whole reason those are pinned. No pre-extraction stall
/// has been observed here; this shape exists so the classifier cannot announce
/// one as #2000 the first time it happens.
///
/// A cargo warning is used as the body because it is the realistic thing to
/// find in an archive run's log ahead of nextest's own first line. What the
/// fixture actually pins is weaker and does not depend on that choice: a
/// non-empty log carrying neither marker is `pre-list`.
const PRE_EXTRACTION_LINE: &str = "warning: unused manifest key: profile.soak.strip-debuginfo";

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Join `\`-continued shell lines, so a flag on a continuation line is read as
/// part of the command that owns it rather than as a command of its own.
fn fold_continuations(text: &str) -> Vec<String> {
    let mut commands = Vec::new();
    let mut pending = String::new();
    for line in text.lines() {
        let trimmed = line.trim();
        if let Some(head) = trimmed.strip_suffix('\\') {
            pending.push_str(head.trim_end());
            pending.push(' ');
        } else {
            pending.push_str(trimmed);
            commands.push(std::mem::take(&mut pending));
        }
    }
    if !pending.is_empty() {
        commands.push(pending);
    }
    commands
}

fn read(path: &str) -> String {
    let full = repo_root().join(path);
    std::fs::read_to_string(&full).unwrap_or_else(|e| panic!("{} is readable: {e}", full.display()))
}

/// Run the shipped wrapper's classifier over a log holding `body`.
fn classify(body: &str) -> String {
    let dir = tempfile::tempdir().expect("tempdir");
    let log = dir.path().join("run.log");
    std::fs::write(&log, body).expect("write log");

    let output = Command::new("bash")
        .arg(repo_root().join(WRAPPER))
        .arg("--classify")
        .arg(&log)
        .output()
        .expect("the wrapper runs");

    assert!(
        output.status.success(),
        "`{WRAPPER} --classify` exited {:?}: {}",
        output.status.code(),
        String::from_utf8_lossy(&output.stderr),
    );
    String::from_utf8_lossy(&output.stdout).trim().to_string()
}

/// The whole point, and it is a THREE-way split rather than a two-way one.
///
/// A healthy archive run prints two markers, 0.101s apart, and they delimit
/// three phases:
///
/// ```text
/// (archive read + unpack)                  <- pre-list
/// Extracted 66 files ... in 1.09s
/// (each archived binary run with --list)   <- list-phase
/// Starting 1794 tests across 24 binaries
/// (test bodies)                            <- execution-phase
/// ```
///
/// Extraction happens INSIDE the wrapped command — `--extract-to` is a
/// `cargo nextest run` flag — so a wedged archive read produces a log with
/// *neither* marker. A two-way test on `Starting` alone calls that `list-phase`
/// and the wrapper then announces it as "the #2000 shape", asks the reader to
/// file the run URL on #2000, and tells them a per-test timeout cannot see it.
/// None of that is true of a stall in extraction.
///
/// Requiring `Extracted` for the `list-phase` verdict is also what makes that
/// verdict match the incidents it names: both recorded stalls reached it.
#[test]
fn the_classifier_separates_the_three_phases_on_recorded_bytes() {
    // A healthy run prints both markers, in this order.
    assert_eq!(
        classify(&format!("{HUNG_LAST_LINE}\n{HEALTHY_STARTING_LINE}")),
        "execution-phase",
        "a run that announced its test count had started running tests",
    );
    assert_eq!(
        classify(HUNG_LAST_LINE),
        "list-phase",
        "both hung runs reached `Extracted` and never reached `Starting` — \
         that, and only that, is the #2000 shape",
    );
    assert_eq!(
        classify(PRE_EXTRACTION_LINE),
        "pre-list",
        "a log with neither marker stopped during the archive read, which is a \
         different failure from #2000 and must not be reported as one",
    );

    // The `Starting` line alone still reads as execution-phase: the phases are
    // ordered, so reaching the later marker settles it whatever precedes.
    assert_eq!(classify(HEALTHY_STARTING_LINE), "execution-phase");
}

/// `pre-list` and `unknown` are different verdicts, and neither may claim #2000.
///
/// They are adjacent — both mean "no marker was seen" — and collapsing them
/// would be tempting, but they call for different first moves: `unknown` says
/// the capture failed and the log is the thing to check, `pre-list` says the
/// capture worked and the archive read is the thing to check.
#[test]
fn an_empty_log_and_a_marker_less_log_are_told_apart() {
    assert_eq!(
        classify(""),
        "unknown",
        "an empty capture determines nothing"
    );
    assert_eq!(
        classify(PRE_EXTRACTION_LINE),
        "pre-list",
        "a capture that worked and holds no marker determines the phase",
    );
}

/// Read one arm's body out of the wrapper's `case` block.
fn annotation_for(script: &str, verdict: &str) -> String {
    let mut body = String::new();
    let mut inside = false;
    for line in script.lines() {
        let trimmed = line.trim();
        if trimmed == format!("{verdict})") {
            inside = true;
            continue;
        }
        if inside {
            if trimmed == ";;" {
                return body;
            }
            body.push_str(trimmed);
            body.push('\n');
        }
    }
    panic!("the wrapper has no `{verdict})` arm — its verdicts and its annotations must agree");
}

/// Only the `list-phase` arm may claim the #2000 shape.
///
/// The classifier returning three strings buys nothing on its own: the
/// annotation IS the deliverable, so the property is that the confident text —
/// "This is the #2000 shape", and the instruction to add the run URL to #2000 —
/// appears on the one verdict whose evidence supports it and on no other. The
/// PR argues exactly this for the `unknown` verdict; it holds equally for the
/// branch that fires most often.
#[test]
fn only_the_list_phase_annotation_claims_the_2000_shape() {
    let script = read(WRAPPER);

    let list_phase = annotation_for(&script, "list-phase");
    assert!(
        list_phase.contains("This is the #2000 shape"),
        "the list-phase arm is the one verdict whose evidence supports the \
         claim, and it must still make it: {list_phase}",
    );
    assert!(
        list_phase.contains("add the run URL to #2000"),
        "the list-phase arm must still ask for the run URL: {list_phase}",
    );
    assert!(
        list_phase.contains("Extracted"),
        "the list-phase arm should say the run got as far as extraction — that \
         is what distinguishes the two real incidents: {list_phase}",
    );

    for verdict in ["pre-list", "execution-phase", "unknown"] {
        let arm = annotation_for(&script, verdict);
        assert!(
            !arm.contains("This is the #2000 shape"),
            "the `{verdict}` arm must not announce the #2000 shape — its \
             evidence does not support it: {arm}",
        );
        assert!(
            !arm.contains("add the run URL to #2000"),
            "the `{verdict}` arm must not ask the reader to file on #2000, \
             which would pollute the issue this script exists to serve: {arm}",
        );
    }
}

/// An unreadable or empty log is `unknown`, never `list-phase`.
///
/// Without a third verdict a missing capture grep-misses and so reads as "never
/// reached the `Starting` line" — so a `mktemp` or `tee` failure would be
/// reported as the exact shape this file exists to identify, fabricating a #2000
/// with no run behind it. The annotation is the whole deliverable of the change,
/// so it must never be more confident than its evidence.
#[test]
fn a_missing_or_empty_log_is_undetermined_rather_than_a_false_positive() {
    assert_eq!(
        classify(""),
        "unknown",
        "an empty capture determines nothing, and must not read as the #2000 shape",
    );

    // A log that does not exist at all — the `mktemp` failure mode.
    let dir = tempfile::tempdir().expect("tempdir");
    let absent = dir.path().join("never-written.log");
    let output = Command::new("bash")
        .arg(repo_root().join(WRAPPER))
        .arg("--classify")
        .arg(&absent)
        .output()
        .expect("the wrapper runs");
    assert_eq!(
        String::from_utf8_lossy(&output.stdout).trim(),
        "unknown",
        "a capture that was never written determines nothing either",
    );
}

/// The wrapper runs under bash 3.2, which is `/bin/bash` on macOS.
///
/// The trap this forbids: expanding an **empty** array under `set -u` is an
/// `unbound variable` error on 3.2. The natural way to make the `timeout(1)`
/// fallback optional is `RUNNER=(timeout …)` / `RUNNER=()` plus
/// `"${RUNNER[@]}"` — and the empty branch would then die on precisely the
/// platform the fallback exists for. `shellcheck` does not flag it. Measured on
/// `/bin/bash` 3.2.57 before this guard was written.
#[test]
fn the_wrapper_avoids_array_expansion_for_bash_3_2() {
    let script = read(WRAPPER);
    let offenders: Vec<&str> = script
        .lines()
        .filter(|line| !line.trim_start().starts_with('#'))
        .filter(|line| line.contains("[@]") || line.contains("[*]"))
        .collect();
    assert!(
        offenders.is_empty(),
        "expanding a possibly-empty array under `set -u` is an unbound-variable \
         error on bash 3.2 (`/bin/bash` on macOS) — use a function or a direct \
         branch instead: {offenders:#?}",
    );

    // `$@` itself is always safe (it is empty-expandable in every bash), and the
    // script must still forward its arguments — so this guard cannot be
    // satisfied by dropping argument forwarding altogether.
    assert!(
        script.contains("cargo nextest run \"$@\""),
        "the wrapper must still forward its arguments to nextest",
    );
}

/// `--classify` with no path is a usage error, not a raw `set -u` diagnostic.
#[test]
fn classify_without_a_path_reports_usage() {
    let output = Command::new("bash")
        .arg(repo_root().join(WRAPPER))
        .arg("--classify")
        .output()
        .expect("the wrapper runs");

    assert_eq!(output.status.code(), Some(2), "usage errors exit 2");
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("usage:"),
        "a missing argument must print usage, not `$2: unbound variable`",
    );
}

/// The colour codes are the trap, so pin that the classifier survives them
/// rather than merely that it works on the two lines above.
///
/// `CARGO_TERM_COLOR: always` is set at `ci.yml`'s top-level `env:`, so every
/// archive run's log is coloured. A classifier that read the raw bytes would
/// call the healthy line `list-phase`, which is a false #2000 report on every
/// genuine in-test hang.
#[test]
fn the_classifier_is_not_defeated_by_colour_codes() {
    assert!(
        HEALTHY_STARTING_LINE.contains('\u{1b}'),
        "the fixture must carry the colour codes it exists to test",
    );
    assert!(
        !HEALTHY_STARTING_LINE.contains("Starting 1794 tests"),
        "the fixture must be the interleaved form — if the colour codes ever stop \
         splitting `Starting` from its count, this guard is testing nothing",
    );
    assert!(
        read(CI_YML).contains("CARGO_TERM_COLOR: always"),
        "the colour codes above are a consequence of this setting; if it goes, \
         re-derive the fixtures rather than deleting this assertion",
    );

    // A plain-text run must classify the same way, so the ANSI handling cannot
    // be the only thing making the test above pass.
    assert_eq!(
        classify("    Starting 10 tests across 3 binaries"),
        "execution-phase",
    );
}

/// nextest pluralises its own announcement, so a one-test run prints
/// `Starting 1 test across 1 binary`.
///
/// This is not hypothetical tidiness: the first version of the classifier
/// required the literal `tests across`, and an end-to-end run of a single-test
/// selection was reported as a list-phase stall — a false #2000 on a run that had
/// plainly started. The singular form is what a bisected or `-E`-narrowed re-run
/// produces, which is precisely when someone is reading this output to decide
/// whether they are looking at #2000 or at a stuck test.
#[test]
fn the_classifier_reads_nextests_singular_announcement() {
    for (line, expected) in [
        (
            "    Starting 1 test across 1 binary (6538 tests skipped)",
            "execution-phase",
        ),
        (
            "    Starting 4 tests across 1 binary (6535 tests skipped)",
            "execution-phase",
        ),
        // Both observed verbatim from local runs of this very wrapper.
        (
            "   Extracted 66 files to /home/runner/work in 0.82s",
            "list-phase",
        ),
    ] {
        assert_eq!(classify(line), expected, "classifying: {line}");
    }
}

/// The temp-file template's `X`s must be the last characters of the name.
///
/// `mktemp .../nextest-run.XXXXXX.log` fails outright on GNU and BusyBox
/// (`mktemp: : Invalid argument`), so every archive step would have died on this
/// line — and it is invisible locally, because BSD `mktemp` instead returns the
/// template *verbatim*, handing concurrent runs one shared fixed filename. A
/// green macOS run therefore says nothing about the platform CI uses.
#[test]
fn the_temp_template_ends_in_its_placeholder() {
    let script = read(WRAPPER);
    // The comment above the command also says "mktemp", so match the command:
    // a non-comment line carrying both the call and a placeholder.
    let template = script
        .lines()
        .find(|line| {
            let code = line.trim_start();
            !code.starts_with('#') && code.contains("mktemp") && code.contains("XXX")
        })
        .expect("the wrapper makes a temp log");
    let start = template.find("XXX").expect("the template has placeholders");
    let tail = &template[start..];
    let placeholder_end = tail.find(|c| c != 'X').unwrap_or(tail.len());
    assert!(
        tail[placeholder_end..].starts_with('"') || tail[placeholder_end..].starts_with('\''),
        "the X placeholder must end the template — a suffix after it is an error \
         on GNU/BusyBox and a silently fixed filename on BSD: {template}",
    );
}

/// Every workflow file in the repository.
///
/// The population is walked rather than listed for the same reason the archive
/// runs inside it are: a workflow added later must be covered without anyone
/// remembering to name it here.
fn workflow_files() -> Vec<(String, String)> {
    let dir = repo_root().join(".github/workflows");
    let mut files: Vec<(String, String)> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("{} is readable: {e}", dir.display()))
        .map(|entry| entry.expect("dir entry").path())
        .filter(|path| {
            path.extension()
                .is_some_and(|ext| ext == "yml" || ext == "yaml")
        })
        .map(|path| {
            let name = path
                .file_name()
                .expect("workflow has a file name")
                .to_string_lossy()
                .into_owned();
            let body = std::fs::read_to_string(&path)
                .unwrap_or_else(|e| panic!("{} is readable: {e}", path.display()));
            (name, body)
        })
        .collect();
    files.sort();
    assert!(
        !files.is_empty(),
        "no workflow files found under .github/workflows — every guard below \
         would pass vacuously",
    );
    files
}

/// Every archive run in **every workflow** goes through the wrapper.
///
/// A bare `cargo nextest run --archive-file` is one that can stall unbounded and
/// unattributed — the state both incidents were in. The list phase is identical
/// in every archive-consuming job, so scoping this to the `test` job would leave
/// `test-oracle`, `soak`, `sweeps` and `censuses` able to reproduce it silently.
///
/// # The set is derived, never restated
///
/// An earlier revision asserted `wrapped >= 6`, which is the shape `CLAUDE.md`'s
/// "Assert the property. Measure the count. Never let a count BE the property"
/// names as a change detector for the literal rather than a guard — and its own
/// failure message told the reader to edit the number, which that section calls
/// "the defect, not a mitigation".
///
/// So the population is computed instead. `--archive-file` appears on two kinds
/// of line, and the split is what makes the property statable: `nextest archive`
/// **creates** an archive and cannot stall in a list phase it never runs, while
/// every other occurrence **consumes** one to run tests and must be wrapped.
/// Adding or removing a job moves the derived set with it.
///
/// # …and the file set is derived too
///
/// This read `ci.yml` alone until `nightly-soak.yml`'s soak step was found
/// unwrapped by review — a seventh archive-consuming run, in the worse place to
/// have one, since nobody watches a green nightly and that job carries no
/// `timeout-minutes` at all. The guard's own doc boasted that "the population is
/// computed instead", and it was: over one file. Computing it over one file is
/// how the seventh got in, so the walk is now over the directory.
#[test]
fn no_archive_run_bypasses_the_wrapper() {
    // Fold `\`-continued shell lines into one command first. `--archive-file`
    // is a continuation of `cargo nextest archive` in the soak-build step, so a
    // line-based read sees a bare `--archive-file nextest-soak.tar.zst` with no
    // verb on it and files an archive *creation* as an unwrapped run. That was a
    // live false positive of this guard, not a hypothetical one.
    let mut consuming: Vec<String> = Vec::new();
    for (name, body) in workflow_files() {
        for cmd in fold_continuations(&body) {
            if cmd.starts_with('#') || !cmd.contains("--archive-file") {
                continue;
            }
            if cmd.contains("nextest archive") {
                continue;
            }
            consuming.push(format!("{name}: {cmd}"));
        }
    }

    assert!(
        !consuming.is_empty(),
        "no archive-consuming nextest run found in any workflow — the guard would \
         pass vacuously. Either the jobs were removed (update this test in the \
         same change) or the line shape moved and this filter no longer sees them",
    );

    let bypassing: Vec<&String> = consuming
        .iter()
        .filter(|cmd| !cmd.contains(WRAPPER))
        .collect();
    assert!(
        bypassing.is_empty(),
        "{} of {} archive-consuming runs call nextest directly, so a stall in them \
         is unbounded and unattributed — route them through {WRAPPER}: {bypassing:#?}",
        bypassing.len(),
        consuming.len(),
    );
}

/// GitHub's default job bound when a job declares no `timeout-minutes`.
const GITHUB_DEFAULT_JOB_MINUTES: u64 = 360;

/// The wrapper's default stall bound, in minutes, read out of the wrapper.
fn wrapper_default_minutes(script: &str) -> u64 {
    let marker = "${NEXTEST_STALL_TIMEOUT:-";
    let start = script.find(marker).unwrap_or_else(|| {
        panic!("{WRAPPER} must default NEXTEST_STALL_TIMEOUT with `{marker}…}}`")
    }) + marker.len();
    let value = script[start..]
        .split('}')
        .next()
        .expect("the default is brace-terminated")
        .trim();

    let (digits, unit) = value.split_at(
        value
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(value.len()),
    );
    let n: u64 = digits
        .parse()
        .unwrap_or_else(|e| panic!("{WRAPPER}'s default `{value}` is not a number: {e}"));
    match unit {
        "m" => n,
        "h" => n * 60,
        "s" => {
            assert_eq!(
                n % 60,
                0,
                "a sub-minute default `{value}` cannot be compared here"
            );
            n / 60
        }
        other => panic!("{WRAPPER}'s default `{value}` has an unhandled unit `{other}`"),
    }
}

/// Every job carrying a wrapped run, with the job bound it runs under.
///
/// Returns `(workflow, job, effective_bound_minutes, declared)`.
fn wrapped_jobs() -> Vec<(String, String, u64, bool)> {
    let mut found = Vec::new();
    for (file, body) in workflow_files() {
        let mut in_jobs = false;
        let mut job: Option<String> = None;
        let mut declared: Option<u64> = None;
        let mut wrapped = false;

        // Close the job we were accumulating, recording it if it was wrapped.
        let mut flush =
            |job: &mut Option<String>, declared: &mut Option<u64>, wrapped: &mut bool| {
                if let Some(name) = job.take() {
                    if *wrapped {
                        found.push((
                            file.clone(),
                            name,
                            declared.unwrap_or(GITHUB_DEFAULT_JOB_MINUTES),
                            declared.is_some(),
                        ));
                    }
                }
                *declared = None;
                *wrapped = false;
            };

        for line in body.lines() {
            if line == "jobs:" {
                in_jobs = true;
                continue;
            }
            if !in_jobs {
                continue;
            }
            // A job key is the only thing at indent 2 with nothing after its colon.
            let is_job_key = line.starts_with("  ")
                && !line.starts_with("   ")
                && line.trim_end().ends_with(':')
                && line
                    .trim()
                    .trim_end_matches(':')
                    .chars()
                    .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-');
            if is_job_key {
                flush(&mut job, &mut declared, &mut wrapped);
                job = Some(line.trim().trim_end_matches(':').to_string());
                continue;
            }
            if let Some(rest) = line.trim().strip_prefix("timeout-minutes:") {
                declared = Some(
                    rest.trim()
                        .parse()
                        .unwrap_or_else(|e| panic!("{file}: unparseable timeout-minutes: {e}")),
                );
            }
            if line.contains(WRAPPER) {
                wrapped = true;
            }
        }
        flush(&mut job, &mut declared, &mut wrapped);
    }
    found
}

/// The step bound must fire before the job bound, in every wrapped job.
///
/// # Why this is asserted rather than noted
///
/// The whole deliverable is the annotation, and the annotation only ever prints
/// if the **step** bound is what stops the run. If a job bound fires first the
/// job dies with a bare cancellation and no annotation — which is exactly the
/// state #2000 opened in, restored with nothing going red.
///
/// The two numbers live in two files and nothing related them. This PR's own
/// description reasons about #1994 as a future event ("if both land, the step
/// bound fires first"); #1994 has since merged and set every wrapped `ci.yml`
/// job to `timeout-minutes: 30`, so 25 < 30 holds today — by coincidence of
/// authorship, not by construction. Drop any wrapped job to 20 for an unrelated
/// reason and #2000's blindness comes back silently.
///
/// A job that declares no bound inherits GitHub's 360-minute default, which is
/// the real bound in force, so that is what it is compared against —
/// `nightly-soak.yml`'s soak job is in that position deliberately.
#[test]
fn the_step_bound_fires_before_every_job_bound() {
    let default = wrapper_default_minutes(&read(WRAPPER));
    let jobs = wrapped_jobs();

    assert!(
        !jobs.is_empty(),
        "no job containing a wrapped run was found — this guard would pass \
         vacuously. Either the wrapper is no longer referenced from any \
         workflow, or the job-key shape moved and the parser no longer sees it",
    );

    // Without at least one declared bound the comparison degenerates to
    // `default < 360`, which is true of almost any value and would keep passing
    // while the property it stands for had quietly stopped being tested.
    assert!(
        jobs.iter().any(|(_, _, _, declared)| *declared),
        "no wrapped job declares a `timeout-minutes`, so every comparison below \
         is against GitHub's {GITHUB_DEFAULT_JOB_MINUTES}-minute default and this \
         guard has stopped measuring anything: {jobs:#?}",
    );

    let violations: Vec<String> = jobs
        .iter()
        .filter(|(_, _, bound, _)| default >= *bound)
        .map(|(file, job, bound, declared)| {
            let source = if *declared {
                "declared"
            } else {
                "GitHub default"
            };
            format!("{file}:{job} bound {bound}m ({source}) <= step bound {default}m")
        })
        .collect();

    assert!(
        violations.is_empty(),
        "NEXTEST_STALL_TIMEOUT defaults to {default}m, which must be strictly less \
         than every wrapped job's bound — otherwise the job dies first, with a bare \
         cancellation and no phase annotation, which is the state #2000 opened in. \
         Either raise the job bound or lower {WRAPPER}'s default: {violations:#?}",
        violations = violations,
    );
}

/// The wrapper bounds the run, and the bound is nowhere near the baseline.
///
/// Measured baseline for `Test (2/6)` over the 16 most recent completed runs:
/// 70–94 s. The two stalls ran 21 and 73 minutes. Anything in that gap ends a
/// stall in minutes while leaving a legitimate slow run untouched; what must not
/// happen is falling back on GitHub's 6-hour default, which is what held the
/// merge queue.
///
/// # Comment lines are filtered, and the assertions are on constructs
///
/// An earlier revision asserted `script.contains("timeout")`, `"124"` and
/// `"137"` over the **whole file**, comments included — and this script carries
/// a comment block explaining each of those three things. Deleting the
/// `timeout …` invocation *and* the entire `if [ "$status" -eq 124 ] …` block —
/// i.e. both halves of what this test is named after — left all three
/// assertions passing on the surviving prose. Verified by doing it.
///
/// That is the direction-of-containment trap: the assertion was true of a
/// superset that includes the documentation of the behaviour, so it detected
/// the comment and not the code. Both siblings in this file already filter
/// comments; this one now does too, and asserts on the constructs rather than
/// on bare digits that any prose could supply.
#[test]
fn the_wrapper_bounds_the_run() {
    let script = read(WRAPPER);
    let code: Vec<&str> = script
        .lines()
        .filter(|line| !line.trim_start().starts_with('#'))
        .collect();

    assert!(
        code.iter()
            .any(|line| line.contains("timeout ") && line.contains("--kill-after")),
        "{WRAPPER} must bound the run with a real `timeout … --kill-after …` \
         invocation on a non-comment line — an unbounded run falls back on \
         GitHub's 6-hour default, which is what blocked the merge queue in #2000",
    );

    // GNU `timeout` reports 124; with `--kill-after` a hard kill reports 137.
    // Missing the second reads a killed stall as an ordinary test failure.
    for status in ["124", "137"] {
        assert!(
            code.iter()
                .any(|line| line.contains("$status") && line.contains(status)),
            "{WRAPPER} must test `$status` against {status} on a non-comment \
             line — without it a timed-out run is reported as an ordinary test \
             failure and the phase annotation never prints",
        );
    }
}

/// The wrapper is executable, or every step calling it dies on `Permission
/// denied` — a failure that arrives as a red CI run on an unrelated PR.
#[test]
fn the_wrapper_is_executable() {
    let path = repo_root().join(WRAPPER);
    assert!(path.is_file(), "{WRAPPER} exists");

    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mode = std::fs::metadata(&path)
            .expect("metadata")
            .permissions()
            .mode();
        assert!(
            mode & 0o111 != 0,
            "{WRAPPER} must be executable (mode {mode:o}); `git update-index --chmod=+x` it",
        );
    }
}
