//! CLI surface test for `ferro normalize --direction` (#1863).
//!
//! `--direction` was declared with no `value_parser`
//! (`src/bin/ferro.rs:284-286` on `64fc76df`), and the function behind it,
//! `parse_shuffle_direction` (`src/cli/parse.rs:50-55`), ended in
//! `_ => ShuffleDirection::ThreePrime`. Between them, `ferro normalize
//! --direction banana` ran a **3'** normalization and exited **0** — the user
//! got a plausible answer to a question they did not ask.
//!
//! This is the CLI half of #1016. #1017 fixed the Python boundary only: it made
//! `python_helpers::parse_direction` fallible and raised `PyValueError`, and its
//! description notes that "`parse_direction` has no non-Python (pure-Rust/CLI)
//! callers, so no core path needed changes". That is true of `parse_direction`
//! itself, and is exactly why the CLI's independent near-duplicate was never
//! brought into scope.
//!
//! **These tests must drive the binary, not the library.** The defect lives in
//! the clap declaration — the seam between the flag and the parse function — so
//! a library-level test calling `parse_shuffle_direction` directly would pass
//! against the bug and keep passing after the fix. Only the binary witnesses the
//! flag.
//!
//! No reference is needed. Argument validation happens before any provider is
//! constructed, and `ferro normalize` already runs without `--reference` (see
//! `cli_normalize_tsv.rs`), so both arms are cheap and hermetic.
//!
//! Scope note: this pins the CLI flag only. `pub fn parse_shuffle_direction`
//! keeps its lenient fallback and its signature — it is public, re-exported at
//! `src/cli/mod.rs:18`, and doctested, so tightening it is a SemVer decision
//! rather than a bug fix, and it is entangled with the open question in #1857
//! of whether the 5' direction survives at all. With clap validating upstream,
//! the fallback is no longer reachable from the CLI.

use std::process::Command;

/// Every spelling the shipped parser resolves to a direction, and must keep
/// resolving after the fix. The 3' forms are the load-bearing ones: on
/// `origin/main` they are **not** matched explicitly — they resolve only by
/// falling into the same catch-all that swallowed `banana` — so a fix that
/// merely removed the catch-all without listing them would break the default.
const VALID_DIRECTIONS: &[&str] = &[
    "3prime", "3'", "3", "5prime", "5'", "5", // case-insensitive: `parse_shuffle_direction`
    // lowercases its input, so these work today and must keep working.
    "5PRIME", "5Prime", "3PRIME", "3Prime",
];

/// Values that must be rejected. `banana` is the reported case; the rest are the
/// near-misses that motivated #1016 (a typo produces wrong-but-plausible
/// output), plus the empty string.
const INVALID_DIRECTIONS: &[&str] = &["banana", "5prim", "3prine", "five", "threeprime", ""];

/// A variant that `ferro normalize` handles with no reference directory, taken
/// from `cli_normalize_tsv.rs`.
const VARIANT: &str = "NM_000088.3:c.16del";

struct Run {
    stdout: String,
    stderr: String,
    success: bool,
}

/// Run `ferro normalize --direction <direction> <VARIANT>`.
fn normalize_with_direction(direction: &str) -> Run {
    let out = Command::new(env!("CARGO_BIN_EXE_ferro"))
        .args(["normalize", "--direction", direction, VARIANT])
        .output()
        .expect("run ferro normalize");
    Run {
        stdout: String::from_utf8_lossy(&out.stdout).into_owned(),
        stderr: String::from_utf8_lossy(&out.stderr).into_owned(),
        success: out.status.success(),
    }
}

/// The reported defect: an unrecognized `--direction` must be a hard failure,
/// not a silent 3' run.
#[test]
fn an_unrecognized_direction_is_rejected() {
    for bad in INVALID_DIRECTIONS {
        let run = normalize_with_direction(bad);
        assert!(
            !run.success,
            "`--direction {bad:?}` exited 0 — this is the #1863 bug: the flag has \
             no value_parser and `parse_shuffle_direction` falls back to ThreePrime, \
             so a bogus direction silently 3'-shifts. stdout: {:?} stderr: {:?}",
            run.stdout, run.stderr
        );
    }
}

/// The rejection must name the accepted values, so the user can fix it without
/// reading the source. This is what `--format` and `--error-mode` already do,
/// and using clap's own `value_parser` is what makes the message consistent
/// with them rather than hand-rolled.
#[test]
fn the_rejection_names_the_accepted_values() {
    let run = normalize_with_direction("banana");
    let combined = format!("{}{}", run.stdout, run.stderr);
    assert!(
        !run.success,
        "`--direction banana` must fail; got: {combined:?}"
    );
    for expected in ["3prime", "5prime"] {
        assert!(
            combined.contains(expected),
            "the rejection must name the accepted value {expected:?} so the user can \
             correct the invocation; got: {combined:?}"
        );
    }
}

/// Every spelling the shipped parser accepted must keep working. Enumerated
/// from `parse_shuffle_direction`'s match arms rather than from the `--help`
/// text, because the 3' aliases are absent from the arms and survive only via
/// the catch-all — guessing from the doc comment ("3prime or 5prime") would
/// have silently dropped `3`, `3'`, and every mixed-case form.
#[test]
fn every_documented_spelling_is_still_accepted() {
    for good in VALID_DIRECTIONS {
        let run = normalize_with_direction(good);
        assert!(
            run.success,
            "`--direction {good:?}` must still be accepted; stdout: {:?} stderr: {:?}",
            run.stdout, run.stderr
        );
    }
}

/// The sharp shape of the regression, stated as a difference rather than as an
/// exit code: before the fix, `--direction banana` produced output
/// **byte-identical** to `--direction 3prime`. That equality is the defect — it
/// is what "silently 3'-shifts" means — so it is asserted directly. An exit-code
/// check alone would keep passing if a future change made garbage fail for some
/// unrelated reason while still normalizing.
#[test]
fn an_unrecognized_direction_does_not_masquerade_as_three_prime() {
    let three_prime = normalize_with_direction("3prime");
    assert!(
        three_prime.success,
        "the control arm must succeed; stderr: {:?}",
        three_prime.stderr
    );

    let garbage = normalize_with_direction("banana");
    assert!(
        !(garbage.success && garbage.stdout == three_prime.stdout),
        "`--direction banana` produced the same successful output as \
         `--direction 3prime` ({:?}) — the user asked for a direction that does \
         not exist and was silently given 3'.",
        three_prime.stdout
    );
}
