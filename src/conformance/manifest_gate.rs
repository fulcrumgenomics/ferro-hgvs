//! One place that decides what an ABSENT reference manifest means.
//!
//! Dozens of reference-aware guards need a ferro-prepared manifest — real
//! contigs, real transcripts, real cdot — which is multi-gigabyte and takes
//! ~25 min to build, so it is not provisioned for a pull request. Every one of
//! those guards therefore ends its gate the same way:
//!
//! ```ignore
//! let Some(provider) = provider() else {
//!     eprintln!("<module>: skipping — no manifest");
//!     return;
//! };
//! ```
//!
//! That is the right behaviour for a developer without a prepared reference,
//! and exactly the wrong behaviour in the one job that *does* prepare one: a
//! `ferro prepare` that produced nothing, a manifest written to a path the test
//! step cannot see, or a module that quietly stopped being selected all read as
//! a green pass. **The runtime is the only tell** — these guards take 10–70 s
//! armed and 0.02–0.08 s skipped — and no exit status, summary line or check
//! row distinguishes the two.
//!
//! That is not hypothetical. `multi_member_cis_axis::axis_normalized` was
//! provably red on `main` for at least two days across 13 merged commits while
//! every PR reported green, because the guard took 0.04 s and returned.
//!
//! So the skip is made conditional here rather than at ~50 call sites.
//! `FERRO_REQUIRE_MANIFEST=1` turns every such skip into a panic naming the
//! guard and the fix. The nightly reference-aware job — the only job that
//! prepares a manifest — sets it, so absence can no longer masquerade as a
//! pass; it stays unset locally and in PR CI, where skipping is the useful
//! default.
//!
//! This is deliberately the same shape as `tests/it/common/bulk_fixtures.rs`,
//! which solved the identical problem for the four release-asset corpora. Two
//! conventions for "a missing input must not read as a pass" would be one too
//! many.
//!
//! It lives in the library, rather than beside `bulk_fixtures` in the test
//! tree, for the same reason the rest of this module does: the guards that need
//! it are not all in one crate. Most are in the `it` integration binary, and
//! two are `#[cfg(test)]` units inside `src/` itself
//! (`accession_inventory::tests::corpus_missing_ng_against_real_reference_when_available`
//! and `ng_placement_builder::tests::derive_for_real_reference_when_available`),
//! which cannot see `tests/it/common/`. A second copy for those two is exactly
//! the drift this module exists to prevent.

/// The environment variable that promotes a manifest-absence skip to a failure.
pub const REQUIRE_ENV: &str = "FERRO_REQUIRE_MANIFEST";

/// The environment variable naming the prepared manifest itself.
pub const MANIFEST_ENV: &str = "FERRO_MANIFEST";

/// Whether an absent reference manifest must fail the guard rather than skip it.
///
/// True when `FERRO_REQUIRE_MANIFEST` is set to anything other than the empty
/// string or `0`. Read on every call rather than cached in a `OnceLock`: these
/// guards are few and slow, the read is trivially cheap beside them, and a
/// cache would make the variable unsettable from within a test — which is what
/// this module's own coverage needs.
pub fn manifest_is_required() -> bool {
    value_requires_manifest(std::env::var(REQUIRE_ENV).ok().as_deref())
}

/// The predicate itself, over an explicit value so it is testable without
/// touching the process environment.
fn value_requires_manifest(value: Option<&str>) -> bool {
    matches!(value, Some(v) if !v.is_empty() && v != "0")
}

/// Called when a reference manifest is not available: panics under
/// `FERRO_REQUIRE_MANIFEST`, otherwise reports the skip on stderr.
///
/// `context` identifies the guard that is standing down — a module or test
/// name — so the failure text says which coverage was about to be dropped.
///
/// Returns `()` rather than `!` deliberately, and for the same reason
/// `bulk_fixtures::absent` does: the call sites keep their own `return`, so
/// what a reader sees at each one is still an ordinary early exit rather than
/// control flow that depends on an environment variable.
pub fn absent(context: &str) {
    assert!(
        !manifest_is_required(),
        "{}",
        missing_manifest_message(context)
    );
    eprintln!(
        "{context}: skipping — no reference manifest \
         (set {MANIFEST_ENV}=<prepared-dir>/manifest.json)"
    );
}

/// The failure text, factored out so it can be asserted on rather than
/// re-typed in a test.
fn missing_manifest_message(context: &str) -> String {
    format!(
        "{REQUIRE_ENV} is set, but `{context}` found no reference manifest.\n\
         \n\
         This guard needs a ferro-prepared reference, which it locates through\n\
         {MANIFEST_ENV} (or `benchmark-output/manifest.json`). Without one the guard\n\
         would SKIP GREEN, silently dropping its coverage — and a skip is\n\
         indistinguishable from a pass in every signal CI reports — so under\n\
         {REQUIRE_ENV} an absent manifest is a failure instead.\n\
         \n\
         Prepare a reference with:\n\
         \n    ferro prepare --output-dir benchmark-output --genome grch38 --ensembl \\\
         \n        --derive-ng-placements tests/fixtures/mutalyzer-normalize/ng_accessions.txt\n\
         \n\
         then point {MANIFEST_ENV} at `benchmark-output/manifest.json`.\n\
         \n\
         In CI this means the nightly's `ferro prepare` step did not leave a manifest\n\
         where the test step reads it."
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The message must name the guard that stood down and the remedy. A guard
    /// whose failure text does not say what to run gets worked around rather
    /// than fixed.
    #[test]
    fn the_missing_manifest_message_names_the_guard_and_the_remedy() {
        let message = missing_manifest_message("multi_member_cis_axis::axis_normalized");
        assert!(message.contains("multi_member_cis_axis::axis_normalized"));
        assert!(message.contains("ferro prepare"));
        assert!(message.contains(REQUIRE_ENV));
        assert!(message.contains(MANIFEST_ENV));
    }

    /// The remedy must survive being **pasted**, which "contains `ferro
    /// prepare`" does not check. The command spans two lines joined by a shell
    /// `\`, and a blank line after that continuation ends the command early:
    /// the shell then runs `ferro prepare` *without* `--derive-ng-placements`
    /// — a 25-minute build producing a manifest with `ng_hosted_transcripts`
    /// unset, i.e. exactly the reference the guards under it go on to
    /// disagree with — and executes the flag line as a command of its own.
    ///
    /// So this asserts the shape a shell reads rather than the substring a
    /// reader recognises: every continued line is followed immediately by the
    /// line that continues it.
    #[test]
    fn the_remedy_command_survives_being_pasted_into_a_shell() {
        let message = missing_manifest_message("some_module::some_guard");
        let lines: Vec<&str> = message.lines().collect();
        for (index, line) in lines.iter().enumerate() {
            if !line.trim_end().ends_with('\\') {
                continue;
            }
            let next = lines
                .get(index + 1)
                .unwrap_or_else(|| panic!("line {index} continues onto nothing: {line:?}"));
            assert!(
                !next.trim().is_empty(),
                "line {index} ends with a shell continuation but line {} is blank, \
                 so a pasted command stops there:\n{message}",
                index + 1
            );
        }
        // Anti-vacuity: the loop above says nothing if the message carries no
        // continuation at all.
        assert!(
            lines.iter().any(|line| line.trim_end().ends_with('\\')),
            "no continued line found — the remedy no longer spans two lines, so this \
             guard is vacuous and should be re-aimed:\n{message}"
        );
    }

    /// The predicate's contract, pinned against the spellings CI and a
    /// developer shell actually produce. `0` and the empty string are the two
    /// ways a variable is "set but off", and reading either as *on* would make
    /// every manifest-less run in the project fail at once.
    ///
    /// Uses the private helper over an explicit value rather than the process
    /// environment: mutating the environment is a process-global side effect
    /// and tests run in threads, so a test that set the real variable could
    /// change what a concurrently running sibling sees.
    #[test]
    fn only_a_non_empty_non_zero_value_requires_the_manifest() {
        assert!(!value_requires_manifest(None));
        assert!(!value_requires_manifest(Some("")));
        assert!(!value_requires_manifest(Some("0")));
        assert!(value_requires_manifest(Some("1")));
        assert!(value_requires_manifest(Some("true")));
    }

    /// The promotion is off by default, so a developer with no prepared
    /// reference keeps the skip. If this ever reads `true` in an unset
    /// environment, every manifest-less run in the project turns red at once.
    #[test]
    fn a_manifest_is_not_required_unless_the_variable_asks_for_it() {
        if std::env::var_os(REQUIRE_ENV).is_none() {
            assert!(!manifest_is_required());
        }
    }
}
