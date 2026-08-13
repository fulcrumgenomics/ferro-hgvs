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
    absent_gated(manifest_is_required(), context);
}

/// The gate itself, over an explicit `required` so the promotion can be
/// asserted without touching the process environment — the same reason
/// [`value_requires_manifest`] takes a value rather than reading the variable.
///
/// This split is what puts the *promotion* under test rather than only its
/// message text. Asserting it through the real variable is not an option:
/// mutating the environment is a process-global side effect and tests run in
/// threads, so a test that set it could change what a concurrent sibling sees.
fn absent_gated(required: bool, context: &str) {
    assert!(!required, "{}", missing_manifest_message(context));
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

/// Called when a reference manifest IS available but does not serve a resource
/// the guard names — a transcript the prepared reference does not carry, say.
/// Panics under `FERRO_REQUIRE_MANIFEST`, otherwise reports the skip on stderr.
///
/// This is the same question as [`absent`] one level down, and it is promoted
/// by the same variable for that reason. The guard's coverage is gone either
/// way, and the promotion is about the coverage rather than about the cause:
/// there is no third state in which the manifest is present, the resource is
/// missing, and the guard still checked anything.
///
/// It is a **separate function** only because the remedy differs — "prepare a
/// reference" against "this reference does not carry `X`" — and a failure whose
/// text names the wrong fix gets worked around rather than fixed.
///
/// It matters more than it looks, because the instrument that found the guards
/// [`absent`] protects cannot see this one. Those were censused by duration:
/// 0.02–0.08 s skipped against 9–70 s armed. Here the provider load happens
/// *before* the lookup, so a stood-down guard costs the same wall clock as an
/// armed one — measured 13.17 s either way on the blessed reference. Neither
/// the exit status, the summary line, the check row, nor the runtime
/// distinguishes them.
///
/// `resource` is the thing that was not served; `reason` is the provider's own
/// error, kept verbatim so a version mismatch reads differently from an
/// accession the reference never had.
pub fn unserved(context: &str, resource: &str, reason: &str) {
    unserved_gated(manifest_is_required(), context, resource, reason);
}

/// The gate itself, over an explicit `required` — see [`absent_gated`] for why
/// the promotion is factored out rather than asserted through the variable.
fn unserved_gated(required: bool, context: &str, resource: &str, reason: &str) {
    assert!(
        !required,
        "{}",
        unserved_resource_message(context, resource, reason)
    );
    eprintln!("{context}: skipping — the reference manifest does not serve `{resource}`: {reason}");
}

/// The failure text for [`unserved`], factored out so it can be asserted on
/// rather than re-typed in a test.
fn unserved_resource_message(context: &str, resource: &str, reason: &str) -> String {
    format!(
        "{REQUIRE_ENV} is set and a reference manifest was found, but `{context}` could\n\
         not get `{resource}` from it: {reason}\n\
         \n\
         The manifest is present, so this is not the absence {REQUIRE_ENV} usually\n\
         reports — the prepared reference simply does not carry what this guard names.\n\
         The guard would SKIP GREEN, silently dropping its coverage, and unlike an\n\
         absent manifest this one does not even show up as a faster run: the provider\n\
         is loaded before the lookup, so a stood-down guard costs the same wall clock\n\
         as an armed one.\n\
         \n\
         Either the `ferro prepare` inputs no longer supply `{resource}` (an accession\n\
         withdrawn or re-versioned upstream is the common case), or the guard is\n\
         naming an accession it should no longer name. Confirm with:\n\
         \n    ferro check --reference <prepared-dir>\n\
         \n\
         Fix the reference, or re-point the guard at an accession the reference has\n\
         and say in the commit why it moved. Do not restore the skip."
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

    /// The unserved-resource message must name the guard, the resource and the
    /// provider's own reason, and must not read as the absent-manifest case —
    /// the manifest was found, so pointing at `ferro prepare --output-dir` is
    /// the wrong remedy and sends the reader to check a file that is fine.
    #[test]
    fn the_unserved_resource_message_names_the_resource_and_a_different_remedy() {
        let message = unserved_resource_message(
            "real_data_normalization_tests",
            "NM_001408491.1",
            "Reference not found: NM_001408491.1",
        );
        assert!(message.contains("real_data_normalization_tests"));
        assert!(message.contains("NM_001408491.1"));
        assert!(message.contains("Reference not found: NM_001408491.1"));
        assert!(message.contains(REQUIRE_ENV));
        assert!(message.contains("ferro check"));
        assert!(
            message.contains("does not carry"),
            "the message must say the reference lacks the resource, not that it is absent"
        );
        assert!(
            !message.contains("found no reference manifest"),
            "must not be confusable with the absent-manifest message"
        );
    }

    /// The two messages are different text for different causes. Pinned because
    /// the cheap implementation of `unserved` is to call `absent` with a
    /// concatenated context, which would silently give both causes the same
    /// remedy.
    #[test]
    fn the_two_messages_are_not_the_same_message() {
        assert_ne!(
            missing_manifest_message("m"),
            unserved_resource_message("m", "NM_1.1", "why")
        );
    }

    /// The promotion itself, which is the whole point of this module and which
    /// asserting on message *text* does not check: with the manifest required,
    /// a stand-down must panic rather than report and continue.
    ///
    /// Without this, deleting the `assert!` from `unserved` leaves every other
    /// test in this file green — the exact shape (a guard that passes without
    /// checking the thing it names) that this module exists to eliminate.
    #[test]
    #[should_panic(expected = "does not carry")]
    fn an_unserved_resource_panics_when_the_manifest_is_required() {
        unserved_gated(true, "some_module", "NM_1.1", "Reference not found: NM_1.1");
    }

    /// The other half: with the manifest not required the same call must return
    /// normally, because that is the local and PR-CI default. A `should_panic`
    /// test alone would be satisfied by a function that always panics.
    #[test]
    fn an_unserved_resource_only_reports_when_the_manifest_is_not_required() {
        unserved_gated(
            false,
            "some_module",
            "NM_1.1",
            "Reference not found: NM_1.1",
        );
    }

    /// The same pair for [`absent`]. Same defect class, same gap — pinned here
    /// so the two stand-down paths cannot drift apart.
    #[test]
    #[should_panic(expected = "found no reference manifest")]
    fn an_absent_manifest_panics_when_the_manifest_is_required() {
        absent_gated(true, "some_module");
    }

    #[test]
    fn an_absent_manifest_only_reports_when_the_manifest_is_not_required() {
        absent_gated(false, "some_module");
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
