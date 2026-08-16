//! Ties the tag release-plz will actually mint to the tag prefixes
//! `release-wheels.yml` will actually accept.
//!
//! # The failure this exists for
//!
//! v0.14.0 published to crates.io and never reached PyPI. `release-wheels.yml`
//! gated all four of its building jobs on
//! `startsWith(github.event.release.tag_name, 'v')`, every release before it had
//! been tagged `v0.13.1`, `v0.13.0`, `v0.12.0` … and v0.14.0 was tagged
//! **`ferro-hgvs-v0.14.0`**. All ten jobs skipped — the four wheel builds,
//! `check-metadata`, the two smoke matrices, `Upload to GitHub Release`,
//! `Publish to PyPI`, and `verify-pypi`, whose entire purpose is to catch a
//! publish that did not happen. `publish.yml` reported `success` throughout, so
//! the only visible signal said the release had gone fine.
//!
//! Nothing renamed the tag on purpose. release-plz derives the tag format from
//! the **workspace**: `v{version}` for a single-package workspace,
//! `{package}-v{version}` for a multi-package one. #1977 added `tests-soak` — a
//! `publish = false` test driver — as a workspace member to buy a CI
//! compile-time saving, and that made the workspace multi-package. (release-plz
//! counts a package by whether it carries release metadata, not by whether it is
//! publishable, which is why a `publish = false` member was enough.) A PR about
//! build times renamed the release tag, in a diff where no line mentioned tags.
//!
//! So the gate's expected prefix and the tag actually minted were two facts that
//! had to agree, living in two files, with nothing relating them. This test is
//! that relation.
//!
//! # What it does
//!
//! [`minted_tag`] derives the tag from the release configuration — from
//! `release-plz.toml`'s `git_tag_name` when it is pinned (it is, to
//! `v{{ version }}`), and otherwise by modelling release-plz's package-count
//! default off `Cargo.toml`. [`gate_prefixes_for`] reads the `startsWith(…,
//! '<prefix>')` literals out of **one** job's `if:` condition. The assertions
//! are that every gated job accepts the minted tag, that no gated job accepts
//! the `test-fixtures-v1` data release, that all four gated jobs carry the
//! *same* gate, and that the two jobs which can publish still require a real
//! `release` event.
//!
//! The gate is deliberately narrow — one accepted prefix, matching the pinned
//! format and nothing else. That is what makes this guard able to see the pin
//! drift: a gate widened to also accept a format the pin cannot mint stays green
//! when the pin is reverted or deleted, which is a guard that has been
//! configured not to ring.
//!
//! # Everything is read per-job, and counted
//!
//! Two misses that a whole-workflow reading has, both of which reproduce real
//! failure modes:
//!
//! - Unioning every job's prefixes into one set answers "does *some* job accept
//!   this tag", when the question is "does *each* job accept it". These four are
//!   siblings, not a chain — nothing downstream re-checks the tag, it only
//!   `needs:` these four — so a gate dropped from one of four near-identical
//!   blocks is a hole, and it is the v0.14.0 skip-chain verbatim.
//! - Reading only the `startsWith` literals ignores whatever *else* the
//!   condition says. Someone repairing a too-narrow gate does not delete the arm
//!   that works, they add breadth: `… || contains(tag_name, 'v')` is valid GHA,
//!   passes `actionlint`, admits `test-fixtures-v1`, and leaves both extracted
//!   prefixes exactly as they were. [`gate_prefixes_for`] therefore also asserts
//!   that every mention of `tag_name` in the condition is inside a `startsWith`
//!   this guard read, so an arm it cannot see cannot be added silently.
//!
//! # What it does NOT catch
//!
//! - **It does not run release-plz.** The package-count fallback is a *model* of
//!   release-plz's default, written from the behaviour observed on
//!   `ferro-hgvs-v0.14.0`. If upstream changes that default the model goes stale
//!   and this test stays green while the real tag drifts. The
//!   `git_tag_name` pin is what takes the model out of the loop — with it
//!   present the tag is read as a literal and nothing is inferred — which is a
//!   reason to keep the pin, not merely a convenience. Note the corollary: with
//!   the pin present, a future workspace change cannot fail this guard, because
//!   the package-count model is never consulted.
//! - **It does not check that the release was created with that tag.** A human
//!   tagging by hand, or a release-plz run with a different config, is outside
//!   what any test of these two files can see.
//! - **`startsWith(tag, 'v')` still admits any `v`-prefixed name.** The pinned
//!   format is all this repo will *mint*, but releases are also published here by
//!   hand — `test-fixtures-v1` is one — under no naming rule, and a future
//!   `validation-corpus-v1` would fire the whole matrix. GitHub expressions have
//!   no regex, so requiring a digit after the `v` means either ten `'v0'…'v9'`
//!   arms or hoisting the decision into a one-job `bash` regex the four jobs
//!   consume through `needs:`. The second also removes the four-way duplication
//!   above; both are follow-ups, deliberately not folded into a release-blocking
//!   fix.
//! - **It does not check the rest of the pipeline.** A release can still skip
//!   silently through a broken `needs:` chain, a downstream job-level `if:`, an
//!   unapproved deployment environment, or a missing trusted-publisher config.
//!   This guard covers exactly one seam.
//! - **It reads the gate textually.** A rewrite that expressed the same
//!   condition without a literal `startsWith(github.event.release.tag_name,
//!   '…')` would make [`gate_prefixes_for`] read nothing — which is why every
//!   extraction below asserts it found something, rather than letting an empty
//!   read pass as agreement.

use std::collections::BTreeSet;
use std::path::PathBuf;

/// The workflow whose gate decides whether wheels get built and published.
const WHEELS_YML: &str = ".github/workflows/release-wheels.yml";

/// The release automation whose configuration decides what the tag is called.
const RELEASE_PLZ_TOML: &str = "release-plz.toml";

/// The four jobs that carry the gate. Everything else in the workflow reaches
/// them through `needs:` and skips with them, which is why a wrong gate takes
/// out the whole pipeline rather than one job.
const GATED_JOBS: &[&str] = &["linux", "macos", "windows", "sdist"];

/// The two jobs that can ship something.
const PUBLISHING_JOBS: &[&str] = &["upload", "publish-pypi"];

/// The data release that hosts the bulk test corpora. It is published through
/// the same `release: published` event, carries no code, and must never trigger
/// a wheel build — it did once, and had to be cancelled by hand.
const DATA_RELEASE_TAG: &str = "test-fixtures-v1";

/// The `startsWith` call whose literal argument is a prefix this gate accepts.
const GATE_CALL: &str = "startsWith(github.event.release.tag_name, '";

/// The context the gate reads. Every mention of it in a gated job's condition
/// must be inside a [`GATE_CALL`], or the guard is reading less than the gate
/// says.
const TAG_CONTEXT: &str = "tag_name";

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn read(relative: &str) -> String {
    let path = repo_root().join(relative);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

/// Where the expected tag came from, so a failure can say which half to fix.
#[derive(Debug, PartialEq, Eq)]
enum TagSource {
    /// `release-plz.toml` pins `git_tag_name`; the tag is that template
    /// rendered, and nothing about the workspace enters into it.
    Pinned(String),
    /// No pin, so release-plz's own default applies and we have to model it from
    /// the package count.
    DerivedFromPackageCount(usize),
}

/// The tag release-plz will mint for the next release of this workspace.
fn minted_tag() -> (String, TagSource) {
    let manifest: toml::Value =
        toml::from_str(&read("Cargo.toml")).expect("Cargo.toml parses as TOML");
    let package = manifest
        .get("package")
        .expect("Cargo.toml declares a root [package]");
    let name = package
        .get("name")
        .and_then(toml::Value::as_str)
        .expect("the root package is named");
    let version = package
        .get("version")
        .and_then(toml::Value::as_str)
        .expect("the root package is versioned");

    let config: toml::Value =
        toml::from_str(&read(RELEASE_PLZ_TOML)).expect("release-plz.toml parses as TOML");

    if let Some(template) = config
        .get("workspace")
        .and_then(|w| w.get("git_tag_name"))
        .and_then(toml::Value::as_str)
    {
        // Check what the message says: the template must interpolate the
        // *version*, specifically. `{{ package }}` alone also renders to
        // something other than itself, and would give every release one tag.
        assert!(
            template.contains("{{ version }}") || template.contains("{{version}}"),
            "{RELEASE_PLZ_TOML}'s git_tag_name = {template:?} does not reference \
             {{{{ version }}}}, so every release would reuse one tag."
        );
        return (
            render(template, name, version),
            TagSource::Pinned(template.to_string()),
        );
    }

    // No pin: release-plz's default is keyed on how many packages the workspace
    // holds — `v{version}` for one, `{package}-v{version}` for more. Modelled,
    // not observed; see the module docs.
    let members = manifest
        .get("workspace")
        .and_then(|w| w.get("members"))
        .and_then(toml::Value::as_array)
        .map_or(0, Vec::len);
    let packages = 1 + members;
    let template = if packages > 1 {
        "{{ package }}-v{{ version }}"
    } else {
        "v{{ version }}"
    };
    (
        render(template, name, version),
        TagSource::DerivedFromPackageCount(packages),
    )
}

/// Render a release-plz tag template. Only the two placeholders that can appear
/// in a tag name are substituted; an unrecognised one is left in place, where it
/// shows up in the failure message rather than being silently dropped.
fn render(template: &str, package: &str, version: &str) -> String {
    template
        .replace("{{ package }}", package)
        .replace("{{package}}", package)
        .replace("{{ version }}", version)
        .replace("{{version}}", version)
}

/// Every job-level `if:` in `release-wheels.yml`, by job name.
fn job_conditions() -> Vec<(String, String)> {
    let workflow: serde_yaml::Value =
        serde_yaml::from_str(&read(WHEELS_YML)).expect("release-wheels.yml parses as YAML");
    let jobs = workflow
        .get("jobs")
        .and_then(serde_yaml::Value::as_mapping)
        .expect("release-wheels.yml declares jobs");
    jobs.iter()
        .filter_map(|(name, body)| {
            let name = name.as_str()?.to_string();
            let condition = body.get("if")?.as_str()?.to_string();
            Some((name, condition))
        })
        .collect()
}

/// The `if:` condition of one named job, or a failure naming the job.
///
/// A building job with no condition at all is the whole defect in one line: it
/// fires on every release, including a data release.
fn condition_of(job: &str) -> String {
    job_conditions()
        .into_iter()
        .find(|(name, _)| name == job)
        .map(|(_, condition)| condition)
        .unwrap_or_else(|| {
            panic!(
                "{WHEELS_YML} has no `{job}` job, or it carries no `if:` at all. Nothing \
                 downstream re-checks the release tag — it only `needs:` the building jobs — \
                 so an ungated one fires the whole matrix on any release, and a gate present \
                 on its three siblings does not cover it."
            )
        })
}

/// The tag prefixes **one** job will accept, read out of the
/// `startsWith(github.event.release.tag_name, '<prefix>')` calls in its
/// job-level condition.
///
/// Per-job on purpose: unioning the four answers "does some job accept this
/// tag", and a gate dropped from one of four near-identical blocks is exactly
/// the skip-chain that cost v0.14.0 its wheels.
fn gate_prefixes_for(job: &str, condition: &str) -> BTreeSet<String> {
    let mut prefixes = BTreeSet::new();
    let mut rest = condition;
    while let Some(at) = rest.find(GATE_CALL) {
        let after = &rest[at + GATE_CALL.len()..];
        let end = after
            .find('\'')
            .expect("a startsWith() literal closes its quote");
        prefixes.insert(after[..end].to_string());
        rest = &after[end..];
    }

    assert!(
        !prefixes.is_empty(),
        "{WHEELS_YML}'s `{job}` job does not gate on \
         `startsWith(github.event.release.tag_name, …)`: {condition:?}. Either the tag gate \
         was removed — in which case the `{DATA_RELEASE_TAG}` data release fires the whole \
         wheel matrix again — or it was rewritten in a form this guard cannot read, in which \
         case every assertion here is vacuous. Both need a human."
    );

    // Every arm this guard cannot read is an arm it cannot judge, and the arm
    // someone actually adds when loosening a gate is an *extra* one beside the
    // ones that work. Counting mentions of the context catches that, where
    // reading the literals alone does not.
    let read_calls = condition.matches(GATE_CALL).count();
    let mentions = condition.matches(TAG_CONTEXT).count();
    assert_eq!(
        mentions, read_calls,
        "{WHEELS_YML}'s `{job}` job mentions `{TAG_CONTEXT}` {mentions} times but only \
         {read_calls} of those are `{GATE_CALL}…')` calls this guard can read: {condition:?}.\n\
         An arm this guard cannot read is an arm it cannot judge — a `contains(…, 'v')` or a \
         `startsWith(…,'…')` written without the space added *beside* the working arms leaves \
         the extracted prefixes untouched, so every assertion below would pass while the gate \
         admits something new. Write the extra arm in the same form, or widen this guard \
         deliberately."
    );

    prefixes
}

/// Each gated job's accepted prefixes, in `GATED_JOBS` order.
fn gated_job_prefixes() -> Vec<(&'static str, BTreeSet<String>)> {
    GATED_JOBS
        .iter()
        .map(|job| {
            let condition = condition_of(job);
            (*job, gate_prefixes_for(job, &condition))
        })
        .collect()
}

/// Every gated job must accept the tag the release automation is actually going
/// to mint.
///
/// This is the assertion that would have failed on the v0.14.0 tree: the gate
/// held `v`, the workspace had gone multi-package, and the minted tag was
/// `ferro-hgvs-v0.14.0`.
#[test]
fn the_wheels_gate_accepts_the_tag_release_plz_will_mint() {
    let (tag, source) = minted_tag();
    let declining: Vec<(&str, BTreeSet<String>)> = gated_job_prefixes()
        .into_iter()
        .filter(|(_, prefixes)| !prefixes.iter().any(|p| tag.starts_with(p)))
        .collect();
    assert!(
        declining.is_empty(),
        "release-plz will tag the next release `{tag}` ({source:?}), and these jobs in \
         {WHEELS_YML} do not accept it: {declining:?}.\n\
         Each of those jobs would skip, and everything that `needs:` it skips with it — \
         including `verify-pypi`, which exists to catch exactly this. `publish.yml` would \
         still report success, so the release would look green while the wheels never \
         reached PyPI. That is what happened to v0.14.0.\n\
         Fix whichever half is wrong: `git_tag_name` in {RELEASE_PLZ_TOML} is the pin, and \
         the gate must match it."
    );
}

/// …and must still reject the data release, which is why the gate exists.
///
/// The obvious repair for the failure above is a looser predicate. A
/// `contains('v')` or a bare suffix match would accept `{DATA_RELEASE_TAG}` and
/// re-break the thing the gate was added for, so the two halves are asserted
/// together and neither can be satisfied alone.
#[test]
fn the_wheels_gate_still_rejects_the_data_release() {
    let accepting: Vec<(&str, Vec<String>)> = gated_job_prefixes()
        .into_iter()
        .filter_map(|(job, prefixes)| {
            let matched: Vec<String> = prefixes
                .into_iter()
                .filter(|p| DATA_RELEASE_TAG.starts_with(p))
                .collect();
            (!matched.is_empty()).then_some((job, matched))
        })
        .collect();
    assert!(
        accepting.is_empty(),
        "{WHEELS_YML} accepts the `{DATA_RELEASE_TAG}` data release: {accepting:?}. That \
         release carries the bulk test corpora and no code; firing the wheel matrix on it \
         uploads wheels onto a data release and attempts a PyPI publish. It happened once \
         and had to be cancelled by hand."
    );
}

/// All four building jobs carry the *same* gate.
///
/// They are siblings, not a chain — nothing downstream re-checks the tag, it
/// only `needs:` these four — so a gate dropped from one job, or narrowed in one
/// job, is a hole rather than a redundancy. Requiring the four prefix sets to be
/// equal is what makes hand-editing three of four near-identical blocks fail
/// here instead of at the next release: `condition_of` fails on a job that lost
/// its `if:` entirely, `gate_prefixes_for` fails on one whose gate no longer
/// reads, and this fails on one whose gate merely disagrees with its siblings.
#[test]
fn every_building_job_carries_the_tag_gate() {
    let per_job = gated_job_prefixes();
    let (first_job, expected) = per_job.first().expect("GATED_JOBS is not empty");
    let disagreeing: Vec<&(&str, BTreeSet<String>)> = per_job
        .iter()
        .filter(|(_, prefixes)| prefixes != expected)
        .collect();
    assert!(
        disagreeing.is_empty(),
        "the building jobs in {WHEELS_YML} do not all gate on the same tag prefixes. \
         `{first_job}` accepts {expected:?}, and these disagree: {disagreeing:?}.\n\
         Nothing downstream re-checks the tag, so the narrowest gate decides which jobs run \
         and the widest decides which releases reach them — a release accepted by three of \
         the four builds wheels, skips one, and takes `check-metadata` and every publishing \
         job down with the one it skipped."
    );
}

/// Widening the tag gate must not let a `workflow_dispatch` dry-run publish.
///
/// `workflow_dispatch` is deliberately outside the tag gate — the dry-run form
/// takes a blank tag — so the only thing keeping a dry run from uploading to a
/// GitHub Release or to PyPI is these two jobs' own conditions. A change to the
/// gate above cannot weaken them, but it is the change most likely to be
/// accompanied by a well-meant edit here.
///
/// A floor, and only a floor: this requires each job to *mention* a real release
/// event, not that it fires on nothing else. `upload` deliberately also runs on a
/// dispatch carrying an explicit non-empty tag — that is how the v0.14.0 rescue
/// re-attached its wheels — so the assertion cannot be stronger than a substring
/// test without pinning that second arm too.
#[test]
fn only_a_real_release_can_publish() {
    for job in PUBLISHING_JOBS {
        let condition = condition_of(job);
        assert!(
            condition.contains("github.event_name == 'release'"),
            "{WHEELS_YML}'s `{job}` job can ship artifacts and its condition does not \
             require a real release event: {condition:?}"
        );
    }
}
