//! Ties the tag release-plz will actually mint to the tag prefixes **every**
//! release-triggered workflow will actually accept.
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
//! # The second failure, which is why the workflow set is DERIVED
//!
//! The first version of this file named its subject: `const WHEELS_YML: &str =
//! ".github/workflows/release-wheels.yml"`, with the gated and publishing job
//! names hardcoded beside it. That is a roster, and a roster has the same drift
//! problem one level up from the one this guard was written to close — a second
//! release-triggered workflow is simply invisible to it.
//!
//! It did not stay hypothetical for long. #2055 added
//! `release-binaries.yml`, carrying the same `release: published` trigger and a
//! copy of the same gate, and **deleting that gate outright left all four of the
//! old assertions passing, EXIT=0**. The workflow set is therefore derived from
//! `.github/workflows/` rather than named, following the precedent
//! `spec_fixture_setup_filter.rs` (#2033) sets: derive the set, assert equality
//! in both directions, and floor it so a matcher that has gone blind cannot pass
//! as "nothing to check".
//!
//! [`release_workflows`] keeps every `.github/workflows/*.{yml,yaml}` whose `on:`
//! block declares a `release` trigger whose `types:` includes `published` — or
//! declares no `types:` at all, which in GitHub Actions means every activity
//! type. `publish.yml` triggers on `push: main` and is correctly out of scope.
//!
//! # What it does
//!
//! [`minted_tag`] derives the tag from the release configuration — from
//! `release-plz.toml`'s `git_tag_name` when it is pinned (it is, to
//! `v{{ version }}`), and otherwise by modelling release-plz's package-count
//! default off `Cargo.toml`. [`gate_prefixes_for`] reads the `startsWith(…,
//! '<prefix>')` literals out of **one** job's `if:` condition.
//!
//! The assertions are that every job a release reaches *directly* is gated, that
//! every gate accepts the minted tag, that no gate accepts the
//! `test-fixtures-v1` data release, that every gate in every release workflow
//! agrees on one prefix set, and that every job which can ship something still
//! requires a real `release` event.
//!
//! The gate is deliberately narrow — one accepted prefix, matching the pinned
//! format and nothing else. That is what makes this guard able to see the pin
//! drift: a gate widened to also accept a format the pin cannot mint stays green
//! when the pin is reverted or deleted, which is a guard that has been
//! configured not to ring.
//!
//! # A deleted gate is invisible to any check that iterates over gates
//!
//! This is the measurement that decides the shape of the file, and it is worth
//! stating in one line because it is not obvious: **both** of the assertions
//! that read every gate and check it — the minted-tag/data-release pair and the
//! cross-workflow agreement — passed on a tree where the gate had been *deleted*
//! from `release-binaries.yml`. They iterate over the gates they can find, and a
//! deleted gate is not one of them.
//!
//! Only [`every_release_reachable_root_job_is_tag_gated`] fires, because it
//! iterates over **jobs** and asks each one whether it is gated. It derives the
//! root set — every job with no `needs:`, i.e. every job a release reaches
//! without passing through another — which on `release-wheels.yml` is exactly
//! `{linux, macos, windows, sdist}` and on `release-binaries.yml` is exactly
//! `{build}`. That derivation is what lets the hardcoded `GATED_JOBS` roster be
//! deleted rather than gain a sibling.
//!
//! # Everything is read per-job, and counted
//!
//! Two misses that a whole-workflow reading has, both of which reproduce real
//! failure modes:
//!
//! - Unioning every job's prefixes into one set answers "does *some* job accept
//!   this tag", when the question is "does *each* job accept it". The four wheel
//!   builds are siblings, not a chain — nothing downstream re-checks the tag, it
//!   only `needs:` those four — so a gate dropped from one of four
//!   near-identical blocks is a hole, and it is the v0.14.0 skip-chain verbatim.
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
//!   what any test of these files can see.
//! - **`startsWith(tag, 'v')` still admits any `v`-prefixed name.** The pinned
//!   format is all this repo will *mint*, but releases are also published here by
//!   hand — `test-fixtures-v1` is one — under no naming rule, and a future
//!   `validation-corpus-v1` would fire every matrix. GitHub expressions have
//!   no regex, so requiring a digit after the `v` means either ten `'v0'…'v9'`
//!   arms or hoisting the decision into a one-job `bash` regex the building jobs
//!   consume through `needs:`. The second also removes the duplication across
//!   sibling jobs; both are follow-ups, deliberately not folded into a
//!   release-blocking fix.
//! - **It does not check the rest of the pipeline.** A release can still skip
//!   silently through a broken `needs:` chain, a downstream job-level `if:`, an
//!   unapproved deployment environment, or a missing trusted-publisher config.
//!   This guard covers exactly one seam.
//! - **It reads the gate textually.** A rewrite that expressed the same
//!   condition without a literal `startsWith(github.event.release.tag_name,
//!   '…')` would make [`gate_prefixes_for`] read nothing — which is why every
//!   extraction below asserts it found something, rather than letting an empty
//!   read pass as agreement.
//! - **[`RELEASE_UNREACHABLE_CONDITIONS`] is a closed list, not an evaluator.**
//!   A root job whose `if:` cannot be true on a release needs no tag gate, but
//!   deciding that in general means evaluating GitHub expressions. The list
//!   holds the two spellings that mean it unambiguously; anything else must
//!   carry the gate, so the failure direction is "a legitimate exemption is
//!   rejected until somebody adds it", not "an ungated job passes".

use std::collections::BTreeSet;
use std::path::PathBuf;

use serde_yaml::Value;

/// Where the workflows live. The release-triggered subset is derived from this
/// directory rather than named; see the module docs for the roster that used to
/// sit here and what it missed.
const WORKFLOW_DIR: &str = ".github/workflows";

/// The release automation whose configuration decides what the tag is called.
const RELEASE_PLZ_TOML: &str = "release-plz.toml";

/// The `release` activity type this repository's release workflows trigger on.
/// A `release:` trigger declaring no `types:` at all fires on *every* activity
/// type and so includes this one; [`release_trigger_types`] represents that as
/// an empty set rather than as an absent trigger.
const PUBLISHED: &str = "published";

/// A floor on the derived release-workflow set, so a parser that has silently
/// stopped matching fails loudly rather than passing vacuously over an empty
/// set — the direct analogue of `spec_fixture_setup_filter.rs`'s
/// `MINIMUM_CONSUMERS`.
///
/// **One, not two, and deliberately.** `release-wheels.yml` is the only
/// release-triggered workflow on `main` today; `release-binaries.yml` (#2055) is
/// the second and has not landed. A floor of 2 would make this file red on the
/// branch it is written for, which would make it un-landable except stacked on
/// #2055 — and the point of deriving the set is that the two changes are
/// independent. Raise it to 2 when `release-binaries.yml` lands; the floor's
/// whole job (a blind derivation cannot read as "nothing to check") is done at
/// 1 either way.
const MINIMUM_RELEASE_WORKFLOWS: usize = 1;

/// A floor on the derived publishing-job set, for the same reason. Today
/// `release-wheels.yml`'s `upload` and `publish-pypi` supply both;
/// `release-binaries.yml` adds a third.
const MINIMUM_PUBLISHING_JOBS: usize = 2;

/// The data release that hosts the bulk test corpora. It is published through
/// the same `release: published` event, carries no code, and must never trigger
/// a build matrix. It did once, and had to be cancelled by hand.
const DATA_RELEASE_TAG: &str = "test-fixtures-v1";

/// The `startsWith` call whose literal argument is a prefix a gate accepts.
const GATE_CALL: &str = "startsWith(github.event.release.tag_name, '";

/// The context a gate reads. Every mention of it in a gated job's condition
/// must be inside a [`GATE_CALL`], or the guard is reading less than the gate
/// says.
const TAG_CONTEXT: &str = "tag_name";

/// The condition a job that can ship something must carry.
const RELEASE_EVENT: &str = "github.event_name == 'release'";

/// The step text that marks a job as able to ship something, in `run:` scripts
/// and in `uses:` action references alike.
///
/// Derived from what the jobs *do* rather than from a roster of job names, which
/// is the same argument the workflow set is derived on: `PUBLISHING_JOBS` used
/// to name `["upload", "publish-pypi"]`, which said nothing at all about a
/// second workflow's `upload`.
const PUBLISHING_MARKERS: &[&str] = &[
    "gh release upload",
    "cargo publish",
    "pypa/gh-action-pypi-publish",
];

/// Job-level `if:` conditions that cannot be true on a `release` event, and so
/// exempt a root job from carrying a tag gate.
///
/// A closed, exact-match list rather than an expression evaluator — see the
/// module docs. Nothing uses it today: every root job in every release workflow
/// currently carries a real gate.
const RELEASE_UNREACHABLE_CONDITIONS: &[&str] = &[
    "github.event_name != 'release'",
    "github.event_name == 'workflow_dispatch'",
];

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn read(relative: &str) -> String {
    let path = repo_root().join(relative);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

/// A workflow's `on:` block.
///
/// Looked up as **both** the string `"on"` and the boolean `true`: YAML 1.1
/// resolves a bare `on` key to a boolean, YAML 1.2 leaves it a string, and which
/// one a given serde_yaml release produces is not something this guard's
/// correctness should depend on. Getting that wrong yields an empty workflow set
/// rather than an error, which is why [`MINIMUM_RELEASE_WORKFLOWS`] exists.
fn on_block(workflow: &Value) -> Option<&Value> {
    workflow.as_mapping()?.iter().find_map(|(key, value)| {
        let is_on = match key {
            Value::String(name) => name == "on",
            Value::Bool(resolved) => *resolved,
            _ => false,
        };
        is_on.then_some(value)
    })
}

/// The `release` activity types a workflow triggers on, or `None` if it declares
/// no `release` trigger at all.
///
/// An **empty** set is not "no types" — it is what GitHub Actions means by a
/// `release:` trigger with no `types:` key, namely every activity type. The
/// scalar (`on: release`) and sequence (`on: [release, push]`) shorthands are
/// the same case: neither can carry a `types:` filter.
fn release_trigger_types(workflow: &Value) -> Option<BTreeSet<String>> {
    match on_block(workflow)? {
        Value::String(event) => (event == "release").then(BTreeSet::new),
        Value::Sequence(events) => events
            .iter()
            .any(|event| event.as_str() == Some("release"))
            .then(BTreeSet::new),
        Value::Mapping(events) => {
            let release = events
                .iter()
                .find_map(|(event, body)| (event.as_str() == Some("release")).then_some(body))?;
            Some(match release.get("types").and_then(Value::as_sequence) {
                None => BTreeSet::new(),
                Some(types) => types
                    .iter()
                    .filter_map(Value::as_str)
                    .map(str::to_string)
                    .collect(),
            })
        }
        _ => None,
    }
}

/// Every workflow a published release fires, as `(repo-relative path, parsed
/// document)` pairs sorted by path.
fn release_workflows() -> Vec<(String, Value)> {
    let dir = repo_root().join(WORKFLOW_DIR);
    let entries = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("read {}: {e}", dir.display()))
        .filter_map(Result::ok);

    let mut workflows = Vec::new();
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
        let Some(types) = release_trigger_types(&document) else {
            continue;
        };
        if types.is_empty() || types.contains(PUBLISHED) {
            workflows.push((format!("{WORKFLOW_DIR}/{name}"), document));
        }
    }
    workflows.sort_by(|(left, _), (right, _)| left.cmp(right));

    assert!(
        workflows.len() >= MINIMUM_RELEASE_WORKFLOWS,
        "found {} release-triggered workflow(s) in {WORKFLOW_DIR} ({:?}); at least \
         {MINIMUM_RELEASE_WORKFLOWS} trigger on `release: {PUBLISHED}` today, so the `on:` \
         parse has gone blind rather than the repository having shrunk. Every assertion in \
         this file iterates over that set, so an empty one passes them all.",
        workflows.len(),
        workflows.iter().map(|(file, _)| file).collect::<Vec<_>>(),
    );
    workflows
}

/// A workflow's jobs, as `(name, body)` pairs.
fn jobs<'a>(file: &str, workflow: &'a Value) -> Vec<(String, &'a Value)> {
    let jobs = workflow
        .get("jobs")
        .and_then(Value::as_mapping)
        .unwrap_or_else(|| panic!("{file} declares no jobs mapping"));
    let jobs: Vec<(String, &Value)> = jobs
        .iter()
        .filter_map(|(name, body)| Some((name.as_str()?.to_string(), body)))
        .collect();
    assert!(
        !jobs.is_empty(),
        "{file} triggers on a release and declares no jobs this guard could read"
    );
    jobs
}

/// The jobs a release reaches **directly** — those with no `needs:`.
///
/// Everything else is downstream and skips with whatever it depends on, which is
/// why gating the roots is sufficient and why the roots are where a missing gate
/// is a hole rather than a redundancy.
fn root_jobs<'a>(file: &str, workflow: &'a Value) -> Vec<(String, &'a Value)> {
    let roots: Vec<(String, &Value)> = jobs(file, workflow)
        .into_iter()
        .filter(|(_, body)| body.get("needs").is_none())
        .collect();
    assert!(
        !roots.is_empty(),
        "{file} triggers on a release and every one of its jobs declares `needs:`. That is \
         either a cycle or a `needs:` shape this guard cannot read; either way it leaves no \
         job for the tag gate to sit on."
    );
    roots
}

/// A condition with its `${{ }}` wrapper, quoting style and whitespace
/// normalised, so [`RELEASE_UNREACHABLE_CONDITIONS`] can be matched exactly.
fn normalize_condition(condition: &str) -> String {
    let trimmed = condition.trim();
    let inner = trimmed
        .strip_prefix("${{")
        .and_then(|rest| rest.strip_suffix("}}"))
        .unwrap_or(trimmed);
    inner
        .replace('"', "'")
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

fn is_unreachable_on_a_release(condition: &str) -> bool {
    let normalized = normalize_condition(condition);
    RELEASE_UNREACHABLE_CONDITIONS
        .iter()
        .any(|exempt| normalize_condition(exempt) == normalized)
}

/// Whether a job's steps can ship something, read from the step text rather than
/// from the job's name.
fn job_publishes(job: &Value) -> bool {
    let steps = job
        .get("steps")
        .and_then(Value::as_sequence)
        .map(Vec::as_slice)
        .unwrap_or(&[]);
    steps.iter().any(|step| {
        ["run", "uses"].iter().any(|key| {
            step.get(key).and_then(Value::as_str).is_some_and(|text| {
                PUBLISHING_MARKERS
                    .iter()
                    .any(|marker| text.contains(marker))
            })
        })
    })
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

/// The tag prefixes **one** job will accept, read out of the
/// `startsWith(github.event.release.tag_name, '<prefix>')` calls in its
/// job-level condition.
///
/// Per-job on purpose: unioning a workflow's jobs answers "does some job accept
/// this tag", and a gate dropped from one of four near-identical blocks is
/// exactly the skip-chain that cost v0.14.0 its wheels.
fn gate_prefixes_for(file: &str, job: &str, condition: &str) -> BTreeSet<String> {
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
        "{file}'s `{job}` job does not gate on \
         `startsWith(github.event.release.tag_name, …)`: {condition:?}. Either the tag gate \
         was removed — in which case the `{DATA_RELEASE_TAG}` data release fires the whole \
         matrix again — or it was rewritten in a form this guard cannot read, in which \
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
        "{file}'s `{job}` job mentions `{TAG_CONTEXT}` {mentions} times but only \
         {read_calls} of those are `{GATE_CALL}…')` calls this guard can read: {condition:?}.\n\
         An arm this guard cannot read is an arm it cannot judge — a `contains(…, 'v')` or a \
         `startsWith(…,'…')` written without the space added *beside* the working arms leaves \
         the extracted prefixes untouched, so every assertion below would pass while the gate \
         admits something new. Write the extra arm in the same form, or widen this guard \
         deliberately."
    );

    prefixes
}

/// One job's tag gate, located.
#[derive(Debug)]
struct Gate {
    /// Repo-relative path of the workflow the gate lives in.
    file: String,
    /// The job carrying it.
    job: String,
    /// The tag prefixes it accepts.
    prefixes: BTreeSet<String>,
}

/// Every tag gate in every release-triggered workflow.
///
/// Read from **all** jobs rather than only the root ones: a downstream job may
/// carry a gate of its own, and if it does it is subject to the same agreement
/// and same accept/reject assertions as the roots.
fn all_release_gates() -> Vec<Gate> {
    let mut gates = Vec::new();
    for (file, workflow) in &release_workflows() {
        for (job, body) in jobs(file, workflow) {
            let Some(condition) = body.get("if").and_then(Value::as_str) else {
                continue;
            };
            if !condition.contains(GATE_CALL) {
                continue;
            }
            gates.push(Gate {
                file: file.clone(),
                job: job.clone(),
                prefixes: gate_prefixes_for(file, &job, condition),
            });
        }
    }
    assert!(
        !gates.is_empty(),
        "no `{GATE_CALL}…')` call was found in any release-triggered workflow. Either every \
         gate has been deleted, or this guard's extraction has stopped matching — and a \
         guard that reads no gates passes every assertion about gates."
    );
    gates
}

/// Every job a release reaches directly carries a tag gate.
///
/// **This is the assertion that catches a gate being deleted**, and it is the
/// only one that can: the two below iterate over the gates they find, so a gate
/// that is no longer there is not among them. Measured — deleting the gate line
/// from `release-binaries.yml` leaves both of them green.
///
/// It also subsumes the old `every_building_job_carries_the_tag_gate`'s
/// per-workflow half without naming a single job: the root set is derived from
/// `needs:`, so a fifth wheel build or a second binaries matrix is covered the
/// day it is added.
#[test]
fn every_release_reachable_root_job_is_tag_gated() {
    let mut checked = 0usize;
    for (file, workflow) in &release_workflows() {
        for (job, body) in root_jobs(file, workflow) {
            checked += 1;
            let condition = body.get("if").and_then(Value::as_str).unwrap_or_else(|| {
                panic!(
                    "{file}'s `{job}` job declares no `needs:` and carries no `if:`, so a \
                     release reaches it directly and ungated. That is the whole defect in one \
                     line: the `{DATA_RELEASE_TAG}` data release fires it, and every release \
                     does, whatever the tag says."
                )
            });
            if is_unreachable_on_a_release(condition) {
                continue;
            }
            // Asserts the gate is present and readable, and that no arm of the
            // condition touches the tag outside a gate this guard read.
            gate_prefixes_for(file, &job, condition);
        }
    }
    eprintln!("{checked} release-reachable root job(s) checked");
}

/// Every gate accepts the tag release-plz will actually mint, and none of them
/// accepts the data release.
///
/// The first half is the assertion that would have failed on the v0.14.0 tree:
/// the gate held `v`, the workspace had gone multi-package, and the minted tag
/// was `ferro-hgvs-v0.14.0`.
///
/// The second half is asserted with it rather than beside it because the obvious
/// repair for the first is a looser predicate: a `contains('v')` or a bare
/// suffix match would accept `test-fixtures-v1` and re-break the thing the gate
/// was added for. Neither half can be satisfied alone.
#[test]
fn all_release_gates_accept_the_minted_tag_and_reject_the_data_release() {
    let (tag, source) = minted_tag();
    let gates = all_release_gates();

    let declining: Vec<&Gate> = gates
        .iter()
        .filter(|gate| !gate.prefixes.iter().any(|prefix| tag.starts_with(prefix)))
        .collect();
    assert!(
        declining.is_empty(),
        "release-plz will tag the next release `{tag}` ({source:?}), and these gates do not \
         accept it: {declining:?}.\n\
         Each of those jobs would skip, and everything that `needs:` it skips with it — \
         including `verify-pypi`, which exists to catch exactly this. `publish.yml` would \
         still report success, so the release would look green while nothing shipped. That is \
         what happened to v0.14.0.\n\
         Fix whichever half is wrong: `git_tag_name` in {RELEASE_PLZ_TOML} is the pin, and \
         every gate must match it."
    );

    let accepting: Vec<(&Gate, Vec<&String>)> = gates
        .iter()
        .filter_map(|gate| {
            let matched: Vec<&String> = gate
                .prefixes
                .iter()
                .filter(|prefix| DATA_RELEASE_TAG.starts_with(prefix.as_str()))
                .collect();
            (!matched.is_empty()).then_some((gate, matched))
        })
        .collect();
    assert!(
        accepting.is_empty(),
        "these gates accept the `{DATA_RELEASE_TAG}` data release: {accepting:?}. That \
         release carries the bulk test corpora and no code; firing a build matrix on it \
         uploads artifacts onto a data release and attempts a publish. It happened once and \
         had to be cancelled by hand."
    );

    eprintln!(
        "{} gates checked, all accepting {:?}",
        gates.len(),
        gates[0].prefixes
    );
}

/// Every release workflow gates on the *same* tag prefixes.
///
/// This is the cross-file strengthening of the old
/// `every_building_job_carries_the_tag_gate`, which compared four jobs inside
/// one file. Nothing downstream re-checks the tag, so the narrowest gate decides
/// which jobs run and the widest decides which releases reach them — and once
/// there are two release workflows, a release accepted by one and declined by
/// the other ships half a release, which is materially the v0.14.0 outcome with
/// a different half missing.
#[test]
fn every_release_workflow_gates_on_the_same_tag_prefixes() {
    let gates = all_release_gates();
    let (first, rest) = gates.split_first().expect("all_release_gates is non-empty");
    for gate in rest {
        assert_eq!(
            gate.prefixes,
            first.prefixes,
            "the release tag gates disagree. `{}`'s `{}` job accepts {:?}, `{}`'s `{}` job \
             accepts {:?}.\n  \
             accepted by `{}` only: {:?}\n  \
             accepted by `{}` only: {:?}\n\
             All release workflows must agree on one gate: nothing downstream re-checks the \
             tag, so a release the wheels build and the binaries decline is a release that \
             looks green and ships half of what it should.",
            first.file,
            first.job,
            first.prefixes,
            gate.file,
            gate.job,
            gate.prefixes,
            first.job,
            first
                .prefixes
                .difference(&gate.prefixes)
                .collect::<Vec<_>>(),
            gate.job,
            gate.prefixes
                .difference(&first.prefixes)
                .collect::<Vec<_>>(),
        );
    }
}

/// Widening the tag gate must not let a `workflow_dispatch` dry-run publish.
///
/// `workflow_dispatch` is deliberately outside the tag gate — the dry-run form
/// takes a blank tag — so the only thing keeping a dry run from uploading to a
/// GitHub Release or to PyPI is the publishing jobs' own conditions. A change to
/// the gate above cannot weaken them, but it is the change most likely to be
/// accompanied by a well-meant edit here.
///
/// The publisher set is **derived from what the steps do** — `gh release
/// upload`, `cargo publish`, `pypa/gh-action-pypi-publish` — rather than from a
/// list of job names. The list it replaces held `["upload", "publish-pypi"]`,
/// which named `release-wheels.yml`'s two jobs and could say nothing whatever
/// about a second workflow's `upload`.
///
/// A floor, and only a floor: this requires each job to *mention* a real release
/// event, not that it fires on nothing else. `upload` deliberately also runs on a
/// dispatch carrying an explicit non-empty tag — that is how the v0.14.0 rescue
/// re-attached its wheels — so the assertion cannot be stronger than a substring
/// test without pinning that second arm too.
#[test]
fn every_publishing_job_requires_a_release() {
    let mut publishers = Vec::new();
    for (file, workflow) in &release_workflows() {
        for (job, body) in jobs(file, workflow) {
            if !job_publishes(body) {
                continue;
            }
            publishers.push(format!("{file}:{job}"));
            let condition = body.get("if").and_then(Value::as_str).unwrap_or_else(|| {
                panic!(
                    "{file}'s `{job}` job can ship artifacts and carries no `if:` at all, so a \
                     `workflow_dispatch` dry-run publishes."
                )
            });
            assert!(
                condition.contains(RELEASE_EVENT),
                "{file}'s `{job}` job can ship artifacts and its condition does not require a \
                 real release event: {condition:?}"
            );
        }
    }

    assert!(
        publishers.len() >= MINIMUM_PUBLISHING_JOBS,
        "found {} publishing job(s) ({publishers:?}) across the release workflows; at least \
         {MINIMUM_PUBLISHING_JOBS} ship something today, so the scan for {PUBLISHING_MARKERS:?} \
         has gone blind rather than the pipeline having shrunk — and a scan that finds no \
         publisher asserts nothing about publishing.",
        publishers.len(),
    );
}
