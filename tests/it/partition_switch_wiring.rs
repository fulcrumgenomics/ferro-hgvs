//! The `FERRO_PARTITION` refusal must be delivered by every entry point.
//!
//! `partition_rule()` used to panic on a value naming no arm, in every build,
//! from inside `canonicalize_from_sequence` — infallible, under every public
//! normalization entry point, and reached across the FFI boundary by the Python
//! bindings. A development-only bake-off switch must not be able to abort a
//! process that merely happens to have the variable set, so a release build of
//! the library now falls safe to the shipped `live` rule and retains the refusal
//! in `partition_switch_startup_error()` for a caller to report.
//!
//! That trade only holds if the callers actually report it. Falling safe without
//! reporting is the *original* defect this switch was given a refusal to prevent:
//! a bake-off served the shipped rule under a candidate's name, coming back with
//! a clean empty diff that reads as "the candidate changes nothing". So the
//! wiring is the contract, and the first two tests pin both halves of it:
//!
//! - the behaviour, end to end, through the `ferro` binary — a misspelling exits
//!   non-zero with the message, a real arm still runs, and `--help` keeps working;
//! - the completeness, over `Cargo.toml`'s own `[[bin]]` **and** `[[example]]`
//!   tables — so a target added later cannot skip the call and leave
//!   `README.md`'s "every binary and example in this repository reports it"
//!   quietly false.
//!
//! The second is a source scan, so it proves the call is *written*, not that it
//! is reached. The first is what proves it is reached, and it is written against
//! the binary a user actually runs.
//!
//! # A third test, about a different call site the arms are the only way to reach
//!
//! `canonicalize_from_sequence` calls `coalesce_compensating_gap_split` behind
//! `compensating_gap_coalesce_applies(partition_rule(), kind)` — the arm **and**
//! the axis (#1711). `merge.rs`'s own unit tests pin that predicate exhaustively,
//! and they can only pin the predicate: `partition_rule()` caches into a
//! process-global `OnceLock`, so an in-process test cannot select the `canonical`
//! / `canonical-coalesced` arms the pass runs on. That reads as "the call site
//! cannot be tested", and it is not — a **subprocess** gets its own `OnceLock`,
//! which is exactly what the two tests above already exploit, against the mock
//! provider and so with no prepared reference.
//!
//! Measured before this test existed, and the reason it now does: replacing
//! that call site's axis argument with a hardcoded `CisKind::Cds` — #1711
//! restored verbatim — left the whole suite passing and clippy exiting 0. A
//! well-guarded predicate wired to nothing is indistinguishable from a fix.
//!
//! It lives here rather than beside the predicate because what it asks is
//! whether the call site consults the axis at all. That is a wiring question,
//! and the answer to it survives either reading of what `DNA/delins.md:47`
//! reaches.

use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

/// A value no build has an arm for. Deliberately the near-miss `README.md` uses,
/// since a plausible misspelling is the input the refusal exists for.
const MISSPELLED: &str = "canonicl";

/// The repository root, from this test target's manifest directory.
fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// A deletion inside the `CCCC` homopolymer at c.16-19 of the mock `NM_000088.3`,
/// so the normalizer 3'-shifts it and the positive control below exercises the
/// partitioner rather than only the parser. Within the mock CDS (60 bases), which
/// the more obvious `c.459A>G` is not — that one is rejected as past the CDS end
/// under the default strict mode.
const SHIFTED: &str = "NM_000088.3:c.16del";

/// Run `ferro normalize` over `inputs` with `FERRO_PARTITION` set to `value`.
///
/// The mock provider serves every accession used here, so no `--reference` and
/// no prepared reference directory is needed — which is what makes driving the
/// real binary cheap enough to be the default answer whenever a `OnceLock` puts
/// a switch out of an in-process test's reach.
fn normalize(inputs: &[&str], value: &str) -> std::process::Output {
    let mut input = tempfile::Builder::new()
        .suffix(".txt")
        .tempfile()
        .expect("create input file");
    for line in inputs {
        writeln!(input, "{line}").expect("write input");
    }
    input.flush().expect("flush input");

    Command::new(env!("CARGO_BIN_EXE_ferro"))
        .arg("normalize")
        .arg("-i")
        .arg(input.path())
        .env("FERRO_PARTITION", value)
        .output()
        .expect("run ferro normalize")
}

/// Run `ferro normalize` over one variant with `FERRO_PARTITION` set to `value`.
fn normalize_one(value: &str) -> std::process::Output {
    normalize(&[SHIFTED], value)
}

/// A misspelled arm name stops the CLI before it reads input, and a real one does
/// not.
///
/// Both directions are asserted in one test on purpose: a binary that refused
/// *every* value would pass the first assertion while making the switch useless,
/// and that failure would look identical in a log to the fix working.
#[test]
fn the_ferro_cli_refuses_a_misspelled_partition_switch_and_still_serves_a_real_arm() {
    let refused = normalize_one(MISSPELLED);
    assert!(
        !refused.status.success(),
        "a value naming no arm must fail the process; stdout was:\n{}",
        String::from_utf8_lossy(&refused.stdout)
    );
    let stderr = String::from_utf8_lossy(&refused.stderr);
    assert!(
        stderr.contains(MISSPELLED) && stderr.contains("is not a partitioner this build has"),
        "the refusal must name the value given and say what it is; stderr was:\n{stderr}"
    );
    assert!(
        stderr.contains("canonical-coalesced"),
        "the refusal must name THIS build's arms, so a value that exists on \
         another branch is reported as absent here; stderr was:\n{stderr}"
    );
    assert!(
        refused.stdout.is_empty(),
        "the refusal must precede any output that could be mistaken for a \
         measurement; stdout was:\n{}",
        String::from_utf8_lossy(&refused.stdout)
    );

    let accepted = normalize_one("canonical");
    assert!(
        accepted.status.success(),
        "a value that IS an arm must still run; stderr was:\n{}",
        String::from_utf8_lossy(&accepted.stderr)
    );
}

/// `--help` still answers under a stale `FERRO_PARTITION`.
///
/// Reading the help is how somebody finds out what the arms are, so refusing it
/// would leave them holding a rejected value and no way to look up the right one.
/// This is why the check sits after `Cli::parse()` rather than before it — clap
/// answers `--help` by exiting from inside `parse`, and nothing has been read or
/// normalized at that point either way.
#[test]
fn help_still_works_when_the_partition_switch_names_no_arm() {
    let out = Command::new(env!("CARGO_BIN_EXE_ferro"))
        .arg("--help")
        .env("FERRO_PARTITION", MISSPELLED)
        .output()
        .expect("run ferro --help");
    assert!(
        out.status.success(),
        "`--help` must not be refused; stderr was:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        !out.stdout.is_empty(),
        "`--help` must still print the help text"
    );
}

/// Every target declared in `Cargo.toml` — `[[bin]]` **and** `[[example]]` —
/// whose source calls `partition_switch_startup_error`.
///
/// Parsed as TOML rather than substring-matched, for the reason
/// `generator_completeness.rs` gives: a `path` inside a comment must not count,
/// and the two tables have to be told apart.
///
/// # Why the denominator is the manifest, and not a list kept here
///
/// This test used to read the `[[bin]]` table and then **append one example by
/// hand** — `dump_normalized_corpus`, the harness somebody had remembered to
/// wire. Every other bake-off example was outside its denominator, so the
/// completeness it asserted was completeness over a hand-written list.
///
/// Measured on `examples/dump_confluence_divergences.rs`, which was outside it:
/// a **release** build run with `FERRO_PARTITION=canonicl` printed a census
/// byte-identical to the `live` one, on empty stderr, at exit 0. A misspelled
/// arm and a correct run of the shipped arm were the same observation. A
/// **debug** build was no better as an instrument — the library's own refusal
/// fires there, but `catch_unwind` in the harness swallows it per row, so the
/// run still exits 0 and reports `divergent: 0`, which reads as a candidate
/// that converges everything.
///
/// So the fix is not a more central refusal — `partition_rule_from_env` already
/// refuses centrally, and in a debug build it did refuse, 47 377 times, without
/// the run failing. What was wrong is *where the refusal is delivered*: at a
/// process boundary, which is per-entry-point by construction. Making the
/// denominator the manifest's own tables is what removes the "targets somebody
/// remembered to add" step, since a target that is not in `Cargo.toml` does not
/// build at all.
///
/// **That last clause is a property of `Cargo.toml`, not of cargo.** Cargo
/// auto-discovers `src/bin/*.rs`, `examples/*.rs` and `benches/*.rs` unless
/// told not to, and an auto-discovered target would build, run, and sit outside
/// this denominator — the same hole under a different name. `autoexamples`,
/// `autobins` and `autobenches` are therefore all `false` in `[package]`, which
/// is what makes "not in `Cargo.toml`" and "does not build" the same statement.
/// If any of the three is ever turned back on, this test stops being complete
/// even though it still passes.
///
/// The obligation is uniform rather than scoped to the targets that normalize
/// today. A target that does not normalize loses nothing by refusing a value
/// naming no arm — the caller asked for something that does not exist — and a
/// uniform rule has no per-target judgement to get wrong when an example starts
/// normalizing later.
///
/// # The one table deliberately left out: `[[bench]]`
///
/// Three of the four declared benches normalize (`benches/benchmarks.rs`,
/// `benches/baseline_head.rs`, `benches/seqfirst_align.rs`), so the argument
/// above reaches them and they are **not** covered. The obstacle is mechanical
/// rather than principled: each ends in `criterion_main!`, which *generates* the
/// `fn main`, so there is no hand-written entry point to put the call in without
/// hand-expanding that macro and taking on its argument handling. A bake-off run
/// through `cargo bench` under a misspelled arm therefore still times the
/// shipped rule under a candidate's name. Wiring them means replacing
/// `criterion_main!` in all four; until someone does, the exclusion is stated
/// here rather than left to be rediscovered from the absence of a table read.
///
/// Both denominators are asserted non-empty, so "0 of 0 targets are missing it"
/// cannot pass as a result — the failure mode this repository has been bitten by
/// often enough to make it a house rule.
#[test]
fn every_binary_and_example_target_reports_a_refused_partition_switch() {
    let manifest_text =
        std::fs::read_to_string(repo_root().join("Cargo.toml")).expect("read Cargo.toml");
    let manifest: toml::Value = manifest_text.parse().expect("parse Cargo.toml as TOML");

    /// The targets declared under one manifest table, as `(name, source path)`.
    ///
    /// `default_dir` supplies cargo's own path convention for the table, since
    /// `[[example]]` entries here declare only a name — the file is
    /// `examples/<name>.rs` — while `[[bin]]` entries all declare a path
    /// explicitly and are required to keep doing so.
    fn targets(
        manifest: &toml::Value,
        table: &str,
        default_dir: Option<&str>,
    ) -> BTreeSet<(String, String)> {
        let mut found = BTreeSet::new();
        for entry in manifest
            .get(table)
            .and_then(toml::Value::as_array)
            .unwrap_or_else(|| panic!("Cargo.toml declares at least one [[{table}]] target"))
        {
            let name = entry
                .get("name")
                .and_then(toml::Value::as_str)
                .unwrap_or_else(|| panic!("every [[{table}]] has a name"));
            let path = match entry.get("path").and_then(toml::Value::as_str) {
                Some(path) => path.to_string(),
                None => match default_dir {
                    Some(dir) => format!("{dir}/{name}.rs"),
                    None => panic!("[[{table}]] {name} declares no path"),
                },
            };
            found.insert((name.to_string(), path));
        }
        found
    }

    let bins = targets(&manifest, "bin", None);
    assert!(
        bins.len() >= 5,
        "expected the five known binaries (ferro, ferro-web, ferro-benchmark and \
         the two spec generators); found {bins:?}"
    );
    let examples = targets(&manifest, "example", Some("examples"));
    assert!(
        examples.len() >= 20,
        "expected the repository's ~23 declared examples — the bake-off \
         harnesses, the window extractors and the artifact generators; found \
         {examples:?}"
    );

    let missing: Vec<&str> = bins
        .iter()
        .chain(examples.iter())
        .filter(|(_, path)| {
            let text = std::fs::read_to_string(repo_root().join(path))
                .unwrap_or_else(|e| panic!("read {path}: {e}"));
            !text.contains("partition_switch_startup_error()")
        })
        .map(|(name, _)| name.as_str())
        .collect();

    assert!(
        missing.is_empty(),
        "these entry points run without reporting a refused FERRO_PARTITION, so \
         a misspelled arm serves the shipped rule under a candidate's name and \
         the run cannot be told apart from a correct `live` one: {missing:?}. \
         Call `ferro_hgvs::normalize::partition_switch_startup_error()` after \
         argument parsing and fail on `Some`."
    );
}

/// The two `partition_block_canonical` arms — the ones
/// `PartitionRule::cuts_with_canonical` names, and the only ones on which
/// `coalesce_compensating_gap_split` runs at all.
const CANONICAL_ARMS: [&str; 2] = ["canonical", "canonical-coalesced"];

/// A cis allele on the **non-coding** `n.` axis whose re-derivation carries a
/// compensating gap, so the pass fires on it once the axis test is bypassed.
///
/// The mock `NR_000123.1` is `ACGTACGTACGT` with no CDS, so `n.` is
/// transcript-relative and the axis is `CisKind::Tx` — `false` under both
/// `AxisFrame::is_coding_dna` and `AxisFrame::carries_translated_frame`, which
/// is what makes it a test of the *kind* rather than of the frame flag.
///
/// The shape is the one #1711 was found on: an `inv` and a `del` one unchanged
/// base apart, the same family as the spec conformance corpus rows
/// (`s01-n3{m,p}-pair-del-inv-p4-sep1`) whose merge that corpus's own negative
/// guard for the **rejected** SVD-WG010 proposal counts.
const NONCODING_COMPENSATING_GAP: &str = "NR_000123.1:n.[2_4inv;6del]";

/// What the pass makes of [`NONCODING_COMPENSATING_GAP`] when the axis test is
/// bypassed: one `delins` spanning the members *and the unchanged bases between
/// them* — two variants merged across a gap `general.md:34` keeps individual.
///
/// Measured on the mutated builds rather than reasoned about, on both arms of
/// [`CANONICAL_ARMS`].
const NONCODING_MERGED: &str = "NR_000123.1:n.2_6delinsACGA";

/// The coding-axis control, in the same shape, on the mock `NM_000088.3`.
///
/// It is here so this test cannot pass in a build where the pass is simply never
/// called — which is the one mutation the negative half above is blind to, and
/// which would look exactly like the fix. Removing the call site entirely splits
/// this row into `c.[2dup;4_5inv;8_9del]` instead (measured), so the assertion
/// below fails.
const CODING_COMPENSATING_GAP: &str = "NM_000088.3:c.[3_7inv;9del]";

/// The single spanning `delins` [`CODING_COMPENSATING_GAP`] must still reach on
/// the coding axis, where `DNA/delins.md:47`'s carve-out does apply.
const CODING_MERGED: &str = "NM_000088.3:c.3_9delinsTGGGCA";

/// Each input's normalized output, read out of `ferro normalize`'s stdout.
///
/// A row the normalizer left alone is printed as the input alone, with no arrow,
/// so it maps to itself; anything else would silently drop exactly the rows a
/// caller most wants to see unchanged.
fn normalized_outputs(inputs: &[&str], arm: &str) -> BTreeMap<String, String> {
    let out = normalize(inputs, arm);
    assert!(
        out.status.success(),
        "`ferro normalize` failed under FERRO_PARTITION={arm}; stderr was:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout).into_owned();
    let mut parsed = BTreeMap::new();
    for line in stdout.lines() {
        let (input, output) = line.split_once(" -> ").unwrap_or((line, line));
        if inputs.contains(&input) {
            parsed.insert(input.to_string(), output.to_string());
        }
    }
    for input in inputs {
        assert!(
            parsed.contains_key(*input),
            "no output for {input} under FERRO_PARTITION={arm}; stdout was:\n{stdout}"
        );
    }
    parsed
}

/// The compensating-gap coalesce consults the **axis** at its call site, not
/// only in its predicate. #1711.
///
/// Both directions in one test, for the reason the misspelling test above gives:
/// a build that declined the pass on *every* axis would satisfy the negative
/// half while making the pass dead, and in a log that is indistinguishable from
/// the fix working.
///
/// The negative half deliberately does **not** pin the split spelling. What
/// `general.md:34` governs is that the members stay individual; which individual
/// members the re-derivation lands on is a separate question, and pinning it
/// here would redden this guard for changes that leave its subject untouched.
#[test]
fn the_cli_declines_the_compensating_gap_coalesce_off_the_coding_dna_axis() {
    for arm in CANONICAL_ARMS {
        let outputs =
            normalized_outputs(&[NONCODING_COMPENSATING_GAP, CODING_COMPENSATING_GAP], arm);

        let noncoding = &outputs[NONCODING_COMPENSATING_GAP];
        assert_ne!(
            noncoding, NONCODING_MERGED,
            "on FERRO_PARTITION={arm} the `n.` row merged across the unchanged \
             base between its members. That is the **rejected** SVD-WG010 \
             proposal, and `DNA/delins.md:47` — the only clause that licenses \
             merging across unchanged bases — reaches `c.` and nothing else \
             (`delins-payload-coincidence-carve-out-is-coding-dna-scoped`). The \
             call site must pass the description's own axis to \
             `compensating_gap_coalesce_applies`, not a constant (#1711)."
        );
        assert!(
            noncoding.contains(';'),
            "on FERRO_PARTITION={arm} the `n.` row collapsed to one member \
             ({noncoding}); off the coding DNA axis `general.md:34` governs and \
             the members stay individual (#1711)"
        );

        let coding = &outputs[CODING_COMPENSATING_GAP];
        assert_eq!(
            coding, CODING_MERGED,
            "on FERRO_PARTITION={arm} the `c.` row did not reach the spanning \
             `delins`, so this test is no longer measuring anything: the pass is \
             not running on the arm it is scoped to, and the assertion above \
             would hold in a build that had simply deleted the call. Check the \
             call site before re-blessing this string."
        );
    }
}
