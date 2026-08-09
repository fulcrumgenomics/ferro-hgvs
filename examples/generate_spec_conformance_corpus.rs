//! Report the exhaustive spec-derived conformance corpus, and check its coverage
//! claims against the spec checkout.
//!
//! # What this is, and what it is not
//!
//! It is **not** the way the corpus reaches the test suite. The corpus is
//! enumerated at test time from `ferro_hgvs::conformance::spec_corpus`, so every
//! case runs on every `cargo test` with no artifact to stage — see
//! `tests/it/spec_conformance_axis.rs`. This binary exists for the two jobs a
//! library module cannot do:
//!
//! 1. **Check the rule inventory against the spec checkout.** Every clause the
//!    inventory cites carries a verbatim quote, and every file it declares a
//!    blanket classification for carries a pinned clause-unit count. Both are
//!    re-derived from `assets/hgvs-nomenclature` here, so a submodule bump that
//!    moves a clause, or adds one under a blanket exemption, **fails the build**
//!    instead of leaving a citation pointing at unrelated prose. This is the same
//!    contract `generate_spec_fixture` enforces for the `rulings` section.
//! 2. **Print the coverage table**, per spec file, over that inventory.
//!
//! # Coverage claims are checked, not asserted
//!
//! An inventory entry classified `generatable` names the rule tags that exercise
//! it. If no row in the corpus carries one of those tags, this refuses: a claim
//! of coverage the corpus does not back is exactly the defect the corpus exists
//! to remove.
//!
//! # Usage
//!
//! ```text
//! cargo run --features dev --example generate_spec_conformance_corpus
//! cargo run --features dev --example generate_spec_conformance_corpus -- --seeds 4
//! cargo run --features dev --example generate_spec_conformance_corpus -- --extended-scale
//! cargo run --features dev --example generate_spec_conformance_corpus -- --dump <path>
//! ```

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use serde::{Deserialize, Serialize};

use ferro_hgvs::conformance::completeness::{Allowance, CaptureCounts, CaptureLedger};
use ferro_hgvs::conformance::spec_corpus::{
    enumerate, CorpusBounds, RowKind, SpecCorpus, BLOCK_LADDER, EXTENDED_BLOCK_LADDER,
    SEPARATION_LADDER,
};

/// Where the committed rule inventory lives, relative to the repo root.
const INVENTORY_PATH: &str = "tests/fixtures/spec-corpus/spec_rule_inventory.json";

/// Where the spec checkout lives, relative to the repo root.
const SPEC_ROOT: &str = "assets/hgvs-nomenclature";

/// Designs the enumeration is allowed to drop.
///
/// Deliberately a ceiling rather than an exact figure: every drop reason is a
/// legitimate outcome of an adversarial enumeration (a region too narrow for a
/// design, a palindromic inversion, a design whose members cancel), and the
/// per-reason census is printed in full. What the ceiling buys is that a change
/// which starts dropping *most* of the corpus refuses instead of quietly
/// shrinking it — the failure mode `CaptureLedger` exists for.
const MAX_DROPPED_FRACTION: f64 = 0.55;

#[derive(Parser, Debug)]
#[command(about = "Report the exhaustive spec-derived conformance corpus and its coverage")]
struct Cli {
    /// Reference-sequence seeds. Two cores per seed. Prefix-stable.
    #[arg(long, default_value_t = 1)]
    seeds: u32,
    /// Include the extended scale ladder, which crosses `MAX_SHIFT_TRACT`
    /// (32768). Off by default: see `EXTENDED_BLOCK_LADDER`.
    #[arg(long)]
    extended_scale: bool,
    /// Write the corpus rows here as JSON, for offline inspection. The suite does
    /// not read this — it enumerates the corpus itself.
    #[arg(long)]
    dump: Option<PathBuf>,
    /// Write the coverage report here instead of to stdout.
    #[arg(long)]
    report: Option<PathBuf>,
    /// Skip the spec-checkout verification. For a checkout without the submodule;
    /// the coverage table is still printed, and says it was unverified.
    #[arg(long)]
    skip_spec_check: bool,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(report) => {
            match &cli.report {
                Some(path) => {
                    if let Err(error) = write_to(path, &report) {
                        eprintln!("error: {error}");
                        return ExitCode::FAILURE;
                    }
                    println!("wrote the coverage report to {}", path.display());
                }
                None => print!("{report}"),
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: &Cli) -> Result<String, String> {
    let bounds = CorpusBounds {
        seeds: cli.seeds,
        extended_scale: cli.extended_scale,
    };

    // Route the population through the ledger. This is the LAST fallible step
    // before anything is written or reported, per `completeness`'s contract: a
    // design either becomes a row or is accounted for by reason.
    let attempts = enumerate(&bounds);
    let mut ledger = CaptureLedger::new("spec conformance corpus designs");
    for attempt in &attempts {
        match attempt {
            Ok(row) => {
                let _ = ledger.record::<(), String>(row.id.as_str(), Ok(()));
            }
            Err((id, reason)) => ledger.record_drop(id.as_str(), reason.label()),
        }
    }
    let attempted = ledger.attempted();
    #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
    let ceiling = (attempted as f64 * MAX_DROPPED_FRACTION) as usize;
    let counts = ledger
        .finish_with(Allowance::at_most(
            ceiling,
            "every drop reason is a legitimate outcome of an adversarial enumeration \
             (a region too narrow for a design, a palindromic inversion, a design whose \
             members cancel); the ceiling exists so a change that drops most of the \
             corpus refuses rather than quietly shrinking it",
        ))
        .map_err(|shortfall| shortfall.to_string())?;

    let corpus = SpecCorpus::from_attempts(bounds.clone(), attempts);

    let inventory = Inventory::load(Path::new(INVENTORY_PATH))?;
    let spec_root = (!cli.skip_spec_check).then(|| PathBuf::from(SPEC_ROOT));
    inventory.verify(spec_root.as_deref(), &corpus)?;

    if let Some(path) = &cli.dump {
        let dump = Dump::of(&corpus, &counts);
        let mut rendered = serde_json::to_string_pretty(&dump)
            .map_err(|e| format!("serializing the dump: {e}"))?;
        rendered.push('\n');
        write_to(path, &rendered)?;
        println!("wrote {} rows to {}", corpus.rows.len(), path.display());
    }

    Ok(render_report(
        &corpus,
        &counts,
        &inventory,
        spec_root.is_some(),
    ))
}

fn write_to(path: &Path, contents: &str) -> Result<(), String> {
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent)
            .map_err(|e| format!("creating {}: {e}", parent.display()))?;
    }
    std::fs::write(path, contents).map_err(|e| format!("writing {}: {e}", path.display()))
}

// ---------------------------------------------------------------------------
// The dump
// ---------------------------------------------------------------------------

/// One row, flattened for JSON.
#[derive(Serialize)]
struct DumpRow {
    id: String,
    kind: String,
    stratum: &'static str,
    shape: &'static str,
    mechanism: &'static str,
    members: usize,
    multi_member: bool,
    separation: usize,
    block_len: usize,
    geometry: &'static str,
    region: &'static str,
    scale_bands: Vec<&'static str>,
    negative_guards: Vec<&'static str>,
    prohibition: Option<String>,
    rules: Vec<&'static str>,
    spellings: Vec<String>,
}

/// The corpus as an inspectable artifact, carrying its own completeness claim.
#[derive(Serialize)]
struct Dump {
    description: &'static str,
    generator: &'static str,
    seeds: u32,
    extended_scale: bool,
    designs_considered: usize,
    /// The ledger's claim, stamped in so a reader of the dump sees how much the
    /// pass dropped without re-running the generator.
    completeness: CaptureCounts,
    rows: Vec<DumpRow>,
}

impl Dump {
    fn of(corpus: &SpecCorpus, counts: &CaptureCounts) -> Self {
        Self {
            description: "Exhaustive spec-derived conformance corpus. Enumerated at test time by \
                          `ferro_hgvs::conformance::spec_corpus`; this dump is for offline \
                          inspection only and no test reads it.",
            generator: "cargo run --features dev --example generate_spec_conformance_corpus",
            seeds: corpus.bounds.seeds,
            extended_scale: corpus.bounds.extended_scale,
            designs_considered: corpus.designs_considered,
            completeness: counts.clone(),
            rows: corpus
                .rows
                .iter()
                .map(|row| DumpRow {
                    id: row.id.clone(),
                    kind: row.kind.to_string(),
                    stratum: row.stratum,
                    shape: row.shape.label(),
                    mechanism: row.mechanism.label(),
                    members: row.members,
                    multi_member: row.is_multi_member(),
                    separation: row.separation,
                    block_len: row.block_len,
                    geometry: row.geometry.label(),
                    region: row.region.label(),
                    scale_bands: row.scale_bands.clone(),
                    negative_guards: row.negative_guards.clone(),
                    prohibition: row
                        .prohibition
                        .map(|(rule, strength)| format!("{rule} ({})", strength.label())),
                    rules: row.rules.clone(),
                    spellings: row.spellings.clone(),
                })
                .collect(),
        }
    }
}

// ---------------------------------------------------------------------------
// The rule inventory
// ---------------------------------------------------------------------------

/// How a clause is classified against this instrument.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "kebab-case")]
enum Classification {
    /// Inputs exercising it can be synthesised, and the corpus does.
    Generatable,
    /// The spec gives a worked answer but the rule is not otherwise
    /// mechanisable, so the existing spec-fixture guards
    /// (`ferro_produces_the_form_the_spec_states`) are what cover it.
    PublishedExampleOnly,
    /// It cannot be exercised by synthesis, and the entry says why.
    #[serde(rename = "not-generatable")]
    NotGenerated,
}

impl Classification {
    fn label(self) -> &'static str {
        match self {
            Self::Generatable => "generatable",
            Self::PublishedExampleOnly => "published-example-only",
            Self::NotGenerated => "not-generatable",
        }
    }
}

/// One clause, cited exactly and quoted verbatim.
#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct ClauseEntry {
    /// `docs/recommendations/general.md:34`, relative to the spec checkout root.
    clause: String,
    /// A verbatim substring of that line. Checked against the checkout, so a
    /// submodule bump that moves the clause fails the build.
    quote: String,
    classification: Classification,
    /// Rule tags in the corpus that exercise it. Required and non-empty for
    /// `generatable`; must be empty otherwise.
    #[serde(default)]
    exercised_by: Vec<String>,
    /// Why it is not generatable. Required for `not-generatable`.
    #[serde(default)]
    reason: Option<String>,
}

/// A blanket classification for every clause unit in one file that no
/// [`ClauseEntry`] names.
///
/// The pinned `clause_units` count is what makes a blanket honest: a submodule
/// bump that adds a clause under one of these fails the build rather than
/// enlarging an exemption in silence.
#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct FileScope {
    /// `docs/recommendations/protein/substitution.md`.
    file: String,
    /// Clause units the extractor finds in it, pinned.
    clause_units: usize,
    classification: Classification,
    #[serde(default)]
    reason: Option<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Inventory {
    description: String,
    /// How a clause unit is recognised. Documented in the artifact so the
    /// denominator is not a property of whoever last read the generator.
    clause_unit_definition: String,
    clauses: Vec<ClauseEntry>,
    files: Vec<FileScope>,
}

impl Inventory {
    fn load(path: &Path) -> Result<Self, String> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| format!("reading the rule inventory {}: {e}", path.display()))?;
        serde_json::from_str(&text)
            .map_err(|e| format!("parsing the rule inventory {}: {e}", path.display()))
    }

    /// Check the inventory against the spec checkout and against the corpus.
    ///
    /// Four refusals, each of which is a way the inventory could otherwise read
    /// as stronger than it is:
    ///
    /// - a cited clause whose quote is no longer on that line;
    /// - a file whose pinned clause-unit count has moved;
    /// - a `generatable` clause no corpus row exercises;
    /// - a `not-generatable` clause — or blanket file scope — with no stated
    ///   reason.
    fn verify(&self, spec_root: Option<&Path>, corpus: &SpecCorpus) -> Result<(), String> {
        let mut problems = Vec::new();

        for entry in &self.clauses {
            match entry.classification {
                Classification::Generatable => {
                    if entry.exercised_by.is_empty() {
                        problems.push(format!(
                            "{}: classified generatable but names no rule tag",
                            entry.clause
                        ));
                    }
                }
                _ => {
                    if !entry.exercised_by.is_empty() {
                        problems.push(format!(
                            "{}: classified {} but names rule tags; only a generatable \
                             clause may claim coverage",
                            entry.clause,
                            entry.classification.label()
                        ));
                    }
                }
            }
            if entry.classification == Classification::NotGenerated
                && entry.reason.as_deref().unwrap_or("").trim().is_empty()
            {
                problems.push(format!(
                    "{}: classified not-generatable with no reason. A rule that cannot be \
                     generated is a legitimate answer; a rule silently skipped is the defect.",
                    entry.clause
                ));
            }
        }

        // A blanket file scope exempts every clause unit in that file, so it is
        // the larger exemption, not the smaller one. It must state a reason on
        // exactly the same terms as a single clause.
        for scope in &self.files {
            if scope.classification == Classification::NotGenerated
                && scope.reason.as_deref().unwrap_or("").trim().is_empty()
            {
                problems.push(format!(
                    "{}: the whole file is classified not-generatable with no reason. A rule \
                     that cannot be generated is a legitimate answer; a rule silently skipped \
                     is the defect.",
                    scope.file
                ));
            }
        }

        // Coverage claims must be backed by rows that exist.
        let tags: BTreeSet<&str> = corpus
            .rows
            .iter()
            .flat_map(|row| row.rules.iter().copied())
            .collect();
        for entry in &self.clauses {
            if entry.classification != Classification::Generatable {
                continue;
            }
            for tag in &entry.exercised_by {
                if !tags.contains(tag.as_str()) {
                    problems.push(format!(
                        "{}: claims coverage by `{tag}`, which no corpus row carries",
                        entry.clause
                    ));
                }
            }
        }

        // Every rule tag the corpus carries must be claimed by some clause, or
        // the corpus is generating shapes nobody can trace to a rule.
        let claimed: BTreeSet<&str> = self
            .clauses
            .iter()
            .flat_map(|entry| entry.exercised_by.iter().map(String::as_str))
            .collect();
        for tag in &tags {
            if !claimed.contains(tag) {
                problems.push(format!(
                    "rule tag `{tag}` is carried by corpus rows but claimed by no inventory \
                     clause; either cite the clause or drop the tag"
                ));
            }
        }

        if let Some(root) = spec_root {
            for entry in &self.clauses {
                if let Err(error) = check_quote(root, &entry.clause, &entry.quote) {
                    problems.push(error);
                }
            }
            for scope in &self.files {
                let path = root.join(&scope.file);
                match clause_units(&path) {
                    Ok(found) if found == scope.clause_units => {}
                    Ok(found) => problems.push(format!(
                        "{}: pinned at {} clause units, found {found}. A clause added under a \
                         blanket classification must be classified, not absorbed.",
                        scope.file, scope.clause_units
                    )),
                    Err(error) => problems.push(error),
                }
            }
            // Every recommendation file must have a scope row, so "every clause"
            // has a mechanical denominator.
            match recommendation_files(root) {
                Ok(found) => {
                    let declared: BTreeSet<&str> =
                        self.files.iter().map(|s| s.file.as_str()).collect();
                    for file in &found {
                        if !declared.contains(file.as_str()) {
                            problems.push(format!(
                                "{file} is in the spec checkout but has no inventory scope row"
                            ));
                        }
                    }
                    for file in &declared {
                        if !found.iter().any(|f| f == file) {
                            problems.push(format!(
                                "{file} has an inventory scope row but is not in the spec checkout"
                            ));
                        }
                    }
                }
                Err(error) => problems.push(error),
            }
        }

        if problems.is_empty() {
            return Ok(());
        }
        Err(format!(
            "the rule inventory does not hold ({} problems):\n  {}",
            problems.len(),
            problems.join("\n  ")
        ))
    }
}

/// Check that `quote` really is a substring of the cited `file:line`.
fn check_quote(root: &Path, clause: &str, quote: &str) -> Result<(), String> {
    let (file, line_text) = clause
        .rsplit_once(':')
        .ok_or_else(|| format!("{clause}: not a `file:line` citation"))?;
    let line: usize = line_text
        .parse()
        .map_err(|_| format!("{clause}: `{line_text}` is not a line number"))?;
    let path = root.join(file);
    let text = std::fs::read_to_string(&path)
        .map_err(|e| format!("{clause}: reading {}: {e}", path.display()))?;
    let actual = text
        .lines()
        .nth(line.saturating_sub(1))
        .ok_or_else(|| format!("{clause}: the file has fewer than {line} lines"))?;
    if actual.contains(quote) {
        return Ok(());
    }
    Err(format!(
        "{clause}: the quote is no longer on that line.\n      quoted: {quote}\n      line:   {}",
        actual.trim()
    ))
}

/// Clause units in one markdown file.
///
/// A **clause unit** is a line that opens a bullet (`- ` / `* `), a numbered item
/// (`1. `), or an admonition block (`!!! `). That is the definition the inventory
/// records and the denominator every coverage figure is over. It is deliberately
/// syntactic: a semantic definition of "a clause" would be a judgement, and a
/// judgement cannot be a denominator a submodule bump is checked against.
fn clause_units(path: &Path) -> Result<usize, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("counting clause units in {}: {e}", path.display()))?;
    Ok(text.lines().filter(|line| is_clause_unit(line)).count())
}

fn is_clause_unit(line: &str) -> bool {
    let trimmed = line.trim_start();
    if line.starts_with("!!! ") {
        return true;
    }
    if trimmed.starts_with("- ") || trimmed.starts_with("* ") {
        return true;
    }
    // A numbered item: digits, then `.`, then whitespace.
    let digits: String = trimmed.chars().take_while(char::is_ascii_digit).collect();
    !digits.is_empty() && trimmed[digits.len()..].starts_with(". ")
}

/// Every markdown file under `docs/recommendations`, as repo-relative paths.
fn recommendation_files(root: &Path) -> Result<Vec<String>, String> {
    let base = root.join("docs/recommendations");
    let mut files = Vec::new();
    let mut stack = vec![base.clone()];
    while let Some(dir) = stack.pop() {
        let entries = std::fs::read_dir(&dir)
            .map_err(|e| format!("reading the spec checkout {}: {e}", dir.display()))?;
        for entry in entries {
            let path = entry
                .map_err(|e| format!("reading {}: {e}", dir.display()))?
                .path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().and_then(|e| e.to_str()) == Some("md") {
                let relative = path
                    .strip_prefix(root)
                    .map_err(|_| format!("{} is outside the spec root", path.display()))?;
                files.push(relative.to_string_lossy().replace('\\', "/"));
            }
        }
    }
    files.sort();
    Ok(files)
}

// ---------------------------------------------------------------------------
// The report
// ---------------------------------------------------------------------------

fn render_report(
    corpus: &SpecCorpus,
    counts: &CaptureCounts,
    inventory: &Inventory,
    spec_verified: bool,
) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "# Exhaustive spec-derived conformance corpus\n");
    let _ = writeln!(
        out,
        "bounds: seeds={} extended_scale={}",
        corpus.bounds.seeds, corpus.bounds.extended_scale
    );
    let _ = writeln!(
        out,
        "designs considered: {}  ->  rows: {}  ->  spellings: {}",
        corpus.designs_considered,
        corpus.rows.len(),
        corpus.spellings()
    );
    let _ = writeln!(
        out,
        "completeness: attempted={} succeeded={} dropped={}",
        counts.attempted, counts.succeeded, counts.dropped
    );
    for (reason, count) in &corpus.drops {
        let _ = writeln!(out, "  dropped ({reason}): {count}");
    }
    let _ = writeln!(
        out,
        "\nmulti-member rows: {} of {} ({:.1}%). Real corpora are 592 of 9,949,738 (0.006%).",
        corpus.multi_member_rows(),
        corpus.rows.len(),
        100.0 * corpus.multi_member_rows() as f64 / corpus.rows.len().max(1) as f64
    );

    let table = |title: &str, rows: Vec<(String, String)>, out: &mut String| {
        let _ = writeln!(out, "\n## {title}\n");
        for (key, value) in rows {
            let _ = writeln!(out, "  {key:<34} {value}");
        }
    };

    table(
        "Rows by kind",
        corpus
            .by_kind()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
        &mut out,
    );
    table(
        "Rows by stratum",
        corpus
            .by_stratum()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
        &mut out,
    );
    table(
        "Rows by combining mechanism (never inferred from brackets)",
        corpus
            .by_mechanism()
            .into_iter()
            .map(|(k, v)| {
                (
                    k.to_string(),
                    format!(
                        "{} rows, {} counted multi-member",
                        v.rows, v.multi_member_rows
                    ),
                )
            })
            .collect(),
        &mut out,
    );
    table(
        "Rows by member geometry (#1456)",
        corpus
            .by_geometry()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
        &mut out,
    );
    table(
        "Rows by scale band (#1460)",
        corpus
            .by_scale_band()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
        &mut out,
    );
    table(
        "Rows by transcript region (#1478)",
        corpus
            .by_region()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
        &mut out,
    );
    table(
        "Rows by reference shape",
        corpus
            .by_shape()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
        &mut out,
    );
    table(
        "Prohibited rows by clause strength",
        corpus
            .by_prohibition_strength()
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
        &mut out,
    );
    table(
        "Negative guards (behaviour that would implement a REJECTED proposal)",
        corpus
            .by_negative_guard()
            .into_iter()
            .map(|(k, v)| (k.to_string(), format!("{v} rows")))
            .collect(),
        &mut out,
    );
    table(
        "Rule tags, and whether each is reached in a multi-member row",
        corpus
            .by_rule()
            .into_iter()
            .map(|(k, v)| {
                (
                    k.to_string(),
                    format!("{} rows, {} multi-member", v.rows, v.multi_member_rows),
                )
            })
            .collect(),
        &mut out,
    );

    let _ = writeln!(out, "\n## Declared scale bounds\n");
    let _ = writeln!(out, "  block ladder:      {BLOCK_LADDER:?}");
    let _ = writeln!(out, "  separation ladder: {SEPARATION_LADDER:?}");
    let _ = writeln!(
        out,
        "  extended ladder:   {EXTENDED_BLOCK_LADDER:?} ({})",
        if corpus.bounds.extended_scale {
            "included"
        } else {
            "DECLARED OUT OF BOUNDS — pass --extended-scale"
        }
    );

    // Coverage over the inventory, per file.
    let _ = writeln!(out, "\n## Rule inventory coverage, per spec file\n");
    let _ = writeln!(out, "  clause unit = {}", inventory.clause_unit_definition);
    let _ = writeln!(
        out,
        "  spec checkout verification: {}",
        if spec_verified {
            "quotes and clause-unit counts CHECKED against assets/hgvs-nomenclature"
        } else {
            "SKIPPED (--skip-spec-check) — the coverage table below is unverified"
        }
    );
    let _ = writeln!(
        out,
        "\n  {:<48} {:>5} {:>5} {:>5} {:>5}  blanket classification",
        "file", "units", "cited", "gen", "pub"
    );
    let mut cited_per_file: BTreeMap<&str, (usize, usize, usize)> = BTreeMap::new();
    for entry in &inventory.clauses {
        let file = entry.clause.rsplit_once(':').map_or("?", |(f, _)| f);
        let slot = cited_per_file.entry(file).or_default();
        slot.0 += 1;
        match entry.classification {
            Classification::Generatable => slot.1 += 1,
            Classification::PublishedExampleOnly => slot.2 += 1,
            Classification::NotGenerated => {}
        }
    }
    let mut totals = (0usize, 0usize, 0usize, 0usize);
    for scope in &inventory.files {
        let (cited, generatable, published) = cited_per_file
            .get(scope.file.as_str())
            .copied()
            .unwrap_or_default();
        totals.0 += scope.clause_units;
        totals.1 += cited;
        totals.2 += generatable;
        totals.3 += published;
        let _ = writeln!(
            out,
            "  {:<48} {:>5} {:>5} {:>5} {:>5}  {}",
            scope.file.trim_start_matches("docs/recommendations/"),
            scope.clause_units,
            cited,
            generatable,
            published,
            scope.classification.label()
        );
    }
    let _ = writeln!(
        out,
        "  {:<48} {:>5} {:>5} {:>5} {:>5}",
        "TOTAL", totals.0, totals.1, totals.2, totals.3
    );
    let _ = writeln!(
        out,
        "\n  Every clause unit is accounted for: one of {} explicit clause entries, or the \n  \
         blanket classification of its file with a pinned unit count. An unclassified clause \n  \
         is impossible — the generator refuses.",
        inventory.clauses.len()
    );

    let _ = writeln!(out, "\n## Not generatable, with the reason\n");
    for entry in &inventory.clauses {
        if entry.classification != Classification::NotGenerated {
            continue;
        }
        let _ = writeln!(
            out,
            "  {}\n      {}",
            entry.clause,
            entry.reason.as_deref().unwrap_or("(no reason — refused)")
        );
    }

    let _ = writeln!(out, "\n## What the suite measures over this corpus\n");
    let _ = writeln!(
        out,
        "  tests/it/spec_conformance_axis.rs enumerates the corpus at test time and runs the\n  \
         four properties over it. It is not driven by any artifact this binary writes."
    );
    let by_kind = corpus.by_kind();
    let _ = writeln!(
        out,
        "  family rows (all four properties):     {}",
        by_kind.get(&RowKind::Family).copied().unwrap_or(0)
    );
    let _ = writeln!(
        out,
        "  single rows (validity + idempotency):  {}",
        by_kind.get(&RowKind::Single).copied().unwrap_or(0)
    );
    let _ = writeln!(
        out,
        "  conflict rows (must be refused):       {}",
        by_kind.get(&RowKind::Conflict).copied().unwrap_or(0)
    );
    let _ = writeln!(
        out,
        "  prohibited rows (must be refused):     {}",
        by_kind.get(&RowKind::Prohibited).copied().unwrap_or(0)
    );

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The clause-unit definition the inventory records must be the one the code
    /// implements, or the pinned per-file counts are measuring something else.
    #[test]
    fn clause_units_are_bullets_numbered_items_and_admonitions() {
        assert!(is_clause_unit("- a bullet"));
        assert!(is_clause_unit("    - an indented bullet"));
        assert!(is_clause_unit("* a star bullet"));
        assert!(is_clause_unit("1.  a numbered item"));
        assert!(is_clause_unit("13. a two-digit numbered item"));
        assert!(is_clause_unit("!!! note \"a question\""));
        assert!(!is_clause_unit("# a heading"));
        assert!(!is_clause_unit("ordinary prose"));
        assert!(!is_clause_unit(
            "    !!! an indented admonition is a continuation"
        ));
        assert!(!is_clause_unit("1.2 a decimal is not an item"));
        assert!(!is_clause_unit(""));
    }

    /// A quote that has moved must be reported, and the report must show the
    /// line it found — a mismatch whose message does not say what is there
    /// leaves the reader to guess whether the line or the citation moved.
    #[test]
    fn a_moved_quote_is_reported_with_the_line_it_found() {
        let dir = std::env::temp_dir().join(format!("ferro-spec-corpus-{}", std::process::id()));
        let file = dir.join("docs/x.md");
        std::fs::create_dir_all(file.parent().expect("parent")).expect("mkdir");
        std::fs::write(&file, "first\nsecond line here\nthird\n").expect("write");

        assert!(check_quote(&dir, "docs/x.md:2", "second line").is_ok());
        let error = check_quote(&dir, "docs/x.md:2", "absent")
            .expect_err("a quote that is not on the line must be refused");
        assert!(error.contains("no longer on that line"), "{error}");
        assert!(error.contains("second line here"), "{error}");
        assert!(check_quote(&dir, "docs/x.md:99", "anything").is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }
}
