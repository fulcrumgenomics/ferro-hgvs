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
    enumerate, CorpusBounds, RowKind, SpecCorpus, BLOCK_LADDER, DNA_SYMBOLS, EXTENDED_BLOCK_LADDER,
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
    // Refuse a `FERRO_PARTITION` naming no arm before any work, so a
    // misspelling cannot be served by the shipped rule under a candidate's
    // name. See `tests/it/partition_switch_wiring.rs`.
    if let Some(message) = ferro_hgvs::normalize::partition_switch_startup_error() {
        eprintln!("error: {message}");
        return ExitCode::FAILURE;
    }
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
            // The alphabet the corpus draws from is a hand-transcribed copy of
            // `background/standards.md`'s DNA table, and the clause-unit
            // definition cannot see a markdown table. So pin it directly.
            if let Err(error) = check_dna_symbol_table(root) {
                problems.push(error);
            }
            // `standards.md` is the only scoped file that currently carries a
            // table (verified 2026-08-09, see `FILES_KNOWN_TO_CARRY_TABLES`).
            // A submodule bump could add one to any other file with every
            // pinned `clause_units` count staying exactly the same, since a
            // clause unit is a bullet/numbered item/admonition and never a
            // table row. Pin the absence directly, so that silently enlarges
            // into a loud refusal instead.
            if let Err(error) = check_no_undeclared_tables(root, &self.files) {
                problems.push(error);
            }
            // Every inventoried file must have a scope row, so "every clause"
            // has a mechanical denominator.
            match inventoried_files(root) {
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

/// The DNA nucleotide symbols `background/standards.md` tabulates, split into
/// the ones a description may state and the ones its `†` footnote marks as
/// "used in alignment only".
///
/// Returned in table order, because [`spec_corpus::DNA_SYMBOLS`] is indexed
/// positionally (`&DNA_SYMBOLS[4..]` skips the four unambiguous bases), so a
/// reordering of the table is a change to what the ambiguity stratum draws.
fn dna_symbol_table(root: &Path) -> Result<(Vec<char>, Vec<char>), String> {
    let path = root.join("docs/background/standards.md");
    let text = std::fs::read_to_string(&path)
        .map_err(|e| format!("reading the symbol table {}: {e}", path.display()))?;
    let mut lines = text.lines().skip_while(|line| line.trim() != "### DNA");
    lines.next();
    let (mut usable, mut alignment_only) = (Vec::new(), Vec::new());
    let mut seen_table = false;
    let mut rows_seen = 0usize;
    for line in lines {
        let trimmed = line.trim();
        if !trimmed.starts_with('|') {
            if seen_table {
                break;
            }
            continue;
        }
        let cell = trimmed
            .trim_matches('|')
            .split('|')
            .next()
            .unwrap_or_default()
            .trim();
        // Find the table by its HEADER, not by "the first `|`-led line after
        // `### DNA`". Any pipe-bearing line in between — a footnote, a nested
        // list continuation — otherwise becomes row 0, and the error below then
        // names that line as a malformed header rather than saying the table was
        // not found. The behaviour stayed loud either way; the diagnosis was
        // wrong, which is worse than a bare failure.
        //
        // Locating the header this way does NOT weaken the positional rule the
        // comment below defends: `Symbol` is a column name and never a DNA
        // symbol, so no real row can be mistaken for it. Once the header is
        // found, the underline and every symbol row are still recognised by
        // POSITION relative to it, which is the part that must not become a
        // content test.
        if !seen_table {
            if cell != "Symbol" {
                continue;
            }
            seen_table = true;
        }
        let row = rows_seen;
        rows_seen += 1;
        // The header row and the alignment underline directly beneath it are
        // not symbols. Both are recognised by POSITION, never by content: `-`
        // is itself a DNA symbol in this table, and it survives a content test
        // only because it currently carries a `†`. Drop that footnote upstream
        // and a content test would discard a real symbol silently — `usable`
        // would shrink with nothing failing, which is precisely the silent
        // shrink this whole check exists to prevent.
        if row == 0 {
            if cell != "Symbol" {
                return Err(format!(
                    "docs/background/standards.md: the DNA table under `### DNA` does not open \
                     with a `Symbol` header (found `{cell}`) — the table was located by that \
                     header, so this means it moved mid-parse"
                ));
            }
            continue;
        }
        if row == 1 {
            if cell.is_empty() || !cell.chars().all(|c| c == ':' || c == '-') {
                return Err(format!(
                    "docs/background/standards.md: expected an alignment underline directly \
                     under the DNA table header, found `{cell}`"
                ));
            }
            continue;
        }
        let footnoted = cell.ends_with('†');
        let symbol = cell.trim_end_matches('†');
        let mut chars = symbol.chars();
        let (Some(first), None) = (chars.next(), chars.next()) else {
            return Err(format!(
                "docs/background/standards.md: `{symbol}` is not a single-character DNA symbol"
            ));
        };
        if footnoted {
            alignment_only.push(first);
        } else {
            usable.push(first);
        }
    }
    if usable.is_empty() {
        return Err(
            "docs/background/standards.md: found no DNA symbol table under `### DNA`".to_string(),
        );
    }
    Ok((usable, alignment_only))
}

/// Check that the corpus's alphabet is the spec's, and that the symbols it
/// withholds are exactly the ones the spec footnotes as alignment-only.
///
/// This exists because the clause-unit denominator cannot reach it. A clause
/// unit is a bullet, a numbered item or an admonition, so `standards.md`'s seven
/// counted units are its table of contents — the symbol tables themselves are
/// markdown tables and are invisible to that count. Pinning the file's unit
/// count therefore says nothing about the alphabet, while the corpus reads its
/// alphabet from precisely there.
fn check_dna_symbol_table(root: &Path) -> Result<(), String> {
    let (usable, alignment_only) = dna_symbol_table(root)?;
    if usable != DNA_SYMBOLS {
        return Err(format!(
            "the corpus alphabet is not the spec's DNA symbol table.\n      \
             DNA_SYMBOLS: {DNA_SYMBOLS:?}\n      standards.md: {usable:?}"
        ));
    }
    if !alignment_only.contains(&'X') {
        return Err(format!(
            "docs/background/standards.md no longer footnotes `X` as alignment-only, but the \
             corpus still emits it as a prohibited base (`standards.md:39-alignment-only-symbols`). \
             Footnoted: {alignment_only:?}"
        ));
    }
    Ok(())
}

/// Scoped files already known to carry a markdown or HTML table.
///
/// `docs/background/standards.md` is the only one today: [`check_dna_symbol_table`]
/// reads its DNA table directly, and the corpus has no dependency on its RNA
/// table, Genetic Code table or Amino Acid Descriptions table, or on the ISCN
/// cytogenetic band tables (`git grep` for `RNA_SYMBOLS`, `GENETIC_CODE`,
/// `cytoBand` etc. across `src/` and `examples/` turns up nothing — those are
/// the declared protein-axis and chromosome-scale gaps recorded on `protein/`
/// and `DNA/complex.md`'s scope rows, not a corpus dependency).
///
/// Every other file in the denominator was hand-verified table-free against
/// the pinned submodule checkout on 2026-08-09 (`find docs/background
/// docs/recommendations -name '*.md' | xargs grep -lE '^\s*\|.*\|'` and the
/// same for `<table` matched only this file). Adding a file to this list is
/// itself a change that must be reviewed on its own: it means a submodule
/// bump added tabular content the clause-unit denominator cannot see, and
/// whoever adds the row must first decide whether the corpus needs a
/// table-reading guard for it — the way `check_dna_symbol_table` reads this
/// file — before declaring it accounted for.
const FILES_KNOWN_TO_CARRY_TABLES: &[&str] = &["docs/background/standards.md"];

/// Fail if any scoped file **outside** [`FILES_KNOWN_TO_CARRY_TABLES`] carries
/// markdown or HTML table syntax.
///
/// This is the mechanized half of the 2026-08-09 table-blindness audit: a
/// clause unit is a bullet, a numbered item or an admonition
/// ([`is_clause_unit`]), so a pinned `clause_units` count cannot detect a
/// table being added to (or growing inside) a file that does not yet have
/// one — the count would stay exactly the same while a whole new rule surface
/// went uncovered, the way `standards.md`'s DNA alphabet did before
/// `check_dna_symbol_table` existed. Rather than re-run that audit by hand on
/// every submodule bump, pin its negative result: today's table-free files
/// must stay table-free, or the build refuses and names the file and line so
/// a human re-does the audit on exactly the file that changed.
fn check_no_undeclared_tables(root: &Path, files: &[FileScope]) -> Result<(), String> {
    let mut problems = Vec::new();
    for scope in files {
        if FILES_KNOWN_TO_CARRY_TABLES.contains(&scope.file.as_str()) {
            continue;
        }
        let path = root.join(&scope.file);
        // Accumulated, not `?`-propagated. This function's whole shape is an
        // audit that visits every scope and reports the problems together, and
        // a `?` here abandoned that on the first unreadable file — discarding
        // the problems already found and never visiting the remaining scopes.
        // The result was an audit that reported "reading X: …" while any number
        // of undeclared tables sat unlooked-at behind it.
        let text = match std::fs::read_to_string(&path) {
            Ok(text) => text,
            Err(e) => {
                problems.push(format!("reading {}: {e}", path.display()));
                continue;
            }
        };
        if let Some(line) = first_table_line(&text) {
            problems.push(format!(
                "{}:{line} now carries table syntax (a markdown `|` row or an HTML `<table>`), \
                 which the clause-unit count cannot see. Audit it for a rule-bearing table the \
                 way `check_dna_symbol_table` reads `standards.md`'s DNA table, then add this \
                 file to `FILES_KNOWN_TO_CARRY_TABLES` once the table either needs no guard or \
                 has one.",
                scope.file
            ));
        }
    }
    if problems.is_empty() {
        Ok(())
    } else {
        Err(problems.join("\n  "))
    }
}

/// The 1-based line number of the first table-like line in `text`, if any.
///
/// A markdown table row is recognised the same way [`dna_symbol_table`] reads
/// one: a line whose content, after leading whitespace, starts with `|`. That
/// is deliberately loose (it does not require a following header-separator
/// row), because the audit this backs is "does anything here look like table
/// syntax", not "is this specifically a well-formed table" — a `|`-led line
/// that is not a genuine table is itself worth a human's attention on a file
/// this inventory says has none. An HTML `<table` tag (case-insensitive, as
/// `background/standards.md`'s Genetic Code table opens) is the other form
/// this repository's spec checkout uses for tabular content.
fn first_table_line(text: &str) -> Option<usize> {
    text.lines().enumerate().find_map(|(index, line)| {
        let trimmed = line.trim_start();
        let is_table_row =
            trimmed.starts_with('|') || trimmed.to_ascii_lowercase().starts_with("<table");
        is_table_row.then_some(index + 1)
    })
}

/// Every markdown file the inventory is a denominator over, as repo-relative
/// paths: `docs/recommendations` plus `docs/background`.
///
/// `docs/background` is included because the corpus does not merely cite it, it
/// is DERIVED from it — the alphabet from `standards.md`, the `c.` coordinate
/// zones the region axis walks from `numbering.md`. `docs/consultation` is
/// deliberately not here; see the inventory's own `description`.
fn inventoried_files(root: &Path) -> Result<Vec<String>, String> {
    let mut files = Vec::new();
    let mut stack = vec![
        root.join("docs/recommendations"),
        root.join("docs/background"),
    ];
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
            scope
                .file
                .trim_start_matches("docs/recommendations/")
                .trim_start_matches("docs/"),
            scope.clause_units,
            cited,
            generatable,
            published,
            scope.classification.label()
        );
    }
    // Files that are CITED but carry no scope row — `consultation/` has none by
    // design. Without this the loop above computes their counts and then drops
    // them, so a clause is verified and then vanishes from the accounting and
    // `TOTAL` under-reports. The quote was checked either way; what was wrong is
    // that the table reads as a complete accounting and was not one.
    //
    // `units` is `-` rather than `0`: these files have no pinned unit count
    // because they have no scope row, and printing `0` would claim they contain
    // no clause units, which is a different and false statement.
    let scoped: BTreeSet<&str> = inventory.files.iter().map(|s| s.file.as_str()).collect();
    for (file, (cited, generatable, published)) in &cited_per_file {
        if scoped.contains(file) {
            continue;
        }
        totals.1 += cited;
        totals.2 += generatable;
        totals.3 += published;
        let _ = writeln!(
            out,
            "  {:<48} {:>5} {:>5} {:>5} {:>5}  cited only, no scope row",
            file.trim_start_matches("docs/recommendations/")
                .trim_start_matches("docs/"),
            "-",
            cited,
            generatable,
            published
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

    /// The symbol table is read out of a markdown table, which the clause-unit
    /// denominator cannot see — so the parser that replaces it has to be pinned
    /// on its own. Both halves matter: which symbols a description may state,
    /// and which the `†` footnote withholds.
    #[test]
    fn the_dna_symbol_table_is_read_with_its_alignment_only_footnote() {
        // A `TempDir`, not a PID-named path under `temp_dir()`: the manual
        // `remove_dir_all` those used is skipped on a panicking assertion, so a
        // failing run leaked its scratch tree and the next one inherited it.
        // Drop-based cleanup runs on the unwind too.
        let scratch = tempfile::TempDir::new().expect("scratch dir");
        let dir = scratch.path();
        let file = dir.join("docs/background/standards.md");
        std::fs::create_dir_all(file.parent().expect("parent")).expect("mkdir");
        std::fs::write(
            &file,
            "### RNA\n\n| Symbol |\n|:--:|\n|   a    |\n\n### DNA\n\n\
             | Symbol |   Meaning    |\n|:------:|:------------:|\n\
             |   A    |      A       |\n|   N    | A, C, G or T |\n\
             |   X†   | A, C, G or T |\n|   -†   |     none     |\n\n\
             † used in alignment only.\n\n## Genetic Code\n",
        )
        .expect("write");

        let (usable, alignment_only) = dna_symbol_table(dir).expect("the DNA table");
        // Table order, not sorted: `DNA_SYMBOLS` is indexed positionally.
        assert_eq!(usable, vec!['A', 'N']);
        assert_eq!(alignment_only, vec!['X', '-']);
        // The preceding `### RNA` table must not leak in, and the prose after
        // the table must end it.
        assert!(!usable.contains(&'a'), "{usable:?}");

        // A checkout with no such table is a refusal, not an empty answer.
        std::fs::write(&file, "### DNA\n\nno table here\n").expect("rewrite");
        assert!(dna_symbol_table(dir).is_err());
    }

    /// A quote that has moved must be reported, and the report must show the
    /// line it found — a mismatch whose message does not say what is there
    /// leaves the reader to guess whether the line or the citation moved.
    #[test]
    fn a_moved_quote_is_reported_with_the_line_it_found() {
        let scratch = tempfile::TempDir::new().expect("scratch dir");
        let dir = scratch.path();
        let file = dir.join("docs/x.md");
        std::fs::create_dir_all(file.parent().expect("parent")).expect("mkdir");
        std::fs::write(&file, "first\nsecond line here\nthird\n").expect("write");

        assert!(check_quote(dir, "docs/x.md:2", "second line").is_ok());
        let error = check_quote(dir, "docs/x.md:2", "absent")
            .expect_err("a quote that is not on the line must be refused");
        assert!(error.contains("no longer on that line"), "{error}");
        assert!(error.contains("second line here"), "{error}");
        assert!(check_quote(dir, "docs/x.md:99", "anything").is_err());
    }

    /// [`first_table_line`] must find both table syntaxes this spec checkout
    /// actually uses — a markdown `|`-led row and an HTML `<table>` tag — and
    /// must not fire on the literal `|` character HGVS uses for methylation
    /// state changes (`general.md:94`, `DNA/other.md:42-47`) or on the `|`
    /// EBNF alternation operator (`grammar.md:16-21`), neither of which opens
    /// a line.
    #[test]
    fn first_table_line_finds_markdown_and_html_tables_but_not_a_bare_pipe_character() {
        assert_eq!(first_table_line("prose\n| a | b |\n|---|---|\n"), Some(2));
        assert_eq!(first_table_line("prose\n  <table class=\"gc\">\n"), Some(2));
        assert_eq!(
            first_table_line("prose\n  <TABLE>\n"),
            Some(2),
            "the check must be case-insensitive"
        );
        assert_eq!(
            first_table_line(
                "- `|` (pipe) is used to indicate a change of state.\n\
                 - A \"pipe\" (`|`) separates alternatives, e.g. `\"A\" | \"B\"`.\n"
            ),
            None,
            "a `|` that does not open the line is not a table row"
        );
        assert_eq!(first_table_line("prose only\nno tables here\n"), None);
    }

    /// The guard this audit adds: a file the inventory does not already list
    /// in `FILES_KNOWN_TO_CARRY_TABLES` must stay table-free, and a table that
    /// appears in one must fail the build by name and line rather than
    /// passing silently under an unchanged `clause_units` count.
    #[test]
    fn check_no_undeclared_tables_flags_a_table_in_a_file_not_on_the_allow_list() {
        let scratch = tempfile::TempDir::new().expect("scratch dir");
        let dir = scratch.path();
        let background = dir.join("docs/background");
        std::fs::create_dir_all(&background).expect("mkdir");
        std::fs::write(
            background.join("standards.md"),
            "the allow-listed file may carry a table\n\n| a | b |\n|---|---|\n",
        )
        .expect("write standards.md");
        std::fs::write(
            background.join("numbering.md"),
            "- a bullet is fine\n- so is another\n",
        )
        .expect("write numbering.md");

        let files = vec![
            FileScope {
                file: "docs/background/standards.md".to_string(),
                clause_units: 7,
                classification: Classification::NotGenerated,
                reason: None,
            },
            FileScope {
                file: "docs/background/numbering.md".to_string(),
                clause_units: 2,
                classification: Classification::NotGenerated,
                reason: None,
            },
        ];
        // Both table-free by this file set's lights: `standards.md` is on the
        // allow list (so its table is skipped, not read), and `numbering.md`
        // carries no table syntax at all.
        assert!(check_no_undeclared_tables(dir, &files).is_ok());

        // A submodule bump adds a table to a file the inventory does not
        // expect one in.
        std::fs::write(
            background.join("numbering.md"),
            "- a bullet is fine\n\n| newly | added |\n|---|---|\n",
        )
        .expect("rewrite numbering.md");
        let error = check_no_undeclared_tables(dir, &files)
            .expect_err("a table outside the allow list must be refused");
        assert!(error.contains("docs/background/numbering.md:3"), "{error}");
        // The allow-listed file's own table must not be (mis)reported as a
        // problem — only the boilerplate guidance, which names
        // `standards.md` as the model to follow, may mention it.
        assert!(!error.contains("docs/background/standards.md:"), "{error}");
    }
}
