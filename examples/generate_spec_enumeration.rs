//! Generator for the exhaustive, non-redundant HGVS spec **test enumeration**.
//!
//! `generate_spec_fixture` covers exactly one slice of the spec corpus:
//! *harvested string × `normalize()` × the default error mode*. This generator
//! enumerates the cross-product that slice does **not** touch, and emits only
//! the rows that are not already asserted somewhere else:
//!
//! | dimension          | what it adds                                                   |
//! |--------------------|----------------------------------------------------------------|
//! | `negative-reject`  | spec-stated negatives (`<code class="invalid">` + prose) as *reject* assertions |
//! | `negative-repair`  | the same negatives as *repair-to-canonical* assertions          |
//! | `error-mode`       | `parse` under strict / lenient / silent, kept only where the modes actually differ |
//! | `grammar-form`     | `syntax.yaml` form × legal coordinate system                    |
//! | `output-invariant` | the MUST-level output rules checked over every emitted output   |
//!
//! Every row carries its provenance (spec `file:line@SHA`, the dimension that
//! produced it, and — for negatives — the RFC 2119 normative level), and a
//! `dedup_note` recording why it is not a restatement of existing coverage.
//!
//! Like `hgvs_spec_normalization.json`, the output is a **generated build
//! artifact**: gitignored, regenerated on demand. Only the hand-curated
//! overrides, this generator, and the test driver are committed.
//!
//! Run:   `cargo run --features dev --example generate_spec_enumeration`
//! Check: `cargo run --features dev --example generate_spec_enumeration -- --check`
//! Census: `cargo run --features dev --example generate_spec_enumeration -- --census`

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser as ClapParser;
use serde::{Deserialize, Serialize};

#[path = "common/spec_harvest.rs"]
mod spec_harvest;

use spec_harvest::{negatives, prefix, sources};

#[derive(ClapParser, Debug)]
#[command(about = "Generate the exhaustive HGVS spec test enumeration")]
struct Cli {
    /// Path to the vendored spec checkout.
    #[arg(long, default_value = "assets/hgvs-nomenclature")]
    spec_dir: PathBuf,

    /// Hand-curated overrides (repair targets, per-row pins).
    #[arg(
        long,
        default_value = "tests/fixtures/grammar/hgvs_spec_enumeration_overrides.json"
    )]
    overrides: PathBuf,

    /// Existing normalization fixture, read for deduplication.
    #[arg(
        long,
        default_value = "tests/fixtures/grammar/hgvs_spec_normalization.json"
    )]
    normalization_fixture: PathBuf,

    /// Output path for the generated enumeration.
    #[arg(
        long,
        default_value = "tests/fixtures/grammar/hgvs_spec_enumeration.json"
    )]
    output: PathBuf,

    /// Regenerate in memory and exit non-zero if the on-disk file differs.
    #[arg(long)]
    check: bool,

    /// Print the per-dimension gross / duplicate / net census to stderr.
    #[arg(long)]
    census: bool,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: &Cli) -> anyhow::Result<()> {
    let candidates = sources::discover(&cli.spec_dir)?;
    let prose = negatives::discover(&cli.spec_dir)?;
    let overrides = overrides::load(&cli.overrides)?;
    let existing = dedup::ExistingCoverage::load(&cli.normalization_fixture)?;
    let sha = spec_sha(&cli.spec_dir)?;

    let build = builder::build(
        &cli.spec_dir,
        &candidates,
        &prose,
        &overrides,
        &existing,
        &sha,
    )?;
    let rendered = render::render(&build, &sha)?;

    if cli.census {
        build.census.report(&mut std::io::stderr())?;
    }

    if cli.check {
        let on_disk = std::fs::read_to_string(&cli.output)
            .map_err(|e| anyhow::anyhow!("read {}: {e}", cli.output.display()))?;
        if on_disk != rendered {
            eprintln!(
                "enumeration {} is out of date; rerun: \
                 cargo run --features dev --example generate_spec_enumeration",
                cli.output.display()
            );
            return Err(anyhow::anyhow!("enumeration out of date"));
        }
    } else {
        std::fs::write(&cli.output, rendered)?;
    }
    Ok(())
}

/// Resolve the pinned spec submodule commit, short form.
///
/// `GIT_*` is cleared for the same reason as in `generate_spec_fixture`: git
/// exports `GIT_DIR`/`GIT_WORK_TREE` pointing at the OUTER repo when it runs a
/// hook, which would otherwise win over `-C spec_dir` (#1046).
fn spec_sha(spec_dir: &Path) -> anyhow::Result<String> {
    let out = std::process::Command::new("git")
        .arg("-C")
        .arg(spec_dir)
        .args(["rev-parse", "--short", "HEAD"])
        .env_remove("GIT_DIR")
        .env_remove("GIT_WORK_TREE")
        .env_remove("GIT_INDEX_FILE")
        .output()
        .map_err(|e| anyhow::anyhow!("git rev-parse: {e}"))?;
    if !out.status.success() {
        anyhow::bail!(
            "git rev-parse failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    Ok(String::from_utf8_lossy(&out.stdout).trim().to_string())
}

// ---------- row ----------

/// One enumerated assertion.
///
/// The driver (`tests/it/spec_enumeration_tests.rs`) replays each row: it
/// recomputes `observed` and requires it to match the pinned value, and
/// separately requires that no *new* MUST-level violation appears beyond the
/// pinned set. That keeps the suite green while still recording divergence.
#[derive(Debug, Serialize)]
pub struct Row {
    /// Stable, unique identity: `<dimension>/<discriminator>`.
    pub id: String,
    pub dimension: String,
    /// The ferro operation under test: `parse` / `normalize` / `preprocess` /
    /// `display-roundtrip` / `invariant-check`.
    pub operation: String,
    /// `strict` / `lenient` / `silent`, or `default` when the row does not vary
    /// the error mode.
    pub error_mode: String,
    /// Verbatim spec text (or the syntax.yaml example / derived form).
    pub input: String,
    /// What ferro is actually run against — `input` with a default accession
    /// prepended when the spec wrote a bare illustrative fragment.
    pub target: String,
    /// `spec-mandated` when the spec states the expected value outright;
    /// `pinned-baseline` when it does not and the row records current
    /// behaviour instead. **Never** invent a spec expectation.
    pub expectation: String,
    /// RFC 2119 level. Only `must` rows may ever hard-fail.
    pub normative_level: String,
    /// The spec-stated expected value, when there is one.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expected: Option<String>,
    /// Ferro's behaviour at generation time.
    pub observed: String,
    pub status: String,
    /// `docs/recommendations/<file>.md:<line>@<sha>` — or the artifact the row
    /// was derived from when the spec has no single line to cite.
    pub spec_citation: String,
    /// Why this row is not a restatement of an existing assertion.
    pub dedup_note: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
}

// ---------- dedup ----------

/// Deduplication against coverage the repo already has.
///
/// The rule, stated once and applied uniformly: a generated row is a
/// **duplicate** iff the tuple `(operation, error_mode, target)` is already
/// asserted by an existing test. Concretely, the 1033 rows of
/// `hgvs_spec_normalization.json` already assert, for every one of their
/// targets:
///
/// * `("normalize", "default", target)` — `pinned_v21_normalization_behavior`
///   pins `normalize()`'s exact output and warning codes;
/// * `("reject", "default", target)` for rows whose `spec_expected` is `null` —
///   the same test pins the rejection;
/// * `("display-roundtrip", "default", target)` —
///   `idempotency_tests::test_parsing_idempotency_spec_examples` parses every
///   fixture target, re-renders and re-parses it;
/// * `("idempotency", "default", target)` —
///   `idempotency_tests::test_normalization_idempotency_reference_free_corpus`
///   normalizes every fixture target twice.
///
/// Those four are therefore *saturated* over the harvested corpus and this
/// generator emits nothing for them. What remains genuinely open is recorded in
/// the census.
mod dedup {
    use super::*;

    #[derive(Debug, Deserialize)]
    struct FixtureDoc {
        rows: Vec<FixtureRow>,
    }

    #[derive(Debug, Deserialize)]
    struct FixtureRow {
        input: String,
        #[serde(default)]
        input_prefixed: Option<String>,
        current: String,
        spec_expected: Option<String>,
        status: String,
    }

    #[derive(Debug, Default)]
    pub struct ExistingCoverage {
        /// `(operation, error_mode, target)` tuples already asserted.
        keys: BTreeSet<(String, String, String)>,
        /// target → ferro's pinned `normalize()` rendering.
        pub current_by_target: BTreeMap<String, String>,
        /// target → pinned status.
        pub status_by_target: BTreeMap<String, String>,
        /// Every raw spec input the normalization fixture retained.
        pub inputs: BTreeSet<String>,
    }

    /// Operations the normalization fixture + idempotency suite already
    /// saturate over every harvested target.
    pub const SATURATED_OPS: &[&str] = &["normalize", "display-roundtrip", "idempotency"];

    impl ExistingCoverage {
        pub fn load(path: &Path) -> anyhow::Result<Self> {
            let text = std::fs::read_to_string(path).map_err(|e| {
                anyhow::anyhow!(
                    "read {}: {e}\nHint: regenerate it first with \
                     `cargo run --features dev --example generate_spec_fixture`",
                    path.display()
                )
            })?;
            let doc: FixtureDoc = serde_json::from_str(&text)
                .map_err(|e| anyhow::anyhow!("parse {}: {e}", path.display()))?;
            let mut cov = ExistingCoverage::default();
            for r in &doc.rows {
                let target = r.input_prefixed.clone().unwrap_or_else(|| r.input.clone());
                for op in SATURATED_OPS {
                    cov.keys
                        .insert(((*op).to_string(), "default".to_string(), target.clone()));
                }
                if r.spec_expected.is_none() {
                    cov.keys
                        .insert(("reject".to_string(), "default".to_string(), target.clone()));
                }
                cov.current_by_target
                    .insert(target.clone(), r.current.clone());
                cov.status_by_target.insert(target, r.status.clone());
                cov.inputs.insert(r.input.clone());
            }
            Ok(cov)
        }

        pub fn covers(&self, operation: &str, error_mode: &str, target: &str) -> bool {
            self.keys.contains(&(
                operation.to_string(),
                error_mode.to_string(),
                target.to_string(),
            ))
        }

        /// Every HGVS string ferro actually **emits** for the harvested corpus,
        /// as `(target, rendered, pinned_status)`.
        ///
        /// This is deliberately wider than the blessed `preserved` set: a
        /// `false-acceptance` row is precisely a case where ferro took a string
        /// the spec forbids and rendered *something* — and what it rendered is
        /// exactly what an output invariant must judge. Error strings are
        /// excluded (there is no output to check).
        pub fn emitted_outputs(&self) -> Vec<(String, String, String)> {
            self.current_by_target
                .iter()
                .filter(|(_, c)| !c.starts_with("parse error") && !c.starts_with("normalize error"))
                .map(|(t, c)| {
                    let status = self
                        .status_by_target
                        .get(t)
                        .cloned()
                        .unwrap_or_else(|| "unknown".to_string());
                    (t.clone(), c.clone(), status)
                })
                .collect()
        }
    }
}

// ---------- overrides ----------

mod overrides {
    use super::*;

    #[derive(Debug, Default, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct Overrides {
        /// Prose explaining the file's contract to the next auditor. Ignored by
        /// the generator; present so `deny_unknown_fields` does not reject it.
        #[serde(default, rename = "_readme")]
        #[allow(dead_code)]
        pub readme: Vec<String>,
        /// Hand-curated repair targets: a spec-stated negative mapped to the
        /// canonical form the **spec itself** gives on the same line. Curated
        /// by hand precisely because automatic extraction is unreliable — the
        /// spec routinely cites a bare coordinate on the same line as the
        /// corrected variant.
        #[serde(default)]
        pub repairs: BTreeMap<String, Repair>,
        /// Per-row pins keyed by `Row::id`.
        #[serde(default)]
        pub by_id: BTreeMap<String, RowOverride>,
    }

    #[derive(Debug, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct Repair {
        /// The canonical form the spec states. `null` means the spec forbids
        /// the input but states no replacement, so no repair row is emitted.
        pub repaired: Option<String>,
        /// `docs/…md:<line>` for the sentence that states the pair.
        pub citation: String,
        /// RFC 2119 level of the sentence: `must` or `should`.
        pub normative_level: String,
        /// Set when evaluating the repair needs real reference bases (e.g. any
        /// 3'-rule example). Such rows are emitted but marked, and the driver
        /// does not assert them under the hermetic mock provider.
        #[serde(default)]
        pub requires_reference: bool,
        #[serde(default)]
        pub note: Option<String>,
    }

    #[derive(Debug, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct RowOverride {
        #[serde(default)]
        pub status: Option<String>,
        #[serde(default)]
        pub note: Option<String>,
    }

    pub fn load(path: &Path) -> anyhow::Result<Overrides> {
        if !path.exists() {
            return Ok(Overrides::default());
        }
        let text = std::fs::read_to_string(path)?;
        serde_json::from_str(&text).map_err(|e| anyhow::anyhow!("parse {}: {e}", path.display()))
    }
}

// ---------- census ----------

mod census {
    use super::*;

    /// Per-dimension gross / duplicate / not-applicable / net accounting.
    #[derive(Debug, Default, Serialize)]
    pub struct DimensionCensus {
        /// Every assertion the dimension can enumerate before any filtering.
        pub gross: usize,
        /// Rows dropped because `(operation, error_mode, target)` is already
        /// asserted by existing tests.
        pub duplicate: usize,
        /// Rows dropped because the assertion is not legally applicable to the
        /// construct (documented per dimension).
        pub not_applicable: usize,
        /// `gross - duplicate - not_applicable`.
        pub net: usize,
    }

    #[derive(Debug, Default, Serialize)]
    pub struct Census {
        pub by_dimension: BTreeMap<String, DimensionCensus>,
    }

    impl Census {
        pub fn tally(&mut self, dimension: &str, gross: usize, dup: usize, na: usize) {
            let e = self.by_dimension.entry(dimension.to_string()).or_default();
            e.gross += gross;
            e.duplicate += dup;
            e.not_applicable += na;
            e.net = e.gross.saturating_sub(e.duplicate + e.not_applicable);
        }

        pub fn report(&self, w: &mut impl std::io::Write) -> anyhow::Result<()> {
            writeln!(
                w,
                "{:<20} {:>8} {:>8} {:>8} {:>8}",
                "dimension", "gross", "dup", "n/a", "net"
            )?;
            let (mut g, mut d, mut n, mut net) = (0, 0, 0, 0);
            for (k, v) in &self.by_dimension {
                writeln!(
                    w,
                    "{:<20} {:>8} {:>8} {:>8} {:>8}",
                    k, v.gross, v.duplicate, v.not_applicable, v.net
                )?;
                g += v.gross;
                d += v.duplicate;
                n += v.not_applicable;
                net += v.net;
            }
            writeln!(w, "{:<20} {:>8} {:>8} {:>8} {:>8}", "TOTAL", g, d, n, net)?;
            Ok(())
        }
    }
}

// ---------- invariants ----------

/// MUST- and SHOULD-level **output** invariants, checked over rendered HGVS
/// strings.
///
/// Scope discipline, stated up front: this is a *well-formedness* tripwire, not
/// a correctness oracle. It answers "is this legal HGVS?", never "is this the
/// right answer". Only rules decidable from the rendered string without
/// reference bases are implemented; anything needing sequence context (3'-rule
/// shiftedness, dup-over-ins prioritisation, frameshift `Ter` counts) is
/// deliberately absent rather than approximated.
mod invariants {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct Violation {
        pub rule_id: &'static str,
        pub level: &'static str,
        pub message: String,
        pub citation: &'static str,
    }

    /// A single member of a variant: its location text and its edit.
    #[derive(Debug)]
    struct Member {
        axis: char,
        loc: String,
        kw: String,
        arg: String,
    }

    const EDIT_KEYWORDS: &[&str] = &["delins", "del", "dup", "ins", "inv", "ext", "fs"];

    /// Split a rendered variant into per-member `(axis, location, edit)` triples.
    /// Members are the `;`/`,`-separated entries of an allele list, or the whole
    /// body when there is no list. Returns an empty vector for shapes the
    /// splitter does not understand — the checker then reports nothing, which is
    /// the safe direction (no false positives).
    fn members(rendered: &str) -> Vec<Member> {
        let mut out = Vec::new();
        // Split off any accession. The axis letter is the char before the first
        // `.` that follows a `:`; bare forms start with `<axis>.`.
        let body_start = match rendered.find(':') {
            Some(i) => i + 1,
            None => 0,
        };
        let body = &rendered[body_start..];
        let Some(dot) = body.find('.') else {
            return out;
        };
        let axis = body[..dot].chars().next_back().unwrap_or('?');
        if !matches!(axis, 'c' | 'g' | 'n' | 'm' | 'o' | 'r' | 'p') {
            return out;
        }
        let rest = &body[dot + 1..];
        // Nested references (`ins[NC_…:g.…]`) and uncertainty parens defeat the
        // naive splitter; skip those shapes wholesale.
        if rest.contains(':') {
            return out;
        }
        for raw in rest.split([';', ',']) {
            let m = raw.trim_matches(|c| c == '[' || c == ']' || c == '(' || c == ')');
            if m.is_empty() {
                continue;
            }
            let Some((idx, kw)) = EDIT_KEYWORDS
                .iter()
                .filter_map(|k| m.find(*k).map(|i| (i, *k)))
                .min_by_key(|(i, _)| *i)
            else {
                // No edit keyword: substitution (`>`), identity (`=`) or a
                // protein sub. Handled by the string-level rules below.
                out.push(Member {
                    axis,
                    loc: m.to_string(),
                    kw: String::new(),
                    arg: String::new(),
                });
                continue;
            };
            out.push(Member {
                axis,
                loc: m[..idx].to_string(),
                kw: kw.to_string(),
                arg: m[idx + kw.len()..].to_string(),
            });
        }
        out
    }

    /// Both endpoints of `loc` as plain integers, or `None` when either carries
    /// an offset / UTR marker / uncertainty (which this checker never judges).
    fn plain_range(loc: &str) -> Option<(u64, u64)> {
        let (a, b) = loc.split_once('_')?;
        if a.is_empty() || b.is_empty() {
            return None;
        }
        Some((a.parse().ok()?, b.parse().ok()?))
    }

    fn is_plain_position(loc: &str) -> bool {
        !loc.is_empty() && loc.bytes().all(|b| b.is_ascii_digit())
    }

    /// Check one rendered variant string against the catalog.
    pub fn check(rendered: &str) -> Vec<Violation> {
        let mut v = Vec::new();
        let ms = members(rendered);

        // B7 — no composite `dupinv` / `dupins`.
        if rendered.contains("dupinv") || rendered.contains("dupins") {
            v.push(Violation {
                rule_id: "B7-no-composite-dup",
                level: "must",
                message: format!("composite dup edit in {rendered}"),
                citation: "docs/recommendations/DNA/duplication.md:92",
            });
        }

        // B4 — substitution is exactly one residue for one residue.
        for cap in substitution_pairs(rendered) {
            if cap.0.chars().count() != 1 || cap.1.chars().count() != 1 {
                v.push(Violation {
                    rule_id: "B4-substitution-one-to-one",
                    level: "must",
                    message: format!("`{}>{}` is not a 1-to-1 substitution", cap.0, cap.1),
                    citation: "docs/recommendations/DNA/delins.md:73",
                });
            } else if cap.0 == cap.1 {
                // B13 — an unchanged residue is `=`, never X>X.
                v.push(Violation {
                    rule_id: "B13-silent-is-equals",
                    level: "must",
                    message: format!("unchanged residue written as `{}>{}`", cap.0, cap.1),
                    citation: "docs/recommendations/checklist.md:64",
                });
            }
        }

        // S6 — the shortest possible frameshift is `fsTer2`; `fsTer1` means the
        // codon itself is the stop, which is a nonsense substitution. The `1`
        // must be the whole count: `fsTer17` is a perfectly legal frameshift.
        if has_ter_count_one(rendered) {
            v.push(Violation {
                rule_id: "S6-no-fsTer1",
                level: "must",
                message: format!("`fsTer1` in {rendered}; a nonsense change is `Ter`/`*`"),
                citation: "docs/recommendations/protein/frameshift.md:22",
            });
        }

        // B11 — start-loss is never `p.Met1<aa>`.
        if let Some(rest) = protein_body(rendered) {
            let stripped = rest.trim_start_matches('(');
            if let Some(tail) = stripped.strip_prefix("Met1") {
                let next = tail.trim_start_matches(')');
                if !next.is_empty()
                    && !next.starts_with('?')
                    && !next.starts_with('=')
                    && next.chars().next().is_some_and(|c| c.is_ascii_uppercase())
                {
                    v.push(Violation {
                        rule_id: "B11-start-loss-not-met1xxx",
                        level: "must",
                        message: format!(
                            "start-loss written as {rendered}; use p.0/p.0?/p.(Met1?)"
                        ),
                        citation: "docs/recommendations/protein/substitution.md:49",
                    });
                }
            }
            // B10 — `X` never denotes a stop codon. `Xaa` is the spec's own
            // code for an *unknown* amino acid and is entirely legal, so only a
            // bare `X` not followed by `aa` is a violation.
            if rest.replace("Xaa", "").contains('X') {
                v.push(Violation {
                    rule_id: "B10-no-X-for-stop",
                    level: "must",
                    message: format!("`X` used as an amino-acid code in {rendered}"),
                    citation: "docs/recommendations/checklist.md:63",
                });
            }
        }

        // B15 — an allele must not mix reference types.
        let axes: BTreeSet<char> = ms.iter().map(|m| m.axis).collect();
        if axes.len() > 1 {
            v.push(Violation {
                rule_id: "B15-no-mixed-reference-types",
                level: "must",
                message: format!("mixed coordinate types in one allele: {rendered}"),
                citation: "docs/recommendations/DNA/alleles.md:23",
            });
        }

        // S1 — two members of a cis list must not describe the same position.
        // `(;)` is the *unknown-phase* operator and `];[` is trans: neither is a
        // cis list, so neither can violate the one-position-per-cis-list rule.
        if rendered.contains(';') && !rendered.contains("];[") && !rendered.contains("(;)") {
            let mut seen: BTreeMap<&str, usize> = BTreeMap::new();
            for m in &ms {
                if m.loc.is_empty() {
                    continue;
                }
                *seen.entry(m.loc.as_str()).or_default() += 1;
            }
            for (loc, n) in seen {
                if n > 1 {
                    v.push(Violation {
                        rule_id: "S1-no-conflicting-cis-members",
                        level: "must",
                        message: format!("{n} cis members at position `{loc}` in {rendered}"),
                        citation: "docs/recommendations/DNA/delins.md:19",
                    });
                }
            }
        }

        for m in &ms {
            match m.kw.as_str() {
                // B5 — an inversion is at least two residues.
                "inv" if is_plain_position(&m.loc) => v.push(Violation {
                    rule_id: "B5-inversion-min-length",
                    level: "must",
                    message: format!("single-position inversion `{}inv`", m.loc),
                    citation: "docs/recommendations/DNA/inversion.md:16",
                }),
                // S5 / B16 — an insertion is between two *flanking* positions.
                "ins" => {
                    if is_plain_position(&m.loc) {
                        v.push(Violation {
                            rule_id: "S5-insertion-spans-two-positions",
                            level: "must",
                            message: format!("insertion at one position `{}ins`", m.loc),
                            citation: "docs/recommendations/protein/insertion.md:18",
                        });
                    } else if let Some((a, b)) = plain_range(&m.loc) {
                        if b != a + 1 {
                            v.push(Violation {
                                rule_id: "B16-insertion-flanks-adjacent",
                                level: "must",
                                message: format!("insertion flanks {a}_{b} are not adjacent"),
                                citation: "docs/recommendations/protein/insertion.md:17",
                            });
                        }
                    }
                }
                _ => {}
            }

            // B6 — a del/dup range is two *different* positions, 5'→3'.
            // `o.`/`m.` are circular: an end < start is legal there.
            if matches!(m.kw.as_str(), "del" | "dup" | "delins" | "inv") {
                if let Some((a, b)) = plain_range(&m.loc) {
                    if a == b {
                        v.push(Violation {
                            rule_id: "B6-range-two-distinct-positions",
                            level: "must",
                            message: format!("range `{a}_{b}` repeats one position"),
                            citation: "docs/recommendations/DNA/deletion.md:15",
                        });
                    } else if a > b && !matches!(m.axis, 'o' | 'm') {
                        v.push(Violation {
                            rule_id: "B6-range-two-distinct-positions",
                            level: "must",
                            message: format!("range `{a}_{b}` is not ordered 5'->3'"),
                            citation: "docs/recommendations/DNA/deletion.md:17",
                        });
                    }
                }
            }

            // B8 — no size-number forms (`g.123del3`, `g.123dup6`, `c.1_2ins6`).
            if matches!(m.kw.as_str(), "del" | "dup" | "ins")
                && !m.arg.is_empty()
                && m.arg.bytes().all(|b| b.is_ascii_digit())
            {
                v.push(Violation {
                    rule_id: "B8-no-size-number-forms",
                    level: "must",
                    message: format!("size-number form `{}{}`", m.kw, m.arg),
                    citation: "docs/recommendations/checklist.md:49",
                });
            }

            // B9 — amino acids after a stop are not listed.
            if m.axis == 'p' && matches!(m.kw.as_str(), "ins" | "delins") {
                if let Some(i) = m.arg.find("Ter") {
                    let tail = m.arg[i + 3..].trim_end_matches(')');
                    // Digits after `Ter` are a repeat count, not residues.
                    if tail.chars().next().is_some_and(|c| c.is_ascii_alphabetic()) {
                        v.push(Violation {
                            rule_id: "B9-nothing-after-stop",
                            level: "must",
                            message: format!("residues listed after `Ter` in `{}`", m.arg),
                            citation: "docs/recommendations/protein/delins.md:45",
                        });
                    }
                }
            }

            // Advisory: the long forms `delA` / `dupT` / `delNNNins` are legal
            // input but the spec recommends against them in output.
            if matches!(m.kw.as_str(), "del" | "dup")
                && !m.arg.is_empty()
                && m.arg.bytes().all(|b| b.is_ascii_alphabetic())
            {
                v.push(Violation {
                    rule_id: "A1-no-longform-del-dup",
                    level: "should",
                    message: format!("long form `{}{}`", m.kw, m.arg),
                    citation: "docs/recommendations/DNA/deletion.md:30",
                });
            }
        }

        // Advisory: a lone single-member `[X]` is "misleading" — but a
        // single-member bracket is *mandatory* in a trans genotype, so this
        // fires only when there is exactly one bracket group.
        if is_lone_single_member_allele(rendered) {
            v.push(Violation {
                rule_id: "A2-lone-single-member-bracket",
                level: "should",
                message: format!("lone single-member allele bracket in {rendered}"),
                citation: "docs/recommendations/DNA/alleles.md:99",
            });
        }

        v.sort_by(|a, b| a.rule_id.cmp(b.rule_id).then(a.message.cmp(&b.message)));
        v.dedup_by(|a, b| a.rule_id == b.rule_id && a.message == b.message);
        v
    }

    /// True when the rendered string carries a frameshift whose new-stop count
    /// is exactly 1 (`fsTer1` / `fs*1`), i.e. not `fsTer17`.
    fn has_ter_count_one(rendered: &str) -> bool {
        for marker in ["fsTer", "fs*"] {
            let mut from = 0usize;
            while let Some(i) = rendered[from..].find(marker) {
                let at = from + i + marker.len();
                let digits: String = rendered[at..]
                    .chars()
                    .take_while(|c| c.is_ascii_digit())
                    .collect();
                if digits == "1" {
                    return true;
                }
                from = at.max(from + i + 1);
            }
        }
        false
    }

    /// True when the body is a *lone* single-member allele bracket — the shape
    /// `<axis>.[X]` with no second group and no separator. A single-member
    /// bracket is **mandatory** in a trans genotype (`c.[A];[B]`), and `[n]`
    /// after a sequence is a repeat count, so both must be excluded.
    fn is_lone_single_member_allele(rendered: &str) -> bool {
        let body = match rendered.find(':') {
            Some(i) => &rendered[i + 1..],
            None => rendered,
        };
        let Some(dot) = body.find('.') else {
            return false;
        };
        let rest = &body[dot + 1..];
        rest.starts_with('[')
            && rest.ends_with(']')
            && rest.matches('[').count() == 1
            && !rest.contains(';')
            && !rest.contains(',')
    }

    /// The protein body of a rendered variant (`…:p.<body>`), if any.
    fn protein_body(rendered: &str) -> Option<&str> {
        let body = match rendered.find(':') {
            Some(i) => &rendered[i + 1..],
            None => rendered,
        };
        body.strip_prefix("p.")
    }

    /// Every `<ref>><alt>` residue pair in a rendered DNA/RNA variant.
    fn substitution_pairs(rendered: &str) -> Vec<(String, String)> {
        let mut out = Vec::new();
        let bytes: Vec<char> = rendered.chars().collect();
        for (i, c) in bytes.iter().enumerate() {
            if *c != '>' {
                continue;
            }
            let mut a = i;
            while a > 0 && is_base(bytes[a - 1]) {
                a -= 1;
            }
            let mut b = i + 1;
            while b < bytes.len() && is_base(bytes[b]) {
                b += 1;
            }
            let lhs: String = bytes[a..i].iter().collect();
            let rhs: String = bytes[i + 1..b].iter().collect();
            if !lhs.is_empty() && !rhs.is_empty() {
                out.push((lhs, rhs));
            }
        }
        out
    }

    fn is_base(c: char) -> bool {
        matches!(
            c,
            'A' | 'C' | 'G' | 'T' | 'U' | 'a' | 'c' | 'g' | 't' | 'u' | 'N' | 'n'
        )
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn ids(s: &str) -> Vec<&'static str> {
            check(s).into_iter().map(|v| v.rule_id).collect()
        }

        #[test]
        fn flags_the_documented_must_violations() {
            assert!(ids("NC_000023.11:g.4GC>TG").contains(&"B4-substitution-one-to-one"));
            assert!(ids("NC_000023.11:g.234inv").contains(&"B5-inversion-min-length"));
            assert!(ids("NG_012232.1:g.123del6").contains(&"B8-no-size-number-forms"));
            assert!(ids("NP_003997.1:p.Leu54Leu").is_empty()); // protein X>X is not a `>` form
            assert!(ids("NC_000023.11:g.123A>A").contains(&"B13-silent-is-equals"));
            assert!(ids("NM_004006.2:c.52insT").contains(&"S5-insertion-spans-two-positions"));
            assert!(ids("NM_004006.2:c.51_54insT").contains(&"B16-insertion-flanks-adjacent"));
            assert!(ids("NP_003997.1:p.Met1Leu").contains(&"B11-start-loss-not-met1xxx"));
        }

        /// The acceptance gate depends on these staying quiet: every one of
        /// these is a form the repo already blesses.
        #[test]
        fn does_not_flag_legal_blessed_output() {
            for legal in [
                "NC_000023.11:g.33038255C>A",
                "NC_000023.11:g.1234_2345del",
                "NC_000023.11:g.1234_1235insACGT",
                "NC_000023.11:g.1234_2345inv",
                "NM_004006.2:c.20_23dup",
                "NP_003997.1:p.Trp24Cys",
                "NP_003997.1:p.Leu54=",
                "LRG_199t1:c.[2376G>C];[3103del]",
                "NM_004006.2:c.(4071_5145)del",
            ] {
                let musts: Vec<&str> = check(legal)
                    .into_iter()
                    .filter(|v| v.level == "must")
                    .map(|v| v.rule_id)
                    .collect();
                assert!(musts.is_empty(), "{legal} wrongly flagged: {musts:?}");
            }
        }

        /// A single-member bracket is *mandatory* in a trans genotype — it must
        /// never be a MUST violation (the §5 rule-2 class of error).
        #[test]
        fn trans_genotype_single_member_brackets_are_legal() {
            let musts: Vec<&str> = check("NM_004006.2:c.[2376G>C];[3103del]")
                .into_iter()
                .filter(|v| v.level == "must")
                .map(|v| v.rule_id)
                .collect();
            assert!(musts.is_empty(), "{musts:?}");
        }
    }
}

// ---------- builder ----------

mod builder {
    use super::census::Census;
    use super::*;
    use ferro_hgvs::error_handling::ErrorConfig;
    use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
    use ferro_hgvs::reference::mock::MockProvider;
    use ferro_hgvs::{parse_hgvs, Normalizer};

    pub struct Build {
        pub rows: Vec<Row>,
        pub census: Census,
    }

    /// Resolve the string ferro is actually run against.
    fn target_for(input: &str) -> String {
        prefix::default_prefixed(input).unwrap_or_else(|| input.to_string())
    }

    fn cite(path: &str, line: Option<usize>, sha: &str) -> String {
        match line {
            Some(l) => format!("{path}:{l}@{sha}"),
            None => format!("{path}@{sha}"),
        }
    }

    pub fn build(
        spec_dir: &Path,
        candidates: &[sources::Candidate],
        prose: &[negatives::ProseNegative],
        ov: &overrides::Overrides,
        existing: &dedup::ExistingCoverage,
        sha: &str,
    ) -> anyhow::Result<Build> {
        let normalizer = Normalizer::new(MockProvider::new());
        let mut rows: Vec<Row> = Vec::new();
        let mut cens = Census::default();

        let raw_negatives = raw_spec_negatives(spec_dir)?;
        negative_reject(
            &raw_negatives,
            prose,
            existing,
            sha,
            &normalizer,
            &mut rows,
            &mut cens,
        );
        negative_repair(
            &raw_negatives,
            prose,
            ov,
            sha,
            &normalizer,
            &mut rows,
            &mut cens,
        );
        error_mode(candidates, existing, sha, &mut rows, &mut cens);
        grammar_form(candidates, existing, sha, &mut rows, &mut cens);
        output_invariant(existing, sha, &mut rows, &mut cens);

        // Apply per-row overrides, then check none is stale.
        let by_id: BTreeSet<String> = rows.iter().map(|r| r.id.clone()).collect();
        let stale: Vec<&str> = ov
            .by_id
            .keys()
            .filter(|k| !by_id.contains(k.as_str()))
            .map(String::as_str)
            .collect();
        if !stale.is_empty() {
            anyhow::bail!("overrides.by_id references unknown row ids: {stale:?}");
        }
        for r in &mut rows {
            if let Some(o) = ov.by_id.get(&r.id) {
                if let Some(s) = &o.status {
                    r.status = s.clone();
                }
                if let Some(n) = &o.note {
                    r.note = Some(n.clone());
                }
            }
        }

        rows.sort_by(|a, b| a.dimension.cmp(&b.dimension).then(a.id.cmp(&b.id)));
        let mut seen = BTreeSet::new();
        for r in &rows {
            anyhow::ensure!(seen.insert(r.id.clone()), "duplicate row id: {}", r.id);
        }
        Ok(Build { rows, census: cens })
    }

    /// Every HGVS-shaped string the spec marks with `<code class="invalid">`,
    /// with markdown escapes undone.
    ///
    /// This deliberately does **not** go through `sources::discover`: that path
    /// feeds the normalization fixture, whose `canonicalize` drops any span
    /// containing a backslash. The spec escapes underscores inside HTML code
    /// spans (`c.(?\_-244)\_(31+1\_32-1)del`), so five genuine negatives are
    /// silently lost there. Unescaping here recovers them without perturbing
    /// the normalization fixture, which must stay byte-stable.
    fn raw_spec_negatives(spec_dir: &Path) -> anyhow::Result<Vec<(String, String, usize)>> {
        let mut out = Vec::new();
        let recs = spec_dir.join("docs/recommendations");
        for path in sources::walkdir_md_pub(&recs) {
            let text = std::fs::read_to_string(&path)?;
            let rel = path
                .strip_prefix(spec_dir)
                .unwrap_or(&path)
                .to_string_lossy()
                .replace('\\', "/");
            for (span, offset) in sources::extract_class_invalid_codespans(&text) {
                let unescaped = unescape_markdown(&span);
                if let Some(input) = sources::canonicalize_pub(&unescaped) {
                    out.push((input, rel.clone(), sources::line_of(&text, offset)));
                } else if let Some(input) = canonicalize_bracket_allele(&unescaped) {
                    out.push((input, rel.clone(), sources::line_of(&text, offset)));
                }
            }
        }
        Ok(out)
    }

    /// Undo markdown backslash escapes (`\_` -> `_`, `\*` -> `*`, ...).
    fn unescape_markdown(s: &str) -> String {
        let mut out = String::with_capacity(s.len());
        let mut chars = s.chars();
        while let Some(c) = chars.next() {
            if c == '\\' {
                if let Some(next) = chars.next() {
                    out.push(next);
                    continue;
                }
            }
            out.push(c);
        }
        out
    }

    /// Accept a bracketed allele list whose members carry the coordinate prefix
    /// (`[c.76A>C;c.83G>C]`). `sources::canonicalize` requires the prefix at the
    /// *head* of the string, so these spec negatives are dropped there.
    fn canonicalize_bracket_allele(s: &str) -> Option<String> {
        let s = s.trim();
        if !(s.starts_with('[') && s.ends_with(']')) {
            return None;
        }
        if s.contains(char::is_whitespace) || s.contains('<') || s.contains('"') {
            return None;
        }
        const PREFIXES: &[&str] = &["c.", "g.", "m.", "n.", "o.", "p.", "r."];
        PREFIXES
            .iter()
            .any(|p| s.contains(p))
            .then(|| s.to_string())
    }

    /// Dimension 1 — every spec-stated negative, as a *reject* assertion.
    ///
    /// Gross = every distinct HGVS-shaped string the spec marks invalid, from
    /// both the `<code class="invalid">` markers and the prose negatives.
    /// Duplicates = those the normalization fixture already retains (it pins
    /// their rejection under `normalize`). Net = the ones the existing
    /// harvester drops — escaped-underscore forms and bare bracket allele lists.
    fn negative_reject(
        raw_negatives: &[(String, String, usize)],
        prose: &[negatives::ProseNegative],
        existing: &dedup::ExistingCoverage,
        sha: &str,
        normalizer: &Normalizer<MockProvider>,
        rows: &mut Vec<Row>,
        cens: &mut Census,
    ) {
        let mut all: BTreeMap<String, (String, Option<usize>, &'static str)> = BTreeMap::new();
        for (input, path, line) in raw_negatives {
            all.entry(input.clone())
                .or_insert((path.clone(), Some(*line), "must"));
        }
        for p in prose {
            let level = match p.normative_level {
                negatives::NormativeLevel::Must => "must",
                negatives::NormativeLevel::Should => "should",
            };
            for s in &p.invalid_spans {
                all.entry(s.clone())
                    .or_insert((p.source_path.clone(), Some(p.source_line), level));
            }
        }

        let (mut dup, mut net) = (0usize, 0usize);
        let gross = all.len();
        for (input, (path, line, level)) in all {
            let target = target_for(&input);
            if existing.covers("reject", "default", &target) {
                dup += 1;
                continue;
            }
            let observed = observe(normalizer, &target);
            let status = if observed.starts_with("parse error") {
                "correctly-rejected"
            } else {
                "false-acceptance"
            };
            net += 1;
            rows.push(Row {
                id: format!("negative-reject/{input}"),
                dimension: "negative-reject".to_string(),
                operation: "reject".to_string(),
                error_mode: "default".to_string(),
                input,
                target,
                expectation: "spec-mandated".to_string(),
                normative_level: level.to_string(),
                expected: None,
                observed,
                status: status.to_string(),
                spec_citation: cite(&path, line, sha),
                dedup_note: "spec-marked negative the normalization fixture's harvester drops \
                             (escaped underscores / bare bracket allele list), so its rejection \
                             is asserted nowhere else"
                    .to_string(),
                note: None,
            });
        }
        cens.tally("negative-reject", gross, dup, 0);
        debug_assert_eq!(net + dup, gross);
    }

    /// Dimension 2 — the same negatives as *repair-to-canonical* assertions.
    ///
    /// This is the genuinely new half of the negative corpus: the normalization
    /// fixture records only that ferro rejects a bad string, never that
    /// `normalize(bad)` yields the good string the spec names on the same line.
    /// Repair targets are hand-curated (`overrides.repairs`) because automatic
    /// extraction is unreliable — the spec routinely cites a bare coordinate on
    /// the same line as the corrected form.
    ///
    /// Not-applicable = negatives the spec forbids without naming a replacement
    /// (`repaired: null`), and repairs that need real reference bases.
    fn negative_repair(
        raw_negatives: &[(String, String, usize)],
        prose: &[negatives::ProseNegative],
        ov: &overrides::Overrides,
        sha: &str,
        normalizer: &Normalizer<MockProvider>,
        rows: &mut Vec<Row>,
        cens: &mut Census,
    ) {
        // Gross = every distinct spec-stated negative (a repair assertion is
        // conceivable for each); n/a = those with no curated target.
        let mut all: BTreeSet<String> = raw_negatives
            .iter()
            .map(|(input, _, _)| input.clone())
            .collect();
        for p in prose {
            all.extend(p.invalid_spans.iter().cloned());
        }
        let gross = all.len();
        let mut na = 0usize;

        for input in &all {
            let Some(rep) = ov.repairs.get(input) else {
                na += 1;
                continue;
            };
            let Some(repaired) = rep.repaired.clone() else {
                na += 1;
                continue;
            };
            let target = target_for(input);
            let observed = observe(normalizer, &target);
            // Distinguish the two very different ways a repair can not happen.
            // Rejecting the bad string outright satisfies the spec — the spec
            // forbids the form, it never obliges an implementation to repair it.
            // Accepting it and rendering something other than the canonical form
            // the spec names is the genuine violation.
            let status = if rep.requires_reference {
                "requires-reference"
            } else if observed == repaired {
                "repaired"
            } else if observed.starts_with("parse error") || observed.starts_with("normalize error")
            {
                "rejected-not-repaired"
            } else {
                "repair-diverges"
            };
            rows.push(Row {
                id: format!("negative-repair/{input}"),
                dimension: "negative-repair".to_string(),
                operation: "normalize".to_string(),
                error_mode: "default".to_string(),
                input: input.clone(),
                target,
                expectation: "spec-mandated".to_string(),
                normative_level: rep.normative_level.clone(),
                expected: Some(repaired),
                observed,
                status: status.to_string(),
                spec_citation: format!("{}@{sha}", rep.citation),
                dedup_note: "the normalization fixture records only that ferro rejects this \
                             string; the canonical form the spec names on the same line is \
                             asserted nowhere"
                    .to_string(),
                note: rep.note.clone(),
            });
        }
        cens.tally("negative-repair", gross, 0, na);
    }

    /// Dimension 3 — `parse` under strict / lenient / silent.
    ///
    /// `ErrorConfig` acts on a *preprocessing* pass over the raw string
    /// (`src/error_handling/preprocessor.rs`): it repairs non-standard input
    /// such as en-dashes, lowercase amino acids and stray whitespace. For a
    /// spec-clean string all three modes are therefore bit-identical, and a
    /// blanket 3x cross-product would be almost entirely waste. This dimension
    /// emits a row **only** where the three modes actually disagree — which is
    /// both the dedup rule and the entire point.
    ///
    /// Not-applicable = inputs where all three modes agree.
    fn error_mode(
        candidates: &[sources::Candidate],
        existing: &dedup::ExistingCoverage,
        sha: &str,
        rows: &mut Vec<Row>,
        cens: &mut Census,
    ) {
        let mut by_input: BTreeMap<String, (String, Option<usize>)> = BTreeMap::new();
        for c in candidates {
            by_input
                .entry(c.input.clone())
                .or_insert((c.source_path.clone(), c.source_line));
        }
        let gross = by_input.len() * 3;
        let (mut na, mut dup) = (0usize, 0usize);

        for (input, (path, line)) in &by_input {
            let target = target_for(input);
            let modes = [
                ("strict", ErrorConfig::strict()),
                ("lenient", ErrorConfig::lenient()),
                ("silent", ErrorConfig::silent()),
            ];
            let outcomes: Vec<(&str, String)> = modes
                .into_iter()
                .map(|(name, cfg)| (name, observe_parse_with_config(&target, cfg)))
                .collect();
            let distinct: BTreeSet<&str> = outcomes.iter().map(|(_, o)| o.as_str()).collect();
            if distinct.len() == 1 {
                // All three modes agree — the single existing `normalize` row
                // already pins this behaviour.
                na += 3;
                continue;
            }
            for (mode, observed) in outcomes {
                if existing.covers("parse", mode, &target) {
                    dup += 1;
                    continue;
                }
                rows.push(Row {
                    id: format!("error-mode/{mode}/{input}"),
                    dimension: "error-mode".to_string(),
                    operation: "parse".to_string(),
                    error_mode: mode.to_string(),
                    input: input.clone(),
                    target: target.clone(),
                    // The spec says nothing about ferro's error modes: this is
                    // ferro policy, pinned, never a spec expectation.
                    expectation: "pinned-baseline".to_string(),
                    normative_level: "n/a".to_string(),
                    expected: None,
                    observed,
                    status: "mode-divergence-pinned".to_string(),
                    spec_citation: cite(path, *line, sha),
                    dedup_note: "the three error modes disagree on this input; the normalization \
                                 fixture pins only the default mode"
                        .to_string(),
                    note: None,
                });
            }
        }
        cens.tally("error-mode", gross, dup, na);
    }

    /// Dimension 4 — `syntax.yaml` form x legal coordinate system.
    ///
    /// The grammar file declares 43 forms across 3 axes (`aa`/`dna`/`rna`) with
    /// 56 canonical examples. The normalization fixture feeds those examples
    /// through `normalize`, but nothing asserts the *form-level* contract: that
    /// the example parses into the coordinate axis its grammar entry declares.
    ///
    /// The axis-substitution cells (a `dna` form re-expressed in each of the
    /// other DNA coordinate types) are enumerated too, since `coordinate_type`
    /// is a grammar element of every DNA form. Those are **pinned baselines**,
    /// not spec expectations: the spec gives no example for them.
    ///
    /// Not-applicable = cells the axis does not admit (a protein form has no
    /// genomic coordinate type; an RNA form has no protein coordinate type).
    fn grammar_form(
        candidates: &[sources::Candidate],
        existing: &dedup::ExistingCoverage,
        sha: &str,
        rows: &mut Vec<Row>,
        cens: &mut Census,
    ) {
        // The syntax.yaml examples, recovered from the harvested candidates.
        // Deduplicated by string: syntax.yaml reuses an example across forms.
        let examples: Vec<&sources::Candidate> = {
            let mut seen = BTreeSet::new();
            candidates
                .iter()
                .filter(|c| c.source_kind == sources::SourceKind::SyntaxYaml)
                .filter(|c| seen.insert(c.input.clone()))
                .collect()
        };

        // Which coordinate types each grammar axis admits. `aa` -> `p` only;
        // the DNA coordinate types are interchangeable grammar-wise; `rna` is
        // `r` only.
        const DNA_AXES: &[char] = &['g', 'c', 'n', 'm', 'o'];
        let gross = examples.len() * 7;
        let (mut na, mut dup) = (0usize, 0usize);

        for c in &examples {
            let target = target_for(&c.input);
            let axis = declared_axis(&c.input);
            let legal: Vec<char> = match axis {
                'p' => vec!['p'],
                'r' => vec!['r'],
                _ => DNA_AXES.to_vec(),
            };
            na += 7 - legal.len();

            for a in legal {
                let derived = if a == axis {
                    target.clone()
                } else {
                    match retarget_axis(&c.input, a) {
                        Some(s) => s,
                        None => {
                            na += 1;
                            continue;
                        }
                    }
                };
                // The declared-axis cell is exactly the existing fixture row.
                if a == axis && existing.covers("normalize", "default", &derived) {
                    // Still emit: the assertion here is the *axis* contract, not
                    // the normalized rendering. But only when the axis check
                    // adds information, i.e. the string parses.
                    if parse_hgvs(&derived).is_err() {
                        dup += 1;
                        continue;
                    }
                }
                let observed = match parse_hgvs(&derived) {
                    Ok(v) => match v.coordinate_axis() {
                        Some(ax) => format!("axis={}", ax.code()),
                        None => "axis=none".to_string(),
                    },
                    Err(e) => format!("parse error: {e}"),
                };
                let expected = format!("axis={a}");
                let spec_stated = a == axis;
                let status = if observed == expected {
                    "form-axis-ok"
                } else if spec_stated {
                    "form-axis-diverges"
                } else {
                    "form-axis-pinned"
                };
                rows.push(Row {
                    id: format!("grammar-form/{a}/{}", c.input),
                    dimension: "grammar-form".to_string(),
                    operation: "parse".to_string(),
                    error_mode: "default".to_string(),
                    input: c.input.clone(),
                    target: derived,
                    expectation: if spec_stated {
                        "spec-mandated".to_string()
                    } else {
                        "pinned-baseline".to_string()
                    },
                    normative_level: if spec_stated { "must" } else { "n/a" }.to_string(),
                    expected: Some(expected),
                    observed,
                    status: status.to_string(),
                    spec_citation: format!("docs/syntax.yaml@{sha}"),
                    dedup_note: "the normalization fixture pins this example's normalized \
                                 rendering; nothing asserts which coordinate axis it parses into"
                        .to_string(),
                    note: None,
                });
            }
        }
        cens.tally("grammar-form", gross, dup, na);
    }

    /// The coordinate-system letter an example declares.
    fn declared_axis(input: &str) -> char {
        for a in ['c', 'g', 'm', 'n', 'o', 'p', 'r'] {
            if input.starts_with(&format!("{a}.")) || input.contains(&format!(":{a}.")) {
                return a;
            }
        }
        '?'
    }

    /// Re-express an example under a different DNA coordinate type, swapping in
    /// the axis-appropriate default accession. Returns `None` when the shape is
    /// not a plain `<accession>:<axis>.<body>` (gene-selector forms, nested
    /// references), since rewriting those would invent syntax the spec never
    /// wrote.
    fn retarget_axis(input: &str, axis: char) -> Option<String> {
        let (_, body) = input.split_once(':')?;
        let rest = body.strip_prefix(&format!("{}.", declared_axis(input)))?;
        if rest.contains(':') {
            return None;
        }
        let acc = prefix::DEFAULTS.iter().find(|(c, _)| *c == axis)?.1;
        Some(format!("{acc}:{axis}.{rest}"))
    }

    /// Dimension 5 — the MUST-level output invariants over every emitted output.
    ///
    /// This is the only dimension that asserts something about what ferro
    /// *produces* rather than what it accepts. It runs the invariant catalog
    /// over every `preserved` output in the normalization fixture — the set the
    /// repo already blesses, which doubles as the checker's zero-false-positive
    /// acceptance gate.
    ///
    /// One row per (output, rule) pair that reports a violation, plus one
    /// summary row per output that is clean. Not-applicable = rules the checker
    /// cannot decide without reference bases; those are absent from the catalog
    /// by design, so they never appear here.
    fn output_invariant(
        existing: &dedup::ExistingCoverage,
        sha: &str,
        rows: &mut Vec<Row>,
        cens: &mut Census,
    ) {
        let outputs = existing.emitted_outputs();
        let gross = outputs.len();
        let mut emitted = 0usize;
        for (target, rendered, pinned_status) in &outputs {
            let vs = invariants::check(rendered);
            if vs.is_empty() {
                continue;
            }
            for v in vs {
                emitted += 1;
                rows.push(Row {
                    id: format!("output-invariant/{}/{}", v.rule_id, target),
                    dimension: "output-invariant".to_string(),
                    operation: "invariant-check".to_string(),
                    error_mode: "default".to_string(),
                    input: target.clone(),
                    target: rendered.clone(),
                    expectation: "spec-mandated".to_string(),
                    normative_level: v.level.to_string(),
                    expected: Some("no violation".to_string()),
                    observed: v.message.clone(),
                    status: if v.level == "must" {
                        "invariant-violation-must".to_string()
                    } else {
                        "invariant-violation-should".to_string()
                    },
                    spec_citation: format!("{}@{sha}", v.citation),
                    dedup_note: "no existing test asserts any property of the string ferro \
                                 emits; the whole suite is input-side"
                        .to_string(),
                    note: Some(format!(
                        "hgvs_spec_normalization.json row status: {pinned_status}"
                    )),
                });
            }
        }
        // The unit of assertion is (emitted output x catalog): the driver
        // re-runs the whole catalog over every preserved output, so all `gross`
        // of them are asserted and none is a restatement — the existing suite is
        // entirely input-side. Clean outputs are deliberately not materialised
        // as rows (700 rows of "no violation" would be noise); only violations
        // are, and `emitted` counts those.
        let _ = emitted;
        cens.tally("output-invariant", gross, 0, 0);
    }

    fn observe(normalizer: &Normalizer<MockProvider>, target: &str) -> String {
        match parse_hgvs(target) {
            Err(e) => format!("parse error: {e}"),
            Ok(v) => match normalizer.normalize_with_diagnostics(&v) {
                Err(e) => format!("normalize error: {e}"),
                Ok(n) => format!("{}", n.result),
            },
        }
    }

    fn observe_parse_with_config(target: &str, cfg: ErrorConfig) -> String {
        match parse_hgvs_with_config(target, cfg) {
            Err(e) => format!("parse error: {e}"),
            Ok(r) => {
                let mut codes: Vec<String> = r
                    .warnings
                    .iter()
                    .map(|w| format!("{:?}", w.error_type))
                    .collect();
                codes.sort();
                codes.dedup();
                if codes.is_empty() {
                    format!("{}", r.result)
                } else {
                    format!("{} warnings={}", r.result, codes.join(","))
                }
            }
        }
    }
}

// ---------- render ----------

mod render {
    use super::*;

    #[derive(Serialize)]
    struct Document<'a> {
        description: &'a str,
        spec: SpecBlock<'a>,
        dedup_rule: &'a str,
        census: &'a census::Census,
        summary: Summary,
        rows: &'a [Row],
    }

    #[derive(Serialize)]
    struct SpecBlock<'a> {
        source: &'a str,
        commit_sha: &'a str,
    }

    #[derive(Serialize, Default)]
    struct Summary {
        total: usize,
        by_dimension: BTreeMap<String, usize>,
        by_status: BTreeMap<String, usize>,
        by_expectation: BTreeMap<String, usize>,
        by_normative_level: BTreeMap<String, usize>,
    }

    pub fn render(build: &builder::Build, sha: &str) -> anyhow::Result<String> {
        let mut summary = Summary {
            total: build.rows.len(),
            ..Default::default()
        };
        for r in &build.rows {
            *summary.by_dimension.entry(r.dimension.clone()).or_default() += 1;
            *summary.by_status.entry(r.status.clone()).or_default() += 1;
            *summary
                .by_expectation
                .entry(r.expectation.clone())
                .or_default() += 1;
            *summary
                .by_normative_level
                .entry(r.normative_level.clone())
                .or_default() += 1;
        }
        let doc = Document {
            description:
                "Exhaustive, non-redundant HGVS spec test enumeration. Generated artifact - \
                 gitignored, regenerate with `cargo run --features dev --example \
                 generate_spec_enumeration`. Complements (never restates) \
                 hgvs_spec_normalization.json.",
            spec: SpecBlock {
                source: "https://github.com/HGVSnomenclature/hgvs-nomenclature",
                commit_sha: sha,
            },
            dedup_rule:
                "A row is a duplicate iff (operation, error_mode, target) is already asserted by \
                 an existing test. hgvs_spec_normalization.json + idempotency_tests saturate \
                 (normalize|display-roundtrip|idempotency|reject, default, <every harvested \
                 target>); those tuples are never re-emitted here.",
            census: &build.census,
            summary,
            rows: &build.rows,
        };
        let mut out = serde_json::to_string_pretty(&doc)?;
        out.push('\n');
        Ok(out)
    }
}
