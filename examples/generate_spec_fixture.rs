//! Generator for the HGVS v21.0 spec normalization fixture.
//!
//! Walks the vendored `assets/hgvs-nomenclature/` markdown + syntax.yaml,
//! extracts every HGVS string, runs each through ferro's `Normalizer`, merges
//! hand overrides, and writes `tests/fixtures/grammar/hgvs_spec_normalization.json`.
//!
//! Two modes, answering two different questions:
//!   - default:  regenerate the fixture in-place. This is the run that
//!     **validates the committed inputs** — it harvests the spec checkout and
//!     resolves every curated override against it — which is why CI and the
//!     pre-push hook use it rather than `--check`.
//!   - --check:  is the *local* artifact current? Exits non-zero only when an
//!     existing fixture differs from a fresh render. An absent fixture is a cold
//!     cache, not drift (the file is gitignored, so there is no committed
//!     baseline), and is simply generated.
//!
//! Run: `cargo run --features dev --bin generate_spec_fixture`
//! Check: `cargo run --features dev --bin generate_spec_fixture -- --check`
//!
//! See: docs/superpowers/specs/2026-05-03-issue-84-hgvs-v21-normalization-fixture-design.md

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use serde::{Deserialize, Deserializer, Serialize};

#[derive(Parser, Debug)]
#[command(about = "Generate the HGVS v21.0 spec normalization fixture")]
struct Cli {
    /// Path to the vendored spec checkout (must be at tag 21.0.0).
    #[arg(long, default_value = "assets/hgvs-nomenclature")]
    spec_dir: PathBuf,

    /// Path to the hand-curated overrides file.
    #[arg(
        long,
        default_value = "tests/fixtures/grammar/hgvs_spec_normalization_overrides.json"
    )]
    overrides: PathBuf,

    /// Path to the output fixture file.
    #[arg(
        long,
        default_value = "tests/fixtures/grammar/hgvs_spec_normalization.json"
    )]
    output: PathBuf,

    /// Is the local fixture current? Exits 1 when an existing fixture differs
    /// from a fresh render, leaving it in place for inspection. An absent
    /// fixture is generated, not treated as a failure: the file is gitignored,
    /// so "never built" is a cold cache rather than drift.
    #[arg(long)]
    check: bool,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    // Every row in this fixture is a `Normalizer` output, so a `FERRO_PARTITION`
    // value naming no arm would record the shipped rule under a candidate's
    // name. The library falls safe to `live` in a release build rather than
    // aborting, and leaves the refusal to its entry point.
    if let Some(message) = ferro_hgvs::normalize::partition_switch_startup_error() {
        eprintln!("error: {message}");
        return ExitCode::FAILURE;
    }
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
    let overrides = overrides::load(&cli.overrides)?;
    let rows = runner::build_rows(&candidates, &overrides)?;
    let decisions = decisions::resolve(&rows, &overrides, &cli.spec_dir)?;
    let rendered = render::render(&rows, &decisions, &cli.spec_dir)?;

    if cli.check {
        check_up_to_date(&cli.output, &rendered, "fixture", "generate_spec_fixture")?;
    } else {
        std::fs::write(&cli.output, rendered)?;
    }
    Ok(())
}

#[path = "common/spec_harvest.rs"]
mod spec_harvest;

#[path = "common/generated_artifact.rs"]
mod generated_artifact;

use generated_artifact::check_up_to_date;
use spec_harvest::{classify, prefix, sources};

// ---------- overrides ----------

mod overrides {
    use super::*;

    #[derive(Debug, Default, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct Overrides {
        #[serde(default)]
        pub by_input: BTreeMap<String, OverrideEntry>,
        /// Curated declarations that two or more harvested rows denote the
        /// **same variant**. See [`EquivalenceClass`].
        #[serde(default)]
        pub equivalence_classes: Vec<EquivalenceClass>,
        /// Curated records of how the project read the spec where its own
        /// clauses conflict. See [`Ruling`].
        #[serde(default)]
        pub rulings: Vec<Ruling>,
    }

    /// A pointer into the vendored spec checkout, carrying the text it points at.
    ///
    /// `clause` is `<path-relative-to-spec-dir>:<line>` or `…:<start>-<end>` —
    /// the same path form the harvester records in `source_paths`, so a citation
    /// and a row's provenance are directly comparable.
    ///
    /// `quote` is what makes the citation self-verifying: the generator asserts
    /// the text appears within the cited lines, so bumping the spec submodule in
    /// a way that moves a clause fails the build instead of silently leaving a
    /// citation pointing at unrelated prose. A bare line number would still
    /// "resolve" against any file long enough to have that line.
    #[derive(Debug, Clone, Deserialize, Serialize)]
    #[serde(deny_unknown_fields)]
    pub struct Citation {
        pub clause: String,
        pub quote: String,
    }

    /// Two or more harvested inputs that denote **one** variant.
    ///
    /// The classes cannot be derived: deciding that two descriptions denote one
    /// variant means applying both to a reference sequence, which many rows
    /// cannot do (the `needs-reference` bucket). So they are curated from the
    /// spec's own equivalence statements — "alternatively", "can also be
    /// described as", "giving an alternative description like" — exactly as the
    /// `by_input` overrides are curated from its validity statements.
    #[derive(Debug, Clone, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct EquivalenceClass {
        /// Stable identifier. This is what the test consumer's pinned
        /// known-non-confluent list keys on, so it must not be renamed casually.
        pub id: String,
        /// Row inputs, verbatim as harvested. Every one must exist as a row.
        pub members: Vec<String>,
        /// Accession all members are evaluated under, replacing whatever
        /// accession each member's row would otherwise carry.
        ///
        /// Needed because the spec routinely writes one member of a pair as a
        /// bare fragment (`c.[235A>T;237G>T]`), which the harvester prefixes
        /// with the per-coordinate-system default rather than with the
        /// accession of the worked example it sits under — and, twice, writes
        /// the two members under *different* accession versions
        /// (`NM_004006.2` / `NM_004006.1`, DNA/insertion.md:74). Comparing
        /// members across different references would report a difference of
        /// reference sequence as a difference of representation.
        ///
        /// Deliberately declared on the class rather than fixed by editing each
        /// row's `input_prefixed`: that would move rows' `current` values and
        /// the census pinned against them, for no gain outside this comparison.
        #[serde(default)]
        pub reference: Option<String>,
        /// Where the spec says these are one variant. Must be non-empty.
        pub citations: Vec<Citation>,
        /// Why this is a class, in prose.
        pub note: String,
    }

    /// Whether the project has taken a position on a contested reading.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
    #[serde(rename_all = "kebab-case")]
    pub enum RulingStatus {
        /// The project has ruled. `governing` names the clause held to govern.
        Decided,
        /// The project has **not** ruled. `governing` must be absent — an
        /// undecided record states the conflict and stops there, so nobody can
        /// read a ruling out of a record that never made one.
        Undecided,
    }

    /// The `README.md` rule a clause-free house choice is made under.
    ///
    /// Only two rules can carry one, and naming which is the whole point of the
    /// field: rule 5's *silent* limb says in terms that ferro "decides under
    /// rule 6 and violates nothing", and rule 6 is where the maintainers choose
    /// among conformant forms. Neither is a claim about what the recommendations
    /// require. A free-text field here would let a record write "the spec" into
    /// the slot reserved for "our own judgement", which is the exact
    /// substitution [`HouseChoice`] exists to make unrepresentable.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
    #[serde(rename_all = "kebab-case")]
    pub enum HouseRule {
        /// `README.md` rule 5, **silent** limb — the recommendations are merely
        /// incomplete here, so ferro decides under rule 6 and violates nothing.
        /// Filing upstream is a feature request, not a bug report.
        RuleFiveSilent,
        /// `README.md` rule 6 — among multiple conformant forms, the maintainers
        /// choose. Use this when every competing form is already conformant and
        /// no clause ranks them.
        RuleSix,
    }

    /// A ruling that is **ours**, not the spec's.
    ///
    /// # Why the ledger needed a shape it did not have
    ///
    /// Before this field the ledger could record two things: a settled reading
    /// of a clause (`decided`, with a `governing` clause) and an unsettled
    /// conflict (`undecided`, with two clauses and no verdict). It could not
    /// record the third and commonest state — *the recommendations are silent
    /// here and we picked X* — because [`validate_ruling_shape`] requires a
    /// `decided` record to name a clause it holds to govern.
    ///
    /// The consequence was measured rather than predicted: the rule that a
    /// reference base is unchanged iff **every** minimum-edit-distance alignment
    /// matches it had no home here and instead lived as prose in six separate
    /// places in `src/` and `tests/`, which then drifted apart. An
    /// HGVS-committee-perspective review of the same rule called it "an
    /// implementer's choice the recommendations neither require nor forbid" and
    /// objected specifically to its being presented as *compliance*.
    ///
    /// # What it must state, and why each part
    ///
    /// * [`Self::under`] — which `README.md` rule the choice is made under. A
    ///   house choice that cannot name its own authority is not a house choice,
    ///   it is an unrecorded preference.
    /// * [`Self::considered_and_rejected`] — what was weighed and put aside.
    ///   The best quality of the existing records is that they record
    ///   *refutations*: `MIN_SEPARATION_NO_FRAME`, `apply_coding_codon_exception`
    ///   and half the ledger's rationales exist because a plausible belief was
    ///   killed by measurement and the belief would otherwise recur. A house
    ///   choice with no rejected alternative is a conclusion with its reasoning
    ///   deleted, and the next reader re-derives the alternative from scratch.
    ///
    /// # What it may NOT state
    ///
    /// A `governing` clause, or a `deviates_from` clause. Both are refused by
    /// [`validate_ruling_shape`], symmetrically with the way an `undecided`
    /// record is refused for naming either. The asymmetry to notice is that the
    /// `undecided` refusal stops a record claiming *more certainty* than the
    /// project has, and this one stops a record claiming *more authority* than
    /// the spec gave it — overclaiming spec authority being the disease this
    /// field was added to treat. Making it unrepresentable is worth more than
    /// discouraging it in prose, because the prose is what drifted.
    ///
    /// A house choice still cites at least one clause, like every other record.
    /// That is deliberate and is **not** a contradiction of "clause-free": the
    /// honest way to assert that the recommendations are silent on a point is to
    /// name the clauses you read and say why none of them reaches, which is
    /// exactly what [`Self::considered_and_rejected`] is for. A record allowed
    /// to assert silence while naming nothing it examined would be an
    /// unfalsifiable claim about the spec.
    #[derive(Debug, Clone, Deserialize, Serialize)]
    #[serde(deny_unknown_fields)]
    pub struct HouseChoice {
        /// The `README.md` rule the choice is made under.
        pub under: HouseRule,
        /// What was considered and rejected in reaching it, in prose. Must not
        /// be blank.
        pub considered_and_rejected: String,
    }

    /// A record of how the project reads the spec where the spec conflicts with
    /// itself, and of what it has *not* decided.
    ///
    /// `style.md:9` binds the spec to RFC 2119, and items 3 and 4 there make
    /// "RECOMMENDED" equal to SHOULD and cover SHOULD NOT. Two SHOULD-strength
    /// clauses can therefore reach the same description and point opposite ways,
    /// with nothing in the text to break the tie — RFC 2119 instead licenses
    /// deviation once the implications are weighed. This record is where that
    /// weighing is written down.
    #[derive(Debug, Clone, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct Ruling {
        /// Stable identifier; pinned by the test consumer.
        pub id: String,
        /// The contested question, in one sentence.
        pub question: String,
        pub status: RulingStatus,
        /// Every clause the record is about, cited and quote-verified. A
        /// conflict record names both sides; a settled carve-out may name one.
        pub clauses: Vec<Citation>,
        /// The `clause` string, drawn from [`Self::clauses`], held to govern.
        /// Required when `decided`, forbidden when `undecided`.
        #[serde(default)]
        pub governing: Option<String>,
        /// The `clause` strings, drawn from [`Self::clauses`], deviated from in
        /// reaching the ruling. Necessarily empty when `undecided`: a record
        /// that has not chosen a side has not deviated from one either.
        #[serde(default)]
        pub deviates_from: Vec<String>,
        /// Present when the ruling is **ours** rather than the spec's: the
        /// recommendations are silent on the point and the project chose under
        /// `README.md` rule 5's silent limb or rule 6. See [`HouseChoice`].
        ///
        /// Mutually exclusive with [`Self::governing`] and
        /// [`Self::deviates_from`], and refused on an `undecided` record —
        /// enforced by [`validate_ruling_shape`], not by convention.
        #[serde(default)]
        pub house_choice: Option<HouseChoice>,
        /// Why. For an undecided record: what the project has and has not
        /// established, and what would settle it.
        pub rationale: String,
        /// A one-sentence, end-user-facing statement of the ruling, transcluded
        /// into the shadow spec's "Why" blocks so the reasoning is rendered for a
        /// human reader without a second copy of it (the full [`Self::rationale`]
        /// stays for agents). Optional: a record only needs one once a shadow-spec
        /// page cites it. The consuming checker
        /// (`tests/it/shadow_spec.rs`) fails a cited ruling that has none.
        // Deserialized only so `deny_unknown_fields` accepts the ledger; the
        // shadow-spec checker reads it, this generator does not.
        #[allow(dead_code)]
        #[serde(default)]
        pub summary: Option<String>,
        /// Row inputs the ruling bears on. Every one must exist as a row.
        #[serde(default)]
        pub applies_to: Vec<String>,
        /// Equivalence-class ids the ruling bears on, if any.
        #[serde(default)]
        pub equivalence_classes: Vec<String>,
        /// What enforces the ruling, or the stated reason nothing does.
        ///
        /// **`Option<Guard>`, not `Guard`.** The type is what carries this, and
        /// it carries two separate things.
        ///
        /// It makes the field **optional at deserialize**. As a bare `Guard`, a
        /// record lacking the key failed to deserialize — and this struct is the
        /// generator's entry to the whole ledger, so the failure was not "one
        /// record is unanswered" but "the ledger cannot be read".
        /// `generate_spec_fixture` exited 1, the generated fixture never existed,
        /// and every test downstream of it failed for a reason mentioning neither
        /// guards nor the record at fault. Measured on this branch rebased onto
        /// `origin/main` at `426a944b`, for two records missing one field each:
        /// **20** failures from this alone, on top of the **28** the ledger
        /// reader's own panic produced.
        ///
        /// And it keeps absence **distinguishable**. As a bare `Guard`, an absent
        /// key and a written `"guard": {}` would both arrive as
        /// `Guard { tests: [], none: None }`, and [`validate_guard_shape`] would
        /// have to give one message for two different mistakes. `None` is
        /// reachable only by the key being absent.
        ///
        /// The `#[serde(default)]` below is **redundant and deliberate**. Serde
        /// already treats a missing `Option` field as `None`, so removing it
        /// changes nothing today — verified, not assumed: with the attribute
        /// deleted the generator still exits 0 on a ledger carrying a guardless
        /// record. It is kept because the *intent* is "absent is allowed here",
        /// and a future change of this field's type away from `Option` would
        /// otherwise silently restore the cascade above. Do not read it as the
        /// thing that makes absence legal; the `Option` is.
        ///
        /// **This is a change of stage, not a relaxation.** A record with no
        /// guard still cannot merge: it is refused by
        /// `tests/it/ruling_guard_field.rs::every_record_declares_a_guard`, in
        /// one targeted failure that names the record and the fix. What the
        /// generator keeps refusing is a guard that is *present and malformed*,
        /// because that is a broken record rather than an unanswered question.
        #[serde(default)]
        pub guard: Option<Guard>,
    }

    /// See [`Ruling::guard`]. Exactly one of the two keys may be set;
    /// `validate_ruling_shape` enforces that, because serde cannot express
    /// "one of these, not both, not neither" on its own.
    #[derive(Debug, Clone, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct Guard {
        /// `<path>.rs::<function>` citations of the enforcing tests.
        #[serde(default)]
        pub tests: Vec<String>,
        /// A stated, reasoned declaration that no test enforces this record.
        #[serde(default)]
        pub none: Option<String>,
    }

    /// `spec_expected` uses the doubly-optional pattern so the override file
    /// can distinguish three states:
    ///   * field absent     → `None`             (use the auto-default)
    ///   * field set to null → `Some(None)`       (force "spec rejects this")
    ///   * field set to "X" → `Some(Some("X"))`  (force a canonical output)
    #[derive(Debug, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct OverrideEntry {
        #[serde(default)]
        pub status: Option<String>,
        #[serde(default, deserialize_with = "deserialize_some")]
        pub spec_expected: Option<Option<String>>,
        /// Free-form note attached to the row in the emitted fixture.
        /// Used to capture nuances where the auto-classified `status`
        /// doesn't fully convey ferro's behavior — e.g. rows whose
        /// `current` is the spec's recommended canonical form even
        /// though the original `input` is `<code class="invalid">`.
        /// Closes-after: #353.
        #[serde(default)]
        pub note: Option<String>,
        /// Force a specific prefixed form for bare fragments. Wins over the
        /// auto-default in `prefix::default_prefixed`.
        #[serde(default)]
        pub input_prefixed: Option<String>,
        /// Mark this row as requiring a real reference sequence to evaluate
        /// (e.g. §2.1 3'-rule shifting examples). The test consumer skips
        /// `Some(true)` rows until #82's `from_manifest` loader lands. There
        /// is no auto-detection — auditors flip this field by hand during
        /// the #83 sweep when they see a row whose correctness depends on
        /// reference data.
        #[serde(default)]
        pub requires_reference: Option<bool>,
    }

    fn deserialize_some<'de, T, D>(d: D) -> Result<Option<T>, D::Error>
    where
        T: Deserialize<'de>,
        D: Deserializer<'de>,
    {
        T::deserialize(d).map(Some)
    }

    pub fn load(path: &Path) -> anyhow::Result<Overrides> {
        if !path.exists() {
            return Ok(Overrides::default());
        }
        let text = std::fs::read_to_string(path)?;
        let parsed: Overrides = serde_json::from_str(&text)
            .map_err(|e| anyhow::anyhow!("parse {}: {e}", path.display()))?;
        Ok(parsed)
    }
}

// ---------- runner ----------

mod runner {
    use super::*;
    use ferro_hgvs::reference::mock::MockProvider;
    use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

    /// What ferro did with one harvested input.
    ///
    /// `current` is the canonical rendering — ferro's default `Display`, and the
    /// string recorded on the row. `equivalent_renderings` holds the *other*
    /// spellings the spec admits for that same normalized variant; conformance
    /// means the spec's stated form matches any of them. The rule lives in
    /// `HgvsVariant::spec_equivalent_renderings` so this generator and the
    /// `ferro_produces_the_form_the_spec_states` oracle cannot drift apart.
    struct FerroOutcome {
        current: String,
        equivalent_renderings: Vec<String>,
        parse_ok: bool,
        normalize_ok: bool,
        expected_warnings: Vec<String>,
    }

    #[derive(Debug, Serialize)]
    pub struct Row {
        pub input: String,
        /// Default-prefixed form for bare illustrative fragments
        /// (e.g. `c.1083A>C` → `NM_004006.2:c.1083A>C`). When `Some`, this is
        /// the string ferro was actually run against, and the test consumer
        /// asserts against this form. Verbatim spec text stays in `input`.
        #[serde(skip_serializing_if = "Option::is_none")]
        pub input_prefixed: Option<String>,
        pub current: String,
        /// `None` (serialized as JSON `null`) means the spec rejects this input
        /// — either via an explicit `<code class="invalid">…</code>` marker in
        /// the spec text or via an override entry with `spec_expected: null`.
        pub spec_expected: Option<String>,
        pub status: String,
        pub coordinate_system: String,
        pub source_kind: String,
        pub source_paths: Vec<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        pub working_group: Option<String>,
        /// Free-form clarifying note. Surfaced from the override file
        /// for rows where the `status` bucket doesn't fully capture
        /// ferro's behavior (e.g. canonical form emitted from a
        /// spec-rejected input). Closes-after: #353.
        #[serde(skip_serializing_if = "Option::is_none")]
        pub note: Option<String>,
        /// Mirrors the override flag — `Some(true)` means the row needs a
        /// real reference sequence to evaluate and is skipped by the test
        /// consumer. Surfaced on the row so consumers don't have to read
        /// the override file.
        #[serde(skip_serializing_if = "Option::is_none")]
        pub requires_reference: Option<bool>,
        /// Warning codes ferro currently emits for this row, sorted alphabetically.
        /// Empty for rows with no warnings.
        #[serde(skip_serializing_if = "Vec::is_empty")]
        pub expected_warnings: Vec<String>,
    }

    fn source_kind_priority(k: sources::SourceKind) -> u8 {
        use sources::SourceKind::*;
        match k {
            Recommendation => 5,
            SyntaxYaml => 4,
            Background => 3,
            Consultation => 2,
            PortedLegacyProbe => 1,
        }
    }

    fn source_kind_str(k: sources::SourceKind) -> &'static str {
        use sources::SourceKind::*;
        match k {
            Recommendation => "recommendation",
            SyntaxYaml => "syntax_yaml",
            Background => "background",
            Consultation => "consultation",
            PortedLegacyProbe => "ported_legacy_probe",
        }
    }

    fn coord_system(input: &str) -> String {
        // Find the first occurrence of `:<sys>.` and return the single-letter
        // system. Strings without a colon may be bare-prefix forms (e.g.
        // `c.123del`) used in some spec examples — handle that too.
        for prefix in ["c.", "g.", "m.", "n.", "o.", "p.", "r."] {
            if input.contains(&format!(":{prefix}")) || input.starts_with(prefix) {
                return prefix[..1].to_string();
            }
        }
        "?".to_string()
    }

    /// Chromosome landmarks the spec uses in place of a numeric position.
    /// They carry no edit meaning, so a coordinate built only from them (e.g.
    /// `g.pter`) is still edit-less.
    const POSITION_LANDMARKS: [&str; 3] = ["pter", "qter", "cen"];

    /// True when `input` carries no edit operator: a coordinate-system prefix
    /// (`c.`/`g.`/`m.`/`n.`/`o.`/`p.`/`r.`) followed only by position tokens —
    /// position characters (digits, `+ - * ?`) or a chromosome landmark
    /// (`pter`/`qter`/`cen`) — joined by `_`.
    ///
    /// The spec formats position references in prose (e.g. "the A-stretch
    /// running from position `c.5690` to `c.5697`") with the same inline-code
    /// markup as variant examples, so the harvester scrapes them; a bare
    /// coordinate is not a describable variant. Edit-less strings are not all
    /// prose references, though: the spec's whole-variant markers (`p.?`, `p.0`,
    /// `r.?`, …) and its edit-less negatives (`c.123-65_-50`, an incomplete
    /// range the spec requires be *rejected*) share this shape. So the caller
    /// pairs this test with evidence that the string is only a position — see
    /// [`is_position_reference`].
    fn is_edit_less_coordinate(input: &str) -> bool {
        // Strip a leading accession: the coordinate follows the *first* colon.
        // Splitting on the last colon instead would look past the junctions of a
        // chimeric description (`NC_1:g.pter_100::NC_2:g.200_cen_qter`) and judge
        // only its tail, wrongly calling the whole string edit-less.
        let core = input.split_once(':').map_or(input, |(_, rest)| rest);
        let rest = match core.get(0..2) {
            Some("c." | "g." | "m." | "n." | "o." | "p." | "r.") => &core[2..],
            _ => return false,
        };
        // A bare prefix (`c.`) is edit-less; but past that, an empty token from
        // `_`-splitting means a malformed range (`c.1_`, `c.1__2`) — not a clean
        // position, so it is not treated as an edit-less coordinate.
        rest.is_empty()
            || rest.split('_').all(|token| {
                !token.is_empty()
                    && (POSITION_LANDMARKS.contains(&token)
                        || token
                            .chars()
                            .all(|ch| matches!(ch, '0'..='9' | '+' | '-' | '*' | '?')))
            })
    }

    /// True when a harvested string is a position reference lifted out of spec
    /// prose rather than a variant, and so should not become a fixture row.
    ///
    /// Two signals, both requiring [`is_edit_less_coordinate`] (no edit
    /// operator) to hold:
    ///
    /// * ferro cannot parse it at all — a position cited in a sentence, e.g.
    ///   `c.` or a truncated fragment; or
    /// * ferro normalizes it by merely appending `=` (e.g. `g.12345678` →
    ///   `NC_000023.11:g.12345678=`), i.e. it invented the "identity" edit the
    ///   string never had. Scoring that as a divergence is noise.
    ///
    /// Edit-less strings that ferro parses *and* leaves alone are kept: the
    /// whole-variant markers (`p.?`, `p.0`) and the spec's edit-less negatives
    /// such as `c.123-65_-50`, which the checklist requires be rejected.
    ///
    /// `parse_ok`/`current` come from `target` (the accession-prefixed form)
    /// while the shape test reads the raw `input`. These are equivalent here:
    /// prefixing only prepends an accession, and `is_edit_less_coordinate`
    /// strips any accession before testing, so both forms share the shape.
    fn is_position_reference(input: &str, target: &str, current: &str, parse_ok: bool) -> bool {
        is_edit_less_coordinate(input) && (!parse_ok || current == format!("{target}="))
    }

    pub fn build_rows(
        candidates: &[sources::Candidate],
        overrides: &overrides::Overrides,
    ) -> anyhow::Result<Vec<Row>> {
        #[derive(Default)]
        struct Agg {
            kinds: Vec<sources::SourceKind>,
            paths: Vec<String>,
            wg: Option<String>,
            spec_rejects: bool,
        }
        let mut by_input: BTreeMap<String, Agg> = BTreeMap::new();
        for c in candidates {
            let a = by_input.entry(c.input.clone()).or_default();
            a.kinds.push(c.source_kind);
            a.paths.push(c.source_path.clone());
            if c.working_group.is_some() {
                a.wg = c.working_group.clone();
            }
            if c.intent == sources::SpanIntent::SpecRejects {
                a.spec_rejects = true;
            }
        }

        let normalizer = Normalizer::with_config(MockProvider::new(), measurement_config());

        let mut rows = Vec::with_capacity(by_input.len());
        for (input, agg) in by_input {
            let kind = agg
                .kinds
                .iter()
                .max_by_key(|k| source_kind_priority(**k))
                .copied()
                .unwrap();
            let mut paths: Vec<String> = agg.paths.clone();
            paths.sort();
            paths.dedup();

            let ov = overrides.by_input.get(&input);

            // input_prefixed resolution. The spec routinely writes bare
            // illustrative fragments like `c.1083A>C` whose accession is
            // implied by the surrounding paragraph. Pin a default accession
            // per coord system so ferro can actually parse them. Override
            // wins for rows that need a non-default prefix.
            let input_prefixed: Option<String> = ov
                .and_then(|o| o.input_prefixed.clone())
                .or_else(|| prefix::default_prefixed(&input));

            // What ferro is run against: the prefixed form when present,
            // else the raw input. This is also the anchor for the
            // spec_expected default (the spec wrote the bare form as
            // shorthand for the prefixed form).
            let target = input_prefixed.as_deref().unwrap_or(&input);
            let FerroOutcome {
                current,
                equivalent_renderings,
                parse_ok,
                normalize_ok,
                expected_warnings,
            } = run_ferro(&normalizer, target);

            // Drop position references harvested from spec prose (e.g.
            // "position `c.5690`"): they are not describable variants, just
            // positions cited in a sentence. An explicit override keeps the row
            // regardless (the auditor deliberately classified it).
            if ov.is_none() && is_position_reference(&input, target, &current, parse_ok) {
                continue;
            }

            // spec_expected resolution. Source order: override > structural
            // extraction > default-to-target.
            //   * override.spec_expected = Some(value) → use value verbatim
            //     (Some(None) forces "spec rejects"; Some(Some(s)) forces a
            //     canonical output)
            //   * spec marks the input via `<code class="invalid">…</code>`
            //     → None (spec rejects)
            //   * otherwise → Some(target) (spec offered the target as
            //     a canonical form; ferro is expected to round-trip it)
            let spec_expected: Option<String> = match ov.and_then(|o| o.spec_expected.clone()) {
                Some(value) => value,
                None => {
                    if agg.spec_rejects {
                        None
                    } else {
                        Some(target.to_string())
                    }
                }
            };

            // Status is derived from the taxonomy table in `classify::classify`.
            // Overrides may pin a non-default status (e.g. for rows the
            // auditor has hand-classified into a sub-bucket).
            let auto_status = classify::classify(
                parse_ok,
                normalize_ok,
                &current,
                &equivalent_renderings,
                spec_expected.as_deref(),
            );
            let status = ov
                .and_then(|o| o.status.clone())
                .unwrap_or_else(|| auto_status.to_string());

            rows.push(Row {
                coordinate_system: coord_system(&input),
                input: input.clone(),
                input_prefixed,
                current,
                spec_expected,
                status,
                source_kind: source_kind_str(kind).to_string(),
                source_paths: paths,
                working_group: agg.wg,
                note: ov.and_then(|o| o.note.clone()),
                requires_reference: ov.and_then(|o| o.requires_reference),
                expected_warnings,
            });
        }

        let row_inputs: std::collections::BTreeSet<&str> =
            rows.iter().map(|r| r.input.as_str()).collect();
        let unknown: Vec<&str> = overrides
            .by_input
            .keys()
            .filter(|k| !row_inputs.contains(k.as_str()))
            .map(String::as_str)
            .collect();
        if !unknown.is_empty() {
            anyhow::bail!(
                "overrides reference unknown inputs (typo or spec drift): {:?}",
                unknown
            );
        }

        rows.sort_by(|a, b| {
            a.coordinate_system
                .cmp(&b.coordinate_system)
                .then(a.status.cmp(&b.status))
                .then(a.input.cmp(&b.input))
        });
        Ok(rows)
    }

    /// The error configuration every row in this fixture is **normalized** under
    /// (#1629).
    ///
    /// Named rather than implied. `Normalizer::new` takes `NormalizeConfig::default()`,
    /// which substitutes **`ErrorConfig::lenient()`** — while `ErrorConfig::default()`
    /// is `strict()` and `ErrorMode`'s own `#[default]` is `Strict`. So "the default"
    /// was never an accurate description of this measurement, and until #1629 it was
    /// the only one the artifact offered. `render` derives the document's
    /// `normalized_under` stamp from this same value, so the artifact's claim and the
    /// normalizer that produced it cannot drift.
    ///
    /// **It covers the normalize half only, and the artifact's field is named for
    /// that.** `run_ferro` reaches the AST through the bare [`parse_hgvs`], which
    /// constructs no `InputPreprocessor` and therefore applies *no* `ErrorConfig` at
    /// all — not lenient's repairs and not strict's refusals (#1632, pinned by
    /// `tests/it/issue_1632_parse_entry_applies_no_mode.rs`). "No mode" is a third
    /// thing, not a synonym for strict, so a document-level stamp reading
    /// `generated_under: lenient` would have claimed a parser-level precondition no
    /// row ever had. Routing the parse through this config instead would change what
    /// the fixture *records*, which is a behaviour change and out of scope here.
    pub fn measurement_config() -> NormalizeConfig {
        NormalizeConfig::default()
    }

    /// A normalizer over the hermetic provider, for callers outside `build_rows`.
    pub fn new_normalizer() -> Normalizer<MockProvider> {
        Normalizer::with_config(MockProvider::new(), measurement_config())
    }

    /// Run one description and report `(rendering, evaluated_cleanly)`.
    ///
    /// `false` means ferro refused it — the rendering is then the parse or
    /// normalize error text, which is not a representation and must not be
    /// compared against one.
    pub fn observe(normalizer: &Normalizer<MockProvider>, target: &str) -> (String, bool) {
        let outcome = run_ferro(normalizer, target);
        let ok = outcome.parse_ok && outcome.normalize_ok;
        (outcome.current, ok)
    }

    fn run_ferro(normalizer: &Normalizer<MockProvider>, input: &str) -> FerroOutcome {
        let failed = |current: String, parse_ok: bool| FerroOutcome {
            current,
            equivalent_renderings: Vec::new(),
            parse_ok,
            normalize_ok: false,
            expected_warnings: Vec::new(),
        };
        match parse_hgvs(input) {
            Err(e) => failed(format!("parse error: {e}"), false),
            Ok(v) => match normalizer.normalize_with_diagnostics(&v) {
                Err(e) => failed(format!("normalize error: {e}"), true),
                Ok(n) => {
                    let mut codes: Vec<String> =
                        n.warnings.iter().map(|w| w.code().to_string()).collect();
                    codes.sort();
                    codes.dedup();
                    FerroOutcome {
                        equivalent_renderings: n.result.spec_equivalent_renderings(),
                        current: format!("{}", n.result),
                        parse_ok: true,
                        normalize_ok: true,
                        expected_warnings: codes,
                    }
                }
            },
        }
    }

    #[cfg(test)]
    mod tests {
        use super::{is_edit_less_coordinate, is_position_reference};

        #[test]
        fn flags_bare_position_references() {
            for bare in [
                "c.5690",
                "c.147",
                "c.1210-33",
                "NM_004006.1:c.5690",
                "AB053210.2:r.1289-365_1289-73",
                "c.", // bare prefix
                "c.*",
                "g.12345678",
                "g.pter", // chromosome landmark, still no edit operator
                "p.?",    // edit-less marker — shape matches; caller's gate keeps it
            ] {
                assert!(
                    is_edit_less_coordinate(bare),
                    "{bare:?} should be an edit-less coordinate"
                );
            }
        }

        #[test]
        fn does_not_flag_real_variants() {
            for variant in [
                "c.5697del",
                "c.76A>T",
                "NM_004006.2:c.1210-33del",
                "g.123_456dup",
                "LRG_199t1:c.11T>G",
                "p.Met1ext-5", // has an edit token (`ext`), not bare
                "g.1234=",
                // Malformed ranges: an empty `_`-token is not a clean position,
                // so these are not classified as edit-less (would otherwise be
                // wrongly dropped as prose if ferro rejects them).
                "c.1_",
                "c.1__2",
                // Chimeric junctions: edit-less-looking tail, but the whole
                // string is a (spec-rejected) description, not a position.
                "NC_000002.12:g.pter_8247756::NC_000011.10:g.15825273_cen_qter",
                "g.[chr11:pter_15825272::chr2:8247757_cen_qter]",
            ] {
                assert!(
                    !is_edit_less_coordinate(variant),
                    "{variant:?} is a real variant, not a bare coordinate"
                );
            }
        }

        #[test]
        fn drops_coordinates_ferro_merely_stamps_with_equals() {
            assert!(is_position_reference(
                "g.12345678",
                "NC_000023.11:g.12345678",
                "NC_000023.11:g.12345678=",
                true,
            ));
            assert!(is_position_reference(
                "g.pter",
                "NC_000023.11:g.pter",
                "NC_000023.11:g.pter=",
                true,
            ));
        }

        #[test]
        fn drops_unparsable_prose_fragments() {
            assert!(is_position_reference(
                "c.",
                "NM_004006.2:c.",
                "parse error",
                false
            ));
        }

        #[test]
        fn keeps_edit_less_spec_negatives_and_markers() {
            // `c.123-65_-50` is an incomplete range the spec requires be
            // rejected: edit-less, but ferro leaves it untouched, so it stays.
            assert!(!is_position_reference(
                "c.123-65_-50",
                "NM_004006.2:c.123-65_-50",
                "NM_004006.2:c.123-65_-50",
                true,
            ));
            // Whole-variant markers parse and round-trip unchanged.
            assert!(!is_position_reference(
                "p.?",
                "NP_000079.2:p.?",
                "NP_000079.2:p.?",
                true
            ));
            // Real variants never match, whatever ferro does with them.
            assert!(!is_position_reference(
                "g.123del",
                "NC_000023.11:g.123del",
                "NC_000023.11:g.123del",
                true,
            ));
        }
    }
}

// ---------- decisions ----------

/// Resolution and validation of the two curated decision-log sections:
/// equivalence classes (Task 1) and rulings (Task 2).
///
/// Everything here is a **build-time gate**. A class naming an input the
/// harvester no longer produces, a citation whose quote has moved in the spec
/// checkout, an `undecided` record that quietly grew a governing clause — each
/// fails generation, which is the same guard the `by_input` overrides already
/// get (`overrides reference unknown inputs`). The test consumer then re-derives
/// the verdicts live; nothing downstream trusts what is written here.
mod decisions {
    use super::*;
    use std::collections::BTreeSet;

    /// One member of a class, as evaluated.
    #[derive(Debug, Serialize)]
    pub struct ResolvedMember {
        /// The harvested input, verbatim.
        pub input: String,
        /// What ferro was actually run against — the row's target, re-anchored
        /// onto the class `reference` when the class declares one.
        pub target: String,
        /// ferro's rendering, or its refusal text when `evaluable` is false.
        pub normalized: String,
        /// False when this member cannot contribute to the comparison: it needs
        /// real reference bases, or ferro refuses it outright. Such a member is
        /// neither evidence of confluence nor of its absence.
        pub evaluable: bool,
        /// Why it is not evaluable. Absent when it is.
        #[serde(skip_serializing_if = "Option::is_none")]
        pub skipped_because: Option<String>,
    }

    #[derive(Debug, Serialize)]
    pub struct ResolvedClass {
        pub id: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        pub reference: Option<String>,
        pub citations: Vec<overrides::Citation>,
        pub note: String,
        pub members: Vec<ResolvedMember>,
        /// The distinct renderings the evaluable members produced. One entry is
        /// confluence; more than one is the defect the class exists to name.
        pub distinct_outputs: Vec<String>,
        pub verdict: &'static str,
    }

    #[derive(Debug, Serialize)]
    pub struct ResolvedRuling {
        pub id: String,
        pub question: String,
        pub status: overrides::RulingStatus,
        pub clauses: Vec<overrides::Citation>,
        #[serde(skip_serializing_if = "Option::is_none")]
        pub governing: Option<String>,
        #[serde(skip_serializing_if = "Vec::is_empty")]
        pub deviates_from: Vec<String>,
        /// Carried through so the generated fixture states, per record, whether
        /// its ruling rests on the recommendations or on the project's own
        /// judgement. Omitted when absent, so no existing record's serialization
        /// moves.
        #[serde(skip_serializing_if = "Option::is_none")]
        pub house_choice: Option<overrides::HouseChoice>,
        pub rationale: String,
        #[serde(skip_serializing_if = "Vec::is_empty")]
        pub applies_to: Vec<String>,
        #[serde(skip_serializing_if = "Vec::is_empty")]
        pub equivalence_classes: Vec<String>,
    }

    #[derive(Debug, Serialize)]
    pub struct Decisions {
        pub equivalence_classes: Vec<ResolvedClass>,
        pub rulings: Vec<ResolvedRuling>,
    }

    /// Verdict strings. Written on the row for readers; the test consumer
    /// re-derives its own rather than reading these.
    const CONFLUENT: &str = "confluent";
    const NON_CONFLUENT: &str = "non-confluent";
    const NOT_EVALUABLE: &str = "not-evaluable";

    /// True when `head` — the text before a target's first colon — is a
    /// reference accession rather than the opening of a bare coordinate
    /// fragment.
    ///
    /// HGVS puts the coordinate system immediately after the reference's colon,
    /// so a fragment carrying no accession *begins* with one of the
    /// single-letter descriptors followed by a dot. Anything else in that
    /// position is the reference — accessions (`NC_000023.11`, `LRG_199t1`,
    /// `NM_004006.2`) never take that shape.
    fn is_accession(head: &str) -> bool {
        !matches!(
            head.as_bytes(),
            [b'g' | b'o' | b'm' | b'c' | b'n' | b'r' | b'p', b'.', ..]
        )
    }

    /// Re-anchor `target` onto `reference`, replacing any accession it carries.
    ///
    /// The accession is everything before the **first** colon, which is where
    /// HGVS puts it; a colon later in the string belongs to an inner reference
    /// (`g.10_11ins[NC_…:g.20_30]`) and must be left alone.
    ///
    /// The first colon only delimits the accession when something
    /// accession-shaped precedes it (#1530). A member that is a bare fragment
    /// *and* carries an inner reference — `g.10_11ins[NC_2.1:g.20_30]` — puts
    /// the inner colon first, and splitting there yields `X.1:g.20_30`: a
    /// different, still well-formed variant, so nothing downstream would catch
    /// it. Latent while every curated member reaches here with a leading
    /// accession via `input_prefixed`; gated now rather than when it goes live.
    ///
    /// This is the same defect `spec_harvest::prefix::default_prefixed` fixed
    /// for the harvester in #955 — an inner `[…:<sys>.…]` read as the variant's
    /// own accession — arrived at from the opposite side: there a bare fragment
    /// was wrongly judged *already* prefixed, here it is wrongly *re*-anchored.
    fn anchor(target: &str, reference: Option<&str>) -> String {
        match reference {
            None => target.to_string(),
            Some(accession) => match target.split_once(':') {
                Some((head, rest)) if is_accession(head) => format!("{accession}:{rest}"),
                _ => format!("{accession}:{target}"),
            },
        }
    }

    /// Fail unless a ruling record is well formed *as a record* — everything
    /// that can be judged from the record alone, with no spec checkout and no
    /// fixture rows.
    ///
    /// Extracted so it can be unit-tested directly. The committed consumer
    /// `ruling_records_are_intact` asserts exactly these properties; a gate that
    /// is weaker than the test it backstops lets a malformed record into the
    /// artifact and only fails later, at the far end of the build.
    fn validate_ruling_shape(owner: &str, ruling: &overrides::Ruling) -> anyhow::Result<()> {
        if ruling.rationale.trim().is_empty() {
            anyhow::bail!("{owner}: a ruling record must say why");
        }
        if ruling.clauses.is_empty() {
            anyhow::bail!("{owner}: a ruling record must cite at least one clause");
        }
        let cited: BTreeSet<&str> = ruling.clauses.iter().map(|c| c.clause.as_str()).collect();
        for role in ruling.governing.iter().chain(ruling.deviates_from.iter()) {
            if !cited.contains(role.as_str()) {
                anyhow::bail!(
                    "{owner}: names {role:?} as governing or deviated-from, but that clause is \
                     not in `clauses` — every role must point at a cited, quote-verified clause"
                );
            }
        }
        // A house choice is checked BEFORE the status match, so that a record
        // which is both `undecided` and a house choice is refused for the reason
        // that names the contradiction rather than falling through to the
        // clause-count arm and being told something true but unhelpful.
        if let Some(house) = &ruling.house_choice {
            if ruling.status != overrides::RulingStatus::Decided {
                anyhow::bail!(
                    "{owner}: a `house_choice` on a non-`decided` record — a choice that has been \
                     made is decided by definition, and a record that has not chosen has no house \
                     choice to record"
                );
            }
            // The refusal this field was added for. A house choice asserts that
            // the recommendations do not reach the point; naming a clause that
            // governs asserts that one does. Both cannot be true, and the
            // failure mode that matters is the one where the second is written
            // down and then cited as compliance.
            if ruling.governing.is_some() {
                anyhow::bail!(
                    "{owner}: a `house_choice` names a governing clause. A house choice is made \
                     under `README.md` rule 5's silent limb or rule 6, which are what the project \
                     does where the recommendations do NOT decide — so it cannot also hold a \
                     clause to govern. Either drop `house_choice` and rule on the clause, or drop \
                     `governing` and record why no clause reaches"
                );
            }
            if !ruling.deviates_from.is_empty() {
                anyhow::bail!(
                    "{owner}: a `house_choice` names a deviated-from clause. There is nothing to \
                     deviate from where the recommendations are silent; a clause that reaches the \
                     point far enough to be departed from is a clause that reaches it"
                );
            }
            if house.considered_and_rejected.trim().is_empty() {
                anyhow::bail!(
                    "{owner}: a `house_choice` must say what was considered and rejected. A \
                     conclusion with its alternatives deleted is re-derived from scratch by the \
                     next reader, which is what these records exist to prevent"
                );
            }
        }
        match ruling.status {
            overrides::RulingStatus::Decided
                if ruling.governing.is_none() && ruling.house_choice.is_none() =>
            {
                anyhow::bail!(
                    "{owner}: a `decided` ruling must name the clause it holds to govern, or \
                     declare itself a `house_choice` under `README.md` rule 5's silent limb or \
                     rule 6"
                )
            }
            // The whole point of the `undecided` state is that nobody has
            // ruled. A governing clause — or a deviated-from one, which
            // implies a side was chosen — is an invented ruling, so it is a
            // build failure rather than a lint.
            overrides::RulingStatus::Undecided
                if ruling.governing.is_some() || !ruling.deviates_from.is_empty() =>
            {
                anyhow::bail!(
                    "{owner}: an `undecided` ruling must not name a governing or deviated-from \
                     clause — record the conflict, not a ruling nobody made"
                )
            }
            // The other half of the same rule, and the half this gate was
            // missing (#1530): `undecided` asserts that clauses *conflict*, so a
            // record citing one clause states no conflict — it states a position
            // while declining to name it as one. Only the `decided` arm may cite
            // a single clause, where it is a settled carve-out naming the clause
            // it applies.
            overrides::RulingStatus::Undecided if ruling.clauses.len() < 2 => anyhow::bail!(
                "{owner}: an `undecided` ruling cites {} clause(s); an unsettled question needs \
                 both sides on the record",
                ruling.clauses.len()
            ),
            _ => {}
        }
        validate_guard_shape(owner, ruling.guard.as_ref())?;
        Ok(())
    }

    /// A `guard`, if the record carries one, must state exactly one of the two
    /// things it can state, and must state it resolvably.
    ///
    /// This is the build-time half of the field, and it exists for the reason
    /// the unknown-equivalence-class check does: a pointer nothing validates is
    /// prose with punctuation. What is checked here is the shape and the
    /// **file**; whether the named function exists, runs and asserts needs the
    /// whole tree walked and is `tests/it/ruling_guard_field.rs`'s job.
    ///
    /// **`None` — no `guard` key at all — is passed through here deliberately**,
    /// and it is the one case this function declines to judge. Refusing it would
    /// put the generator back in the position the `#[serde(default)]` on
    /// [`overrides::Ruling::guard`] exists to get it out of: the generator is a
    /// *build step*, so a refusal from it is not one failure but every failure
    /// downstream of the artifact it did not write. The presence rule is
    /// enforced instead by a test that can name the record and print the fix,
    /// and that test is not optional, so nothing is weakened by declining here.
    ///
    /// Note the asymmetry that makes this safe: a record that says nothing is
    /// still *readable*, so the fixture built from it is correct as far as it
    /// goes. A record whose `guard` contradicts itself is not, which is why
    /// `Some(..)` is still judged in full.
    fn validate_guard_shape(owner: &str, guard: Option<&overrides::Guard>) -> anyhow::Result<()> {
        let Some(guard) = guard else {
            return Ok(());
        };
        match (guard.tests.is_empty(), guard.none.as_deref()) {
            (false, Some(_)) => anyhow::bail!(
                "{owner}: `guard` sets both `tests` and `none` — it has both named an enforcing \
                 test and declared that nothing enforces it"
            ),
            (true, None) => anyhow::bail!(
                "{owner}: `guard` states nothing. Set `tests` to the tests that fail when this \
                 ruling stops holding, or `none` to the reason nothing enforces it. An empty \
                 `guard` is the silence the field exists to reject"
            ),
            (true, Some(reason)) if reason.trim().is_empty() => anyhow::bail!(
                "{owner}: `guard.none` is blank; declining is a first-class answer only when it \
                 says why"
            ),
            (true, Some(_)) => Ok(()),
            (false, None) => {
                for citation in &guard.tests {
                    let (path, _) = split_guard_citation(citation).ok_or_else(|| {
                        anyhow::anyhow!(
                            "{owner}: guard citation {citation:?} is not of the form \
                             `<path>.rs::<function>` — a citation that does not resolve \
                             mechanically is the prose this field replaces"
                        )
                    })?;
                    if !crate_relative(path).is_file() {
                        anyhow::bail!(
                            "{owner}: guard citation {citation:?} names {path:?}, which is not a \
                             file in this tree (typo, or the guard moved)"
                        );
                    }
                }
                Ok(())
            }
        }
    }

    /// Splits `<path>.rs::<function>`; `None` for anything outside that grammar.
    ///
    /// Deliberately duplicated from `tests/it/common/rulings.rs` rather than
    /// shared: this binary and the integration-test crate do not link to each
    /// other, and the two copies are held together by
    /// `the_generator_and_the_ledger_reader_agree_on_the_citation_grammar`
    /// there, which is a check rather than a hope.
    fn split_guard_citation(token: &str) -> Option<(&str, &str)> {
        let (path, function) = token.split_once("::")?;
        if !path.ends_with(".rs") || path.starts_with('/') || path.contains("..") {
            return None;
        }
        if function.is_empty() || function.contains("::") {
            return None;
        }
        if !function.starts_with(|c: char| c.is_ascii_lowercase()) {
            return None;
        }
        if !function
            .chars()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_')
        {
            return None;
        }
        Some((path, function))
    }

    /// A repo-relative path resolved against the crate root, so the generator
    /// answers the same way whatever directory it is invoked from.
    fn crate_relative(path: &str) -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(path)
    }

    /// Split a `path:line` or `path:start-end` clause into its parts.
    fn parse_clause(clause: &str) -> anyhow::Result<(&str, usize, usize)> {
        let (path, lines) = clause
            .rsplit_once(':')
            .ok_or_else(|| anyhow::anyhow!("citation {clause:?} is not `<path>:<line>`"))?;
        let (first, last) = match lines.split_once('-') {
            Some((a, b)) => (a, b),
            None => (lines, lines),
        };
        let parse = |s: &str| {
            s.parse::<usize>()
                .map_err(|_| anyhow::anyhow!("citation {clause:?} has a non-numeric line"))
        };
        let (first, last) = (parse(first)?, parse(last)?);
        if first == 0 || last < first {
            anyhow::bail!("citation {clause:?} has an empty or inverted line range");
        }
        Ok((path, first, last))
    }

    /// Fail unless the citation resolves against the spec checkout **and** the
    /// quoted text is still on the cited lines.
    fn verify_citation(
        spec_dir: &Path,
        owner: &str,
        citation: &overrides::Citation,
    ) -> anyhow::Result<()> {
        let (path, first, last) = parse_clause(&citation.clause)?;
        let full = spec_dir.join(path);
        let text = std::fs::read_to_string(&full).map_err(|e| {
            anyhow::anyhow!(
                "{owner}: citation {:?} names {} which cannot be read: {e}",
                citation.clause,
                full.display()
            )
        })?;
        let lines: Vec<&str> = text.lines().collect();
        if last > lines.len() {
            anyhow::bail!(
                "{owner}: citation {:?} runs past the end of {path} ({} lines)",
                citation.clause,
                lines.len()
            );
        }
        if citation.quote.trim().is_empty() {
            anyhow::bail!("{owner}: citation {:?} has an empty quote", citation.clause);
        }
        // Join on a space so a quote may span the cited lines. Both sides are
        // whitespace-collapsed, so re-wrapping the spec's prose does not
        // invalidate a citation while moving it elsewhere still does.
        let collapse = |s: &str| s.split_whitespace().collect::<Vec<_>>().join(" ");
        let haystack = collapse(&lines[first - 1..last].join(" "));
        if !haystack.contains(&collapse(&citation.quote)) {
            anyhow::bail!(
                "{owner}: citation {:?} no longer quotes the spec. Expected to find\n  {:?}\nwithin those lines, which currently read\n  {:?}\nThe spec submodule moved under this citation — re-point it rather than deleting the quote.",
                citation.clause,
                citation.quote,
                haystack
            );
        }
        Ok(())
    }

    pub fn resolve(
        rows: &[runner::Row],
        overrides: &overrides::Overrides,
        spec_dir: &Path,
    ) -> anyhow::Result<Decisions> {
        let by_input: BTreeMap<&str, &runner::Row> =
            rows.iter().map(|r| (r.input.as_str(), r)).collect();
        let normalizer = runner::new_normalizer();

        let mut seen_ids: BTreeSet<&str> = BTreeSet::new();
        let mut classes = Vec::with_capacity(overrides.equivalence_classes.len());
        for class in &overrides.equivalence_classes {
            let owner = format!("equivalence class {:?}", class.id);
            if !seen_ids.insert(class.id.as_str()) {
                anyhow::bail!("{owner}: duplicate id");
            }
            if class.members.len() < 2 {
                anyhow::bail!("{owner}: a class needs at least two members");
            }
            if class.citations.is_empty() {
                anyhow::bail!(
                    "{owner}: a class must cite where the spec says its members are one variant"
                );
            }
            for citation in &class.citations {
                verify_citation(spec_dir, &owner, citation)?;
            }

            let mut members = Vec::with_capacity(class.members.len());
            let mut distinct: BTreeSet<String> = BTreeSet::new();
            for input in &class.members {
                let row = by_input.get(input.as_str()).ok_or_else(|| {
                    anyhow::anyhow!(
                        "{owner}: member {input:?} is not a fixture row (typo or spec drift)"
                    )
                })?;
                let target = anchor(
                    row.input_prefixed.as_deref().unwrap_or(&row.input),
                    class.reference.as_deref(),
                );
                let (normalized, ok) = runner::observe(&normalizer, &target);
                // Two ways a member cannot contribute. A row the auditor already
                // marked reference-dependent is skipped without running, since
                // the hermetic provider would only produce a lookup error; a row
                // ferro refuses is skipped on the refusal itself.
                let skipped_because =
                    if row.status == "needs-reference" || row.requires_reference == Some(true) {
                        Some("needs real reference bases".to_string())
                    } else if !ok {
                        Some(format!("ferro does not evaluate it: {normalized}"))
                    } else {
                        None
                    };
                if skipped_because.is_none() {
                    distinct.insert(normalized.clone());
                }
                members.push(ResolvedMember {
                    input: input.clone(),
                    target,
                    normalized,
                    evaluable: skipped_because.is_none(),
                    skipped_because,
                });
            }

            let verdict = match distinct.len() {
                0 | 1 if members.iter().filter(|m| m.evaluable).count() < 2 => NOT_EVALUABLE,
                1 => CONFLUENT,
                _ => NON_CONFLUENT,
            };
            classes.push(ResolvedClass {
                id: class.id.clone(),
                reference: class.reference.clone(),
                citations: class.citations.clone(),
                note: class.note.clone(),
                members,
                distinct_outputs: distinct.into_iter().collect(),
                verdict,
            });
        }

        let class_ids: BTreeSet<&str> = classes.iter().map(|c| c.id.as_str()).collect();
        let mut seen_rulings: BTreeSet<&str> = BTreeSet::new();
        let mut rulings = Vec::with_capacity(overrides.rulings.len());
        for ruling in &overrides.rulings {
            let owner = format!("ruling {:?}", ruling.id);
            if !seen_rulings.insert(ruling.id.as_str()) {
                anyhow::bail!("{owner}: duplicate id");
            }
            validate_ruling_shape(&owner, ruling)?;
            for citation in &ruling.clauses {
                verify_citation(spec_dir, &owner, citation)?;
            }
            for input in &ruling.applies_to {
                if !by_input.contains_key(input.as_str()) {
                    anyhow::bail!(
                        "{owner}: applies_to {input:?} is not a fixture row (typo or spec drift)"
                    );
                }
            }
            for id in &ruling.equivalence_classes {
                if !class_ids.contains(id.as_str()) {
                    anyhow::bail!("{owner}: names unknown equivalence class {id:?}");
                }
            }
            rulings.push(ResolvedRuling {
                id: ruling.id.clone(),
                question: ruling.question.clone(),
                status: ruling.status,
                clauses: ruling.clauses.clone(),
                governing: ruling.governing.clone(),
                deviates_from: ruling.deviates_from.clone(),
                house_choice: ruling.house_choice.clone(),
                rationale: ruling.rationale.clone(),
                applies_to: ruling.applies_to.clone(),
                equivalence_classes: ruling.equivalence_classes.clone(),
            });
        }

        Ok(Decisions {
            equivalence_classes: classes,
            rulings,
        })
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn anchor_replaces_only_the_leading_accession() {
            assert_eq!(
                anchor("NM_004006.2:c.[235A>T;237G>T]", Some("LRG_199t1")),
                "LRG_199t1:c.[235A>T;237G>T]"
            );
            // A colon inside an inserted-range payload is not the accession.
            assert_eq!(
                anchor("NC_1.1:g.10_11ins[NC_2.1:g.20_30]", Some("NC_9.9")),
                "NC_9.9:g.10_11ins[NC_2.1:g.20_30]"
            );
            assert_eq!(
                anchor("c.235_237delinsTAT", Some("X.1")),
                "X.1:c.235_237delinsTAT"
            );
            assert_eq!(anchor("NM_1.1:c.1del", None), "NM_1.1:c.1del");
        }

        /// The discriminating case for [`anchor`]: a member with **no** leading
        /// accession whose only colon is an inner reference's. `split_once(':')`
        /// alone splits there and returns a different variant
        /// (`X.1:g.20_30`) — silently, since the result is still well-formed.
        #[test]
        fn anchor_prefixes_a_bare_fragment_that_carries_an_inner_reference() {
            assert_eq!(
                anchor("g.10_11ins[NC_2.1:g.20_30]", Some("X.1")),
                "X.1:g.10_11ins[NC_2.1:g.20_30]"
            );
            // Every coordinate-system descriptor, so the shape test cannot be
            // narrowed to the one letter this case happens to use.
            for prefix in ["g", "o", "m", "c", "n", "r", "p"] {
                assert_eq!(
                    anchor(&format!("{prefix}.1_2ins[NC_2.1:g.20_30]"), Some("X.1")),
                    format!("X.1:{prefix}.1_2ins[NC_2.1:g.20_30]")
                );
            }
            // …and the converse still holds: a *real* leading accession is
            // replaced, not prefixed, even though it too precedes a colon.
            assert_eq!(
                anchor("NC_1.1:g.10_11ins[NC_2.1:g.20_30]", Some("X.1")),
                "X.1:g.10_11ins[NC_2.1:g.20_30]"
            );
        }

        fn citation(clause: &str) -> overrides::Citation {
            overrides::Citation {
                clause: clause.to_string(),
                quote: "quoted".to_string(),
            }
        }

        fn ruling(
            status: overrides::RulingStatus,
            clauses: &[&str],
            governing: Option<&str>,
            deviates_from: &[&str],
        ) -> overrides::Ruling {
            overrides::Ruling {
                id: "r".to_string(),
                question: "q".to_string(),
                status,
                clauses: clauses.iter().map(|c| citation(c)).collect(),
                governing: governing.map(str::to_string),
                deviates_from: deviates_from.iter().map(|c| c.to_string()).collect(),
                house_choice: None,
                rationale: "because".to_string(),
                summary: None,
                applies_to: Vec::new(),
                equivalence_classes: Vec::new(),
                guard: Some(guard_citing(&[
                    "examples/generate_spec_fixture.rs::a_guard",
                ])),
            }
        }

        /// A well-formed [`overrides::HouseChoice`], for the tests below.
        fn house(under: overrides::HouseRule) -> overrides::HouseChoice {
            overrides::HouseChoice {
                under,
                considered_and_rejected: "the position-wise reading, rejected because it is \
                                          not cheaper"
                    .to_string(),
            }
        }

        /// A `guard` naming the given citations. The default fixture cites a
        /// file that really exists, so `validate_guard_shape`'s file check is
        /// exercised positively rather than skipped.
        fn guard_citing(tests: &[&str]) -> overrides::Guard {
            overrides::Guard {
                tests: tests.iter().map(|t| t.to_string()).collect(),
                none: None,
            }
        }

        /// The gate this module exists to be: the *build* must refuse every
        /// malformed ruling that `ruling_records_are_intact` refuses, or the
        /// committed test is stricter than the generator that backstops it.
        #[test]
        fn an_undecided_ruling_must_put_both_sides_on_the_record() {
            use overrides::RulingStatus::{Decided, Undecided};

            // The one the generator used to let through: `undecided` is a claim
            // that two clauses conflict, so a single citation states no conflict.
            assert!(
                validate_ruling_shape("r", &ruling(Undecided, &["a.md:1"], None, &[])).is_err()
            );
            assert!(validate_ruling_shape(
                "r",
                &ruling(Undecided, &["a.md:1", "b.md:2"], None, &[])
            )
            .is_ok());

            // Already gated, re-asserted here so the arm cannot be rewritten in
            // terms of the clause count alone.
            assert!(validate_ruling_shape(
                "r",
                &ruling(Undecided, &["a.md:1", "b.md:2"], Some("a.md:1"), &[])
            )
            .is_err());
            assert!(validate_ruling_shape(
                "r",
                &ruling(Undecided, &["a.md:1", "b.md:2"], None, &["a.md:1"])
            )
            .is_err());

            // A `decided` record *may* cite one clause — a settled carve-out
            // names only the clause it applies. The two-clause rule is the
            // undecided arm's alone.
            assert!(
                validate_ruling_shape("r", &ruling(Decided, &["a.md:1"], Some("a.md:1"), &[]))
                    .is_ok()
            );
            assert!(validate_ruling_shape("r", &ruling(Decided, &["a.md:1"], None, &[])).is_err());
        }

        /// The build must refuse a house choice that claims spec authority.
        ///
        /// This is the mirror of [`an_undecided_ruling_must_put_both_sides_on_the_record`]
        /// and is written the same way — positive control first, then one
        /// mutation per refusal — because the two arms guard opposite
        /// overclaims. An `undecided` record naming a governing clause claims
        /// the project decided something it did not; a `house_choice` naming one
        /// claims the *spec* decided something it did not. The second is the
        /// costlier of the two here: a governing clause is quotable, so once one
        /// is written down the choice starts being cited as compliance.
        #[test]
        fn a_house_choice_may_not_claim_spec_authority() {
            use overrides::HouseRule::{RuleFiveSilent, RuleSix};
            use overrides::RulingStatus::{Decided, Undecided};

            // Positive control: a decided record with NO governing clause is
            // legal exactly when it declares itself a house choice. The second
            // assertion is the one that makes the first mean anything — without
            // the `house_choice` the same record is refused.
            let mut choice = ruling(Decided, &["a.md:1"], None, &[]);
            choice.house_choice = Some(house(RuleSix));
            assert!(
                validate_ruling_shape("r", &choice).is_ok(),
                "positive control"
            );
            let mut bare = choice.clone();
            bare.house_choice = None;
            assert!(validate_ruling_shape("r", &bare).is_err());

            // Rule 5's silent limb is the other admissible authority.
            let mut silent_limb = choice.clone();
            silent_limb.house_choice = Some(house(RuleFiveSilent));
            assert!(validate_ruling_shape("r", &silent_limb).is_ok());

            // The refusal this field exists for.
            let mut governs = choice.clone();
            governs.governing = Some("a.md:1".to_string());
            assert!(validate_ruling_shape("r", &governs).is_err());

            // …and its sibling: a clause reaching far enough to be departed from
            // is a clause that reaches.
            let mut deviates = ruling(Decided, &["a.md:1", "b.md:2"], None, &["b.md:2"]);
            deviates.house_choice = Some(house(RuleSix));
            assert!(validate_ruling_shape("r", &deviates).is_err());

            // A choice that has been made is decided by definition.
            let mut unsettled = ruling(Undecided, &["a.md:1", "b.md:2"], None, &[]);
            unsettled.house_choice = Some(house(RuleSix));
            assert!(validate_ruling_shape("r", &unsettled).is_err());

            // A conclusion with its alternatives deleted.
            let mut mute = choice.clone();
            mute.house_choice = Some(overrides::HouseChoice {
                under: RuleSix,
                considered_and_rejected: "   ".to_string(),
            });
            assert!(validate_ruling_shape("r", &mute).is_err());
        }

        /// `under` is a closed enum, so a record cannot write "the spec" into
        /// the slot reserved for the project's own authority.
        ///
        /// Checked at *deserialize*, which is the only place it can be checked:
        /// by the time [`validate_ruling_shape`] runs the value is already one
        /// of the two variants. A free-text field would have made this
        /// unenforceable, which is why it is not one.
        #[test]
        fn the_house_rule_is_a_closed_set() {
            let with = |under: &str| {
                serde_json::from_value::<overrides::HouseChoice>(serde_json::json!({
                    "under": under,
                    "considered_and_rejected": "something",
                }))
            };
            assert!(with("rule-six").is_ok());
            assert!(with("rule-five-silent").is_ok());
            assert!(with("general.md:34").is_err());
            assert!(with("rule-2").is_err());
        }

        /// An absent `house_choice` deserializes, so every existing record is
        /// unaffected by the field being added.
        #[test]
        fn an_absent_house_choice_deserializes() {
            let parsed: overrides::Ruling = serde_json::from_value(serde_json::json!({
                "id": "a-record",
                "status": "decided",
                "question": "does a record predating the field still parse?",
                "rationale": "it must, or adding the field re-statuses the whole ledger",
                "governing": "a.md:1",
                "clauses": [{ "clause": "a.md:1", "quote": "q" }],
            }))
            .expect("a record with no `house_choice` must deserialize");
            assert!(parsed.house_choice.is_none());
        }

        #[test]
        fn every_named_role_must_be_a_cited_clause() {
            use overrides::RulingStatus::Decided;
            assert!(
                validate_ruling_shape("r", &ruling(Decided, &["a.md:1"], Some("z.md:9"), &[]))
                    .is_err()
            );
            assert!(validate_ruling_shape(
                "r",
                &ruling(Decided, &["a.md:1"], Some("a.md:1"), &["z.md:9"])
            )
            .is_err());
        }

        /// The build refuses a `guard` that is present and says nothing usable.
        ///
        /// The positive control is first, so a `should`-fail below cannot pass
        /// because the fixture is broken for an unrelated reason.
        #[test]
        fn a_ruling_must_state_what_enforces_it() {
            use overrides::RulingStatus::Decided;
            let base = || ruling(Decided, &["a.md:1"], Some("a.md:1"), &[]);

            assert!(
                validate_ruling_shape("r", &base()).is_ok(),
                "positive control"
            );

            // Silence — the state the field exists to abolish. Written out as an
            // empty object, so it is the AUTHOR's blank rather than an absent
            // key; the absent key is `an_absent_guard_is_left_to_the_ledger_test`.
            let mut silent = base();
            silent.guard = Some(overrides::Guard {
                tests: Vec::new(),
                none: None,
            });
            assert!(validate_ruling_shape("r", &silent).is_err());

            // Declining IS an answer, and needs no test to exist.
            let mut declined = base();
            declined.guard = Some(overrides::Guard {
                tests: Vec::new(),
                none: Some("process record; it mandates no output".to_string()),
            });
            assert!(validate_ruling_shape("r", &declined).is_ok());

            // …but only when it says why.
            let mut blank = base();
            blank.guard = Some(overrides::Guard {
                tests: Vec::new(),
                none: Some("  ".to_string()),
            });
            assert!(validate_ruling_shape("r", &blank).is_err());

            // Both at once is a contradiction, not a belt-and-braces.
            let mut both = base();
            both.guard = Some(overrides::Guard {
                tests: vec!["examples/generate_spec_fixture.rs::a_guard".to_string()],
                none: Some("and also nothing".to_string()),
            });
            assert!(validate_ruling_shape("r", &both).is_err());
        }

        /// An **absent** `guard` deserializes and passes the generator, and is
        /// left to `tests/it/ruling_guard_field.rs::every_record_declares_a_guard`.
        ///
        /// Both halves matter and neither implies the other. The deserialize
        /// half is what stops a missing field taking the whole ledger — and
        /// therefore the generated fixture, and everything downstream of it —
        /// with it; the validate half is what stops the generator re-imposing
        /// the same refusal one function later, which would restore the cascade
        /// while looking like a shape check.
        ///
        /// Read alongside [`a_ruling_must_state_what_enforces_it`]: absence is
        /// passed through, a written-out blank is still refused. If those two
        /// ever agree, the distinction the `Option` exists for has been lost.
        #[test]
        fn an_absent_guard_is_left_to_the_ledger_test() {
            use overrides::RulingStatus::Decided;

            // Deserializes: a record with no `guard` key at all.
            let json = serde_json::json!({
                "id": "a-record",
                "status": "decided",
                "question": "does an absent guard deserialize?",
                "rationale": "it must, or one unanswered record breaks the build",
                "governing": "a.md:1",
                "clauses": [{ "clause": "a.md:1", "quote": "q" }],
            });
            let parsed: overrides::Ruling =
                serde_json::from_value(json).expect("a record with no `guard` must deserialize");
            assert!(
                parsed.guard.is_none(),
                "an absent key must arrive as `None`, distinguishable from a written `{{}}`"
            );

            // …and passes validation, which is where the generator stops caring.
            let mut absent = ruling(Decided, &["a.md:1"], Some("a.md:1"), &[]);
            absent.guard = None;
            assert!(
                validate_ruling_shape("r", &absent).is_ok(),
                "the generator must not refuse an absent guard — refusing it here is what \
                 turned one unanswered record into a fixture-wide cascade"
            );
        }

        /// A guard citation must be resolvable, and must name a file that is
        /// really there — the same class of check as "names unknown equivalence
        /// class", which is the only other structured pointer a record carries.
        #[test]
        fn a_guard_citation_must_name_a_file_in_this_tree() {
            use overrides::RulingStatus::Decided;
            let with = |tests: &[&str]| {
                let mut r = ruling(Decided, &["a.md:1"], Some("a.md:1"), &[]);
                r.guard = Some(guard_citing(tests));
                r
            };

            assert!(validate_ruling_shape(
                "r",
                &with(&["examples/generate_spec_fixture.rs::a_guard"])
            )
            .is_ok());

            // Outside the grammar: a module path names no file, and a bare path
            // names no proposition.
            assert!(validate_ruling_shape("r", &with(&["merge::a_guard"])).is_err());
            assert!(
                validate_ruling_shape("r", &with(&["examples/generate_spec_fixture.rs"])).is_err()
            );
            // Well-formed and pointing at nothing.
            assert!(validate_ruling_shape("r", &with(&["tests/it/no_such_file.rs::a"])).is_err());
        }

        /// The grammar itself, since the generator keeps its own copy of it.
        #[test]
        fn the_guard_citation_grammar_admits_a_path_and_nothing_else() {
            assert_eq!(
                split_guard_citation("tests/it/foo.rs::a_guard_name"),
                Some(("tests/it/foo.rs", "a_guard_name"))
            );
            for rejected in [
                "merge::a_guard_name",
                "tests/it/foo.rs",
                "a_guard_name",
                "tests/it/foo.rs::mod::a_guard",
                "tests/it/foo.rs::AGuard",
                "tests/it/foo.rs::",
                "/abs/foo.rs::a_guard",
                "tests/../../foo.rs::a_guard",
            ] {
                assert_eq!(split_guard_citation(rejected), None, "{rejected}");
            }
        }

        #[test]
        fn a_ruling_must_cite_a_clause_and_say_why() {
            use overrides::RulingStatus::Decided;
            assert!(
                validate_ruling_shape("r", &ruling(Decided, &[], Some("a.md:1"), &[])).is_err()
            );
            let mut silent = ruling(Decided, &["a.md:1"], Some("a.md:1"), &[]);
            silent.rationale = "   ".to_string();
            assert!(validate_ruling_shape("r", &silent).is_err());
        }

        #[test]
        fn parses_single_lines_and_ranges() {
            assert_eq!(
                parse_clause("docs/recommendations/DNA/delins.md:47").unwrap(),
                ("docs/recommendations/DNA/delins.md", 47, 47)
            );
            assert_eq!(
                parse_clause("docs/recommendations/DNA/delins.md:44-47").unwrap(),
                ("docs/recommendations/DNA/delins.md", 44, 47)
            );
            assert!(parse_clause("docs/recommendations/DNA/delins.md").is_err());
            assert!(parse_clause("a.md:0").is_err());
            assert!(parse_clause("a.md:9-4").is_err());
        }
    }
}

// ---------- render ----------

mod render {
    use super::*;
    use ferro_hgvs::conformance::error_mode_stamp::ErrorModeStamp;

    #[derive(Serialize)]
    struct Document<'a> {
        description: &'a str,
        spec: SpecBlock<'a>,
        generated_utc: String,
        /// The **normalizer's** error-handling precondition every row below was
        /// measured under (#1629). Derived from `runner::measurement_config`, so
        /// it cannot disagree with the normalizer that produced the rows.
        ///
        /// Named `normalized_under` rather than `generated_under` because that
        /// is its exact reach: the parse half of every row runs through the bare
        /// `parse_hgvs`, which applies no `ErrorConfig` at all (#1632). See
        /// `runner::measurement_config`.
        normalized_under: ErrorModeStamp,
        summary: Summary,
        /// Curated "these spellings are one variant" declarations, resolved
        /// against ferro's current output. See `decisions`.
        equivalence_classes: &'a [decisions::ResolvedClass],
        /// Curated records of contested spec readings. See `decisions`.
        rulings: &'a [decisions::ResolvedRuling],
        rows: &'a [runner::Row],
    }

    #[derive(Serialize)]
    struct SpecBlock<'a> {
        source: &'a str,
        tag: &'a str,
        commit_sha: String,
    }

    #[derive(Serialize, Default)]
    struct Summary {
        total: usize,
        by_status: BTreeMap<String, usize>,
        by_coordinate_system: BTreeMap<String, usize>,
        by_source_kind: BTreeMap<String, usize>,
        /// Equivalence-class verdicts, so the decision log's headline number
        /// sits beside the row census rather than having to be counted by hand.
        equivalence_classes_by_verdict: BTreeMap<String, usize>,
        rulings_by_status: BTreeMap<String, usize>,
    }

    fn read_submodule_commit(spec_dir: &Path) -> anyhow::Result<String> {
        // Clear the ambient GIT_* env before resolving the submodule pin. Git
        // exports GIT_DIR/GIT_WORK_TREE (pointing at the OUTER repo) when it
        // runs a hook, and those override the `-C spec_dir` repo discovery — so
        // under the pre-push `--check` hook this would otherwise resolve the
        // outer branch HEAD instead of the submodule's commit (#1046).
        let out = std::process::Command::new("git")
            .arg("-C")
            .arg(spec_dir)
            .args(["rev-parse", "HEAD"])
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

    pub fn render(
        rows: &[runner::Row],
        decisions: &decisions::Decisions,
        spec_dir: &Path,
    ) -> anyhow::Result<String> {
        let commit_sha = read_submodule_commit(spec_dir)?;
        let mut summary = Summary {
            total: rows.len(),
            ..Default::default()
        };
        for c in &decisions.equivalence_classes {
            *summary
                .equivalence_classes_by_verdict
                .entry(c.verdict.to_string())
                .or_default() += 1;
        }
        for r in &decisions.rulings {
            let key = serde_json::to_value(r.status)?
                .as_str()
                .unwrap_or("unknown")
                .to_string();
            *summary.rulings_by_status.entry(key).or_default() += 1;
        }
        for r in rows {
            *summary.by_status.entry(r.status.clone()).or_default() += 1;
            *summary
                .by_coordinate_system
                .entry(r.coordinate_system.clone())
                .or_default() += 1;
            *summary
                .by_source_kind
                .entry(r.source_kind.clone())
                .or_default() += 1;
        }
        let doc = Document {
            description:
                "HGVS v21.0 spec normalization fixture - pins ferro's current normalize() output \
                 for every variant string in the spec, plus the curated decision log: \
                 `equivalence_classes` (spellings the spec states denote one variant, and whether \
                 ferro converges them) and `rulings` (how the project reads clauses that conflict \
                 at equal RFC 2119 strength, including what it has not decided). \
                 Companion to issue #84.",
            spec: SpecBlock {
                source: "https://github.com/HGVSnomenclature/hgvs-nomenclature",
                tag: "21.0.0",
                commit_sha,
            },
            generated_utc: "fixture-byte-stable".to_string(),
            normalized_under: ErrorModeStamp::of(&runner::measurement_config().error_config),
            summary,
            equivalence_classes: &decisions.equivalence_classes,
            rulings: &decisions.rulings,
            rows,
        };
        let mut out = serde_json::to_string_pretty(&doc)?;
        out.push('\n');
        Ok(out)
    }
}
