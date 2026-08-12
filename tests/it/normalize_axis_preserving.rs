//! Normalization does not move a description onto a different coordinate axis.
//!
//! # What this asserts
//!
//! For every description the normalizer accepts, the coordinate-system letter of
//! the output equals that of the input: `g.` in gives `g.` out, and likewise for
//! `c.`/`n.`/`r.`/`p.`/`m.`/`o.`. Re-expressing one axis on another is the
//! *projector's* job (`src/project/`), and the dependency runs one way —
//! `src/normalize/` names no projector — so today the property holds by
//! construction. That is exactly why it is worth pinning: nothing in the type
//! system stops a future normalization rule from reaching for a projection, and
//! a consumer that stored `c.` strings would get `g.` ones back with no error.
//!
//! # The one sanctioned rewrite, and why it is not a counter-example
//!
//! A `g.` description on a mitochondrial reference (`NC_012920` / `NC_001807`)
//! is re-labelled `m.` (#487, `Normalizer::normalize_dispatch`). The reference
//! molecule is unchanged — only the coordinate-system letter HGVS requires for
//! it is corrected — so this is a label repair rather than a change of axis.
//! It is admitted here as an enumerated exception and each occurrence is
//! re-checked against [`Accession::is_mitochondrial`], so the carve-out cannot
//! launder a `g.` -> `m.` rewrite on any other accession.
//!
//! # The near-miss to keep in mind: the accession may change, the axis may not
//!
//! Legacy-selector resolution rewrites `NG_012337.1(NM_003002.2):c.274G>T` to
//! `NM_003002.2:c.274G>T`. The *accession* changes; the axis letter does not.
//! So the comparison here is on [`CoordinateAxis`] alone and never on the
//! description's prefix — `the_check_reads_the_axis_letter_not_the_prefix`
//! pins that distinction directly, because a prefix comparison would report
//! that rewrite as a violation and the corpora (which run with no reference
//! provider) cannot exhibit it.
//!
//! # Denominators are part of the result
//!
//! "0 violations" is meaningless without the population it was measured over, and
//! a census that compared nothing looks identical to a clean one. Every test here
//! prints its full accounting — rows, unparsed, unnormalizable, axis-less,
//! compared, and the per-axis breakdown of what was compared — and asserts the
//! compared count is non-zero before believing a zero. The first attempt at this
//! measurement reported `compared: 0`, from a mis-guessed JSON field name; the
//! denominator is what caught it.
//!
//! # What this does NOT claim
//!
//! Nothing about whether the normalized *description* is right — only about its
//! axis. Nothing about descriptions the normalizer rejects, or that the parser
//! rejects: both are counted and skipped, not judged. And nothing about the
//! projector, whose whole purpose is to change axis.
//!
//! The bulk fixtures are release assets rather than git objects
//! (`scripts/fetch-test-fixtures.sh`). Following the sibling corpus suites
//! (`clinvar_hgvs_tests`, `cmrg_exhaustive_tests`, `paraphase_exhaustive_tests`),
//! a missing fixture skips green rather than failing — so an absent corpus is
//! reported as a skip, never as a pass. Under `FERRO_REQUIRE_BULK_FIXTURES`,
//! which CI sets, that skip becomes a hard failure instead; see
//! `common::bulk_fixtures`.

use ferro_hgvs::{parse_hgvs, CoordinateAxis, HgvsVariant, MockProvider, Normalizer};
use flate2::read::GzDecoder;
use rayon::prelude::*;
use serde::Deserialize;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::Read;
use std::time::Instant;

/// Slim deserialization shape shared by the four bulk fixtures: `input` is the
/// only field this census reads. `&'a str` borrows into the decompressed buffer
/// — see `cmrg_exhaustive_tests::load_fixture_bytes` for why that matters here.
#[derive(Deserialize)]
struct BulkFixture<'a> {
    #[serde(borrow)]
    test_cases: Vec<BulkCase<'a>>,
}

#[derive(Deserialize)]
struct BulkCase<'a> {
    #[serde(borrow)]
    input: &'a str,
}

impl AsRef<str> for BulkCase<'_> {
    fn as_ref(&self) -> &str {
        self.input
    }
}

/// The generated HGVS spec-normalization fixture. Read for the axes the bulk
/// corpora barely reach — it is the only population here carrying an `o.` row,
/// and it holds roughly nineteen times as many `r.` rows as all four bulk
/// corpora combined (measured: 130 against 7).
#[derive(Deserialize)]
struct SpecFixture {
    rows: Vec<SpecRow>,
}

#[derive(Deserialize)]
struct SpecRow {
    input: String,
    /// Default-prefixed form for the bare illustrative fragments the spec states
    /// without an accession (#84); `input` verbatim when there is none.
    #[serde(default)]
    input_prefixed: Option<String>,
}

fn load_gzipped(path: &str) -> Option<Vec<u8>> {
    // Absent means "skip" locally and "fail" under `FERRO_REQUIRE_BULK_FIXTURES`,
    // which CI sets — see `common::bulk_fixtures`.
    crate::common::bulk_fixtures::present_or_skip(path)?;
    let file = File::open(path).unwrap_or_else(|e| panic!("failed to open {path}: {e}"));
    let mut buf = Vec::new();
    GzDecoder::new(file)
        .read_to_end(&mut buf)
        .unwrap_or_else(|e| panic!("failed to decompress {path}: {e}"));
    Some(buf)
}

/// One corpus's accounting. Every input lands in exactly one of the four
/// disposition counters, so `rows` is recoverable from the census and a silently
/// dropped population is not representable.
#[derive(Default)]
struct AxisCensus {
    /// Rows offered to the census.
    rows: usize,
    /// `parse_hgvs` declined — not this module's subject.
    unparsed: usize,
    /// `normalize` returned `Err` — likewise not this module's subject.
    unnormalizable: usize,
    /// One side has no single coordinate axis: an RNA fusion (`::`, which joins
    /// two transcripts) or a bare `[0]`/`[?]` marker. There is no axis letter to
    /// compare, so there is nothing to assert.
    axis_less: usize,
    /// Rows where both sides had an axis and the two were compared.
    compared: usize,
    /// Of `compared`, those whose axis letter is unchanged.
    preserved: usize,
    /// Of `compared`, the #487 `g.` -> `m.` label repair on a mitochondrial
    /// reference, verified per occurrence.
    mito_relabels: usize,
    /// Input-axis breakdown of `compared`, so a zero can be read per axis.
    per_axis: BTreeMap<&'static str, usize>,
    /// Of `compared`, how many violations occurred — uncapped, unlike
    /// [`Self::violations`], which holds only a printable sample. Reporting
    /// `violations.len()` as the count would describe any regression wider than
    /// [`MAX_REPORTED_VIOLATIONS`] as exactly that many rows, understating it at
    /// precisely the moment the number matters most.
    violations_seen: usize,
    /// A readable sample of the axis changes that are not the sanctioned
    /// relabel, capped so a broad regression prints a sample rather than
    /// millions of lines.
    violations: Vec<String>,
}

/// How the census classifies one `(input, output)` pair.
///
/// Extracted so `the_check_reads_the_axis_letter_not_the_prefix` can exercise
/// the census's *own* decision instead of restating it. A test that re-derives
/// the comparison it means to pin passes no matter what the census does — which
/// is what that test did before this was factored out.
#[derive(Debug, PartialEq, Eq)]
enum Disposition {
    /// The axis letter is unchanged.
    Preserved,
    /// The #487 `g.` -> `m.` label repair on a mitochondrial reference.
    MitochondrialRelabel,
    /// Any other axis change.
    Violation,
}

/// Classify one pair, or `None` when either side has no single coordinate axis.
///
/// The comparison is on [`CoordinateAxis`] and on nothing else — in particular
/// not on the accession or the rendered prefix, both of which normalization is
/// allowed to move (legacy-selector resolution does exactly that).
fn classify(variant: &HgvsVariant, normalized: &HgvsVariant) -> Option<Disposition> {
    let (before, after) = (variant.coordinate_axis()?, normalized.coordinate_axis()?);
    Some(if before == after {
        Disposition::Preserved
    } else if is_mitochondrial_relabel(variant, before, after) {
        Disposition::MitochondrialRelabel
    } else {
        Disposition::Violation
    })
}

const MAX_REPORTED_VIOLATIONS: usize = 20;

impl AxisCensus {
    fn of_one(input: &str) -> Self {
        let mut census = Self {
            rows: 1,
            ..Self::default()
        };
        let Ok(variant) = parse_hgvs(input) else {
            census.unparsed = 1;
            return census;
        };
        // An empty provider, deliberately: the property is about the axis label,
        // which no reference can inform. Normalization degrades leniently
        // without one, so a large majority of every corpus still normalizes.
        let normalizer = Normalizer::new(MockProvider::new());
        let Ok(normalized) = normalizer.normalize(&variant) else {
            census.unnormalizable = 1;
            return census;
        };
        let (Some(before), Some(after)) = (variant.coordinate_axis(), normalized.coordinate_axis())
        else {
            census.axis_less = 1;
            return census;
        };
        census.compared = 1;
        *census.per_axis.entry(before.code()).or_default() += 1;
        match classify(&variant, &normalized).expect("both sides carry an axis here") {
            Disposition::Preserved => census.preserved = 1,
            Disposition::MitochondrialRelabel => census.mito_relabels = 1,
            Disposition::Violation => {
                census.violations_seen = 1;
                census.violations.push(format!(
                    "`{input}` [{}] normalized to `{normalized}` [{}]",
                    before.code(),
                    after.code()
                ));
            }
        }
        census
    }

    fn merge(mut self, other: Self) -> Self {
        self.rows += other.rows;
        self.unparsed += other.unparsed;
        self.unnormalizable += other.unnormalizable;
        self.axis_less += other.axis_less;
        self.compared += other.compared;
        self.preserved += other.preserved;
        self.mito_relabels += other.mito_relabels;
        self.violations_seen += other.violations_seen;
        for (axis, count) in other.per_axis {
            *self.per_axis.entry(axis).or_default() += count;
        }
        for violation in other.violations {
            if self.violations.len() < MAX_REPORTED_VIOLATIONS {
                self.violations.push(violation);
            }
        }
        self
    }

    /// Print the accounting, then assert on it.
    ///
    /// The print is unconditional and comes first: a run that asserts without
    /// reporting its denominator cannot be told apart from one that compared
    /// nothing, which is the failure this module is written against.
    fn report_and_assert(&self, label: &str, elapsed_secs: f64) {
        eprintln!("\n=== axis preservation: {label} ===");
        eprintln!("  rows offered      : {}", self.rows);
        eprintln!("  unparsed          : {}", self.unparsed);
        eprintln!("  unnormalizable    : {}", self.unnormalizable);
        eprintln!("  no single axis    : {}", self.axis_less);
        eprintln!("  COMPARED          : {}", self.compared);
        eprintln!("    axis preserved  : {}", self.preserved);
        eprintln!("    #487 g. -> m.   : {}", self.mito_relabels);
        eprintln!("    violations      : {}", self.violations_seen);
        eprintln!("  compared by input axis:");
        for (axis, count) in &self.per_axis {
            eprintln!("    {axis}. : {count}");
        }
        eprintln!("  elapsed: {elapsed_secs:.2}s\n");

        assert!(
            self.compared > 0,
            "{label}: compared 0 rows of {} offered — a structural zero, not a clean result",
            self.rows
        );
        assert_eq!(
            self.violations_seen,
            0,
            "{label}: normalization changed the coordinate axis on {} of {} compared rows \
             ({} shown):\n  {}",
            self.violations_seen,
            self.compared,
            self.violations.len(),
            self.violations.join("\n  ")
        );
        assert_eq!(
            self.preserved + self.mito_relabels + self.violations_seen,
            self.compared,
            "{label}: the dispositions do not account for every compared row"
        );
        // The doc on `AxisCensus` claims every input lands in exactly one of the
        // four disposition counters. Check it rather than assert it in prose: a
        // population that stopped being counted would otherwise shrink the
        // denominator silently, which is the one failure this module is written
        // against.
        assert_eq!(
            self.unparsed + self.unnormalizable + self.axis_less + self.compared,
            self.rows,
            "{label}: the dispositions do not account for every row offered"
        );
    }
}

/// Whether an axis change is the #487 mitochondrial label repair.
///
/// Checked per occurrence rather than by matching the letters alone: the
/// accession must itself be mitochondrial, so the carve-out admits only the
/// rewrite `normalize_dispatch` actually performs.
fn is_mitochondrial_relabel(
    variant: &HgvsVariant,
    before: CoordinateAxis,
    after: CoordinateAxis,
) -> bool {
    before == CoordinateAxis::Genomic
        && after == CoordinateAxis::Mitochondrial
        && variant
            .accession()
            .is_some_and(|accession| accession.is_mitochondrial())
}

fn census_of<I>(inputs: I) -> AxisCensus
where
    I: IntoParallelIterator,
    I::Item: AsRef<str>,
{
    inputs
        .into_par_iter()
        .map(|input| AxisCensus::of_one(input.as_ref()))
        .reduce(AxisCensus::default, AxisCensus::merge)
}

/// Run the census over one gzipped bulk fixture, skipping green when the
/// release-asset fixture is absent — the idiom the sibling corpus suites use —
/// or failing, under `FERRO_REQUIRE_BULK_FIXTURES`.
fn assert_axis_preserved_over_bulk_fixture(path: &str) {
    // `load_gzipped` has already reported the skip, or panicked under
    // `FERRO_REQUIRE_BULK_FIXTURES`.
    let Some(buf) = load_gzipped(path) else {
        return;
    };
    let fixture: BulkFixture<'_> =
        serde_json::from_slice(&buf).unwrap_or_else(|e| panic!("failed to parse {path}: {e}"));
    let start = Instant::now();
    let census = census_of(&fixture.test_cases[..]);
    census.report_and_assert(path, start.elapsed().as_secs_f64());
}

#[test]
fn axis_is_preserved_over_the_clinvar_500k_corpus() {
    assert_axis_preserved_over_bulk_fixture("tests/fixtures/bulk/clinvar_hgvs_500k.json.gz");
}

#[test]
fn axis_is_preserved_over_the_clinvar_unique_corpus() {
    assert_axis_preserved_over_bulk_fixture("tests/fixtures/bulk/clinvar_hgvs_unique.json.gz");
}

#[test]
fn axis_is_preserved_over_the_cmrg_corpus() {
    assert_axis_preserved_over_bulk_fixture(
        "tests/fixtures/validation/cmrg_genes_exhaustive.json.gz",
    );
}

#[test]
fn axis_is_preserved_over_the_paraphase_corpus() {
    assert_axis_preserved_over_bulk_fixture(
        "tests/fixtures/validation/paraphase_genes_exhaustive.json.gz",
    );
}

/// The spec fixture is small, and it is here for axis *breadth* rather than
/// volume. Measured over the four bulk corpora: `o.` is reached **0** times and
/// `r.` **7** times, against 1 and 130 in this fixture. A census over ten
/// million observed variants that never reaches one of the seven axes at all,
/// and reaches another seven times, is not a census over seven axes.
///
/// Those are the figures the per-axis breakdown prints, so they can be
/// re-derived from any run rather than taken from this comment.
#[test]
fn axis_is_preserved_over_the_spec_normalization_corpus() {
    crate::common::spec_fixture::ensure_spec_fixture();
    let path = crate::common::spec_fixture::spec_fixture_path();
    let content = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
    let fixture: SpecFixture = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("failed to parse {}: {e}", path.display()));
    let inputs: Vec<String> = fixture
        .rows
        .into_iter()
        .map(|row| row.input_prefixed.unwrap_or(row.input))
        .collect();
    let start = Instant::now();
    let census = census_of(inputs);
    census.report_and_assert(
        "hgvs_spec_normalization.json",
        start.elapsed().as_secs_f64(),
    );
}

/// The comparison is on the axis letter, never on the description's prefix.
///
/// Legacy-selector resolution replaces an `NG_`-anchored selector with the bare
/// transcript, so the accession moves while the axis stays put. A check written
/// on the prefix would call that a violation. This pins the distinction
/// hermetically — the two spellings differ in accession and agree in axis —
/// because the corpus censuses run with no reference provider and therefore
/// never trigger the rewrite themselves (measured: 0 accession changes across
/// all four bulk corpora, so those censuses cannot stand in for this test).
///
/// It asserts on [`classify`], the census's own decision, rather than
/// re-deriving the comparison. That distinction is the whole test: an earlier
/// version called only `parse_hgvs` and `coordinate_axis`, so it survived a
/// mutation that made the census compare accessions as well as axes — exactly
/// the check its name forbids. Routed through `classify`, that mutation now
/// turns this red.
#[test]
fn the_check_reads_the_axis_letter_not_the_prefix() {
    let selector_form = parse_hgvs("NG_012337.1(NM_003002.2):c.274G>T").expect("parses");
    let resolved_form = parse_hgvs("NM_003002.2:c.274G>T").expect("parses");

    assert_ne!(
        selector_form.accession().map(|a| a.to_string()),
        resolved_form.accession().map(|a| a.to_string()),
        "the two spellings must differ in accession for this to be the near-miss it claims"
    );
    assert_eq!(
        classify(&selector_form, &resolved_form),
        Some(Disposition::Preserved),
        "the census must read legacy-selector resolution as axis-preserving; a comparison \
         that also looked at the accession or the rendered prefix would call it a violation"
    );
    assert_eq!(
        selector_form.coordinate_axis().map(CoordinateAxis::code),
        Some("c")
    );
}

/// The sanctioned exception is scoped to the mitochondrial accessions, and
/// nothing else can pass through it.
///
/// Without this, a regression that re-expressed arbitrary `g.` descriptions as
/// `m.` would be absorbed by the carve-out and the corpus censuses would stay
/// green.
#[test]
fn the_mitochondrial_carve_out_admits_only_mitochondrial_accessions() {
    for accession in ["NC_012920.1", "NC_001807.4"] {
        let variant = parse_hgvs(&format!("{accession}:g.7623C>T")).expect("parses");
        assert!(
            is_mitochondrial_relabel(
                &variant,
                CoordinateAxis::Genomic,
                CoordinateAxis::Mitochondrial
            ),
            "{accession} is a mitochondrial reference, so #487's relabel applies to it"
        );
    }

    let chromosomal = parse_hgvs("NC_000017.11:g.7623C>T").expect("parses");
    assert!(
        !is_mitochondrial_relabel(
            &chromosomal,
            CoordinateAxis::Genomic,
            CoordinateAxis::Mitochondrial
        ),
        "a chromosomal accession must not be admitted by the mitochondrial carve-out"
    );
    // …and the carve-out is about `g.` -> `m.` specifically, not about any pair
    // of axes on a mitochondrial accession.
    let mito = parse_hgvs("NC_012920.1:g.7623C>T").expect("parses");
    assert!(!is_mitochondrial_relabel(
        &mito,
        CoordinateAxis::Genomic,
        CoordinateAxis::Coding
    ));
}
