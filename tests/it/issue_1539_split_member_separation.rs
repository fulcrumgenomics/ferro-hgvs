//! Issue #1539: a member emitted by a split must itself obey `general.md:34`.
//!
//! Five real ClinVar/CMRG rows (three distinct variants; two pairs are one
//! variant on two transcripts) are lone `delins` descriptions that ferro
//! returns **unchanged and conformant**. A splitting pass has been observed to
//! cut each into an allele whose first member holds a reference base that is
//! unchanged in *every* minimal alignment and sits between two changed ones,
//! across a codon boundary — the shape `general.md:34` says to describe
//! individually and which `general.md:35`'s one-amino-acid exception cannot
//! rescue. Those parents were clean, so the members are regressions introduced
//! by the split rather than pre-existing defects it surfaced.
//!
//! # What this file asserts, and what it deliberately does not
//!
//! It asserts the **property**, never the string. Pinning today's five outputs
//! would make this a change detector, not a guard: the repo's adjudication
//! policy (`CLAUDE.md`, #1531) is explicit that "a test that merely pins today's
//! output is not an adjudication record". The authority here is
//! `assets/hgvs-nomenclature/docs/recommendations/general.md:34-35`, quoted in
//! [`adjudicate_member`], and the ruling is: whatever ferro chooses to emit for
//! these inputs, no member of it may carry a forced-unchanged reference column
//! between two changed ones outside the codon exception. A future change may
//! legitimately split any of these five; it may not split them into that.
//!
//! # Why it lives here rather than in the biocommons corpus
//!
//! `axis_normalized_hermetic` in `biocommons_normalize_tests.rs` is the other
//! per-PR hermetic gate and would otherwise be the natural home (#1539 suggests
//! it). It is the wrong home for these five: that corpus's `cases.json` is
//! regenerated verbatim from a pinned upstream SHA by
//! `scripts/refresh-biocommons-fixtures.py`, and its `reference-windows.json` is
//! regenerated from `cases.json`, so ClinVar/CMRG rows added there are erased by
//! the next refresh. This file reproduces that gate's shape instead — a
//! `WindowProvider` over a committed, generated `reference-windows.json` — so it
//! sees real reference bases, needs no `FERRO_MANIFEST`, and cannot skip.

use std::collections::BTreeSet;
use std::path::Path;
use std::sync::{Arc, Mutex};

use ferro_hgvs::conformance::reference_window::{WindowFixture, WindowProvider};
use ferro_hgvs::convert::mapper::CoordinateMapper;
use ferro_hgvs::hgvs::location::TxPos;
use ferro_hgvs::hgvs::variant::Accession;
use ferro_hgvs::reference::provider::GenomicPlacement;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::spdi::convert::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, FerroError, HgvsVariant, Normalizer, ReferenceProvider};

const INPUTS_TXT: &str = "tests/fixtures/split-member-separation/inputs.txt";
const WINDOWS_JSON: &str = "tests/fixtures/split-member-separation/reference-windows.json";
const REGENERATE: &str =
    "cargo run --features dev --example extract_split_member_separation_windows";

// ---------------------------------------------------------------------------
// The property: forced-unchanged reference columns
// ---------------------------------------------------------------------------

/// Cap on the alignment matrix, so a pathological member cannot hang the gate.
/// Every member these inputs produce is orders of magnitude below it; exceeding
/// it is reported as a deferral, which the gate then fails on rather than
/// quietly counting as clean.
const MAX_ALIGNMENT_CELLS: usize = 4_000_000;

/// Which reference columns of `reference -> payload` are unchanged in **every**
/// minimal alignment.
///
/// A column is *forced unchanged* when no minimum-cost edit path consumes it via
/// a deletion or a mismatch — i.e. every optimal alignment pairs it with an
/// equal payload base. Computed from the forward and backward Levenshtein
/// tables: `d[i][j]` is the cost of aligning `reference[..i]` to `payload[..j]`,
/// `e[i][j]` the cost of aligning `reference[i..]` to `payload[j..]`, and column
/// `i` can change iff some `j` completes a path of total cost `d[n][m]` through
/// a deletion of `reference[i]` or a substitution of it for a *different*
/// `payload[j]`.
///
/// Asking "in every minimal alignment" rather than "in some" is what keeps the
/// verdict sound: a column that survives one alignment but not another is left
/// unflagged. The check can therefore under-report and never over-report, which
/// is the correct bias for a gate that blocks merges.
///
/// Returns `None` when the matrix would exceed [`MAX_ALIGNMENT_CELLS`].
///
/// # Relationship to `common::minimal_alignment`
///
/// This is the closed form of the notion
/// `rulings[unchanged-is-read-over-every-minimal-alignment]` decides on, and it
/// is the detector that record names as where the ruling came from. It reaches
/// the answer in `Θ(n·m)` without ever materialising an alignment, which is
/// what makes it usable as a gate.
///
/// `common::minimal_alignment` reaches the same answer the other way — by
/// enumerating every minimal alignment and intersecting — and is
/// correspondingly exponential. The two are kept separate rather than one
/// rewritten in terms of the other precisely so that each can check the other;
/// `minimal_alignment_enumeration_proptest::the_closed_form_and_the_enumerator_agree`
/// asserts they do. Note this function is hardcoded to substitution cost 1,
/// which is the model that record's own worked example implies but never
/// states.
///
/// `pub(crate)` only so that cross-check can name it.
pub(crate) fn forced_unchanged_columns(reference: &[u8], payload: &[u8]) -> Option<Vec<bool>> {
    let (n, m) = (reference.len(), payload.len());
    if (n + 1).checked_mul(m + 1)? > MAX_ALIGNMENT_CELLS {
        return None;
    }

    // Forward: d[i][j] = edit distance between reference[..i] and payload[..j].
    let mut d = vec![vec![0u32; m + 1]; n + 1];
    for (j, cell) in d[0].iter_mut().enumerate() {
        *cell = j as u32;
    }
    for i in 1..=n {
        d[i][0] = i as u32;
        for j in 1..=m {
            let substitute = d[i - 1][j - 1] + u32::from(reference[i - 1] != payload[j - 1]);
            d[i][j] = substitute.min(d[i - 1][j] + 1).min(d[i][j - 1] + 1);
        }
    }

    // Backward: e[i][j] = edit distance between reference[i..] and payload[j..].
    let mut e = vec![vec![0u32; m + 1]; n + 1];
    for (j, cell) in e[n].iter_mut().enumerate() {
        *cell = (m - j) as u32;
    }
    for i in (0..n).rev() {
        e[i][m] = (n - i) as u32;
        for j in (0..m).rev() {
            let substitute = e[i + 1][j + 1] + u32::from(reference[i] != payload[j]);
            e[i][j] = substitute.min(e[i + 1][j] + 1).min(e[i][j + 1] + 1);
        }
    }

    let total = d[n][m];
    let mut forced = Vec::with_capacity(n);
    for i in 0..n {
        let mut can_change = false;
        for j in 0..=m {
            // reference[i] consumed as a deletion...
            let deleted = d[i][j] + 1 + e[i + 1][j] == total;
            // ...or as a mismatch against payload[j].
            let mismatched =
                j < m && reference[i] != payload[j] && d[i][j] + 1 + e[i + 1][j + 1] == total;
            if deleted || mismatched {
                can_change = true;
                break;
            }
        }
        forced.push(!can_change);
    }
    Some(forced)
}

/// A maximal run of forced-unchanged columns strictly inside the member.
/// `offset` is the index of the run's first column; `len` its width.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct InteriorRun {
    offset: usize,
    len: usize,
}

/// Maximal forced-unchanged runs that are strictly interior to the member.
///
/// A run touching either end is not a separation — it is an untrimmed flank, a
/// different defect. [`adjudicate_member`] declines to judge a member with one
/// rather than reasoning about it here.
fn interior_runs(forced: &[bool]) -> Vec<InteriorRun> {
    let mut runs = Vec::new();
    let mut i = 0;
    while i < forced.len() {
        if !forced[i] {
            i += 1;
            continue;
        }
        let start = i;
        while i < forced.len() && forced[i] {
            i += 1;
        }
        if start > 0 && i < forced.len() {
            runs.push(InteriorRun {
                offset: start,
                len: i - start,
            });
        }
    }
    runs
}

/// The verdict for one emitted member.
#[derive(Debug, Clone, PartialEq, Eq)]
enum MemberVerdict {
    /// The member names no reference column — a pure insertion or duplication.
    /// An insertion consumes no reference position, so its closed reference
    /// interval is empty (`lo = A + 1`, `hi = A`) and the separation rule has
    /// nothing inside it to bite on.
    NoColumns,
    /// Adjudicated and conformant.
    Clean,
    /// One entry per interior run that demands a further split.
    Violation(Vec<String>),
    /// Not adjudicated. Never counted as clean: a deferral is the guard going
    /// blind, not evidence of health, so the gate fails on one.
    Deferred(String),
}

/// Adjudicate one member against `general.md:34-35`.
///
/// > "two variants separated by one or more nucleotides should be described
/// > individually and **not** as a 'delins'."
/// >   — `assets/hgvs-nomenclature/docs/recommendations/general.md:34`
/// >
/// > "**exception**: two variants separated by one nucleotide, together
/// > affecting one amino acid, should be described as a 'delins'."
/// >   — `general.md:35`
///
/// `reference` and `payload` are the member's reference bases and its
/// replacement. `codon_of_column(i)` gives the 0-based codon index of reference
/// column `i`, or `None` when that column has no reading frame (a genomic or
/// non-coding axis, an intronic offset, a UTR position) — without a reading
/// frame there is no amino acid, so `:35` cannot apply and `:34` stands. That is
/// the conclusion `MIN_SEPARATION_NO_FRAME` reaches in `normalize::merge`, from
/// the same clause.
///
/// # Why the verdict is sound
///
/// A run of columns unchanged in every minimal alignment is not by itself a
/// violation: the rule needs *changed* material on both sides of it. That is
/// supplied by the endpoint-trim precondition. When the member is trimmed
/// (`reference[0] != payload[0]` and `reference[n-1] != payload[m-1]`) then in
/// *any* alignment column 0 is either changed or preceded by an insertion — a
/// matched first column must pair with some `payload[j]`, `j > 0`, leaving
/// `payload[..j]` inserted before it — and symmetrically at the far end. So a
/// strictly interior forced-unchanged run always separates two pieces of changed
/// material, in every minimal alignment. A member that is *not* trimmed is
/// deferred rather than judged, which under-reports and never over-reports.
///
/// The codon test is applied to **every** interior run, not just the first. An
/// earlier revision of this analysis tested `runs[0]` alone; where the first gap
/// was codon-internal and a later one was not, it scored the row a correct merge
/// when the correct form still splits — worth 54 missed violations per arm of
/// the measurement behind #1539.
///
/// For a one-column run the pair tested is the run's immediate neighbours. Those
/// are the *closest* two changed columns any alignment can offer, hence the most
/// likely to share a codon: if even they do not, no wider pair can, and `:35` is
/// unavailable under every reading. Testing the closest pair is what makes the
/// codon test conservative in the direction that matters.
fn adjudicate_member(
    reference: &[u8],
    payload: &[u8],
    codon_of_column: &dyn Fn(usize) -> Option<u64>,
) -> MemberVerdict {
    if reference.is_empty() {
        return MemberVerdict::NoColumns;
    }
    let Some(forced) = forced_unchanged_columns(reference, payload) else {
        return MemberVerdict::Deferred(format!(
            "alignment matrix {}x{} exceeds {MAX_ALIGNMENT_CELLS} cells",
            reference.len() + 1,
            payload.len() + 1
        ));
    };
    let runs = interior_runs(&forced);
    if runs.is_empty() {
        return MemberVerdict::Clean;
    }
    // The soundness argument above needs both flanks trimmed. A normalized
    // member always is; anything else is a different defect, not judged here.
    let untrimmed = !payload.is_empty()
        && (reference[0] == payload[0]
            || reference[reference.len() - 1] == payload[payload.len() - 1]);
    if untrimmed {
        return MemberVerdict::Deferred(
            "member is not endpoint-trimmed, so a forced-unchanged interior cannot be shown to \
             separate two changed pieces"
                .to_string(),
        );
    }

    let mut findings = Vec::new();
    for run in runs {
        // `:35` covers a separation of exactly one nucleotide. Two or more is
        // outside the exception under every reading.
        if run.len >= 2 {
            findings.push(format!(
                "{} forced-unchanged reference bases at offset {} separate two changed bases; \
                 general.md:35 covers a one-nucleotide separation only",
                run.len, run.offset
            ));
            continue;
        }
        let before = codon_of_column(run.offset - 1);
        let after = codon_of_column(run.offset + run.len);
        match (before, after) {
            (Some(lhs), Some(rhs)) if lhs == rhs => {
                // Both changes fall in one codon: general.md:35 applies.
            }
            (Some(lhs), Some(rhs)) => findings.push(format!(
                "one forced-unchanged reference base at offset {} separates changes in codons {} \
                 and {}; general.md:35's one-amino-acid exception does not reach across a codon \
                 boundary",
                run.offset, lhs, rhs
            )),
            _ => findings.push(format!(
                "one forced-unchanged reference base at offset {} separates two changed bases on \
                 an axis with no reading frame; general.md:35 needs one amino acid, so \
                 general.md:34 stands",
                run.offset
            )),
        }
    }
    if findings.is_empty() {
        MemberVerdict::Clean
    } else {
        MemberVerdict::Violation(findings)
    }
}

// ---------------------------------------------------------------------------
// Resolving an emitted member to reference columns
// ---------------------------------------------------------------------------

/// One adjudicated member: what it was, what it resolved to, and the verdict.
#[derive(Debug)]
struct MemberReport {
    rendered: String,
    reference: String,
    payload: String,
    verdict: MemberVerdict,
}

/// Flatten a normalized result into the members the rule applies to.
fn members_of(variant: &HgvsVariant) -> Vec<&HgvsVariant> {
    match variant {
        HgvsVariant::Allele(allele) => allele.variants.iter().flat_map(members_of).collect(),
        other => vec![other],
    }
}

/// Resolve `member` to its reference columns and adjudicate it.
///
/// The reference/payload pair comes from the member's SPDI triple, which is a
/// pure coordinate conversion — not a re-normalization — and already models
/// every edit type the normalizer emits: `inv` becomes its reverse complement,
/// `dup` and `ins` become empty-deletion insertions, `sub` a one-column delins.
fn adjudicate<P: ReferenceProvider>(member: &HgvsVariant, provider: &P) -> MemberReport {
    let rendered = member.to_string();
    let spdi = match hgvs_to_spdi(member, provider) {
        Ok(s) => s,
        Err(e) => {
            return MemberReport {
                rendered,
                reference: String::new(),
                payload: String::new(),
                verdict: MemberVerdict::Deferred(format!("cannot resolve reference bases: {e}")),
            }
        }
    };
    // SPDI is 0-based interbase, so the member's closed reference interval is
    // [position + 1, position + deletion.len()] — empty for an insertion, which
    // consumes no reference position.
    let first_base = spdi.position + 1;
    let reference = spdi.deletion.to_ascii_uppercase();
    let payload = spdi.insertion.to_ascii_uppercase();

    let codon_of_column = codon_resolver(member, first_base, provider);
    let verdict = adjudicate_member(
        reference.as_bytes(),
        payload.as_bytes(),
        codon_of_column.as_ref(),
    );
    MemberReport {
        rendered,
        reference,
        payload,
        verdict,
    }
}

/// A column -> codon-index resolver for `member`, whose reference column 0 is
/// 1-based position `first_base` on the member's own sequence.
///
/// Only `c.` and `r.` descriptions carry a reading frame; every other axis is
/// reported frameless. That is not a shortcut but the spec reading — `:35`'s
/// exception is about one amino acid, and a genomic or non-coding sequence has
/// none to speak of.
fn codon_resolver<P: ReferenceProvider>(
    member: &HgvsVariant,
    first_base: u64,
    provider: &P,
) -> Box<dyn Fn(usize) -> Option<u64>> {
    let accession = match member {
        HgvsVariant::Cds(v) => v.accession.transcript_accession(),
        HgvsVariant::Rna(v) => v.accession.transcript_accession(),
        _ => return Box::new(|_| None),
    };
    let Ok(transcript) = provider.get_transcript(&accession) else {
        return Box::new(|_| None);
    };
    Box::new(move |column: usize| {
        let mapper = CoordinateMapper::new(&transcript);
        let tx_pos = TxPos::new(i64::try_from(first_base + column as u64).ok()?);
        let cds = mapper.tx_to_cds(&tx_pos).ok()?;
        // A UTR, intronic or non-positive position has no codon.
        if cds.offset.is_some() || cds.utr3 || cds.base < 1 {
            return None;
        }
        Some((cds.base as u64 - 1) / 3)
    })
}

// ---------------------------------------------------------------------------
// Fixture plumbing
// ---------------------------------------------------------------------------

/// The committed input list: one HGVS string per line, `#` comments and blank
/// lines ignored. `examples/extract_split_member_separation_windows.rs` applies
/// the identical rule when capturing the windows.
fn inputs() -> Vec<String> {
    let content =
        std::fs::read_to_string(INPUTS_TXT).unwrap_or_else(|e| panic!("read {INPUTS_TXT}: {e}"));
    content
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_string)
        .collect()
}

fn window_fixture() -> WindowFixture {
    WindowFixture::from_json_path(Path::new(WINDOWS_JSON))
        .unwrap_or_else(|e| panic!("load {WINDOWS_JSON}: {e}"))
}

fn window_provider() -> WindowProvider {
    window_fixture().to_provider()
}

// ---------------------------------------------------------------------------
// Fidelity: the slice must answer every read the real reference would answer
// ---------------------------------------------------------------------------

/// A [`ReferenceProvider`] that forwards to the hermetic [`WindowProvider`] and
/// records every reference read it **fails**.
///
/// This exists because of a measured false negative, not a hypothetical one. A
/// window captured from a pass is exactly as wide as that pass happened to read,
/// so the fixture serves the recorded pass and errors on any pass that reads
/// further; normalization then bails and returns its input unchanged, which this
/// file's property check scores as perfectly conformant. The gate stayed green
/// with the defect fully present on three of its five rows.
///
/// The lesson generalises past the one fixture bug that caused it: an
/// under-serving slice is indistinguishable, at the output, from a clean
/// normalization. So the gate asserts not only that the members are conformant
/// but that **no reference read failed while producing them** — the hermetic
/// provider must answer everything the pass asks for, or say so loudly.
///
/// # What is audited, and why it is not simply "the `Result`-returning reads"
///
/// Audited: the reads that source reference **bases** from the captured slice —
/// `get_transcript`, `get_transcript_for_variant`, `get_transcript_for_accession`,
/// `get_sequence`, `get_genomic_sequence`. Those are exactly the accesses
/// `examples/common/recording.rs` records and the fixture stores, so a failure
/// among them is a missing window and regenerating is a remedy that works.
///
/// Everything else is delegated verbatim, including two `Result`-returning
/// methods, because auditing them would report a legitimate absence as a fixture
/// defect and send the reader to a regeneration that cannot fix it:
///
/// * `get_protein_sequence` — the fixture has no protein data to lose.
///   [`WindowFixture`] carries no protein field and the recorder captures no
///   protein reads, so [`WindowProvider`] inherits the trait default, which
///   returns `ProteinReferenceNotAvailable` for *every* accession; `has_protein_data`
///   is correspondingly `false`, which is the honest signal a caller checks.
///   Worse, `ReferenceProvider::get_protein_length`'s default implementation
///   binary-searches the length by probing `get_protein_sequence` and reading
///   `Err` as "past the end" — so a *successful* length query is built out of
///   ~17 deliberate failures. Auditing the read would turn intended control flow
///   into a wall of reported defects.
/// * `get_sequence_length` — it sources no bases, and the recorder does not
///   record length lookups, so a fixture can never be regenerated to answer one
///   it does not already cover. Every in-tree caller treats `Err` as "length
///   unknown" and proceeds (`.ok()`, `if let Ok(..)`), and any read that could
///   actually starve normalization goes through `get_sequence` or
///   `get_genomic_sequence`, which *are* audited. So the failure mode this
///   provider exists to catch stays covered without the false positives.
///
/// The `has_*` probes and the placement/selector lookups return no `Result` at
/// all and are delegated for the same reason: an existence check that answers
/// "no" is an answer, not a missing window.
#[derive(Clone)]
struct AuditedProvider {
    inner: WindowProvider,
    failures: Arc<Mutex<Vec<String>>>,
}

impl AuditedProvider {
    fn new(inner: WindowProvider) -> Self {
        Self {
            inner,
            failures: Arc::new(Mutex::new(Vec::new())),
        }
    }

    fn record<T>(&self, what: String, outcome: Result<T, FerroError>) -> Result<T, FerroError> {
        if let Err(e) = &outcome {
            self.failures
                .lock()
                .expect("audit lock")
                .push(format!("{what}: {e}"));
        }
        outcome
    }

    fn failures(&self) -> Vec<String> {
        self.failures.lock().expect("audit lock").clone()
    }
}

impl ReferenceProvider for AuditedProvider {
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
        let outcome = self.inner.get_transcript(id);
        self.record(format!("get_transcript({id})"), outcome)
    }

    fn get_transcript_for_variant(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Arc<Transcript>, FerroError> {
        let outcome = self.inner.get_transcript_for_variant(variant);
        self.record(format!("get_transcript_for_variant({variant})"), outcome)
    }

    fn get_transcript_for_accession(
        &self,
        accession: &Accession,
    ) -> Result<Arc<Transcript>, FerroError> {
        let outcome = self.inner.get_transcript_for_accession(accession);
        self.record(
            format!("get_transcript_for_accession({})", accession.full()),
            outcome,
        )
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        let outcome = self.inner.get_sequence(id, start, end);
        self.record(format!("get_sequence({id}, {start}, {end})"), outcome)
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        let outcome = self.inner.get_genomic_sequence(contig, start, end);
        self.record(
            format!("get_genomic_sequence({contig}, {start}, {end})"),
            outcome,
        )
    }

    // Delegated verbatim — the audit rule is on the type, above. The first two
    // return a `Result` and are still not audited: their `Err` is a legitimate
    // absence no regeneration could remove, not a missing window.
    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        self.inner.get_protein_sequence(accession, start, end)
    }

    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        self.inner.get_sequence_length(id)
    }

    fn has_transcript(&self, id: &str) -> bool {
        self.inner.has_transcript(id)
    }

    fn has_transcript_version_exact(&self, id: &str) -> bool {
        self.inner.has_transcript_version_exact(id)
    }

    fn has_genomic_data(&self) -> bool {
        self.inner.has_genomic_data()
    }

    fn has_protein_data(&self) -> bool {
        self.inner.has_protein_data()
    }

    fn genomic_placement(&self, parent: &Accession) -> Option<GenomicPlacement> {
        self.inner.genomic_placement(parent)
    }

    fn genomic_placement_on_build(
        &self,
        parent: &Accession,
        build: Option<&str>,
    ) -> Option<GenomicPlacement> {
        self.inner.genomic_placement_on_build(parent, build)
    }

    fn resolve_legacy_gene_selector(
        &self,
        selector: &str,
        ng_parent: Option<&Accession>,
    ) -> Option<String> {
        self.inner.resolve_legacy_gene_selector(selector, ng_parent)
    }

    fn sole_hosted_transcript(&self, ng_parent: &Accession) -> Option<String> {
        self.inner.sole_hosted_transcript(ng_parent)
    }
}

/// Adjudicate every member of `variant`.
fn adjudicate_all<P: ReferenceProvider>(variant: &HgvsVariant, provider: &P) -> Vec<MemberReport> {
    members_of(variant)
        .into_iter()
        .map(|member| adjudicate(member, provider))
        .collect()
}

/// Normalize `input` against the hermetic windows and adjudicate the result.
fn normalize_and_adjudicate<P: ReferenceProvider + Clone>(
    input: &str,
    provider: &P,
) -> (String, Vec<MemberReport>) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse `{input}`: {e}"));
    let normalized = Normalizer::new(provider.clone())
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalize `{input}`: {e}"));
    let rendered = normalized.to_string();
    (rendered, adjudicate_all(&normalized, provider))
}

/// Render one row's non-clean members for a failure message.
fn describe(input: &str, output: &str, reports: &[MemberReport]) -> String {
    let mut lines = vec![format!("  {input}\n    -> {output}")];
    for report in reports {
        match &report.verdict {
            MemberVerdict::Violation(findings) => {
                for finding in findings {
                    lines.push(format!(
                        "       VIOLATION {} (ref={} payload={}): {finding}",
                        report.rendered, report.reference, report.payload
                    ));
                }
            }
            MemberVerdict::Deferred(reason) => lines.push(format!(
                "       DEFERRED  {} (ref={} payload={}): {reason}",
                report.rendered, report.reference, report.payload
            )),
            MemberVerdict::Clean | MemberVerdict::NoColumns => {}
        }
    }
    lines.join("\n")
}

// ---------------------------------------------------------------------------
// The gate
// ---------------------------------------------------------------------------

/// The per-PR gate: no member ferro emits for these five inputs may hold a
/// forced-unchanged reference base between two changed ones outside
/// `general.md:35`'s codon exception.
#[test]
fn no_emitted_member_separates_two_changes_with_an_unchanged_base() {
    let provider = AuditedProvider::new(window_provider());
    let inputs = inputs();
    assert_eq!(
        inputs.len(),
        5,
        "{INPUTS_TXT} is the #1539 case set; it should hold the five rows the issue names"
    );

    let mut failures = Vec::new();
    let mut members_with_columns = 0usize;
    for input in &inputs {
        let (output, reports) = normalize_and_adjudicate(input, &provider);
        println!("{input}\n  -> {output}");
        for report in &reports {
            println!(
                "     {} ref={} payload={} {:?}",
                report.rendered, report.reference, report.payload, report.verdict
            );
        }
        members_with_columns += reports
            .iter()
            .filter(|report| !matches!(report.verdict, MemberVerdict::NoColumns))
            .count();
        if reports.iter().any(|report| {
            !matches!(
                report.verdict,
                MemberVerdict::Clean | MemberVerdict::NoColumns
            )
        }) {
            failures.push(describe(input, &output, &reports));
        }
    }

    // A zero is only meaningful if the members were actually examined. Each
    // input is a delins over a multi-base reference span, so every row must
    // contribute at least one member carrying reference columns; fewer means the
    // guard adjudicated nothing and passed vacuously.
    assert!(
        members_with_columns >= inputs.len(),
        "only {members_with_columns} of {} inputs produced a member with reference columns — the \
         guard cannot have examined what it claims to guard",
        inputs.len()
    );

    // ...and it is only meaningful if the slice answered every read. A read the
    // committed windows cannot serve makes normalization bail and return the
    // input unchanged, which is indistinguishable at the output from a clean
    // answer. This assertion is what separates the two, and it is checked before
    // the verdicts so a fixture defect is never reported as conformance.
    let access_failures = provider.failures();
    assert!(
        access_failures.is_empty(),
        "the hermetic slice failed {} reference read(s) during the pass, so any verdict above is \
         about a normalization that ran short of reference data, not about ferro's answer. \
         Regenerate with `{REGENERATE}` (widen the capture if the accession is genuinely new):\n  {}",
        access_failures.len(),
        access_failures.join("\n  ")
    );

    assert!(
        failures.is_empty(),
        "{} of {} #1539 inputs normalized to a form whose members violate general.md:34:\n{}",
        failures.len(),
        inputs.len(),
        failures.join("\n")
    );
}

/// The committed windows must carry every input's transcript. A fixture that
/// lost one would make the gate above panic rather than pass, but naming the
/// regeneration command here turns that into a one-line diagnosis.
#[test]
fn the_committed_windows_cover_every_input() {
    let captured: BTreeSet<String> = window_fixture()
        .transcripts
        .iter()
        .map(|transcript| transcript.id.clone())
        .collect();
    let missing: Vec<String> = inputs()
        .into_iter()
        .filter(|input| {
            let accession = input.split(':').next().unwrap_or_default();
            !captured.contains(accession)
        })
        .collect();
    assert!(
        missing.is_empty(),
        "{WINDOWS_JSON} is missing the transcript for {missing:?} — regenerate with `{REGENERATE}`"
    );
}

/// Every captured transcript must be served **whole**, through both entry
/// points, with bases identical to the stored record.
///
/// This is the hermetic half of the fidelity fix, and it encodes a real defect
/// rather than a tidiness preference. The recorder files the normalizer's
/// `get_genomic_sequence` reads under the transcript's own accession, which
/// registers it in the fixture's `genomic` map, which makes
/// `WindowProvider::is_known_contig` true for it — and that routes
/// `get_sequence` for the transcript through the genomic-window path instead of
/// the whole stored sequence. With the recorder's padded windows, measurement
/// against the manifest showed **every** 200 bp block of all five transcripts
/// failing under `WindowProvider` while succeeding under `MultiFastaProvider`.
///
/// The consequence was not a loud error but a silent pass: a normalization that
/// reads past the captured window bails and returns its input unchanged, and an
/// unchanged lone `delins` is exactly what a conformant answer looks like here.
/// Three of the five rows were scored Clean while carrying the defect in full.
///
/// `serve_transcripts_whole` in the generator widens each transcript's window to
/// the entire sequence so the failure mode cannot recur; this test is what stops
/// a future regeneration from quietly dropping that.
#[test]
fn the_committed_windows_serve_every_transcript_whole() {
    let fixture = window_fixture();
    let provider = fixture.to_provider();
    for transcript in &fixture.transcripts {
        let id = &transcript.id;
        let bases = transcript
            .sequence
            .as_ref()
            .unwrap_or_else(|| panic!("{WINDOWS_JSON}: transcript {id} carries no sequence"));
        let length = bases.len() as u64;
        assert_eq!(
            provider.get_sequence_length(id).ok(),
            Some(length),
            "{id}: get_sequence_length disagrees with the stored sequence"
        );
        for (entry_point, served) in [
            ("get_sequence", provider.get_sequence(id, 0, length)),
            (
                "get_genomic_sequence",
                provider.get_genomic_sequence(id, 0, length),
            ),
        ] {
            let served = served.unwrap_or_else(|e| {
                panic!(
                    "{id}: {entry_point}(0, {length}) failed against the committed windows: {e}\n\
                     The slice does not serve this accession whole, so a pass that reads further \
                     than the recorded one will bail and return its input unchanged — which this \
                     file's property check scores as conformant. Regenerate with `{REGENERATE}`."
                )
            });
            assert_eq!(
                &served, bases,
                "{id}: {entry_point} served different bases than the stored transcript"
            );
        }
    }
}

/// Ground truth for the fidelity property: the hermetic slice must make ferro
/// answer exactly what the full prepared reference makes it answer.
///
/// Manifest-or-skip, in the Layer 2 / Layer 3 relationship this repo already
/// uses (`axis_normalized` vs `axis_normalized_hermetic` in
/// `biocommons_normalize_tests.rs`, `tests/fixtures/CORPUS_LAYOUT.md`): the
/// hermetic gate is what blocks PRs, and this tier is what proves the hermetic
/// gate is measuring the same thing.
///
/// It compares two providers running the *same* code, so it stays correct
/// through any legitimate representation change — which a committed pin of
/// today's five outputs would not: such a pin fires on a normalizer change and
/// misreports it as a fixture defect, and this file's whole point is that it
/// does not pin outputs. The hermetic half of the same property is
/// [`the_committed_windows_serve_every_transcript_whole`] plus the
/// zero-failed-reads assertion in the gate, both of which run in CI and both of
/// which fail on the defect this test was written for.
#[test]
fn hermetic_normalization_matches_the_prepared_reference() {
    // The per-PR gate is
    // `no_emitted_member_separates_two_changes_with_an_unchanged_base`.
    let Some(manifest) = std::env::var_os("FERRO_MANIFEST").map(std::path::PathBuf::from) else {
        crate::common::manifest::absent(
            "issue_1539_split_member_separation::\
             hermetic_normalization_matches_the_prepared_reference",
        );
        return;
    };
    if !manifest.exists() {
        crate::common::manifest::absent(
            "issue_1539_split_member_separation::\
             hermetic_normalization_matches_the_prepared_reference \
             (FERRO_MANIFEST points at a path that does not exist)",
        );
        return;
    }
    let real = Arc::new(
        ferro_hgvs::MultiFastaProvider::from_manifest(&manifest)
            .unwrap_or_else(|e| panic!("from_manifest({}): {e}", manifest.display())),
    );
    let hermetic = window_provider();

    let mut divergences = Vec::new();
    for input in inputs() {
        let parsed = parse_hgvs(&input).unwrap_or_else(|e| panic!("parse `{input}`: {e}"));
        let from_windows = Normalizer::new(hermetic.clone())
            .normalize(&parsed)
            .map(|v| v.to_string());
        let from_manifest = Normalizer::new(Arc::clone(&real))
            .normalize(&parsed)
            .map(|v| v.to_string());
        match (&from_windows, &from_manifest) {
            (Ok(a), Ok(b)) if a == b => {}
            _ => divergences.push(format!(
                "  {input}\n    windows : {from_windows:?}\n    manifest: {from_manifest:?}"
            )),
        }
    }
    assert!(
        divergences.is_empty(),
        "{} of the #1539 inputs normalize differently under the committed windows than under the \
         prepared reference — the hermetic gate is not measuring ferro's real behaviour. \
         Regenerate with `{REGENERATE}`:\n{}",
        divergences.len(),
        divergences.join("\n")
    );
}

// ---------------------------------------------------------------------------
// The guard must be able to fail: the known-bad forms
// ---------------------------------------------------------------------------

/// The split forms observed for these inputs, recorded verbatim in #1539.
///
/// These are *not* an expectation — ferro does not produce them, and this file
/// does not ask it to. They exist so the gate above is demonstrably capable of
/// failing: a property check that cannot reject the known-bad output is
/// worthless, and "the corpus reports zero" means nothing until the corpus is
/// shown able to report non-zero.
///
/// The two sibling rows follow their partner's split at the same offsets, the
/// transcripts differing only by a constant CDS shift (0 for `NM_000089.4`, +3
/// for `NM_001009944.3`).
const KNOWN_BAD_SPLITS: &[&str] = &[
    "NM_000089.3:c.[2148_2152delinsAC;2155_2156inv]",
    "NM_000089.4:c.[2148_2152delinsAC;2155_2156inv]",
    "NM_001277115.2:c.[3545_3548delinsAGAT;3551delinsTGAACTGTAGCTGTAAGA]",
    "NM_000296.4:c.[10317_10330delinsAT;10334T>C]",
    "NM_001009944.3:c.[10320_10333delinsAT;10337T>C]",
];

#[test]
fn the_guard_rejects_the_known_bad_split_forms() {
    let provider = window_provider();
    let mut accepted = Vec::new();
    for form in KNOWN_BAD_SPLITS {
        let parsed = parse_hgvs(form).unwrap_or_else(|e| panic!("parse `{form}`: {e}"));
        let reports = adjudicate_all(&parsed, &provider);
        let deferred: Vec<&MemberReport> = reports
            .iter()
            .filter(|report| matches!(report.verdict, MemberVerdict::Deferred(_)))
            .collect();
        assert!(
            deferred.is_empty(),
            "`{form}` was not adjudicated — the demonstration proves nothing: {deferred:#?}"
        );
        let violations: Vec<&MemberReport> = reports
            .iter()
            .filter(|report| matches!(report.verdict, MemberVerdict::Violation(_)))
            .collect();
        if violations.is_empty() {
            accepted.push(describe(form, form, &reports));
        }
        for report in violations {
            println!(
                "rejected {} (ref={} payload={}): {:?}",
                report.rendered, report.reference, report.payload, report.verdict
            );
        }
    }
    assert!(
        accepted.is_empty(),
        "the separation check accepted {} of {} known-bad split forms from #1539 — it cannot \
         guard what it cannot reject:\n{}",
        accepted.len(),
        KNOWN_BAD_SPLITS.len(),
        accepted.join("\n")
    );
}

// ---------------------------------------------------------------------------
// Unit tests for the property itself
// ---------------------------------------------------------------------------

mod property {
    use super::*;

    /// A reading frame whose CDS coordinates start at `first_cds_base` for
    /// reference column 0.
    fn frame_from(first_cds_base: u64) -> impl Fn(usize) -> Option<u64> {
        move |column: usize| Some((first_cds_base + column as u64 - 1) / 3)
    }

    fn frameless(_: usize) -> Option<u64> {
        None
    }

    #[test]
    fn a_column_matched_by_every_minimal_alignment_is_forced() {
        // GGTCG -> AC costs 4, and every way of reaching 4 keeps the C: column 3
        // is forced unchanged and every other column can change. This is the
        // first member of #1539's worked example.
        let forced = forced_unchanged_columns(b"GGTCG", b"AC").expect("small matrix");
        assert_eq!(forced, vec![false, false, false, true, false]);
    }

    #[test]
    fn a_column_unchanged_in_only_some_alignments_is_not_forced() {
        // GACA -> AGAT aligns at cost 3 two ways: deleting the leading G keeps
        // columns 1 and 3, inserting a leading A keeps columns 0 and 1. Only
        // column 1 survives both, so only column 1 is forced — column 0 is
        // unchanged in one minimal alignment and must not be flagged.
        let forced = forced_unchanged_columns(b"GACA", b"AGAT").expect("small matrix");
        assert_eq!(forced, vec![false, true, false, false]);
    }

    #[test]
    fn runs_touching_either_end_are_not_interior() {
        assert_eq!(interior_runs(&[true, true, false]), vec![]);
        assert_eq!(interior_runs(&[false, true, true]), vec![]);
        assert_eq!(
            interior_runs(&[false, true, true, false]),
            vec![InteriorRun { offset: 1, len: 2 }]
        );
    }

    #[test]
    fn every_interior_run_is_tested_not_only_the_first() {
        // Two one-base gaps in one member: the first is codon-internal, the
        // second straddles a codon boundary. Judging `runs[0]` alone would score
        // this a correct merge — the 54-violations-per-arm error #1539 records.
        let forced = [false, true, false, false, false, true, false];
        assert_eq!(
            interior_runs(&forced),
            vec![
                InteriorRun { offset: 1, len: 1 },
                InteriorRun { offset: 5, len: 1 },
            ]
        );
        let frame = frame_from(1);
        assert_eq!(frame(0), frame(2), "the first gap is codon-internal");
        assert_ne!(frame(4), frame(6), "the second gap straddles a codon");
    }

    #[test]
    fn a_pure_insertion_names_no_reference_column() {
        // An insertion consumes no reference position: its closed interval is
        // empty, so there is nothing for the separation rule to bite on.
        assert_eq!(
            adjudicate_member(b"", b"ACGT", &frameless),
            MemberVerdict::NoColumns
        );
    }

    #[test]
    fn two_changes_in_one_codon_are_a_delins_by_the_exception() {
        // ACG -> TCA at CDS bases 1-3: column 1 is forced unchanged and the
        // flanking changes are both in codon 0, so general.md:35 applies.
        assert_eq!(
            forced_unchanged_columns(b"ACG", b"TCA").expect("small matrix"),
            vec![false, true, false]
        );
        assert_eq!(
            adjudicate_member(b"ACG", b"TCA", &frame_from(1)),
            MemberVerdict::Clean
        );
    }

    #[test]
    fn the_same_two_changes_across_a_codon_boundary_must_split() {
        // The identical block one base later, at CDS bases 2-4: the flanking
        // changes now fall in codons 0 and 1 and :35 cannot reach across.
        match adjudicate_member(b"ACG", b"TCA", &frame_from(2)) {
            MemberVerdict::Violation(findings) => {
                assert_eq!(findings.len(), 1);
                assert!(
                    findings[0].contains("codons 0 and 1"),
                    "unexpected finding: {}",
                    findings[0]
                );
            }
            other => panic!("expected a violation, got {other:?}"),
        }
    }

    #[test]
    fn the_codon_exception_is_unavailable_without_a_reading_frame() {
        match adjudicate_member(b"ACG", b"TCA", &frameless) {
            MemberVerdict::Violation(findings) => {
                assert_eq!(findings.len(), 1);
                assert!(
                    findings[0].contains("no reading frame"),
                    "unexpected finding: {}",
                    findings[0]
                );
            }
            other => panic!("expected a violation, got {other:?}"),
        }
    }

    #[test]
    fn a_two_base_separation_is_outside_the_exception_even_within_one_codon() {
        // ACGT -> TCGA keeps columns 1-2. Every column lies in codon 0, so the
        // one-amino-acid half of :35 is satisfied — and it still does not apply,
        // because :35 covers a separation of exactly one nucleotide.
        assert_eq!(
            forced_unchanged_columns(b"ACGT", b"TCGA").expect("small matrix"),
            vec![false, true, true, false]
        );
        match adjudicate_member(b"ACGT", b"TCGA", &frame_from(1)) {
            MemberVerdict::Violation(findings) => {
                assert_eq!(findings.len(), 1);
                assert!(
                    findings[0].contains("one-nucleotide separation only"),
                    "unexpected finding: {}",
                    findings[0]
                );
            }
            other => panic!("expected a violation, got {other:?}"),
        }
    }

    #[test]
    fn a_fully_changed_block_is_clean() {
        assert_eq!(
            adjudicate_member(b"AAAA", b"TTTT", &frame_from(1)),
            MemberVerdict::Clean
        );
    }

    #[test]
    fn an_untrimmed_member_is_deferred_not_judged() {
        // AGCGA -> ATCTA is two substitutions one base apart, so the rule does
        // apply — but both flanks match, which no normalized member's do, and
        // the soundness argument needs them trimmed. Declining is the documented
        // under-report; claiming would be an over-report on a member shape this
        // check has not reasoned about.
        assert!(matches!(
            adjudicate_member(b"AGCGA", b"ATCTA", &frame_from(2)),
            MemberVerdict::Deferred(_)
        ));
    }
}
