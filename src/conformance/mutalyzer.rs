//! Deserialization schema for the mutalyzer-normalize conformance corpus
//! (`tests/fixtures/mutalyzer-normalize/cases.json`).
//!
//! These types are the **single source of truth** for the corpus schema: both
//! the integration harness (`tests/mutalyzer_normalize_tests.rs`) and the
//! summary generator (`examples/generate_conformance_summary.rs`) deserialize
//! `cases.json` through them, so the generated `failure-patterns.md` cannot
//! drift from the schema the harness enforces. They live in the library (rather
//! than a test crate) only so an `examples/` binary can share them.

use serde::Deserialize;

use super::schema::{validate_cluster_refs, Cluster};
use super::summary::{DispositionKind, MemberRow, SummaryModel};

/// Corpus title used as the generated summary's heading and model title.
pub const CORPUS_TITLE: &str = "mutalyzer-normalize";

/// Top-level corpus document: metadata header plus the case list.
#[derive(Debug, Deserialize)]
#[allow(dead_code)]
pub struct Fixture {
    pub description: String,
    pub source: String,
    pub source_commit: String,
    pub license: String,
    pub refreshed_at: String,
    /// Root-cause cluster taxonomy referenced by the per-disposition `cluster`
    /// field. Optional so a corpus with no taxonomy still parses.
    #[serde(default)]
    pub clusters: Vec<Cluster>,
    pub cases: Vec<Case>,
}

impl Fixture {
    /// Every `(input, cluster_id)` pair referenced by a disposition `cluster`.
    pub fn cluster_refs(&self) -> Vec<(&str, &str)> {
        let mut refs = Vec::new();
        for case in &self.cases {
            let input = case.input.as_str();
            for cluster in [
                case.accepted_divergence
                    .as_ref()
                    .and_then(|d| d.cluster.as_deref()),
                case.known_bug.as_ref().and_then(|d| d.cluster.as_deref()),
                case.improvement.as_ref().and_then(|d| d.cluster.as_deref()),
                case.reference_unavailable
                    .as_ref()
                    .and_then(|d| d.cluster.as_deref()),
                case.spec_citation
                    .as_ref()
                    .and_then(|d| d.cluster.as_deref()),
            ]
            .into_iter()
            .flatten()
            {
                refs.push((input, cluster));
            }
        }
        refs
    }

    /// Validate that every disposition `cluster` ref resolves to a registry
    /// entry (orphan clusters are allowed). See [`validate_cluster_refs`].
    pub fn validate_clusters(&self) -> Result<(), String> {
        validate_cluster_refs(&self.clusters, self.cluster_refs())
    }

    /// Flatten this corpus into the corpus-agnostic summary model. mutalyzer
    /// dispositions carry no `ferro_output`, so that column is always absent.
    pub fn to_summary(&self) -> SummaryModel {
        let mut rows = Vec::new();
        for case in &self.cases {
            let input = case.input.as_str();
            if let Some(d) = &case.accepted_divergence {
                rows.push(MemberRow {
                    cluster: d.cluster.clone(),
                    input: input.to_string(),
                    axis: d.axis.as_str().to_string(),
                    kind: DispositionKind::AcceptedDivergence,
                    ferro_output: None,
                    tracking_issue: None,
                });
            }
            if let Some(d) = &case.known_bug {
                rows.push(MemberRow {
                    cluster: d.cluster.clone(),
                    input: input.to_string(),
                    axis: d.axis.as_str().to_string(),
                    kind: DispositionKind::KnownBug,
                    ferro_output: None,
                    tracking_issue: Some(d.tracking_issue),
                });
            }
            if let Some(d) = &case.improvement {
                rows.push(MemberRow {
                    cluster: d.cluster.clone(),
                    input: input.to_string(),
                    axis: d.axis.as_str().to_string(),
                    kind: DispositionKind::Improvement,
                    ferro_output: None,
                    tracking_issue: Some(d.tracking_issue),
                });
            }
            if let Some(d) = &case.reference_unavailable {
                rows.push(MemberRow {
                    cluster: d.cluster.clone(),
                    input: input.to_string(),
                    axis: d.axis.as_str().to_string(),
                    kind: DispositionKind::ReferenceUnavailable,
                    ferro_output: None,
                    tracking_issue: Some(d.tracking_issue),
                });
            }
            if let Some(d) = &case.spec_citation {
                rows.push(MemberRow {
                    cluster: d.cluster.clone(),
                    input: input.to_string(),
                    axis: d.axis.as_str().to_string(),
                    kind: DispositionKind::SpecCitation,
                    ferro_output: None,
                    tracking_issue: None,
                });
            }
        }
        SummaryModel {
            title: CORPUS_TITLE.to_string(),
            clusters: self.clusters.clone(),
            rows,
        }
    }
}

/// One corpus case: an input plus the expected output on each axis and any
/// per-axis divergence disposition.
#[derive(Debug, Deserialize, Clone)]
#[allow(dead_code)]
pub struct Case {
    #[serde(default)]
    pub keywords: Vec<String>,
    pub input: String,
    #[serde(default)]
    pub normalized: Option<String>,
    #[serde(default)]
    pub genomic: Option<String>,
    #[serde(default)]
    pub coding_protein_descriptions: Option<Vec<Vec<String>>>,
    #[serde(default)]
    pub protein_description: Option<String>,
    #[serde(default)]
    pub rna_description: Option<String>,
    #[serde(default)]
    pub errors: Option<Vec<String>>,
    #[serde(default)]
    pub infos: Option<Vec<String>>,
    #[serde(default)]
    pub noncoding: Option<Vec<String>>,
    #[serde(default = "default_true")]
    pub to_test: bool,
    /// Marks a mismatch on a specific axis as an accepted (non-bug)
    /// divergence — e.g. a ferro policy decision that intentionally
    /// disagrees with mutalyzer while both outputs remain HGVS-spec-allowed.
    /// Tallied into `AxisTally::divergence_accepted` instead of `fail`.
    #[serde(default)]
    pub accepted_divergence: Option<AcceptedDivergence>,
    /// Marks a mismatch on a specific axis as a known ferro bug (xfail):
    /// ferro is wrong and the divergence is tracked by an issue rather than
    /// accepted as policy. Tallied into `AxisTally::known_bug` instead of
    /// `fail`. If ferro starts matching mutalyzer (XPASS), the harness fails
    /// so the annotation and now-fixed row are cleaned up.
    #[serde(default)]
    pub known_bug: Option<KnownBug>,
    /// Marks a mismatch on a specific axis as a tracked *improvement* (not a
    /// bug): ferro's output is valid HGVS but not the spec-*preferred*
    /// canonical form, and `normalize()` should converge on the preferred form
    /// once `tracking_issue` lands. Tallied into `AxisTally::improvement`
    /// instead of `fail`. Distinct from `known_bug` (ferro is *wrong*) and from
    /// `accepted_divergence` (a deliberate, terminal policy where both forms
    /// are equally spec-allowed). Requires a `section` citation so every use is
    /// substantiated by the spec line that makes ferro's form non-preferred.
    /// XPASS-guarded like `known_bug`.
    #[serde(default)]
    pub improvement: Option<Improvement>,
    /// Marks a mismatch on a specific axis as a reference-availability artifact
    /// (the prepared manifest lacks the reference this row needs, so ferro
    /// no-ops), NOT a ferro bug. Tallied into `AxisTally::reference_unavailable`
    /// instead of `fail`. XPASS-guarded: if ferro starts matching, the reference
    /// became available and the annotation must be removed. Distinct from
    /// `known_bug` (ferro is *wrong* on a reference it *can* resolve).
    #[serde(default)]
    pub reference_unavailable: Option<ReferenceUnavailable>,
    /// Documentary annotation for cases where mutalyzer's expected output
    /// has been corrected in `cases.json` to ferro's spec-correct value.
    /// Has NO effect on tally bucketing — the corrected expected string
    /// makes the case PASS by the usual match check; this field only
    /// records the spec citation for review.
    #[serde(default)]
    pub spec_citation: Option<SpecCitation>,
}

/// Closed enum of corpus-runner axes. Replaces the previous free-form
/// `String` field on `AcceptedDivergence` and `SpecCitation` — serde
/// deserialization rejects unknown variants, so a typo in `cases.json`
/// surfaces as a hard parse error rather than silently no-op'ing the
/// annotation (#397 item 2). The variants match the 8 `AxisTally::new`
/// call sites; the lowercase serde rename keeps the JSON wire format
/// stable.
#[derive(Debug, Deserialize, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
#[allow(dead_code)]
pub enum Axis {
    Normalized,
    Genomic,
    Coding,
    ProteinDescription,
    CodingProteinDescriptions,
    RnaDescription,
    Noncoding,
    Errors,
    Infos,
}

impl Axis {
    /// Stable kebab-case string used in summary lines, file paths, and
    /// log messages. Matches the lowercase strings the previous
    /// `&'static str` axis carried.
    pub fn as_str(self) -> &'static str {
        match self {
            Axis::Normalized => "normalized",
            Axis::Genomic => "genomic",
            Axis::Coding => "coding",
            Axis::ProteinDescription => "protein_description",
            Axis::CodingProteinDescriptions => "coding_protein_descriptions",
            Axis::RnaDescription => "rna_description",
            Axis::Noncoding => "noncoding",
            Axis::Errors => "errors",
            Axis::Infos => "infos",
        }
    }
}

impl std::fmt::Display for Axis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Closed enum of accepted-divergence policies. Each policy identifier
/// names a specific reason a mismatch is non-failing (e.g. ferro's spec
/// arbitration on gene-symbol selectors). Adding a new policy requires
/// a code change here, which is intentional: it forces deliberate
/// review of every accepted divergence rather than allowing ad-hoc
/// string typing in `cases.json` (#397 item 2).
#[derive(Debug, Deserialize, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum Policy {
    /// ferro's gene-symbol-selector handling per issue #121 (the
    /// `NM_X(GENE):c.X` rendering and the related dispatcher
    /// arbitration). See PR #146 for the spec arbitration trail.
    #[serde(rename = "ferro-policy-121-gene-symbol-selector")]
    GeneSymbolSelector121,
    /// ferro emits a per-member SHUFFLE_APPLIED for a genuine 3'-shift in a
    /// compound allele; mutalyzer reports only the allele-level
    /// ISORTEDVARIANTS (it sorts/collapses members; ferro preserves order).
    /// Both are spec-defensible. Accepted divergence per #499 (precedent #481).
    #[serde(rename = "ferro-policy-499-shuffle-applied-compound-allele")]
    ShuffleAppliedCompoundAllele499,
    /// ferro's whole-CDS-deletion protein consequence convention. For a deletion
    /// spanning the entire coding sequence (including the initiation codon), the
    /// corpus emits `p.0?` (no protein produced, uncertain) while ferro emits
    /// `p.(Met1?)` (start-loss, uncertain). Both are spec-allowed predicted forms
    /// for a variant that removes the start codon
    /// (`docs/recommendations/protein/`): `p.0?` ("no protein") and `p.(Met1?)`
    /// ("unknown effect on translation initiation") describe the same
    /// undetermined outcome. ferro reports the start-loss form it derives from
    /// the normalized variant; this is a terminal convention difference, not a
    /// bug. (`p.0?` is arguably preferable as a future refinement.)
    #[serde(rename = "ferro-policy-whole-cds-del-met1")]
    WholeCdsDeletionMet1,
    /// ferro renders a deletion inside a genomic homopolymer (mononucleotide
    /// tandem tract) as the spec-valid repeat contraction `N[k]` (e.g.
    /// `g.320802_320804T[1]`) — `deletion_to_repeat` makes the repeat nature
    /// explicit; DNA `repeated.md` L5 allows a repeat unit of "one or more"
    /// nucleotides. mutalyzer keeps a plain `del` for a contraction (while both
    /// use repeat notation for an expansion). Both spec-defensible; ferro's
    /// repeat form is deliberate (`test_deletion_to_repeat_homopolymer_two_removed`).
    /// Accepted divergence per #745.
    #[serde(rename = "ferro-policy-745-homopolymer-repeat-contraction")]
    HomopolymerRepeatContraction745,
    /// #486 errors axis: ferro resolves an Ensembl transcript (`ENST…`)
    /// natively, so a well-formed `ENST…:c.` input normalizes instead of
    /// erroring. Mutalyzer's reference retrieval has no Ensembl backend and
    /// rejects with `ERETR`. ferro's broader reference support is an accepted
    /// divergence, not a defect (#486 lists Ensembl `ERETR` as a "not a bug").
    #[serde(rename = "ferro-policy-486-ensembl-transcript-supported")]
    EnsemblTranscriptSupported486,
    /// #486 errors axis: ferro rejects an HGVS-invalid input at *parse* time
    /// with a structured parse error — single-position insertion
    /// (`EINSERTIONRANGE`, also flanked by `EOVERLAP`), malformed concatenation,
    /// missing axis prefix, or an incomplete edit — whereas mutalyzer rejects
    /// the same input downstream with a semantic code (`EINTRONIC`, `ENOCDS`,
    /// `EREPEATUNSUPPORTED`, `EVARIANTNOTSUPPORTED`). Both reject the input;
    /// only the error taxonomy differs, so the rejection is accepted.
    #[serde(rename = "ferro-policy-486-parse-time-rejection-taxonomy")]
    ParseTimeRejectionTaxonomy486,
    /// #486 errors axis: ferro does not validate that a parenthesised transcript
    /// selector exists (or is version-pinned) in the reference catalog —
    /// selector resolution is a deferred scope boundary (#500-class) — so an
    /// unknown or unversioned selector normalizes rather than erroring.
    /// Mutalyzer rejects with `ENOSELECTORFOUND`. Accepted divergence pending
    /// dedicated selector-resolution work.
    #[serde(rename = "ferro-policy-486-selector-existence-not-validated")]
    SelectorExistenceNotValidated486,
    /// #486 errors axis: ferro is lenient about transcript-relative
    /// coordinate-system constraints that mutalyzer enforces — an `n.` position
    /// on a coding transcript (mutalyzer `ENOCDS`) and an intronic offset on an
    /// `r.` RNA position (mutalyzer `EINTRONICRNA`) both normalize rather than
    /// erroring. Accepted divergence pending stricter coordinate-system
    /// validation.
    #[serde(rename = "ferro-policy-486-transcript-coordinate-lenient")]
    TranscriptCoordinateLenient486,
}

impl Policy {
    /// Stable identifier used in summary lines so reviewers can grep
    /// across runs.
    pub fn as_str(self) -> &'static str {
        match self {
            Policy::GeneSymbolSelector121 => "ferro-policy-121-gene-symbol-selector",
            Policy::ShuffleAppliedCompoundAllele499 => {
                "ferro-policy-499-shuffle-applied-compound-allele"
            }
            Policy::WholeCdsDeletionMet1 => "ferro-policy-whole-cds-del-met1",
            Policy::HomopolymerRepeatContraction745 => {
                "ferro-policy-745-homopolymer-repeat-contraction"
            }
            Policy::EnsemblTranscriptSupported486 => {
                "ferro-policy-486-ensembl-transcript-supported"
            }
            Policy::ParseTimeRejectionTaxonomy486 => {
                "ferro-policy-486-parse-time-rejection-taxonomy"
            }
            Policy::SelectorExistenceNotValidated486 => {
                "ferro-policy-486-selector-existence-not-validated"
            }
            Policy::TranscriptCoordinateLenient486 => {
                "ferro-policy-486-transcript-coordinate-lenient"
            }
        }
    }
}

impl std::fmt::Display for Policy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Closed enum of reasons a corpus row's reference cannot be resolved by the
/// prepared manifest. Like [`Policy`] and [`SpecSection`], adding a variant
/// requires a code change here, which forces deliberate review of every claim
/// that a divergence is a reference-availability artifact rather than a ferro
/// defect. The `#[serde(rename)]` strings are the stable `cases.json` wire
/// format.
#[derive(Debug, Deserialize, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum ReferenceUnavailableReason {
    /// The exact pinned `accession.version` the row references is not carried by
    /// the prepared manifest (a vanilla `ferro prepare` ships only current
    /// RefSeq versions), so ferro cannot resolve the reference and no-ops
    /// (echoes the input). Provisioning the exact version resolves it.
    #[serde(rename = "accession-version-absent")]
    AccessionVersionAbsent,
    /// The row references an Ensembl accession (`ENST`/`ENSG`), which is absent
    /// from the RefSeq-only prepared manifest, so the reference cannot be
    /// resolved and ferro no-ops. Preparing an Ensembl reference resolves it.
    #[serde(rename = "ensembl-absent-from-refseq-manifest")]
    EnsemblAbsentFromRefseqManifest,
}

impl ReferenceUnavailableReason {
    /// Stable identifier used in summary lines so reviewers can grep across
    /// runs. Matches the `#[serde(rename)]` strings byte-for-byte.
    pub fn as_str(self) -> &'static str {
        match self {
            ReferenceUnavailableReason::AccessionVersionAbsent => "accession-version-absent",
            ReferenceUnavailableReason::EnsemblAbsentFromRefseqManifest => {
                "ensembl-absent-from-refseq-manifest"
            }
        }
    }
}

impl std::fmt::Display for ReferenceUnavailableReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Accepted-divergence annotation. When attached to a `Case`, the
/// corpus runner treats a mismatch on `axis` as a tracked-but-non-failing
/// divergence (counted in `divergence_accepted`) rather than a hard FAIL.
///
/// `axis` and `policy` are now closed enums — `serde` rejects unknown
/// variants at fixture-parse time, so misspellings surface as a hard
/// parse error rather than silently no-op'ing the annotation.
#[derive(Debug, Deserialize, Clone)]
#[allow(dead_code)]
pub struct AcceptedDivergence {
    /// Axis this annotation applies to.
    pub axis: Axis,
    /// Closed policy identifier explaining *why* the divergence is
    /// acceptable. Surfaced in the per-axis tally summary so reviewers
    /// can grep for the policy.
    pub policy: Policy,
    /// Optional human-readable note expanding on the policy.
    #[serde(default)]
    pub note: Option<String>,
    /// Root-cause cluster id (see the corpus `clusters` registry), grouping this
    /// divergence under a cross-case pattern in the generated summary.
    #[serde(default)]
    pub cluster: Option<String>,
}

/// Known-bug annotation. When attached to a `Case`, the corpus runner
/// treats a mismatch on `axis` as an expected failure (xfail): ferro is
/// *wrong* on this axis and the divergence is a tracked bug rather than a
/// policy decision. Tallied into `AxisTally::known_bug` instead of `fail`.
///
/// Contrast with [`AcceptedDivergence`], which marks an *intentional*,
/// spec-allowed divergence. A `known_bug` says ferro should eventually
/// match mutalyzer once `tracking_issue` is fixed.
///
/// If ferro *starts* matching mutalyzer on this axis (an XPASS), the
/// harness FAILs loudly: the bug appears fixed, so the annotation and the
/// stale fixture row must be cleaned up (the row demoted to a normal
/// passing case).
#[derive(Debug, Deserialize, Clone)]
#[allow(dead_code)]
pub struct KnownBug {
    /// Axis name this annotation applies to (e.g. `"normalized"`).
    pub axis: Axis,
    /// Issue tracking the fix.
    pub tracking_issue: u64,
    /// Optional human-readable note.
    #[serde(default)]
    pub note: Option<String>,
    /// Root-cause cluster id (see the corpus `clusters` registry), grouping this
    /// divergence under a cross-case pattern in the generated summary.
    #[serde(default)]
    pub cluster: Option<String>,
}

/// Reference-unavailable annotation. When attached to a `Case`, the corpus
/// runner treats a mismatch on `axis` as a **reference-availability artifact**,
/// NOT a ferro defect: the prepared manifest lacks the reference this row needs,
/// so ferro no-ops (echoes the input) and cannot match mutalyzer. Tallied into
/// `AxisTally::reference_unavailable`, distinct from `known_bug`.
///
/// This disposition exists to stop conflating "the oracle has no reference for
/// this input" with "ferro produced the wrong answer." The former is a
/// provisioning gap (tracked by `tracking_issue` — the issue that will make the
/// reference available, e.g. exact-version provisioning or an Ensembl
/// reference), not a normalization bug. Contrast with [`KnownBug`], which
/// asserts ferro is *wrong* on a reference it *can* resolve.
///
/// XPASS-guarded like [`KnownBug`]: if ferro *starts* matching mutalyzer on this
/// axis, the reference became available (the provisioning landed), so the
/// harness FAILs loudly to demote the row and remove the annotation.
#[derive(Debug, Deserialize, Clone)]
#[allow(dead_code)]
pub struct ReferenceUnavailable {
    /// Axis this annotation applies to.
    pub axis: Axis,
    /// Closed reason identifier explaining why the reference is unavailable.
    /// Surfaced in the per-axis tally summary so reviewers can grep for it.
    pub reason: ReferenceUnavailableReason,
    /// Issue tracking the provisioning that will make the reference available.
    pub tracking_issue: u64,
    /// Optional human-readable note expanding on the reason.
    #[serde(default)]
    pub note: Option<String>,
    /// Root-cause cluster id (see the corpus `clusters` registry), grouping this
    /// divergence under a cross-case pattern in the generated summary.
    #[serde(default)]
    pub cluster: Option<String>,
}

/// Tracked-improvement annotation. When attached to a `Case`, the corpus
/// runner treats a mismatch on `axis` as a tracked *enhancement* (xfail):
/// ferro's output is valid HGVS but the spec strongly *prefers* the form
/// mutalyzer emits, and `normalize()` should converge on that preferred form
/// once `tracking_issue` is implemented. Tallied into `AxisTally::improvement`
/// instead of `fail`.
///
/// This is the middle disposition between [`AcceptedDivergence`] and
/// [`KnownBug`]:
/// - `accepted_divergence` — terminal policy; both forms are equally
///   spec-allowed and ferro will *not* change.
/// - `improvement` — ferro's output is *valid but spec-non-preferred*; ferro
///   *will* converge (tracked to `tracking_issue`).
/// - `known_bug` — ferro's output is *wrong / invalid*; ferro *will* be fixed.
///
/// Unlike `known_bug`, an `improvement` also requires a `section`: the spec
/// citation that establishes ferro's current form as non-preferred. This keeps
/// the disposition substantiated rather than a catch-all for deferred work.
///
/// XPASS-guarded like `known_bug`: if ferro *starts* matching mutalyzer on
/// this axis the harness FAILs loudly so the stale annotation and the
/// now-converged fixture row are cleaned up.
#[derive(Debug, Deserialize, Clone)]
#[allow(dead_code)]
pub struct Improvement {
    /// Axis this annotation applies to.
    pub axis: Axis,
    /// Issue tracking the convergence to the spec-preferred form.
    pub tracking_issue: u64,
    /// HGVS spec section establishing ferro's current form as non-preferred.
    pub section: SpecSection,
    /// Optional human-readable note.
    #[serde(default)]
    pub note: Option<String>,
    /// Root-cause cluster id (see the corpus `clusters` registry), grouping this
    /// divergence under a cross-case pattern in the generated summary.
    #[serde(default)]
    pub cluster: Option<String>,
}

/// Closed enum of HGVS spec section identifiers cited from `cases.json`.
/// Mirrors [`Policy`]'s shape: adding a new section requires a code
/// change here, which forces deliberate review of every new citation
/// rather than allowing ad-hoc string typing in `cases.json` (#430,
/// follow-up to #397 item 2). The `#[serde(rename = "...")]` annotation
/// keeps the JSON wire format stable.
#[derive(Debug, Deserialize, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum SpecSection {
    /// HGVS §Prioritization — the spec's edit-prioritization rules (dup
    /// over ins, etc.) that ferro applies during normalization. This is
    /// a ferro-internal label, not a stable spec heading: the upstream
    /// spec submodule (`assets/hgvs-nomenclature/docs/recommendations/`)
    /// has no `Prioritization` anchor today, so binding to a path here
    /// would create a brittle pointer. The human label is the stable
    /// catalog key; the rename string preserves the byte sequence in
    /// `cases.json`.
    #[serde(rename = "HGVS §Prioritization")]
    Prioritization,
    /// HGVS protein reference format — a `p.` description's reference is the
    /// protein sequence accession alone (e.g. `NP_000068.1:p.(…)`); every
    /// example in `docs/recommendations/protein/*` uses this bare form. The
    /// `genomic_context(protein)` wrapper mutalyzer emits
    /// (`NG_007485.1(NP_000068.1):p.…`) has no counterpart in the spec, so
    /// ferro's bare-`NP_` output is the spec-correct value.
    #[serde(rename = "HGVS protein reference (bare NP)")]
    ProteinReference,
    /// HGVS RefSeqGene transcript selection — on a genomic reference
    /// (`NG_`/`NC_`/`LRG_`) the parenthetical selector for a `c.`/`n.`
    /// coordinate must name a transcript (`NM_`) or protein (`NP_`) accession,
    /// not a gene symbol (`assets/hgvs-nomenclature/docs/background/refseq.md`
    /// L38-42: "Gene symbols should not be used as specification"; L138-139:
    /// the transcript used should be indicated). ferro preserves the input's
    /// gene-symbol selector verbatim — valid HGVS, but the transcript-accession
    /// form mutalyzer emits is spec-preferred. (The `NM_(GENE)` policy from #121
    /// is unaffected; this section concerns only the genomic-reference selector
    /// slot, which #121 never examined.)
    #[serde(rename = "HGVS §RefSeqGene transcript selection")]
    RefSeqGeneSelector,
    /// HGVS substitution consequence — a single-base coding substitution
    /// preserves length and reading frame, so its predicted protein
    /// consequence is a missense / nonsense / synonymous / stop-loss extension
    /// / initiation-codon change, **never** a frameshift
    /// (`docs/recommendations/protein/substitution.md:5`; `frameshift.md:18`,
    /// where a frameshift requires an indel whose length is not a multiple of
    /// 3). A `p.…fs…` consequence reported for a substitution can only arise
    /// when the transcript reference is reconstructed from a divergent
    /// genomic embedding; it is invalid against a fixed reference. ferro emits
    /// the in-frame missense against the canonical standalone RefSeq transcript.
    #[serde(rename = "HGVS §Substitution (no frameshift)")]
    SubstitutionConsequence,
    /// HGVS protein initiation codon — a variant affecting the translation
    /// initiation codon (CDS 1–3) has an unpredictable protein consequence,
    /// reported as `p.(Met1?)` (`docs/recommendations/protein/substitution.md:51`,
    /// `deletion.md:62`). ferro derives the consequence from the canonical
    /// normalized variant; when an input's normalized form spans into the
    /// initiation codon (e.g. `c.-1_2dup`), ferro reports `p.(Met1?)`. mutalyzer
    /// computes the consequence from the un-normalized input and so can report
    /// inconsistent results for inputs that share a normalized form, so ferro's
    /// `p.(Met1?)` is the spec-correct value (#504, #512).
    #[serde(rename = "HGVS protein initiation codon (Met1?)")]
    ProteinInitiationCodon,
    /// HGVS transcript flanking — on a genomic-reference `c.` description
    /// (`NG_`/`NC_`/`LRG_`), `pter`/`qter` denote the genomic parent's terminus,
    /// a 5'/3' transcript-flank position. HGVS does **not** permit numbering
    /// flank nucleotides in `c.` coordinates (`background/numbering.md`
    /// transcript-flanking; `background/refseq.md` L46; the flank-numbering
    /// proposal was rejected, `consultation/open-issues.md` — "Use `NC_…:g.…`
    /// **and not** `NC_…(NM_…):c.-N-uM`"). ferro leaves the input verbatim and
    /// refuses the flank `c.` coordinate (#488 Phase 2b); mutalyzer's
    /// `c.-5059del` extrapolates into the flank, which HGVS forbids.
    #[serde(rename = "HGVS §Transcript flanking (not c.-numberable)")]
    TranscriptFlankNotNumberable,
    /// HGVS amino-acid glyph preference — three-letter codes (including `Ter`)
    /// are the preferred canonical form; the single-character `*` is an
    /// accepted alternative only (`assets/hgvs-nomenclature/docs/background/standards.md`).
    /// ferro canonicalizes a C-terminal stop-loss extension to the three-letter
    /// `extTer<n>` form (#224); hgvs-rs carries the legacy `ext*<n>` glyph. Both
    /// are valid HGVS, but ferro's `extTer` is spec-preferred, so the divergence
    /// is a tracked `improvement` rather than a bug.
    #[serde(rename = "HGVS §Standards (three-letter ext Ter preferred over ext*)")]
    ExtensionTerGlyph,
    /// HGVS repeated sequences — coding-DNA codon exception. On a coding DNA
    /// reference (`c.`), repeat notation `unit[N]` may be used **only** for
    /// repeat units whose length is a multiple of 3 (i.e. that cannot affect
    /// the reading frame); any other expansion/contraction must be described
    /// as a `dup`, `ins`, or `del`
    /// (`docs/recommendations/DNA/repeated.md` L21-23: use
    /// `NM_024312.4:c.2692_2693dup` **not** `c.2686A[10]`, and
    /// `NM_024312.4:c.1741_1742insTATATATA` **not** `c.1738TA[6]`). The
    /// restriction applies only to the coding sequence, not introns or UTR.
    /// For a non-codon-aligned unit in coding `c.`, ferro emits the
    /// spec-mandated `dup`/`ins`/`del`; mutalyzer's `unit[N]` form is invalid
    /// HGVS, so ferro's output is the spec-correct value (#487).
    #[serde(rename = "HGVS §Repeated (coding codon exception)")]
    RepeatCodingCodonException,
}

impl SpecSection {
    /// Stable identifier used in summary lines so reviewers can grep
    /// across runs. Matches the `#[serde(rename)]` strings byte-for-byte
    /// to preserve the pre-enum summary-line format.
    pub fn as_str(self) -> &'static str {
        match self {
            SpecSection::Prioritization => "HGVS §Prioritization",
            SpecSection::ProteinReference => "HGVS protein reference (bare NP)",
            SpecSection::RefSeqGeneSelector => "HGVS §RefSeqGene transcript selection",
            SpecSection::SubstitutionConsequence => "HGVS §Substitution (no frameshift)",
            SpecSection::ProteinInitiationCodon => "HGVS protein initiation codon (Met1?)",
            SpecSection::TranscriptFlankNotNumberable => {
                "HGVS §Transcript flanking (not c.-numberable)"
            }
            SpecSection::ExtensionTerGlyph => {
                "HGVS §Standards (three-letter ext Ter preferred over ext*)"
            }
            SpecSection::RepeatCodingCodonException => "HGVS §Repeated (coding codon exception)",
        }
    }
}

impl std::fmt::Display for SpecSection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Spec-citation annotation. Documentary only — does NOT affect tally
/// bucketing. Used to record that the `cases.json` expected output for
/// `axis` has been corrected to ferro's spec-correct value (versus
/// mutalyzer's spec-incorrect output), with a citation to the relevant
/// HGVS spec section.
///
/// `axis` and `section` are now closed enums (see [`Axis`], [`SpecSection`]);
/// serde rejects unknown variants at fixture-parse time, so a misspelled
/// section identifier surfaces as a hard parse error rather than silently
/// losing cataloguability.
#[derive(Debug, Deserialize, Clone)]
#[allow(dead_code)]
pub struct SpecCitation {
    /// Axis this citation applies to.
    pub axis: Axis,
    /// HGVS spec section the citation references (e.g.
    /// `SpecSection::Prioritization`).
    pub section: SpecSection,
    /// Optional human-readable note expanding on the citation.
    #[serde(default)]
    pub note: Option<String>,
    /// Root-cause cluster id (see the corpus `clusters` registry), grouping this
    /// divergence under a cross-case pattern in the generated summary.
    #[serde(default)]
    pub cluster: Option<String>,
}

fn default_true() -> bool {
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    const BASE: &str =
        r#""description":"t","source":"t","source_commit":"t","license":"t","refreshed_at":"t""#;

    fn parse(clusters: &str, cases: &str) -> Fixture {
        let json = format!("{{{BASE},\"clusters\":[{clusters}],\"cases\":[{cases}]}}");
        serde_json::from_str(&json).expect("fixture should deserialize")
    }

    #[test]
    fn cluster_refs_collects_every_disposition_kind() {
        let clusters = r#"
            {"id":"sel","title":"RefSeqGene selector","spec_section":"background/refseq.md"},
            {"id":"np","title":"bare NP","spec_section":"protein"}
        "#;
        let cases = r#"
            {"input":"A","improvement":{"axis":"normalized","tracking_issue":500,
              "section":"HGVS §RefSeqGene transcript selection","cluster":"sel"}},
            {"input":"B","spec_citation":{"axis":"protein_description",
              "section":"HGVS protein reference (bare NP)","cluster":"np"}}
        "#;
        let fixture = parse(clusters, cases);
        let mut refs = fixture.cluster_refs();
        refs.sort();
        assert_eq!(refs, vec![("A", "sel"), ("B", "np")]);
        assert!(fixture.validate_clusters().is_ok());

        let summary = fixture.to_summary();
        assert_eq!(summary.title, "mutalyzer-normalize");
        assert_eq!(summary.rows.len(), 2);
        let imp = summary
            .rows
            .iter()
            .find(|r| r.kind == DispositionKind::Improvement)
            .expect("improvement row");
        assert_eq!(imp.tracking_issue, Some(500));
        assert_eq!(imp.axis, "normalized");
        assert_eq!(imp.cluster.as_deref(), Some("sel"));
        assert!(summary
            .rows
            .iter()
            .any(|r| r.kind == DispositionKind::SpecCitation));
    }

    #[test]
    fn dangling_disposition_cluster_ref_is_rejected() {
        let clusters = r#"{"id":"np","title":"bare NP","spec_section":"protein"}"#;
        let cases = r#"
            {"input":"B","spec_citation":{"axis":"protein_description",
              "section":"HGVS protein reference (bare NP)","cluster":"missing"}}
        "#;
        let err = parse(clusters, cases)
            .validate_clusters()
            .expect_err("a dangling disposition cluster ref must be rejected");
        assert!(err.contains("missing"), "{err}");
    }

    #[test]
    fn dispositions_without_clusters_are_ok() {
        let cases = r#"{"input":"A","improvement":{"axis":"normalized",
            "tracking_issue":500,"section":"HGVS §RefSeqGene transcript selection"}}"#;
        let fixture = parse("", cases);
        assert!(fixture.cluster_refs().is_empty());
        assert!(fixture.validate_clusters().is_ok());
    }
}
