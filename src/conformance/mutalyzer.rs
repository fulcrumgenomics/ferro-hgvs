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
    /// Comparator runtime provenance (#890): the exact comparator versions the
    /// corpus expected values are validated against, plus the ferro
    /// `reference_identity()` hash (#764). Optional so a corpus without it still
    /// parses; recorded so "which mutalyzer/biocommons/hgvs-rs" is unambiguous
    /// (the root-cause provenance gap #890 identified).
    #[serde(default)]
    pub comparator_provenance: Option<ComparatorProvenance>,
    pub cases: Vec<Case>,
}

/// Comparator runtime versions the corpus expected values are validated against
/// (#890). The `mutalyzer-normalize` corpus is *derived from* the mutalyzer
/// normalizer, so `mutalyzer` / `mutalyzer_hgvs_parser` are its authoritative
/// source; the biocommons/hgvs-rs versions are the benchmark suite's other
/// comparators, recorded here for one canonical provenance record. Prefer
/// **deriving** expected values from these pinned comparators (via
/// `ferro-benchmark compare-normalize`) over copying an upstream test suite —
/// the practice that prevents the Mutalyzer2-era drift #882 surfaced.
#[derive(Debug, Deserialize, Default)]
#[allow(dead_code)]
pub struct ComparatorProvenance {
    /// Free-text provenance note (validation scope, derivation practice).
    pub note: String,
    /// Date the versions below were confirmed against the corpus (ISO-8601).
    pub validated_at: String,
    /// ferro prepared-reference identity hash (#764) the corpus was blessed
    /// against; must match `reference_identity_from_manifest`.
    pub reference_identity: String,
    /// mutalyzer normalizer version — the source of this corpus's values.
    pub mutalyzer: String,
    /// mutalyzer-hgvs-parser version (mutalyzer's HGVS parser dependency).
    pub mutalyzer_hgvs_parser: String,
    /// biocommons/hgvs version (cross-check comparator; empty if not recorded).
    #[serde(default)]
    pub biocommons_hgvs: String,
    /// biocommons.seqrepo version (cross-check comparator).
    #[serde(default)]
    pub biocommons_seqrepo: String,
    /// hgvs-rs version/tag (cross-check comparator; projection-fixture source).
    #[serde(default)]
    pub hgvs_rs: String,
}

impl Fixture {
    /// Every `(input, cluster_id)` pair referenced by a disposition `cluster`.
    pub fn cluster_refs(&self) -> Vec<(&str, &str)> {
        let mut refs = Vec::new();
        for case in &self.cases {
            let input = case.input.as_str();
            for cluster in [
                case.improvement.as_ref().and_then(|d| d.cluster.as_deref()),
                case.reference_unavailable
                    .as_ref()
                    .and_then(|d| d.cluster.as_deref()),
            ]
            .into_iter()
            .flatten()
            {
                refs.push((input, cluster));
            }
            // A `Case` may carry multiple per-axis accepted divergences (#870),
            // known bugs (#870), accepted rejections (#870), and spec citations
            // (#827); collect the cluster ref of each.
            for cluster in case
                .accepted_divergences
                .iter()
                .filter_map(|d| d.cluster.as_deref())
                .chain(case.known_bugs.iter().filter_map(|d| d.cluster.as_deref()))
                .chain(
                    case.accepted_rejections
                        .iter()
                        .filter_map(|d| d.cluster.as_deref()),
                )
                .chain(
                    case.spec_citations
                        .iter()
                        .filter_map(|d| d.cluster.as_deref()),
                )
            {
                refs.push((input, cluster));
            }
        }
        refs
    }

    /// Validate that every disposition `cluster` ref resolves to a registry
    /// entry (orphan clusters are allowed, see [`validate_cluster_refs`]) AND
    /// that no `Case` carries more than one `spec_citation` for the same axis.
    /// The matcher keys at most one citation per axis (it `find`s the first
    /// matching one), so a duplicate-axis citation would be silently ignored;
    /// surfacing it here turns an authoring mistake into a hard error at
    /// fixture-validation time (#827).
    pub fn validate_clusters(&self) -> Result<(), String> {
        validate_cluster_refs(&self.clusters, self.cluster_refs())?;
        self.validate_multi_axis_annotation_uniqueness()
    }

    /// Reject any `Case` whose multi-axis annotations (`spec_citations`,
    /// `accepted_divergences`, `known_bugs`, `accepted_rejections`) name the same axis more than once
    /// within a single annotation kind — the matcher `find`s the first match per
    /// axis, so a duplicate-axis entry would be silently ignored. Surfacing it
    /// here turns an authoring mistake into a hard error at fixture-validation
    /// time (#827/#870). See [`validate_clusters`](Self::validate_clusters).
    fn validate_multi_axis_annotation_uniqueness(&self) -> Result<(), String> {
        fn check_unique(
            input: &str,
            kind: &str,
            axes: impl IntoIterator<Item = Axis>,
        ) -> Result<(), String> {
            let mut seen: Vec<Axis> = Vec::new();
            for axis in axes {
                if seen.contains(&axis) {
                    return Err(format!(
                        "case {input:?} has more than one {kind} for axis {:?}; \
                         the matcher honors only one per axis",
                        axis.as_str(),
                    ));
                }
                seen.push(axis);
            }
            Ok(())
        }

        for case in &self.cases {
            check_unique(
                &case.input,
                "spec_citation",
                case.spec_citations.iter().map(|c| c.axis),
            )?;
            check_unique(
                &case.input,
                "accepted_divergence",
                case.accepted_divergences.iter().map(|d| d.axis),
            )?;
            check_unique(
                &case.input,
                "known_bug",
                case.known_bugs.iter().map(|d| d.axis),
            )?;
            check_unique(
                &case.input,
                "accepted_rejection",
                case.accepted_rejections.iter().map(|d| d.axis),
            )?;
        }
        Ok(())
    }

    /// Flatten this corpus into the corpus-agnostic summary model. mutalyzer
    /// dispositions carry no `ferro_output`, so that column is always absent.
    pub fn to_summary(&self) -> SummaryModel {
        let mut rows = Vec::new();
        for case in &self.cases {
            let input = case.input.as_str();
            // One summary row per per-axis accepted divergence (#870): a
            // multi-axis row contributes an `AcceptedDivergence` row per axis.
            for d in &case.accepted_divergences {
                rows.push(MemberRow {
                    cluster: d.cluster.clone(),
                    input: input.to_string(),
                    axis: d.axis.as_str().to_string(),
                    kind: DispositionKind::AcceptedDivergence,
                    ferro_output: None,
                    tracking_issue: None,
                });
            }
            // One summary row per per-axis known bug (#870).
            for d in &case.known_bugs {
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
            // One summary row per per-axis spec citation (#827): a multi-axis
            // row contributes a `SpecCitation` row for each axis it cites.
            for d in &case.spec_citations {
                rows.push(MemberRow {
                    cluster: d.cluster.clone(),
                    input: input.to_string(),
                    axis: d.axis.as_str().to_string(),
                    kind: DispositionKind::SpecCitation,
                    ferro_output: None,
                    tracking_issue: None,
                });
            }
            // `accepted_rejection` (ferro correctly errors; mutalyzer lenient) is
            // an accepted, terminal divergence — aggregate it under the
            // `AcceptedDivergence` summary column. The distinct, reviewable
            // disposition lives in `cases.json`. One summary row per per-axis
            // rejection (#870): a multi-axis row contributes one row per axis.
            for d in &case.accepted_rejections {
                rows.push(MemberRow {
                    cluster: d.cluster.clone(),
                    input: input.to_string(),
                    axis: d.axis.as_str().to_string(),
                    kind: DispositionKind::AcceptedDivergence,
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
    /// Per-axis accepted (non-bug) divergence annotations — e.g. a ferro policy
    /// decision that intentionally disagrees with mutalyzer while both outputs
    /// remain HGVS-spec-allowed. A mismatch on an axis that carries a divergence
    /// for that axis is tallied into `AxisTally::divergence_accepted` instead of
    /// `fail`.
    ///
    /// A `Case` may carry MORE THAN ONE divergence, each scoped to a different
    /// axis (#870) — e.g. a `normalized`-axis divergence AND a `genomic`-axis
    /// divergence on the same row (the copy-range literal-expansion policy
    /// diverges on both axes once the genomic axis routes through the user-facing
    /// `project_variant`). The matcher keys each divergence to its axis, so at
    /// most one per axis takes effect.
    ///
    /// The wire (`cases.json`) key is `accepted_divergence` and accepts either a
    /// single object (the common case, deserialized to a one-element vec) or an
    /// array of objects (a multi-axis row). Absent → empty vec. See
    /// [`de_one_or_many`].
    #[serde(
        default,
        rename = "accepted_divergence",
        deserialize_with = "de_one_or_many"
    )]
    pub accepted_divergences: Vec<AcceptedDivergence>,
    /// Per-axis known-ferro-bug annotations (xfail): ferro is wrong and the
    /// divergence is tracked by an issue rather than accepted as policy. A
    /// mismatch on an axis that carries a known bug for that axis is tallied into
    /// `AxisTally::known_bug` instead of `fail`. If ferro starts matching
    /// mutalyzer (XPASS), the harness fails so the annotation and now-fixed row
    /// are cleaned up.
    ///
    /// Like [`accepted_divergences`](Self::accepted_divergences), a `Case` may
    /// carry MORE THAN ONE known bug, each scoped to a different axis (#870). The
    /// wire key is `known_bug` and accepts a single object or an array. Absent →
    /// empty vec. See [`de_one_or_many`].
    #[serde(default, rename = "known_bug", deserialize_with = "de_one_or_many")]
    pub known_bugs: Vec<KnownBug>,
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
    /// Per-axis spec-citation annotations. A mismatch on an axis that carries a
    /// citation for that axis routes to `AxisTally::spec_overridden` (ferro's
    /// value is the spec-correct one; mutalyzer's expectation is non-preferred).
    ///
    /// A `Case` may carry MORE THAN ONE citation, each scoped to a different
    /// axis (#827) — e.g. a `normalized`-axis citation AND a
    /// `protein_description`-axis citation on the same row. The matcher keys each
    /// citation to its axis, so at most one citation per axis takes effect.
    ///
    /// The wire (`cases.json`) key is `spec_citation` and accepts either a single
    /// citation object (the common case, deserialized to a one-element vec) or an
    /// array of citation objects (a multi-axis row). Absent → empty vec. See
    /// [`de_one_or_many`].
    #[serde(default, rename = "spec_citation", deserialize_with = "de_one_or_many")]
    pub spec_citations: Vec<SpecCitation>,
    /// Marks a non-panic `Err` on a specific axis as ferro *correctly rejecting*
    /// an input mutalyzer leniently reinterprets (the only disposition that
    /// covers an `Err` on a projection axis). Tallied into
    /// `AxisTally::divergence_accepted`. See [`AcceptedRejection`].
    ///
    /// Like [`accepted_divergences`](Self::accepted_divergences) and
    /// [`known_bugs`](Self::known_bugs), a `Case` may carry MORE THAN ONE
    /// rejection, each scoped to a different axis (#870): an input ferro rejects
    /// at parse or projection time fails on *every* projected axis, so the same
    /// rejection is dispositioned on both the `normalized` and `genomic` axes
    /// once the genomic axis routes through the user-facing `project_variant`.
    /// The wire (`cases.json`) key is `accepted_rejection` and accepts either a
    /// single object or an array. Absent → empty vec. See
    /// [`de_one_or_many`].
    #[serde(
        default,
        rename = "accepted_rejection",
        deserialize_with = "de_one_or_many"
    )]
    pub accepted_rejections: Vec<AcceptedRejection>,
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
    /// ferro expands a positional / copy-range insertion payload
    /// (`ins N_M` / `delins N_M`) into the literal copied bases it denotes
    /// (e.g. `ins180_188` -> `insGCGAGGAAA`), whereas mutalyzer preserves the
    /// positional notation (or re-expresses it as a decomposed allele). Both
    /// describe the same edit; the literal inserted sequence is always a
    /// spec-valid HGVS form (`insertion.md` permits an explicit inserted
    /// sequence). Terminal rendering divergence (#654).
    #[serde(rename = "ferro-policy-654-positional-insert-literal-expansion")]
    PositionalInsertLiteralExpansion654,
    /// ferro preserves the input's transcript coordinate system, whereas
    /// mutalyzer re-expresses it: a `c.`-coding transcript addressed in `n.`
    /// (non-coding) coordinates keeps ferro's `n.` numbering rather than being
    /// renumbered to `c.`. ferro is lenient where mutalyzer enforces the
    /// coordinate-system convention; both reference the same position. Mirrors
    /// the errors-axis `transcript-coordinate-lenient` stance (#654).
    #[serde(rename = "ferro-policy-654-transcript-coordinate-preserved")]
    TranscriptCoordinatePreserved654,
    /// mutalyzer emitted no normalized form for the row (empty oracle value)
    /// while ferro normalizes the input to a spec-valid result. There is no
    /// competing normalized form to match, so the divergence is recorded rather
    /// than failing the run (#654).
    #[serde(rename = "ferro-policy-654-mutalyzer-no-normalized-form")]
    MutalyzerNoNormalizedForm654,
    /// For a no-op edit (a reference-identical `delins`, or a bare position with
    /// no operator) ferro renders the explicit positional identity it computed —
    /// `c.1_3=` / `g.274=` — while mutalyzer collapses it to the whole-reference
    /// `=` or a bare position. Both denote "no change" at that locus; ferro's
    /// positional `=` is spec-valid. Terminal rendering divergence (#654).
    #[serde(rename = "ferro-policy-654-explicit-identity-rendering")]
    ExplicitIdentityRendering654,
    /// On `coding_protein_descriptions` ferro enumerates the latest curated
    /// version of every transcript cdot reports overlapping the locus (#710:
    /// collapse superseded versions, prefer curated NM_/NR_ over predicted
    /// XM_/XR_), whereas mutalyzer lists a different overlapping-transcript set
    /// and partitions the input gene's own consequence onto its separate
    /// `protein_description` field. `background/refseq.md` L256/L259 leave "which
    /// transcript / reference sequence should one use" open and L264 endorses
    /// the latest build, so ferro's complete-current-curated-set enumeration is
    /// spec-valid; both representations are correct. Accepted divergence (#763).
    #[serde(rename = "ferro-policy-763-transcript-set-enumeration")]
    TranscriptSetEnumeration763,
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
            Policy::TranscriptSetEnumeration763 => "ferro-policy-763-transcript-set-enumeration",
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
            Policy::PositionalInsertLiteralExpansion654 => {
                "ferro-policy-654-positional-insert-literal-expansion"
            }
            Policy::TranscriptCoordinatePreserved654 => {
                "ferro-policy-654-transcript-coordinate-preserved"
            }
            Policy::MutalyzerNoNormalizedForm654 => "ferro-policy-654-mutalyzer-no-normalized-form",
            Policy::ExplicitIdentityRendering654 => "ferro-policy-654-explicit-identity-rendering",
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

/// Closed enum of reasons ferro *correctly* emits an error for an input that
/// mutalyzer leniently reinterprets. Like [`Policy`], adding a variant is a
/// deliberate code change so every "ferro-correctly-rejects" claim is reviewed.
#[derive(Debug, Deserialize, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum RejectionReason {
    /// ferro rejects a malformed / non-spec input — e.g. a substitution lacking
    /// a reference base (`c.41>CA`), or a concatenated/duplicated accession
    /// (`c.11LRG_199t1:c.11[10]`) — with a structured parse/normalize error,
    /// while mutalyzer leniently reinterprets it and emits a normalized form.
    /// ferro's rejection is the spec-correct behaviour (#654).
    #[serde(rename = "ferro-policy-654-malformed-input-rejected")]
    MalformedInputRejected654,
    /// ferro rejects an input that numbers a transcript's 5'/3' flank in `c.`
    /// (a `pter`/`qter` marker as a coordinate or inside a cross-reference /
    /// delins payload). A coding/non-coding DNA reference does not contain the
    /// gene's flanking sequences and "can not be used to describe variants"
    /// there (HGVS background/refseq.md L45-46), so the flank is not
    /// `c.`-numberable. mutalyzer leniently extrapolates `pter`/`qter` into the
    /// genomic flank and resolves the bases; ferro declines, which is the
    /// spec-correct behaviour (#758, consistent with the `c.pterdel`/`c.qterdel`
    /// `spec_citation` rows and #537).
    #[serde(rename = "ferro-policy-758-transcript-flank-not-numberable-in-c")]
    TranscriptFlankNotNumberableInC758,
    /// ferro rejects a c./n./r. input that names an explicit `accession.version`
    /// the prepared manifest does not carry, where a lenient provider would
    /// silently fall back to a *sibling* version and normalize against its
    /// frame/sequence. A c./r. description is defined only against its named
    /// version (versions differ in length/CDS offset, background/refseq.md), so
    /// substituting a sibling's frame yields a result whose stated reference and
    /// coordinate frame disagree. ferro declines with `TranscriptVersionNotExact`
    /// rather than substitute; mutalyzer (which carries the pinned version)
    /// normalizes it. ferro's decline is the spec-correct behaviour — a clean
    /// reference-unavailable miss rather than a confidently-wrong substitution
    /// (#785). Distinct from the `accession-version-absent` *no-op* disposition,
    /// which covers an absent version with no sibling to substitute (ferro echoes
    /// the input); here a sibling exists and the silent substitution is refused.
    #[serde(rename = "ferro-policy-785-version-substitution-refused")]
    VersionSubstitutionRefused785,
    /// ferro declines genomic projection of a transcript-coordinate variant on an
    /// `NG_` genomic reference whose selector is a gene symbol /
    /// `GENE_v001` / numeric gene-id, or is absent (bare `NG_:c.`). Gene symbols
    /// are not valid HGVS specifications (background/refseq.md L40-41) and a bare
    /// `NG_:c.` gives no transcript specification, so ferro cannot resolve a unique
    /// transcript and declines with `no parent reference`. mutalyzer leniently
    /// resolves the selector and projects anyway; ferro's refusal is the
    /// spec-correct behaviour. Optional convergence with mutalyzer is tracked in
    /// #872 (residual of #858; closed #327 covered explicit-`NM_`-parent inputs).
    #[serde(rename = "ferro-policy-858-nonstandard-or-absent-selector-declined")]
    NonstandardOrAbsentSelectorDeclined858,
    /// ferro declines to project a transcript-coordinate variant whose locus
    /// falls **outside the genomic reference's annotated coverage** of that
    /// transcript. An `NG_` RefSeqGene may annotate only the portion of a
    /// neighboring transcript that lies within its span — e.g. NG_009299.1
    /// (MYH11-centric) annotates NM_017668.3 (NDE1) only partially: `mRNA
    /// complement(135678..>137840)` / `CDS complement(137780..>137840)` (60 bp,
    /// 5'-incomplete), with NDE1 itself `complement(135678..>160896)` extending
    /// beyond the NG_. The authoritative transcript→genome alignment (cdot) puts
    /// the variant ~60 kb outside the NG_'s placed span, so there is no valid NG_
    /// coordinate for it. mutalyzer maps it against the partial annotation stub
    /// (treating the truncation boundary as c.1), yielding a coordinate that
    /// points at a different genomic base; ferro's decline is the spec-correct
    /// behaviour (#853, residual of #480/#655; subsumes #865).
    #[serde(rename = "ferro-policy-853-ng-partial-transcript-coverage-declined")]
    NgPartialTranscriptCoverageDeclined853,
    /// ferro emits **no projection** on this axis for a variant its policy says
    /// has no consequence there — e.g. a variant wholly outside the CDS
    /// (pure-5'UTR / pure-3'UTR) produces no protein prediction (the #857 policy,
    /// tested by `pure_3utr_deletion_has_no_protein_consequence`) — where the
    /// comparator instead emits a value (e.g. mutalyzer `p.(=)`, "analysed, no
    /// change"). The HGVS spec is ambiguous on no-projection vs the "analysed,
    /// unchanged" value for a non-CDS variant, so ferro's decline is a spec-
    /// defensible accepted rejection rather than a bug (#857/#891; #466 tracks the
    /// upstream spec clarification). This is the **only** `RejectionReason` that
    /// opts in to bucketing an empty-projection `Err` (see
    /// [`RejectionReason::disposition_empty_projection`]): the harness represents
    /// "ferro produced nothing" as an `EMPTY_PROJECTION_SENTINEL` `Err`, which is
    /// otherwise a hard FAIL. The empty→output transition is handled honestly:
    /// because `accepted_rejection` buckets only `Err`, if ferro later starts
    /// emitting a value it is NOT masked — an `Ok` that *matches* XPASS-FAILs (the
    /// stale annotation is caught), and an `Ok` that *mismatches* falls through to
    /// a hard FAIL (ferro stopped declining, now produces a differing value). The
    /// `#764` empty-projection count gate additionally catches any *rise* in the
    /// empty count.
    #[serde(rename = "ferro-policy-857-non-cds-no-projection")]
    NonCdsNoProjection857,
}

impl RejectionReason {
    /// Whether this reason opts in to dispositioning an **empty-projection**
    /// `Err` (ferro produced no output on the axis) as an accepted rejection.
    /// `accepted_rejection` normally excludes the `EMPTY_PROJECTION_SENTINEL`
    /// (an empty projection is otherwise a hard FAIL, #651); a reason returning
    /// `true` here relaxes that for its rows only. Closed match (`_ => false`) so
    /// adding a reason is a deliberate, reviewed opt-in — never an ad-hoc flag.
    pub fn disposition_empty_projection(self) -> bool {
        matches!(self, RejectionReason::NonCdsNoProjection857)
    }
    /// Stable identifier used in summary lines so reviewers can grep across runs.
    pub fn as_str(self) -> &'static str {
        match self {
            RejectionReason::MalformedInputRejected654 => {
                "ferro-policy-654-malformed-input-rejected"
            }
            RejectionReason::TranscriptFlankNotNumberableInC758 => {
                "ferro-policy-758-transcript-flank-not-numberable-in-c"
            }
            RejectionReason::VersionSubstitutionRefused785 => {
                "ferro-policy-785-version-substitution-refused"
            }
            RejectionReason::NonstandardOrAbsentSelectorDeclined858 => {
                "ferro-policy-858-nonstandard-or-absent-selector-declined"
            }
            RejectionReason::NgPartialTranscriptCoverageDeclined853 => {
                "ferro-policy-853-ng-partial-transcript-coverage-declined"
            }
            RejectionReason::NonCdsNoProjection857 => "ferro-policy-857-non-cds-no-projection",
        }
    }
}

impl std::fmt::Display for RejectionReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Accepted-**rejection** annotation. When attached to a `Case`, the corpus
/// runner treats a non-panic `Err` on `axis` as a tracked-but-non-failing
/// divergence: ferro *correctly* declines an input that mutalyzer leniently
/// reinterprets. Distinct from [`AcceptedDivergence`] (which covers a string
/// *mismatch* from a successful run) — this is the only disposition that
/// buckets an `Err` on a projection axis. A genuine panic still hard-FAILs, and
/// if ferro *starts* producing output (the rejection is gone) the harness FAILs
/// loudly (XPASS). The empty-projection sentinel also hard-FAILs by default —
/// **except** when `reason.disposition_empty_projection()` is `true` (#903), for
/// a reason like [`RejectionReason::NonCdsNoProjection857`] where ferro's
/// no-output is itself the spec-defensible decline; the `#764` empty-projection
/// count gate remains the backstop against a genuinely-new empty. Tallied into
/// `AxisTally::divergence_accepted` alongside `accepted_divergence`.
#[derive(Debug, Deserialize, Clone)]
#[allow(dead_code)]
pub struct AcceptedRejection {
    /// Axis this annotation applies to.
    pub axis: Axis,
    /// Closed reason identifier explaining why ferro's rejection is correct.
    pub reason: RejectionReason,
    /// Optional human-readable note expanding on the reason.
    #[serde(default)]
    pub note: Option<String>,
    /// Root-cause cluster id (see the corpus `clusters` registry).
    #[serde(default)]
    pub cluster: Option<String>,
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
    /// A deletion inside a DNA homopolymer / short tandem tract: ferro renders
    /// the spec-valid repeat contraction `N[k]` (DNA `repeated.md` L5 allows a
    /// unit of "one or more" nucleotides) while mutalyzer keeps a plain `del`.
    /// Same divergence as the `ferro-policy-745-homopolymer-repeat-contraction`
    /// accepted_divergence; cited here for the secondary axis of a row whose
    /// single `accepted_divergence` slot already documents the other axis (#745).
    #[serde(rename = "HGVS §Repeated (DNA contraction)")]
    RepeatDnaContraction,
    /// HGVS RNA alleles — a predicted (uncertain) cis allele places the `( )`
    /// predicted wrapper **inside** the allele bracket: `r.[(a;b)]`
    /// (`docs/recommendations/RNA/alleles.md`, example
    /// `LRG_199t1:r.[(578c>u;1339a>g;1680del)]`). ferro emits this spec-correct
    /// inside-bracket form; mutalyzer wraps the whole bracket
    /// (`r.([a;b])`), which has no counterpart in the spec. ferro's
    /// `r.[(…)]` is the spec-correct value (#693).
    #[serde(rename = "HGVS §RNA alleles (predicted bracket placement)")]
    RnaAlleleBracketPlacement,
    /// HGVS inserted inversion — an `<range>inv` segment in an `ins`/`delins`
    /// payload inserts the **reverse complement** (inverted copy) of that range's
    /// reference bases (`docs/recommendations/DNA/insertion.md` L52-62, e.g.
    /// `c.940_941ins[903_940inv;851_885inv]` = "an inverted copy of nucleotides";
    /// `docs/recommendations/DNA/inversion.md` L5, L48-49, L53-62 define an
    /// inversion as the reverse complement, explicitly *not* the plain reverse or
    /// plain complement). ferro reverse-complements the `inv` segment and then
    /// normalizes the resulting delins (common-prefix trim); mutalyzer leaves the
    /// segment forward (inversion not applied) and so also skips the trim, which
    /// is invalid HGVS. ferro's output is the spec-correct value (#854).
    #[serde(rename = "HGVS §Inserted inversion (reverse complement)")]
    InsertedInversion,
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
            SpecSection::RepeatDnaContraction => "HGVS §Repeated (DNA contraction)",
            SpecSection::RnaAlleleBracketPlacement => {
                "HGVS §RNA alleles (predicted bracket placement)"
            }
            SpecSection::InsertedInversion => "HGVS §Inserted inversion (reverse complement)",
        }
    }
}

impl std::fmt::Display for SpecSection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Spec-citation annotation. Records that the `cases.json` expected output
/// for `axis` has been corrected to ferro's spec-correct value (versus
/// mutalyzer's spec-incorrect output), with a citation to the relevant
/// HGVS spec section.
///
/// This annotation **does** affect tally bucketing: when a row mismatches
/// mutalyzer on `axis` and carries a matching `spec_citation`, `Tally::record`
/// routes it into the `spec_overridden` bucket (away from `fail`). The
/// *citation* is documentary, but the disposition is not — it is what keeps
/// these spec-overridden rows out of the failure tally.
///
/// Unlike `accepted_divergence`/`known_bug`/`improvement`/
/// `reference_unavailable`/`accepted_rejection`, `spec_citation` has no XPASS
/// guard in `Tally::record`: a row that later converges to mutalyzer (so
/// `actual == expected`) silently demotes to `pass` rather than FAILing to
/// flag the now-stale citation.
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

/// Deserialize a per-axis annotation field, accepting either a single object or
/// an array of objects and normalizing both to a `Vec<T>`. The single-object
/// form (the common case) deserializes to a one-element vec, so existing
/// single-annotation rows need no JSON change; a row that needs the annotation on
/// more than one axis uses the array form. An explicit `null` (like an absent
/// field, and like the prior `Option<T>` shapes this replaced) yields an empty
/// vec. Each element is still validated against `T`'s closed enums, so a typo in
/// either form is a hard parse error.
///
/// Shared by the `spec_citation` (#827), `accepted_divergence` (#870), and
/// `known_bug` (#870) fields — the declared `Vec<T>` type of the annotated field
/// drives the `T` inference, so one helper covers all three (and any future
/// per-axis annotation field gets it for free).
fn de_one_or_many<'de, D, T>(deserializer: D) -> Result<Vec<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    /// Untagged helper: serde tries `One` first, then `Many`, accepting whichever
    /// wire shape the field carries.
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum OneOrMany<T> {
        One(T),
        Many(Vec<T>),
    }

    Ok(match Option::<OneOrMany<T>>::deserialize(deserializer)? {
        None => Vec::new(),
        Some(OneOrMany::One(item)) => vec![item],
        Some(OneOrMany::Many(items)) => items,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn non_cds_no_projection_857_opts_in_to_empty_projection() {
        // #903: exactly one reason opts in to dispositioning an empty projection.
        assert!(RejectionReason::NonCdsNoProjection857.disposition_empty_projection());
        for r in [
            RejectionReason::MalformedInputRejected654,
            RejectionReason::TranscriptFlankNotNumberableInC758,
            RejectionReason::VersionSubstitutionRefused785,
            RejectionReason::NonstandardOrAbsentSelectorDeclined858,
            RejectionReason::NgPartialTranscriptCoverageDeclined853,
        ] {
            assert!(
                !r.disposition_empty_projection(),
                "{} must NOT opt in to empty-projection bucketing",
                r.as_str()
            );
        }
    }

    #[test]
    fn accepted_rejection_857_reason_round_trips() {
        // The new reason deserializes from its serde rename and stringifies back.
        let ar: AcceptedRejection = serde_json::from_str(
            r#"{"axis":"protein_description","reason":"ferro-policy-857-non-cds-no-projection"}"#,
        )
        .expect("accepted_rejection with the #903 reason should deserialize");
        assert_eq!(ar.reason, RejectionReason::NonCdsNoProjection857);
        assert_eq!(ar.reason.as_str(), "ferro-policy-857-non-cds-no-projection");
    }

    const BASE: &str =
        r#""description":"t","source":"t","source_commit":"t","license":"t","refreshed_at":"t""#;

    fn parse(clusters: &str, cases: &str) -> Fixture {
        let json = format!("{{{BASE},\"clusters\":[{clusters}],\"cases\":[{cases}]}}");
        serde_json::from_str(&json).expect("fixture should deserialize")
    }

    #[test]
    fn comparator_provenance_absent_still_parses() {
        // Backward-compat (#890): a corpus header with no `comparator_provenance`
        // block (the pre-#890 shape, and every other corpus) parses with the
        // field as `None`.
        let f = parse("", "");
        assert!(f.comparator_provenance.is_none());
    }

    #[test]
    fn comparator_provenance_parses_and_exposes_versions() {
        // #890: the recorded provenance deserializes into typed fields (not
        // silently dropped as an unknown key), so a future test can assert the
        // recorded reference_identity matches the pinned reference.
        let json = format!(
            "{{{BASE},{},\"clusters\":[],\"cases\":[]}}",
            r#""comparator_provenance":{
                "note":"n","validated_at":"2026-06-30",
                "reference_identity":"c098e18ec6b01df2",
                "mutalyzer":"3.1.1","mutalyzer_hgvs_parser":"0.3.9",
                "biocommons_hgvs":"1.5.6","biocommons_seqrepo":"0.6.11",
                "hgvs_rs":"0.20.2-1-gb6513b3"
            }"#
        );
        let f: Fixture = serde_json::from_str(&json).expect("provenance should deserialize");
        let p = f
            .comparator_provenance
            .expect("comparator_provenance present");
        assert_eq!(p.reference_identity, "c098e18ec6b01df2");
        assert_eq!(p.mutalyzer, "3.1.1");
        assert_eq!(p.mutalyzer_hgvs_parser, "0.3.9");
        assert_eq!(p.biocommons_hgvs, "1.5.6");
        assert_eq!(p.biocommons_seqrepo, "0.6.11");
        assert_eq!(p.hgvs_rs, "0.20.2-1-gb6513b3");
    }

    #[test]
    fn cluster_refs_collects_every_disposition_kind() {
        let clusters = r#"
            {"id":"sel","title":"RefSeqGene selector","spec_section":"background/refseq.md"},
            {"id":"np","title":"bare NP","spec_section":"protein"},
            {"id":"rej","title":"malformed input rejected","spec_section":"DNA"}
        "#;
        let cases = r#"
            {"input":"A","improvement":{"axis":"normalized","tracking_issue":500,
              "section":"HGVS §RefSeqGene transcript selection","cluster":"sel"}},
            {"input":"B","spec_citation":{"axis":"protein_description",
              "section":"HGVS protein reference (bare NP)","cluster":"np"}},
            {"input":"C","accepted_rejection":{"axis":"normalized",
              "reason":"ferro-policy-654-malformed-input-rejected","cluster":"rej"}}
        "#;
        let fixture = parse(clusters, cases);
        let mut refs = fixture.cluster_refs();
        refs.sort();
        assert_eq!(refs, vec![("A", "sel"), ("B", "np"), ("C", "rej")]);
        assert!(fixture.validate_clusters().is_ok());

        let summary = fixture.to_summary();
        assert_eq!(summary.title, "mutalyzer-normalize");
        assert_eq!(summary.rows.len(), 3);
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
        // `accepted_rejection` aggregates under the `AcceptedDivergence` summary
        // column; its `cluster` must still be collected by `cluster_refs`.
        let rej = summary
            .rows
            .iter()
            .find(|r| r.kind == DispositionKind::AcceptedDivergence)
            .expect("accepted_rejection row (mapped to AcceptedDivergence)");
        assert_eq!(rej.input, "C");
        assert_eq!(rej.axis, "normalized");
        assert_eq!(rej.cluster.as_deref(), Some("rej"));
        assert_eq!(rej.tracking_issue, None);
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

    // #827: a single `spec_citation` object deserializes to a one-element vec
    // (backward compat) and an array deserializes to a multi-element vec, with
    // both `cluster_refs` and `to_summary` enumerating every citation.
    #[test]
    fn spec_citation_single_object_and_array_both_parse() {
        // Single-object form (unchanged wire shape) -> one citation.
        let one = r#"{"input":"S","spec_citation":{"axis":"normalized",
            "section":"HGVS §Prioritization"}}"#;
        let fixture = parse("", one);
        assert_eq!(fixture.cases[0].spec_citations.len(), 1);
        assert_eq!(fixture.cases[0].spec_citations[0].axis, Axis::Normalized);

        // Array form -> two citations on two axes, each with its own cluster.
        let clusters = r#"
            {"id":"cre","title":"coding repeat","spec_section":"DNA"},
            {"id":"np","title":"bare NP","spec_section":"protein"}
        "#;
        let many = r#"{"input":"M","spec_citation":[
            {"axis":"normalized","section":"HGVS §Repeated (coding codon exception)","cluster":"cre"},
            {"axis":"protein_description","section":"HGVS protein reference (bare NP)","cluster":"np"}
        ]}"#;
        let fixture = parse(clusters, many);
        let case = &fixture.cases[0];
        assert_eq!(case.spec_citations.len(), 2);
        assert_eq!(case.spec_citations[0].axis, Axis::Normalized);
        assert_eq!(case.spec_citations[1].axis, Axis::ProteinDescription);

        // Both citations contribute a cluster ref and a summary row.
        let mut refs = fixture.cluster_refs();
        refs.sort();
        assert_eq!(refs, vec![("M", "cre"), ("M", "np")]);
        assert!(fixture.validate_clusters().is_ok());
        let spec_rows = fixture
            .to_summary()
            .rows
            .into_iter()
            .filter(|r| r.kind == DispositionKind::SpecCitation)
            .count();
        assert_eq!(spec_rows, 2);
    }

    // #827: an explicit `"spec_citation": null` deserializes to an empty vec,
    // matching both the absent case and the prior `Option<SpecCitation>`
    // null-handling (null is semantically "no citation"), rather than erroring.
    #[test]
    fn spec_citation_explicit_null_is_empty() {
        let fixture = parse("", r#"{"input":"N","spec_citation":null}"#);
        assert!(fixture.cases[0].spec_citations.is_empty());
    }

    #[test]
    fn accepted_rejection_single_object_and_array_both_parse() {
        // Single-object form (unchanged wire shape) -> one rejection.
        let one = r#"{"input":"S","accepted_rejection":{"axis":"normalized",
            "reason":"ferro-policy-654-malformed-input-rejected"}}"#;
        let fixture = parse("", one);
        assert_eq!(fixture.cases[0].accepted_rejections.len(), 1);
        assert_eq!(
            fixture.cases[0].accepted_rejections[0].axis,
            Axis::Normalized
        );

        // Array form (#870) -> the same rejection dispositioned on two axes.
        let many = r#"{"input":"M","accepted_rejection":[
            {"axis":"normalized","reason":"ferro-policy-654-malformed-input-rejected"},
            {"axis":"genomic","reason":"ferro-policy-654-malformed-input-rejected"}
        ]}"#;
        let case = &parse("", many).cases[0];
        assert_eq!(case.accepted_rejections.len(), 2);
        assert_eq!(case.accepted_rejections[0].axis, Axis::Normalized);
        assert_eq!(case.accepted_rejections[1].axis, Axis::Genomic);

        // Absent -> empty vec.
        assert!(parse("", r#"{"input":"N"}"#).cases[0]
            .accepted_rejections
            .is_empty());
    }

    // #870: two rejections for the SAME axis on one case are rejected by
    // `validate_clusters` — the matcher honors only one per axis, so a duplicate
    // is an authoring mistake that must surface loudly (parity with
    // `duplicate_spec_citation_axis_is_rejected`).
    #[test]
    fn duplicate_accepted_rejection_axis_is_rejected() {
        let cases = r#"{"input":"D","accepted_rejection":[
            {"axis":"genomic","reason":"ferro-policy-654-malformed-input-rejected"},
            {"axis":"genomic","reason":"ferro-policy-758-transcript-flank-not-numberable-in-c"}
        ]}"#;
        let err = parse("", cases)
            .validate_clusters()
            .expect_err("two rejections on the same axis must be rejected");
        assert!(err.contains("genomic"), "{err}");
        assert!(err.contains("accepted_rejection"), "{err}");
    }

    // #827: two spec citations for the SAME axis on one case are rejected by
    // `validate_clusters` — the matcher honors only one per axis, so a duplicate
    // is an authoring mistake that must surface loudly.
    #[test]
    fn duplicate_spec_citation_axis_is_rejected() {
        let cases = r#"{"input":"D","spec_citation":[
            {"axis":"normalized","section":"HGVS §Prioritization"},
            {"axis":"normalized","section":"HGVS §Repeated (coding codon exception)"}
        ]}"#;
        let err = parse("", cases)
            .validate_clusters()
            .expect_err("two citations on the same axis must be rejected");
        assert!(err.contains("normalized"), "{err}");

        // Two citations on DIFFERENT axes are fine (the #827 use case).
        let ok = r#"{"input":"O","spec_citation":[
            {"axis":"normalized","section":"HGVS §Prioritization"},
            {"axis":"protein_description","section":"HGVS protein reference (bare NP)"}
        ]}"#;
        assert!(parse("", ok).validate_clusters().is_ok());
    }

    // #827: a typo in any element of the array form is still a hard parse error
    // (the closed `Axis`/`SpecSection` enums apply per element).
    #[test]
    fn spec_citation_array_element_typo_rejected() {
        let json = format!(
            "{{{BASE},\"cases\":[{}]}}",
            r#"{"input":"x","spec_citation":[
                {"axis":"normalized","section":"HGVS §Prioritization"},
                {"axis":"protein_description","section":"not a real section"}
            ]}"#
        );
        let result: Result<Fixture, _> = serde_json::from_str(&json);
        assert!(
            result.is_err(),
            "expected a typo'd section in an array element to be rejected; got {:?}",
            result.ok(),
        );
    }
}
