//! `VariantProjector` orchestrator.

use crate::data::cdot::CdotTranscript;
use crate::data::mapping::MappingInfo;
use crate::data::projection::Projector;
use crate::error::FerroError;
use crate::hgvs::edit::NaEdit;
use crate::hgvs::interval::{CdsInterval, TxInterval};
use crate::hgvs::location::{CdsPos, GenomePos, TxPos};
use crate::hgvs::variant::{
    is_frameshift, Accession, AllelePhase, AlleleVariant, CdsVariant, HgvsVariant, LocEdit,
    TxVariant,
};
use crate::normalize::{NormalizeConfig, Normalizer};
use crate::project::accession::parse_accession;
use crate::project::edit::transform_edit_for_strand;
use crate::project::protein::{
    build_initiator_unknown, build_whole_protein_unknown, cds_has_recognized_start,
    cds_has_valid_start, edit_reaches_initiation_codon, edit_spans_cds_into_3utr,
    predict_indel_protein, predict_stop_region_extension, predict_substitution_protein,
    read_cds_start_codon, whole_exon_deletion_span, RefProteinBundle,
};
use crate::project::result::VariantProjection;
use crate::reference::transcript::Transcript;
use crate::reference::ReferenceProvider;
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

/// Cache key for [`VariantProjector`]'s transcript cache:
/// `(transcript_id, parent_accession)`. The parent component is `None` for
/// non-variant-aware lookups and for inputs without a `genomic_context`; it
/// is `Some` only when the variant-aware path observes an NG/NC parent that
/// could steer the lookup to a non-primary cdot build (issue #332).
type TranscriptCacheKey = (String, Option<String>);

/// Build a bare whole-protein-unknown `p.?` variant for an in-cis frameshift
/// allele (#855), reusing the protein accession and gene symbol of an existing
/// member protein consequence so no transcript fetch is needed. `member` is
/// expected to be an `HgvsVariant::Protein`; any other shape is returned as-is.
fn build_cis_whole_protein_unknown(member: &HgvsVariant) -> HgvsVariant {
    use crate::hgvs::edit::ProteinEdit;
    use crate::hgvs::interval::ProtInterval;
    use crate::hgvs::location::{AminoAcid, ProtPos};
    use crate::hgvs::variant::{LocEdit, ProteinVariant};
    match member {
        HgvsVariant::Protein(pv) => HgvsVariant::Protein(ProteinVariant {
            accession: pv.accession.clone(),
            gene_symbol: pv.gene_symbol.clone(),
            // The location is ignored for a whole-protein unknown; a placeholder
            // `Met1` point is fine. `LocEdit::new` (not `new_predicted`) renders
            // the bare `p.?`, not `p.(?)`.
            loc_edit: LocEdit::new(
                ProtInterval::point(ProtPos::new(AminoAcid::Met, 1)),
                ProteinEdit::whole_protein_unknown(),
            ),
        }),
        other => other.clone(),
    }
}

/// Sort projected genomic allele members into ascending genomic order
/// (recommendations/DNA/alleles.md:118 "Variants should be listed in genomic
/// order"). A minus-strand transcript projects members in descending genomic
/// order, so this re-canonicalizes. Stable (`sort_by_key`); a non-Genome or
/// unresolved-start member sorts to the end without reordering among such
/// members. Scoped to the projection allele path — not the global normalizer —
/// to bound blast radius, and applied by the caller to cis alleles only (#851).
fn sort_genomic_allele_members(members: &mut [HgvsVariant]) {
    members.sort_by_key(|m| match m {
        HgvsVariant::Genome(g) => g
            .loc_edit
            .location
            .start
            .inner()
            .map(|p| p.base)
            .unwrap_or(u64::MAX),
        _ => u64::MAX,
    });
}

/// Combine the per-member protein consequences of an **in-cis** allele into a
/// single description when the members are adjacent single substitutions (#855).
///
/// Per HGVS `delins.md`, adjacent in-cis residue changes are described as one
/// combined delins (`p.(Asp92_Tyr93delinsTyrCys)`), not a bracketed list
/// (`p.[(Asp92Tyr);(Tyr93Cys)]`, explicitly "not correct"); variants separated
/// by one or more unchanged residues are described individually.
///
/// Returns:
/// - `Some(single delins/substitution)` when all members are single
///   substitutions that collapse to one contiguous run;
/// - `Some(allele)` when they form multiple contiguous runs (each rendered, then
///   bracketed) — separated runs stay individual;
/// - `None` (caller keeps the existing bracketed allele) when any member is not a
///   single point substitution, two members hit the *same* residue (ambiguous —
///   needs codon-level combination, not description merging), or there is nothing
///   to merge (every run is a single residue, i.e. the existing bracket is
///   already the correct individual-variant rendering).
fn combine_cis_substitution_proteins(members: &[HgvsVariant]) -> Option<HgvsVariant> {
    use crate::hgvs::edit::{AminoAcidSeq, ProteinEdit};
    use crate::hgvs::interval::ProtInterval;
    use crate::hgvs::location::{AminoAcid, ProtPos};
    use crate::hgvs::variant::{Accession, LocEdit, ProteinVariant};

    struct SubSpan {
        pos: u64,
        ref_aa: AminoAcid,
        alt_aa: AminoAcid,
    }

    let mut spans: Vec<SubSpan> = Vec::with_capacity(members.len());
    let mut accession: Option<Accession> = None;
    let mut gene_symbol: Option<String> = None;
    for m in members {
        let HgvsVariant::Protein(pv) = m else {
            return None;
        };
        let edit = pv.loc_edit.edit.inner()?;
        let ProteinEdit::Substitution {
            reference,
            alternative,
        } = edit
        else {
            return None;
        };
        let start = pv.loc_edit.location.start.inner()?;
        let end = pv.loc_edit.location.end.inner()?;
        if start.number != end.number {
            return None; // not a point substitution
        }
        if accession.is_none() {
            accession = Some(pv.accession.clone());
            gene_symbol = pv.gene_symbol.clone();
        }
        spans.push(SubSpan {
            pos: start.number,
            ref_aa: *reference,
            alt_aa: *alternative,
        });
    }
    if spans.len() < 2 {
        return None;
    }
    spans.sort_by_key(|s| s.pos);

    // Group consecutive residues (no unchanged residue between them).
    let mut groups: Vec<Vec<SubSpan>> = Vec::new();
    for s in spans {
        match groups.last_mut() {
            Some(g) if s.pos <= g.last().expect("non-empty group").pos + 1 => {
                if s.pos == g.last().expect("non-empty group").pos {
                    // Two substitutions on the same residue: this needs codon-level
                    // combination, which per-member description merging cannot do.
                    // Bail to the bracketed allele rather than emit a malformed
                    // single-residue delins.
                    return None;
                }
                g.push(s);
            }
            _ => groups.push(vec![s]),
        }
    }
    // Nothing merged → the existing bracket of individual substitutions is already
    // the spec-correct rendering for separated variants.
    if groups.iter().all(|g| g.len() == 1) {
        return None;
    }

    let accession = accession?;
    let render_group = |g: &[SubSpan]| -> ProteinVariant {
        if g.len() == 1 {
            let s = &g[0];
            ProteinVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new_predicted(
                    ProtInterval::point(ProtPos::new(s.ref_aa, s.pos)),
                    ProteinEdit::Substitution {
                        reference: s.ref_aa,
                        alternative: s.alt_aa,
                    },
                ),
            }
        } else {
            let first = &g[0];
            let last = &g[g.len() - 1];
            let inserted: Vec<AminoAcid> = g.iter().map(|s| s.alt_aa).collect();
            ProteinVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new_predicted(
                    ProtInterval::new(
                        ProtPos::new(first.ref_aa, first.pos),
                        ProtPos::new(last.ref_aa, last.pos),
                    ),
                    ProteinEdit::Delins {
                        sequence: AminoAcidSeq::new(inserted),
                    },
                ),
            }
        }
    };

    if groups.len() == 1 {
        Some(HgvsVariant::Protein(render_group(&groups[0])))
    } else {
        let variants: Vec<HgvsVariant> = groups
            .iter()
            .map(|g| HgvsVariant::Protein(render_group(g)))
            .collect();
        Some(HgvsVariant::Allele(AlleleVariant::new(
            variants,
            AllelePhase::Cis,
        )))
    }
}

/// Map an LRG transcript accession (`LRG_<gnum>t<tnum>`) to its protein
/// accession (`LRG_<gnum>p<tnum>`). LRG defines `p<k>` as the translation of
/// `t<k>`, so the index is shared — a structural fact of the accession, not a
/// cdot inference (the LRG mapping table carries no protein column). Returns
/// `None` for any accession that is not an `LRG_<digits>t<digits>` (#860).
fn lrg_protein_accession(lrg_transcript: &str) -> Option<String> {
    let rest = lrg_transcript.strip_prefix("LRG_")?;
    let (gnum, tnum) = rest.split_once('t')?;
    if gnum.is_empty()
        || tnum.is_empty()
        || !gnum.bytes().all(|b| b.is_ascii_digit())
        || !tnum.bytes().all(|b| b.is_ascii_digit())
    {
        return None;
    }
    Some(format!("LRG_{gnum}p{tnum}"))
}

pub struct VariantProjector<P: ReferenceProvider + Clone> {
    projector: Projector,
    provider: P,
    normalizer: Normalizer<P>,
    /// Cache of fetched transcripts keyed by `(transcript_id, parent_accession)`.
    ///
    /// `project_single_inner` looks up each overlapping transcript via the
    /// provider once per variant; for `project_all` workloads (hundreds of
    /// overlapping transcripts per variant on dense chromosomes) the same
    /// transcript ID is fetched repeatedly. Caching collapses N→1 fetches
    /// per unique key over the projector's lifetime.
    ///
    /// The parent-accession component is non-`None` only on the variant-aware
    /// lookup path (`cached_get_transcript_for_variant`), where the variant's
    /// `genomic_context` (NG/NC parent) can steer
    /// `ReferenceProvider::get_transcript_for_variant` to a different cdot
    /// build for the same `transcript_id` (issue #332). Inputs without a
    /// `genomic_context` and the legacy `cached_get_transcript` entry point
    /// both cache under `(transcript_id, None)`, so they share entries.
    transcript_cache: RwLock<HashMap<TranscriptCacheKey, Arc<Transcript>>>,
    /// Cache of the translated reference protein (ref CDS + ref protein +
    /// ref-with-stop) keyed by the same `(transcript_id, parent_accession)`
    /// identity used by `transcript_cache`. Each entry is stable for the life
    /// of the projector and is consumed by `predict_indel_protein` to avoid
    /// retranslating the full CDS on every indel.
    ///
    /// The parent-aware key matters because one `transcript_id` can resolve to
    /// different transcript records (and therefore different reference CDS
    /// bases) under different NG/NC parents (issue #332). Keying by
    /// `transcript_id` alone would let the first resolved build poison indel
    /// protein predictions for projections on the other build.
    ref_protein_cache: RwLock<HashMap<TranscriptCacheKey, Arc<RefProteinBundle>>>,
    /// Explicit genome-build override for build-agnostic inputs (#715).
    ///
    /// A bare `NG_`/`LRG_` input carries no build, so by default it resolves to
    /// the GRCh38-preferred placement. When set (via [`Self::with_assembly`]),
    /// this build (`"GRCh37"` / `"GRCh38"`) *fills in* for such inputs — it does
    /// **not** override a build the input accession already encodes
    /// (`NC_*.10`/`.11`, `GRCh37(...)`); the accession's build stays
    /// authoritative (see [`Self::build_hint_for_variant`]).
    assembly_override: Option<&'static str>,
}

impl<P: ReferenceProvider + Clone> VariantProjector<P> {
    pub fn new(projector: Projector, provider: P) -> Self {
        let normalizer = Normalizer::with_config(provider.clone(), NormalizeConfig::default());
        Self {
            projector,
            provider,
            normalizer,
            transcript_cache: RwLock::new(HashMap::new()),
            ref_protein_cache: RwLock::new(HashMap::new()),
            assembly_override: None,
        }
    }

    /// Non-variant-aware transcript fetch with caching. Kept for inline tests
    /// that need to populate a `Transcript` for `cached_ref_translation`
    /// without constructing an `HgvsVariant`. Hot-path callers go through
    /// [`Self::cached_get_transcript_for_variant`] so an NG/NC-parented input
    /// resolves to the correct cdot build (issue #332).
    #[cfg(test)]
    fn cached_get_transcript(&self, transcript_id: &str) -> Result<Arc<Transcript>, FerroError> {
        let key = (transcript_id.to_string(), None);
        if let Some(tx) = self
            .transcript_cache
            .read()
            .expect("transcript cache poisoned")
            .get(&key)
        {
            return Ok(Arc::clone(tx));
        }
        let tx = self.provider.get_transcript(transcript_id)?;
        let mut guard = self
            .transcript_cache
            .write()
            .expect("transcript cache poisoned");
        // `tx` is already an `Arc<Transcript>` from the provider; insert directly.
        let entry = guard.entry(key).or_insert_with(|| tx);
        Ok(Arc::clone(entry))
    }

    /// Variant-aware transcript lookup with caching.
    ///
    /// Routes through [`ReferenceProvider::get_transcript_for_variant`] so an
    /// NG/NC-parented input picks the build-correct cdot transcript (issue
    /// #332). The cache key is `(transcript_id, parent_accession)` because
    /// the same `transcript_id` can resolve to different transcripts under
    /// different parents.
    ///
    /// On `ReferenceNotFound` from the variant-aware call, falls back to the
    /// bare provider lookup (matches the prior `tx_for_codon_with_fallback`
    /// contract). Other provider errors propagate.
    /// Build the parent-aware cache key component for a variant.
    ///
    /// Returns `Some(parent.to_string())` when the variant's accession carries
    /// a `genomic_context` (NG/NC parent), `None` otherwise. Shared by
    /// [`Self::cached_get_transcript_for_variant`] and
    /// [`Self::cached_ref_translation`] so both caches use byte-identical keys
    /// for the same variant.
    fn parent_key_for(variant: &HgvsVariant) -> Option<String> {
        variant
            .accession()
            .and_then(|a| a.genomic_context.as_deref())
            .map(|gc| gc.to_string())
    }

    /// The build-bearing genomic accession to stamp as `genomic_context` when
    /// fetching a transcript for [`crate::project::rna::predict_rna`], so the
    /// fetch partitions by genome build the same way the single genome-pivot
    /// path does via its parallel `cache_variant` (#843).
    ///
    /// `predict_rna` reads the transcript *sequence* (e.g. to 3'-shift a
    /// deletion or resolve a range insert), so the served bases must come from
    /// the build the `input` is addressed on — a `g.` allele on `NC_*.10`
    /// (GRCh37) must read the GRCh37 transcript, not the GRCh38 primary. The
    /// raw `cached_get_transcript_for_variant(input, …)` keys on
    /// [`Self::parent_key_for`], which is `None` for a bare `g.` input (its
    /// build-bearing `NC_*` accession lives in `accession`, not
    /// `genomic_context`), collapsing both builds onto the primary transcript.
    ///
    /// Resolution order, mirroring the genome-pivot path's build selection:
    /// 1. If `input`'s accession already carries a `genomic_context` (an
    ///    NG/NC-parented c./n./r. input), reuse it — that context already
    ///    drives the correct build probe, so the cache identity is unchanged
    ///    (a no-op for the existing NG-parented corpus cases).
    /// 2. Else, if `input`'s own accession is a build-bearing genomic
    ///    accession (a `g.` on `NC_*.10`/`.11`), use the accession itself as
    ///    the context — the new build-scoping for bare `g.` inputs.
    /// 3. Else (truly build-agnostic — bare `NM_`/`NG_`/`LRG_`), `None`: the
    ///    primary-build fetch is the only defined answer and no divergence is
    ///    possible, so behavior is unchanged.
    ///
    /// The [`Self::assembly_override`] (`--assembly`) is deliberately *not*
    /// consulted here. It only `fills in` a build for a build-agnostic
    /// `NG_`/`LRG_` parent, and such a parent's transcript bases are
    /// build-independent: RefSeqGene/`LRG_` transcripts are served from the
    /// build-agnostic transcript FASTA, so `predict_rna` (which reads the
    /// transcript *sequence*) produces identical output under either build. The
    /// override changes only the genomic *placement* (the `.genomic` re-anchor),
    /// which the genome-pivot path handles separately via its build-resolved
    /// `cache_variant`. Folding the override in here would stamp a build the
    /// fetch cannot observably honor (no per-build NG transcript bases exist) —
    /// a no-op that no test could exercise — so it is intentionally omitted.
    fn rna_build_context(&self, input: &HgvsVariant) -> Option<Accession> {
        let accession = input.accession()?;
        if let Some(parent) = accession.genomic_context.as_deref() {
            return Some(parent.clone());
        }
        // Reached only when `accession.genomic_context` is None (step 1 returned
        // early otherwise), so the accession is already bare — a build-bearing
        // genomic accession (a `g.` on `NC_*.10`/`.11`) is used directly as the
        // context to stamp.
        if self.provider.infer_genome_build(accession).is_some() {
            return Some(accession.clone());
        }
        None
    }

    /// Fetch the transcript for a `predict_rna` call, build-scoped by `input`
    /// (#843). Builds a minimal `cache_variant` that stamps the build-bearing
    /// genomic accession from [`Self::rna_build_context`] as `genomic_context`
    /// on `transcript_id`, then resolves through the parent-aware
    /// [`Self::cached_get_transcript_for_variant`] so the transcript/ref caches
    /// partition by build — exactly the mechanism the single genome-pivot path
    /// uses (`cache_variant` at the `predict_rna` site). When no build-bearing
    /// context is derivable the cache_variant is bare and the fetch is the
    /// (unchanged) primary-build lookup.
    ///
    /// `predict_rna` walks the transcript *bases*, so unlike the coding/protein
    /// axes (which read build-identical CDS coordinates) it must not run against
    /// another build's sequence. The build-stamped fetch above resolves the
    /// correct build when the provider has it, but the parent-aware cache it
    /// shares with the g.→c./n. fetch can hold a build-agnostic *primary*
    /// transcript that the provider degraded to on a build-keyed miss (the
    /// `ReferenceNotFound` fallback in `cached_get_transcript_for_variant`).
    /// Serving that for a non-primary input would emit a wrong-build `r.` (#843),
    /// so when the input encodes a build, require the resolved transcript's build
    /// to match it; a *definite* mismatch (a known, different build) declines
    /// rather than predicting against the wrong bases. A build-agnostic
    /// transcript (`GenomeBuild::Unknown`, e.g. FASTA-served bases that are
    /// build-independent) is accepted — there is no divergence to get wrong.
    fn predict_rna_transcript(
        &self,
        input: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<Arc<Transcript>, FerroError> {
        let mut accession = parse_accession(transcript_id);
        if let Some(context) = self.rna_build_context(input) {
            accession = accession.with_genomic_context(context);
        }
        // A minimal Cds carrier: `cached_get_transcript_for_variant` only reads
        // `accession()` (for the parent key and the build-aware provider probe),
        // so the location/edit are immaterial.
        let cache_variant = HgvsVariant::Cds(CdsVariant {
            accession,
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: crate::hgvs::edit::Base::A,
                    alternative: crate::hgvs::edit::Base::G,
                },
            ),
        });
        let tx = self.cached_get_transcript_for_variant(&cache_variant, transcript_id)?;
        if let Some(input_build) = self.infer_input_build(input) {
            let tx_build = tx.genome_build;
            if tx_build != crate::reference::transcript::GenomeBuild::Unknown
                && tx_build.to_string() != input_build
            {
                return Err(FerroError::ReferenceNotFound {
                    id: format!(
                        "{transcript_id} on build {input_build}: resolved transcript is \
                         build {tx_build} (refusing to predict r. against wrong-build bases)"
                    ),
                });
            }
        }
        Ok(tx)
    }

    /// Reject parentless c./n./r. fan-out inputs up front.
    ///
    /// `project_normalized_all` seeds a stab query against the contig
    /// derived from cdot's exon table, then projects through every
    /// overlapping transcript via `project_single_inner`. For c./n./r.
    /// inputs `project_single_inner` calls `project_to_genomic`, which
    /// requires an explicit NG/NC parent on `accession.genomic_context`
    /// (see #327). Without the parent every per-transcript projection
    /// fails with `UnsupportedProjection { "no parent reference..." }`,
    /// and the fan-out loop's `log::trace!`-and-drop converts the user
    /// error into `Ok([])` — silently turning a missing-context error
    /// into "no overlaps." Surface the error here, before we touch the
    /// stab-query path (#389 follow-up).
    fn require_parent_for_fanout(
        accession: &crate::hgvs::variant::Accession,
        axis: &str,
    ) -> Result<(), FerroError> {
        if accession.genomic_context.is_none() {
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "project_all requires a parent reference (genomic_context) on {} \
                     inputs to resolve back to g.; accession {} has none",
                    axis,
                    accession.full()
                ),
            });
        }
        Ok(())
    }

    /// Infer the genome build the *input variant itself* encodes.
    ///
    /// Examines the variant's own accession first: a g. variant on
    /// `NC_*.10` carries GRCh37, on `NC_*.11` carries GRCh38, and an
    /// assembly-style `GRCh37(chr1):g.…` carries GRCh37 directly. If the
    /// variant has no inherent build (e.g. a c. variant on a bare NM)
    /// but does have a `genomic_context` (NG/NC parent), the parent is
    /// consulted. Anything else returns `None` (e.g. a bare `NG_`/`LRG_`,
    /// which carries no build), so the caller falls back to the
    /// build-agnostic primary-cdot path (issue #389).
    fn infer_input_build(&self, variant: &HgvsVariant) -> Option<&'static str> {
        let acc = variant.accession()?;
        if let Some(b) = self.provider.infer_genome_build(acc) {
            return Some(b);
        }
        if let Some(parent) = acc.genomic_context.as_deref() {
            return self.provider.infer_genome_build(parent);
        }
        None
    }

    /// The genome build to project `variant` under, honoring the
    /// [`Self::with_assembly`] override (#715).
    ///
    /// The build the input *encodes* is authoritative; the `--assembly` override
    /// only *fills in* when the input is build-agnostic (a bare `NG_`/`LRG_`).
    /// A contradictory override is warned and ignored at the public entry points
    /// (see [`Self::warn_assembly_conflict`]) — it never reinterprets an explicit
    /// `NC_*.10`/`.11` accession, which would emit wrong coordinates under a
    /// right-looking accession (the #655/#702 failure class).
    fn build_hint_for_variant(&self, variant: &HgvsVariant) -> Option<&'static str> {
        self.infer_input_build(variant).or(self.assembly_override)
    }

    /// Emit a warning when an explicit `--assembly` override contradicts the
    /// build the input accession already encodes (#715). The input wins; this
    /// only flags the ignored override. Called once per public projection entry.
    fn warn_assembly_conflict(&self, variant: &HgvsVariant) {
        if let (Some(override_build), Some(input_build)) =
            (self.assembly_override, self.infer_input_build(variant))
        {
            if override_build != input_build {
                log::warn!(
                    "ignoring --assembly {override_build}: input accession already \
                     specifies {input_build}"
                );
            }
        }
    }

    /// Build-aware cdot transcript lookup shared by every projection
    /// entry point: routes through
    /// [`CdotMapper::get_transcript_on_build`] when `build_hint` is
    /// `Some`, falling back to the primary-build
    /// [`CdotMapper::get_transcript`] otherwise. Returns
    /// [`FerroError::ReferenceNotFound`] if the requested build's view
    /// is absent from cdot (#389).
    fn cdot_tx_with_build_hint(
        &self,
        transcript_id: &str,
        build_hint: Option<&'static str>,
    ) -> Result<&CdotTranscript, FerroError> {
        match build_hint {
            Some(b) => self
                .projector
                .mapper()
                .cdot()
                .get_transcript_on_build(transcript_id, b),
            None => self.projector.mapper().cdot().get_transcript(transcript_id),
        }
        .ok_or_else(|| FerroError::ReferenceNotFound {
            id: transcript_id.to_string(),
        })
    }

    /// Extract the contig name and a representative 1-based genomic
    /// position from an already-normalized [`HgvsVariant`], for use as
    /// the input to `Projector::project`'s stab query.
    ///
    /// For [`HgvsVariant::Allele`] the first inner variant is used.
    ///
    /// # Per-axis behavior
    ///
    /// - **Genome**: contig is taken from the variant's accession
    ///   (`chromosome` for assembly-style refs, `full()` otherwise);
    ///   position is the start of the genomic interval.
    /// - **Cds / Tx / Rna**: cdot is consulted to resolve the
    ///   transcript's contig and to project the start position to
    ///   genome — `cds_to_genome_on_build` for c. (handles
    ///   intronic offsets, 5'UTR negative bases, and 3'UTR `utr3`
    ///   natively); `tx_to_genome` for n./r. simple exonic positions.
    ///   For n./r. positions that `tx_to_genome` cannot resolve
    ///   precisely (intronic offsets, downstream/utr3 markers,
    ///   non-positive bases) the function falls back to the
    ///   transcript's first-exon genomic start — the stab query still
    ///   lands inside the transcript, and the per-transcript
    ///   projection downstream re-derives the accurate coordinate. The
    ///   contig is taken from `cdot_tx.contig` (production cdot
    ///   typically keys this by the RefSeq genomic accession, with
    ///   `chr*` aliases handled by `resolve_contig` inside the stab
    ///   query).
    ///
    /// Removes the g.-only gate that pre-#389 defeated PR #379's
    /// fan-out widening. Build hints are inferred from the input's
    /// `genomic_context` parent (or assembly tag) the same way as the
    /// rest of the projector.
    fn extract_contig_and_pos(&self, variant: &HgvsVariant) -> Result<(String, u64), FerroError> {
        let effective = match variant {
            HgvsVariant::Allele(allele) => {
                allele
                    .variants
                    .first()
                    .ok_or_else(|| FerroError::UnsupportedProjection {
                        reason: "cannot project an empty allele to all transcripts".to_string(),
                    })?
            }
            other => other,
        };

        match effective {
            HgvsVariant::Genome(gv) => {
                let contig = if let Some(chr) = &gv.accession.chromosome {
                    chr.to_string()
                } else {
                    gv.accession.full()
                };
                let pos = gv
                    .loc_edit
                    .location
                    .start
                    .inner()
                    .cloned()
                    .ok_or_else(|| FerroError::InvalidCoordinates {
                        msg: "genomic interval start is unknown".to_string(),
                    })?
                    .base;
                Ok((contig, pos))
            }
            HgvsVariant::Cds(c) => {
                Self::require_parent_for_fanout(&c.accession, "c.")?;
                let transcript_id = c.accession.transcript_accession();
                let build_hint = self.build_hint_for_variant(effective);
                let cdot_tx = self.cdot_tx_with_build_hint(&transcript_id, build_hint)?;
                let start_cds = c.loc_edit.location.start.inner().cloned().ok_or_else(|| {
                    FerroError::InvalidCoordinates {
                        msg: "CDS interval start is unknown".to_string(),
                    }
                })?;
                let pos = match self.projector.mapper().cds_to_genome_on_build(
                    &transcript_id,
                    &start_cds,
                    build_hint,
                ) {
                    Ok(r) => r.variant.base,
                    Err(e) => {
                        // Poly-A 3'UTR fallback (#797): a `c.*` position in the
                        // post-transcriptional poly-A tail has no exon, so
                        // `cds_to_genome` declines and `project_variant_all` would
                        // otherwise fail to seed fan-out and miss the transcript.
                        // `try_extend_polya_to_genome` self-gates on a genuine
                        // `c.*` poly-A endpoint; when it confirms one, seed the
                        // stab query from an *in-transcript* anchor — the
                        // 3'-terminal exon's genomic start, which is inside the
                        // exon span on either strand — rather than the walked
                        // downstream coordinate (which is past the exon end by
                        // construction and would land the stab off the transcript).
                        // The per-transcript projection then takes the multi-axis
                        // short-circuit. Propagate the original decline when the
                        // position is not a poly-A endpoint.
                        if self
                            .try_extend_polya_to_genome(cdot_tx, &transcript_id, &start_cds)
                            .is_some()
                        {
                            let last = cdot_tx.exons.last().ok_or(e)?;
                            last[0]
                        } else {
                            return Err(e);
                        }
                    }
                };
                Ok((cdot_tx.contig.clone(), pos))
            }
            HgvsVariant::Tx(t) => {
                Self::require_parent_for_fanout(&t.accession, "n.")?;
                let transcript_id = t.accession.transcript_accession();
                let build_hint = self.build_hint_for_variant(effective);
                let cdot_tx = self.cdot_tx_with_build_hint(&transcript_id, build_hint)?;
                let start_tx = *t.loc_edit.location.start.inner().ok_or_else(|| {
                    FerroError::InvalidCoordinates {
                        msg: "Tx interval start is unknown".to_string(),
                    }
                })?;
                let pos =
                    nr_representative_genome_pos(cdot_tx, start_tx.base, start_tx.offset, false)?;
                Ok((cdot_tx.contig.clone(), pos))
            }
            HgvsVariant::Rna(r) => {
                Self::require_parent_for_fanout(&r.accession, "r.")?;
                let transcript_id = r.accession.transcript_accession();
                let build_hint = self.build_hint_for_variant(effective);
                let cdot_tx = self.cdot_tx_with_build_hint(&transcript_id, build_hint)?;
                let start_rna = *r.loc_edit.location.start.inner().ok_or_else(|| {
                    FerroError::InvalidCoordinates {
                        msg: "RNA interval start is unknown".to_string(),
                    }
                })?;
                let pos = nr_representative_genome_pos(
                    cdot_tx,
                    start_rna.base,
                    start_rna.offset,
                    start_rna.utr3,
                )?;
                Ok((cdot_tx.contig.clone(), pos))
            }
            HgvsVariant::Protein(_) => Err(FerroError::UnsupportedProjection {
                reason: "project_all does not accept protein (p.) inputs".to_string(),
            }),
            HgvsVariant::Mt(_) => Err(FerroError::UnsupportedProjection {
                reason: "project_all does not accept mitochondrial (m.) inputs".to_string(),
            }),
            HgvsVariant::Circular(_) => Err(FerroError::UnsupportedProjection {
                reason: "project_all does not accept circular (o.) inputs".to_string(),
            }),
            HgvsVariant::RnaFusion(_) => Err(FerroError::UnsupportedProjection {
                reason: "project_all does not accept RNA fusion inputs".to_string(),
            }),
            HgvsVariant::GenomeRing(_) => Err(FerroError::UnsupportedProjection {
                reason: "project_all does not accept genome ring inputs".to_string(),
            }),
            HgvsVariant::Supernumerary(_) => Err(FerroError::UnsupportedProjection {
                reason: "project_all does not accept supernumerary (sup) inputs".to_string(),
            }),
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => {
                Err(FerroError::UnsupportedProjection {
                    reason: "project_all does not accept null/unknown allele inputs".to_string(),
                })
            }
            // Allele was already unwrapped to its first inner variant above;
            // if we reach this arm with one it's because the inner variant
            // didn't match any of the supported axes (e.g. p./m./o.).
            HgvsVariant::Allele(_) => Err(FerroError::UnsupportedProjection {
                reason:
                    "project_all does not accept alleles whose first inner variant is not g./c./n./r."
                        .to_string(),
            }),
        }
    }

    fn cached_get_transcript_for_variant(
        &self,
        variant: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<Arc<Transcript>, FerroError> {
        let parent_key = Self::parent_key_for(variant);
        let key = (transcript_id.to_string(), parent_key);
        if let Some(tx) = self
            .transcript_cache
            .read()
            .expect("transcript cache poisoned")
            .get(&key)
        {
            return Ok(Arc::clone(tx));
        }
        let tx = match self.provider.get_transcript_for_variant(variant) {
            Ok(tx) => tx,
            Err(FerroError::ReferenceNotFound { .. }) => {
                self.provider.get_transcript(transcript_id)?
            }
            Err(e) => return Err(e),
        };
        let mut guard = self
            .transcript_cache
            .write()
            .expect("transcript cache poisoned");
        // `tx` is already an `Arc<Transcript>` from the provider; insert directly.
        let entry = guard.entry(key).or_insert_with(|| tx);
        Ok(Arc::clone(entry))
    }

    /// Build (or retrieve) the cached `RefProteinBundle` for a transcript.
    ///
    /// The reference CDS and its translation are stable for a given
    /// `(transcript_id, parent_accession)` pair — `predict_indel_protein`
    /// would otherwise recompute them on every variant. Caching collapses
    /// that to one translation per (id, parent) over the projector's
    /// lifetime.
    ///
    /// The cache is keyed by the same parent-aware identity as
    /// [`Self::cached_get_transcript_for_variant`] so that an NG/NC-parented
    /// input does not reuse the bundle built for a different cdot build of
    /// the same `transcript_id` (issue #332).
    pub(crate) fn cached_ref_translation(
        &self,
        variant: &HgvsVariant,
        transcript_id: &str,
        transcript: &Transcript,
    ) -> Result<Arc<RefProteinBundle>, FerroError> {
        let parent_key = Self::parent_key_for(variant);
        let key = (transcript_id.to_string(), parent_key);
        if let Some(b) = self
            .ref_protein_cache
            .read()
            .expect("ref-protein cache poisoned")
            .get(&key)
        {
            return Ok(Arc::clone(b));
        }
        let bundle = RefProteinBundle::from_transcript(transcript)?;
        let mut guard = self
            .ref_protein_cache
            .write()
            .expect("ref-protein cache poisoned");
        let entry = guard.entry(key).or_insert_with(|| Arc::new(bundle));
        Ok(Arc::clone(entry))
    }

    pub fn with_normalize_config(mut self, config: NormalizeConfig) -> Self {
        self.normalizer = Normalizer::with_config(self.provider.clone(), config);
        self
    }

    /// Set an explicit genome-build override for build-agnostic inputs (#715).
    ///
    /// `build` must be a canonical ferro build string (`"GRCh37"` / `"GRCh38"`);
    /// validate/normalize user-supplied names at the boundary with
    /// [`crate::liftover::aliases::normalize_assembly_name`]. The override only
    /// *fills in* a build for inputs that don't encode one (a bare `NG_`/`LRG_`);
    /// an input whose accession already specifies a build keeps it, and a
    /// contradictory override is warned and ignored (see
    /// [`Self::build_hint_for_variant`]). `None` clears the override.
    pub fn with_assembly(mut self, build: Option<&'static str>) -> Self {
        self.assembly_override = build;
        self
    }

    /// Parse, normalize, and project an HGVS string onto a transcript.
    pub fn project(
        &self,
        hgvs_string: &str,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        let variant = crate::parse_hgvs(hgvs_string)?;
        self.project_variant(&variant, transcript_id)
    }

    /// Normalize and project an already-parsed g. variant onto a transcript.
    ///
    /// The variant is normalized first; for pre-normalized variants use
    /// [`project_normalized`] to skip the redundant normalization step.
    pub fn project_variant(
        &self,
        variant: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        self.warn_assembly_conflict(variant);
        // 1. Normalize the genomic variant. The normalizer is built once at
        // construction time so we don't clone the (potentially heavy) provider
        // on every call.
        let normalized = self.normalizer.normalize(variant)?;
        // 2. Reconcile the requested target against any selector resolution that
        // normalization performed (#783). A legacy LOVD gene-model selector
        // `NG_(GENE_v001):c.` keeps the genomic ref as its accession, so a
        // self-projection caller (the CLI, the conformance harness) derives the
        // target from `transcript_accession()` and gets the genomic ref, not a
        // transcript. Normalization resolves the selector to `NG_(NM_):c.` (#637),
        // so the requested target now names the variant's genomic-context wrapper
        // rather than its transcript — re-point it at the resolved transcript.
        let target = Self::reconcile_self_projection_target(&normalized, transcript_id);
        self.project_variant_inner(&normalized, &target)
    }

    /// Reconcile a self-projection target with a selector/context that
    /// normalization may have rewritten (#783).
    ///
    /// Returns the variant's own (resolved) transcript accession when the
    /// requested target names this variant's genomic-context wrapper rather than
    /// its transcript — i.e. the caller asked to project a transcript-coordinate
    /// variant "against itself" but derived the target from an accession whose
    /// transcript slot was a genomic ref (an unresolved legacy gene-model
    /// selector, now resolved). Otherwise returns the requested target unchanged,
    /// so a genuine transcript_id mismatch (an unrelated target) still errors.
    fn reconcile_self_projection_target(normalized: &HgvsVariant, requested: &str) -> String {
        if let Some(acc) = normalized.accession() {
            let resolved_tx = acc.transcript_accession();
            if resolved_tx != requested {
                if let Some(ctx) = &acc.genomic_context {
                    if ctx.full() == requested {
                        return resolved_tx;
                    }
                }
            }
        }
        requested.to_string()
    }

    /// Project an already-normalized g. variant onto a transcript, skipping the
    /// normalization step.
    ///
    /// Callers that pre-normalize once and then project against many transcripts
    /// should use this method to avoid the cost of re-normalization.
    ///
    /// **Warning**: passing a non-normalized variant will produce coordinates
    /// that are technically valid but may not match other tools' canonical form.
    pub fn project_normalized(
        &self,
        variant: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        self.warn_assembly_conflict(variant);
        // Reconcile the requested target against any selector resolution carried
        // by the (already-normalized) variant, mirroring `project_variant` (#783).
        // A caller that pre-normalizes a legacy `NG_(GENE_v001):c.` once and then
        // projects against the raw-derived `NG_900.1` target would otherwise name
        // the genomic-context wrapper rather than the resolved transcript and trip
        // the mismatch guard. Reconciliation is a no-op for a genuine transcript
        // target, so this keeps the two public entry points behaving identically.
        let target = Self::reconcile_self_projection_target(variant, transcript_id);
        self.project_variant_inner(variant, &target)
    }

    /// Parse, normalize, and project an HGVS string onto the *curated* set of
    /// overlapping transcripts (#656), returning results in clinical priority
    /// order (MANE Select first, then Plus Clinical, then canonical, then
    /// longest CDS).
    ///
    /// The result is the curated enumerated set rather than every
    /// cdot-overlapping record: superseded transcript versions are collapsed
    /// (only the highest version per base accession is kept) and predicted
    /// `XM_`/`XR_` models are dropped when a curated `NM_`/`NR_` transcript
    /// covers the same locus; predicted models are kept only when they are the
    /// sole coverage. See [`select_enumerated_transcript_ids`] for the policy.
    ///
    /// Returns an empty `Vec` when the variant overlaps no known transcripts.
    /// Individual transcript errors are logged at trace level and silently
    /// skipped so that a single bad transcript does not abort the whole call.
    pub fn project_all(&self, hgvs_string: &str) -> Result<Vec<VariantProjection>, FerroError> {
        let variant = crate::parse_hgvs(hgvs_string)?;
        self.project_variant_all(&variant)
    }

    /// Normalize and project an already-parsed g. variant onto the *curated*
    /// set of overlapping transcripts (#656).
    ///
    /// See [`project_all`] for the curated enumeration policy, ordering, and
    /// error-handling semantics.
    pub fn project_variant_all(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Vec<VariantProjection>, FerroError> {
        self.warn_assembly_conflict(variant);
        // 0. An NG_/LRG_ genomic *input* carries parent-relative coordinates that
        //    cdot cannot map; rewrite it into the chromosome (NC_) frame so the
        //    fan-out below can enumerate and project onto the overlapping
        //    transcripts (#480 inverse). Capture the build the *original* input
        //    carries before `deanchor` shadows `variant` — the de-anchored form is
        //    on an `NC_` accession whose build would be circular to read back
        //    (#653/#713); the parent placement must be selected by the input's
        //    build, matching what `deanchor` itself used.
        let input_build = self.build_hint_for_variant(variant);
        let (variant, parent) = self.deanchor_genomic_parent_input(variant);
        // 1. Normalize once in the genome frame, then fan out across the
        //    overlapping transcripts.
        let normalized = self.normalizer.normalize(&variant)?;
        let results = self.project_normalized_all_inner(&normalized)?;
        // 2. The plain (non-parent) path returns the genome-frame result directly.
        let Some(parent) = parent else {
            return Ok(results);
        };
        // 3. Parent path: the fan-out 3'-shifts the variant in the *genome*
        //    frame, so a repeat-region edit can land at a transcript coordinate
        //    two bases off from the HGVS-correct, transcript-frame-normalized
        //    position. Re-project each transcript from its own coding form (which
        //    normalizes in the transcript frame), then re-frame under the parent.
        let placement = self
            .provider
            .genomic_placement_on_build(&parent, input_build);
        let mut framed = Vec::with_capacity(results.len());
        for r in results {
            let r = match r.coding.as_ref() {
                // Re-project from the coding form to pick up transcript-frame
                // normalization. Tolerate a failure the same way the fan-out in
                // `project_normalized_all` tolerates a failing transcript — but
                // log it rather than swallow it silently. The genome-frame `r` we
                // fall back to is itself a valid projection (only a repeat-region
                // edit can be off — by up to the repeat-unit length — from the
                // transcript-frame-optimal position); emitting it is preferable to
                // dropping the transcript outright.
                Some(coding) => match self.project_variant(coding, &r.transcript_id) {
                    Ok(reframed) => reframed,
                    Err(e) => {
                        log::trace!(
                            "project_variant_all: transcript-frame re-projection failed for {}: \
                             {}; keeping genome-frame projection",
                            r.transcript_id,
                            e
                        );
                        r
                    }
                },
                None => r,
            };
            framed.push(self.frame_projection_owned(r, &parent, placement.as_ref()));
        }
        Ok(framed)
    }

    /// Re-frame a [`VariantProjection`] produced from a de-anchored `NC_` variant
    /// back under the original `NG_`/`LRG_` parent (#480): the transcript-relative
    /// descriptions (`coding`/`noncoding`/`protein`) gain the parent as their
    /// `genomic_context` (rendering `NG_(NM_)` / `NG_(NP_)`), and the genomic
    /// description is re-anchored into the parent's own coordinate frame.
    ///
    /// The synthesized gene-symbol selector is dropped: per #121 ferro does not
    /// emit a selector that was not in the input, and the mutalyzer `NG_(NM_)`
    /// form carries none.
    fn frame_projection_owned(
        &self,
        mut result: VariantProjection,
        parent: &Accession,
        placement: Option<&crate::reference::GenomicPlacement>,
    ) -> VariantProjection {
        for field in [
            &mut result.coding,
            &mut result.noncoding,
            &mut result.protein,
        ] {
            if let Some(v) = field.as_ref() {
                *field = Some(Self::relabel_under_parent(v, parent));
            }
        }
        // #860: under a bare LRG genomic parent, echo the input's LRG namespace
        // instead of the resolved NM_/NP_. The transcript index (LRG_<n>t<k>)
        // comes from the resolved NM_ of the coding/noncoding axis via the
        // parent-scoped reverse lookup; the protein is the structural t->p of the
        // same index. Render bare — clear the genomic_context the loop above
        // attached, since the corpus expects `LRG_24t1:c.`, not `LRG_24(LRG_24t1):c.`.
        if parent.is_lrg() {
            let cdot = self.projector.mapper().cdot();
            let base = parent.base();
            let lrg_tx = result
                .coding
                .as_ref()
                .or(result.noncoding.as_ref())
                .and_then(|v| v.accession())
                .map(|a| a.transcript_accession())
                .and_then(|nm| cdot.lrg_transcript_for_parent(&base, &nm));
            if let Some(lrg_tx) = lrg_tx {
                let lrg_tx_acc = parse_accession(&lrg_tx);
                // The rna axis can be an allele needing apply_rna_framing-style
                // recursion; no #860 corpus row exercises rna under an LRG parent,
                // so it is intentionally out of scope here.
                for field in [&mut result.coding, &mut result.noncoding] {
                    if let Some(v) = field.as_mut() {
                        Self::set_bare_accession(v, lrg_tx_acc.clone());
                    }
                }
                if let (Some(protein), Some(lrg_p)) =
                    (result.protein.as_mut(), lrg_protein_accession(&lrg_tx))
                {
                    Self::set_bare_accession(protein, parse_accession(&lrg_p));
                }
            }
        }
        // The RNA axis can be a predicted allele (`r.[...]`), so frame it through
        // `apply_rna_framing` — which recurses into allele members — rather than
        // the single-variant `relabel_under_parent` used above. Drop any
        // synthesized gene symbol (`None`), matching the parent-framing policy the
        // other axes follow here (#121/#693).
        if let Some(rna) = result.rna.as_mut() {
            Self::apply_rna_framing(rna, Some(parent), None);
        }
        if let (Some(HgvsVariant::Genome(gv)), Some(placement)) =
            (result.genomic.as_ref(), placement)
        {
            match Self::reanchor_genome_to_parent(gv.clone(), placement) {
                Ok(mut reanchored) => {
                    reanchored.accession = parent.clone();
                    reanchored.gene_symbol = None;
                    result.genomic = Some(HgvsVariant::Genome(reanchored));
                }
                Err(e) => {
                    // The genomic axis cannot be re-anchored into the parent
                    // frame (endpoint outside the placed span or an uncertain
                    // position). Drop the unframable genomic axis rather than
                    // stamp the parent accession on a chromosome coordinate,
                    // which would be invalid HGVS (#655). The transcript-relative
                    // coding/noncoding/protein forms re-labeled above remain
                    // valid and are preferable to dropping the whole projection
                    // (mirrors the tolerate-and-log fallback in
                    // `project_variant_all`).
                    log::trace!(
                        "frame_projection_owned: cannot re-anchor the genomic axis into the {} \
                         parent frame: {}; reporting genomic = None",
                        parent,
                        e
                    );
                    result.genomic = None;
                }
            }
        }
        result
    }

    /// Return a clone of `variant` whose top-level accession carries `parent` as
    /// its `genomic_context` and whose synthesized gene-symbol selector is
    /// cleared. Used to frame a transcript-relative description under its
    /// `NG_`/`LRG_` parent (#480).
    fn relabel_under_parent(variant: &HgvsVariant, parent: &Accession) -> HgvsVariant {
        let mut v = variant.clone();
        match &mut v {
            HgvsVariant::Cds(c) => {
                c.accession = c.accession.clone().with_genomic_context(parent.clone());
                c.gene_symbol = None;
            }
            HgvsVariant::Tx(t) => {
                t.accession = t.accession.clone().with_genomic_context(parent.clone());
                t.gene_symbol = None;
            }
            HgvsVariant::Rna(r) => {
                r.accession = r.accession.clone().with_genomic_context(parent.clone());
                r.gene_symbol = None;
            }
            HgvsVariant::Protein(p) => {
                p.accession = p.accession.clone().with_genomic_context(parent.clone());
                p.gene_symbol = None;
            }
            _ => {}
        }
        v
    }

    /// Set `variant`'s top-level accession to `acc` with no `genomic_context`
    /// (bare rendering), clearing any synthesized gene symbol. Used by #860 to
    /// echo the input's LRG namespace (`LRG_<n>t<k>`/`LRG_<n>p<k>`) on a fan-out
    /// projection without the `LRG_<n>(inner)` genomic-context wrapper. Only the
    /// single-variant arms are handled (an allele rna axis is out of scope).
    fn set_bare_accession(variant: &mut HgvsVariant, acc: Accession) {
        let acc = acc.without_genomic_context();
        match variant {
            HgvsVariant::Cds(c) => {
                c.accession = acc;
                c.gene_symbol = None;
            }
            HgvsVariant::Tx(t) => {
                t.accession = acc;
                t.gene_symbol = None;
            }
            HgvsVariant::Rna(r) => {
                r.accession = acc;
                r.gene_symbol = None;
            }
            HgvsVariant::Protein(p) => {
                p.accession = acc;
                p.gene_symbol = None;
            }
            _ => {}
        }
    }

    /// Re-frame a predicted RNA (`r.`) consequence so its accession matches the
    /// *input* variant's reference framing (#693).
    ///
    /// [`crate::project::rna::predict_rna`] copies the accession and gene symbol
    /// from the synthesized `coding` form it is handed, which on the
    /// genome-pivot / direct-`n.` paths is the bare transcript accession with a
    /// cdot-looked-up gene symbol (e.g. `NM_058195.3(CDKN2A)`). The single
    /// projection path (`project_single_inner`) — unlike `project_variant_all`'s
    /// `frame_projection_owned` — never re-labels its output under the input's
    /// parent, so the predicted `r.` dropped the input's `NG_(NM_)` genomic
    /// context and substituted a synthesized gene symbol.
    ///
    /// Carry the input's framing onto the `r.` axis: set the input's
    /// `genomic_context` (if any) on the RNA accession, and adopt the input's
    /// gene symbol (`None` when the input carried none) so a symbol that was not
    /// in the input is not fabricated (#121). Only the `r.` axis is touched;
    /// `coding`/`protein` keep their existing (bare) single-path rendering.
    /// A `None` prediction (axis unavailable) passes through unchanged.
    fn reframe_rna_from_input(
        rna: Option<HgvsVariant>,
        input: &HgvsVariant,
    ) -> Option<HgvsVariant> {
        let mut rna = rna?;
        let context = input
            .accession()
            .and_then(|a| a.genomic_context.as_deref().cloned());
        // `gene_symbol()` returns `None` for an `Allele`, so an allele input that
        // carries a user-provided selector on its members (e.g.
        // `NM_x(GENE):c.[...]`) would otherwise lose it. Derive the symbol from
        // the members when they unanimously agree, so the predicted `r.` allele
        // keeps the input's selector instead of stripping it (#693).
        let gene_symbol = match input {
            HgvsVariant::Allele(a) => a.variants.first().and_then(|first| {
                let symbol = first.gene_symbol()?;
                a.variants
                    .iter()
                    .all(|member| member.gene_symbol() == Some(symbol))
                    .then_some(symbol)
            }),
            _ => input.gene_symbol(),
        }
        .map(str::to_string);
        Self::apply_rna_framing(&mut rna, context.as_ref(), gene_symbol.as_deref());
        Some(rna)
    }

    /// Stamp `context`/`gene_symbol` onto every `r.` leaf of `rna` (descending
    /// into a predicted allele so its first-member-derived prefix is re-framed).
    fn apply_rna_framing(
        rna: &mut HgvsVariant,
        context: Option<&Accession>,
        gene_symbol: Option<&str>,
    ) {
        match rna {
            HgvsVariant::Rna(r) => {
                r.accession = match context {
                    Some(parent) => r.accession.clone().with_genomic_context(parent.clone()),
                    // Strip any genomic context the predicted form may carry when
                    // the input had none, mirroring the input's bare framing.
                    None => r.accession.clone().without_genomic_context(),
                };
                r.gene_symbol = gene_symbol.map(str::to_string);
            }
            // A predicted RNA allele renders its accession prefix from the first
            // member (compact form), so every member must carry the input frame.
            HgvsVariant::Allele(a) => {
                for member in &mut a.variants {
                    Self::apply_rna_framing(member, context, gene_symbol);
                }
            }
            _ => {}
        }
    }

    /// Derive the genomic LRG parent (`LRG_<n>`) of an LRG transcript or protein
    /// accession (`LRG_<n>t<m>` / `LRG_<n>p<m>`).
    ///
    /// LRG accessions embed the genomic number plus a `t`/`p` suffix, so the
    /// genomic reference is the leading digits — a structural fact of the
    /// accession, not a cdot inference. Returns `None` for non-LRG accessions
    /// and for a bare genomic LRG (all digits, no suffix), which is already its
    /// own genome reference.
    fn lrg_genomic_parent(accession: &Accession) -> Option<Accession> {
        if !accession.is_lrg() {
            return None;
        }
        let number = &*accession.number;
        let genomic_digits: String = number.chars().take_while(|c| c.is_ascii_digit()).collect();
        if genomic_digits.is_empty() || genomic_digits.len() == number.len() {
            return None;
        }
        Some(Accession::new("LRG", genomic_digits, None))
    }

    /// Project a transcript-coordinate variant (`c.`/`n.`/`r.`) onto its parent
    /// genomic reference and return an [`HgvsVariant::Genome`].
    ///
    /// The output uses the parent accession stored in the input's
    /// `Accession.genomic_context`:
    /// - an `NC_` chromosome parent is already in the genome frame and passes
    ///   through unchanged;
    /// - an `NG_` RefSeqGene or `LRG_` parent whose chromosomal placement is
    ///   known (see [`ReferenceProvider::genomic_placement`]) is re-anchored
    ///   into the parent's own frame (#480);
    /// - an `NG_`/`LRG_` parent with **no** known placement, or whose endpoint
    ///   cannot be re-anchored, is declined (see `# Errors`) rather than
    ///   emitting chromosome (`NC_`) coordinates under the parent accession,
    ///   which would be invalid HGVS (#655).
    ///
    /// Idempotent on `Genome` input.
    ///
    /// This is the genomic-axis *output*. The internal pivot used to derive
    /// coding/protein keeps the chromosome frame — see
    /// [`Self::project_to_genomic_nc`].
    ///
    /// # Errors
    ///
    /// Returns [`FerroError::UnsupportedProjection`] when an `NG_`/`LRG_` parent
    /// has no known chromosomal placement or an endpoint falls outside the placed
    /// span, and [`FerroError::InvalidCoordinates`] when an endpoint has no single
    /// resolved position (uncertain/compound boundary) (#655). See
    /// [`Self::project_to_genomic_nc`] for the rest of the error contract
    /// (`UnsupportedProjection` / `InvalidCoordinates` / `ReferenceNotFound`).
    pub fn project_to_genomic(&self, variant: &HgvsVariant) -> Result<HgvsVariant, FerroError> {
        self.warn_assembly_conflict(variant);
        // A `Genome` input is already in its parent frame: an `NC_` chromosome
        // variant or an `NG_`/`LRG_` genomic variant in its own coordinates.
        // `reanchor_genome_output` assumes an `NC_`-frame pivot produced by
        // `project_to_genomic_nc`, so feeding an already-parent-framed `NG_`/`LRG_`
        // variant through it would either double-transform (mistaking the parent
        // coordinate for an `NC_` coordinate) or decline outright. Pass `Genome`
        // inputs through unchanged to preserve idempotence (#702).
        if matches!(variant, HgvsVariant::Genome(_)) {
            return Ok(variant.clone());
        }

        // #851: project a compound allele by projecting each member into the
        // parent genomic frame and reassembling. Each member carries its own
        // accession + genomic_context, so per-member `project_to_genomic`
        // re-anchors and re-normalizes it (the `project_to_genomic` layer, not
        // `project_to_genomic_nc`). All-or-nothing: a declining member errors the
        // whole allele (mirrors `project_allele_inner`'s `all_have_genomic`).
        if let HgvsVariant::Allele(allele) = variant {
            let mut members = Vec::with_capacity(allele.variants.len());
            for member in &allele.variants {
                members.push(self.project_to_genomic(member)?);
            }
            // Genomic-order sort applies only to a CIS allele — a trans allele
            // `[m1];[m2]` encodes haplotype assignment by member order, so
            // reordering would swap haplotypes (alleles.md:118 is about variants
            // within one bracket/chromosome).
            if allele.phase == AllelePhase::Cis {
                sort_genomic_allele_members(&mut members);
            }
            let mut projected = AlleleVariant::new(members, allele.phase);
            projected.uncertain = allele.uncertain;
            return Ok(HgvsVariant::Allele(projected));
        }

        // #537: a `c.pter`/`c.qter` telomere-flank input resolves to the parent
        // reference's own genomic terminus, not a transcript-mapped coordinate.
        // Handle it before the cdot pivot below, which cannot (and per #534
        // must not) number a transcript-flank marker.
        if let Some(result) = self.project_cds_terminus_to_parent(variant) {
            return result;
        }

        match self.project_to_genomic_nc(variant)? {
            HgvsVariant::Genome(gv) => {
                self.reanchor_and_normalize_genomic(gv, self.build_hint_for_variant(variant))
            }
            other => Ok(other),
        }
    }

    /// Re-anchor an NC-frame pivot into its parent frame, then re-normalize in the
    /// parent's own frame (#737), returning the spec-canonical genomic form.
    ///
    /// On normalize error, emit the un-normalized re-anchored form (the dropped
    /// error is `trace!`-logged so a regression is observable rather than silently
    /// masked — the fallback swallows *any* normalize error, not only missing
    /// bases). The reanchor decline is returned as `Err`; the caller chooses
    /// whether to propagate it (`?`, public `project_to_genomic` contract) or
    /// degrade the genomic axis to `None` (`.ok()`, the multi-axis paths, #655/#702).
    ///
    /// Why renormalize: the pivot was 3'-normalized in transcript space, which is
    /// the genome's 5' end on a minus-strand transcript, so a shiftable del/dup/ins
    /// (or repeat) must be re-shuffled to its genomic-3' position against the output
    /// accession's own sequence.
    fn reanchor_and_normalize_genomic(
        &self,
        gv: crate::hgvs::variant::GenomeVariant,
        build: Option<&str>,
    ) -> Result<HgvsVariant, FerroError> {
        let reanchored = self.reanchor_genome_output(gv, build)?;
        Ok(self.normalizer.normalize(&reanchored).unwrap_or_else(|e| {
            log::trace!(
                "genomic-frame renormalization failed for {reanchored}; \
                 emitting un-normalized re-anchored form: {e}"
            );
            reanchored
        }))
    }

    /// Project a `c.pter`/`c.qter` (telomere-flank marker) coding input onto its
    /// parent reference's own genomic frame (#537).
    ///
    /// `pter` denotes the 5'-most position of the parent reference (`g.1`) and
    /// `qter` the 3'-most (`g.<length>`); a `pter_qter` range spans the whole
    /// reference. Unlike a numeric `c.` coordinate these markers do not number a
    /// transcript position — PR #534 correctly refuses to do so on the `c.`
    /// axis — so they map straight to the parent reference's termini with no
    /// cdot transcript mapping. The emitted `g.1` / `g.<length>` is normalized
    /// downstream (e.g. the 3' rule rolls `g.1del` through a leading homopolymer
    /// run, matching mutalyzer's `g.3del`).
    ///
    /// Returns `None` when `variant` is not a handled terminus case — not a
    /// `Cds` input, an endpoint that is not a `pter`/`qter` marker (numeric,
    /// `?`, `cen`, or a mixed marker/numeric range), or a reversed `qter_pter` —
    /// so the caller falls through to the normal transcript-mapped projection
    /// (which declines or numbers as appropriate). Returns `Some(Err(..))` only
    /// for a *supported* terminus whose parent reference is unresolvable or whose
    /// length is unknown to the provider (e.g. the pinned parent version is
    /// absent — a reference-coverage gap, #645/#672 — not a projection bug); the
    /// marker-pair classification runs *before* the parent-length lookup, so an
    /// unsupported pair always declines with `None` regardless of whether the
    /// parent length could be resolved.
    fn project_cds_terminus_to_parent(
        &self,
        variant: &HgvsVariant,
    ) -> Option<Result<HgvsVariant, FerroError>> {
        use crate::hgvs::interval::GenomeInterval;
        use crate::hgvs::location::{GenomePos, SpecialPosition};
        use crate::hgvs::variant::{GenomeVariant, LocEdit};

        let HgvsVariant::Cds(v) = variant else {
            return None;
        };
        // Both endpoints must be telomere markers. A numeric or `?` endpoint, or
        // a mixed marker/numeric range (e.g. `c.pter_-51`), needs the transcript
        // mapping and is left to the normal path.
        let start = resolve_uncertain_boundary(&v.loc_edit.location.start, "c.", "start").ok()?;
        let end = resolve_uncertain_boundary(&v.loc_edit.location.end, "c.", "end").ok()?;
        let (Some(start_marker), Some(end_marker)) = (start.special, end.special) else {
            return None;
        };

        // Classify the whole-arm terminus BEFORE resolving the parent length, so
        // an unsupported marker pair (`cen`, a reversed `qter_pter`, a mixed
        // `pter_cen`, …) declines here and falls through to the normal rejection
        // path — rather than surfacing an incidental parent-length-lookup error,
        // which is only meaningful for a supported terminus. `pter` → 5'-most
        // (g.1), `qter` → 3'-most (g.<length>), `pter_qter` → the whole reference.
        enum ArmTerminus {
            Pter,
            Qter,
            Whole,
        }
        let terminus = match (start_marker, end_marker) {
            (SpecialPosition::Pter, SpecialPosition::Pter) => ArmTerminus::Pter,
            (SpecialPosition::Qter, SpecialPosition::Qter) => ArmTerminus::Qter,
            (SpecialPosition::Pter, SpecialPosition::Qter) => ArmTerminus::Whole,
            _ => return None,
        };

        // Resolve the parent genomic reference: an explicit `genomic_context`
        // parent, or the structural `LRG_<n>` parent of a bare LRG transcript
        // (#480). No parent → leave it to the normal path, which raises the
        // canonical "no parent reference" error.
        let parent = match v.accession.genomic_context.as_deref().cloned() {
            Some(p) => p,
            None => Self::lrg_genomic_parent(&v.accession)?,
        };
        let length = match self.provider.get_sequence_length(&parent.full()) {
            Ok(n) => n,
            Err(e) => return Some(Err(e)),
        };
        // A zero-length parent is unreachable for any real NG/LRG/NC contig, but
        // guard the endpoint math against it (#526): `qter` would build an
        // invalid 1-based `g.0` and `pter_qter` a reversed `g.1_0`. Decline so
        // the normal path raises the canonical refusal rather than emitting a
        // malformed interval.
        if length == 0 {
            return None;
        }

        let (g_start, g_end) = match terminus {
            ArmTerminus::Pter => (1, 1),
            ArmTerminus::Qter => (length, length),
            ArmTerminus::Whole => (1, length),
        };

        let interval = GenomeInterval::new(GenomePos::new(g_start), GenomePos::new(g_end));
        let gv = GenomeVariant {
            accession: parent,
            gene_symbol: None,
            loc_edit: LocEdit::with_uncertainty(interval, v.loc_edit.edit.clone()),
        };
        Some(Ok(HgvsVariant::Genome(gv)))
    }

    /// Re-anchor a chromosome-frame genomic projection into its parent's own
    /// frame for *output* (#480), or decline.
    ///
    /// `gv` is a genome variant stamped with the parent accession from the
    /// projected transcript's `genomic_context` (as produced by
    /// [`Self::project_to_genomic_nc`]). An `NC_` chromosome parent is already in
    /// the genome frame and is returned unchanged; an `NG_`/`LRG_` parent is
    /// re-anchored into its own frame via
    /// [`ReferenceProvider::genomic_placement_on_build`].
    ///
    /// `build` is the genome build the *original input* carries (from
    /// [`Self::build_hint_for_variant`]); it selects the parent's build-appropriate
    /// placement (#653/#713). `None` (a bare `NG_` with no build signal) keeps the
    /// GRCh38-preferred behavior.
    ///
    /// Shared by the public [`Self::project_to_genomic`] (which propagates the
    /// decline as an error) and the multi-axis single-transcript path in
    /// [`Self::project_single_inner`] (which degrades the genomic axis to `None`
    /// on decline, keeping the still-valid coding/protein axes — #702).
    ///
    /// # Errors
    ///
    /// Returns [`FerroError::UnsupportedProjection`] for an `NG_`/`LRG_` parent
    /// with no known placement or an endpoint outside the placed span, and
    /// [`FerroError::InvalidCoordinates`] for an uncertain/compound endpoint
    /// (#655) — rather than emit a chromosome coordinate under the parent
    /// accession (invalid HGVS).
    fn reanchor_genome_output(
        &self,
        gv: crate::hgvs::variant::GenomeVariant,
        build: Option<&str>,
    ) -> Result<HgvsVariant, FerroError> {
        let reanchored = match self
            .provider
            .genomic_placement_on_build(&gv.accession, build)
        {
            Some(placement) => Self::reanchor_genome_to_parent(gv, &placement)?,
            // No placement: an `NC_` chromosome parent is already in the genome
            // frame and passes through unchanged, but an `NG_`/`LRG_` parent
            // cannot be re-anchored (cdot carries only the transcript's `NC_`
            // alignment). Emitting the `NC_` coordinate under the `NG_`/`LRG_`
            // accession would be invalid HGVS, so decline instead (#480/#655).
            None if &*gv.accession.prefix == "NG" || gv.accession.is_lrg() => {
                return Err(FerroError::UnsupportedProjection {
                    reason: format!(
                        "cannot project to the {} parent frame: no chromosomal \
                         placement is known for this NG_/LRG_ reference, so cdot's \
                         chromosome (NC_) coordinates cannot be re-anchored into the \
                         parent's own frame (#480/#655)",
                        gv.accession.transcript_accession(),
                    ),
                });
            }
            None => gv,
        };
        Ok(HgvsVariant::Genome(reanchored))
    }

    /// Project a transcript-coordinate variant onto its parent genomic reference
    /// in the **chromosome (`NC_`) frame** — the pivot for downstream
    /// genome→CDS/protein derivation, which cdot computes in chromosome
    /// coordinates. The public [`Self::project_to_genomic`] wraps this and
    /// re-anchors the *output* into an `NG_`/`LRG_` parent's own frame (#480).
    ///
    /// The output stamps the parent accession from `Accession.genomic_context`.
    ///
    /// # Errors
    ///
    /// Returns [`FerroError::UnsupportedProjection`] for:
    /// - `p.` / `m.` / `o.` / RNA-fusion / `Allele` / `NullAllele` / `UnknownAllele`
    ///   inputs (see #328 for allele support),
    /// - transcript-coordinate inputs whose `Accession.genomic_context` is
    ///   absent (no parent NG/NC reference) — except a bare LRG transcript,
    ///   whose `LRG_<n>` parent is derived structurally (#480),
    /// - `?` position sentinels in the start or end coordinate.
    ///
    /// Returns [`FerroError::InvalidCoordinates`] when an endpoint is a
    /// compound `(a_b)` `UncertainBoundary::Range` (no single position).
    ///
    /// Returns [`FerroError::ReferenceNotFound`] when the named transcript
    /// is absent from the cdot mapper.
    ///
    /// The output ref-base is **not** validated against the genomic reference;
    /// mismatches surface downstream in `Normalizer`.
    ///
    /// `r.` inputs carrying `Base::U` are translated to DNA on the way to the
    /// g. output (`U`→`T` on the plus strand, `U`→`A` via complement on the
    /// minus strand), so the emitted g. variant is always valid DNA (#395 item 4).
    /// Whether either endpoint of a `c.` input is a confirmed poly-A 3'UTR
    /// position (#797) — i.e. lands in the transcript's post-transcriptional
    /// poly-A tail, where it has no genomic alignment and is mapped only by the
    /// contiguous downstream walk of [`Self::try_extend_polya_to_genome`].
    ///
    /// Used by [`Self::deanchor_genomic_parent_input`] to keep an NG_/LRG_-parented
    /// poly-A `c.*` input in its transcript-coordinate form rather than de-anchoring
    /// it to a synthetic (off-exon) `NC_` `Genome`: the `c.` fan-out seed path
    /// (`extract_contig_and_pos`) handles poly-A endpoints by seeding the stab from
    /// the 3'-terminal exon anchor, whereas the off-exon `Genome` coordinate would
    /// land the stab off the transcript and fail enumeration. Resolves the
    /// transcript via the input's build and returns `false` (decline) when the
    /// transcript or an endpoint cannot be resolved.
    fn cds_input_has_polya_endpoint(&self, variant: &HgvsVariant) -> bool {
        let HgvsVariant::Cds(c) = variant else {
            return false;
        };
        let transcript_id = c.accession.transcript_accession();
        let build_hint = self.build_hint_for_variant(variant);
        let Ok(cdot_tx) = self.cdot_tx_with_build_hint(&transcript_id, build_hint) else {
            return false;
        };
        [
            c.loc_edit.location.start.inner(),
            c.loc_edit.location.end.inner(),
        ]
        .into_iter()
        .flatten()
        .any(|p| {
            self.try_extend_polya_to_genome(cdot_tx, &transcript_id, p)
                .is_some()
        })
    }

    /// Map a `c.*` 3'UTR position that falls in the transcript's poly-A tail to
    /// the contiguous downstream genome coordinate (#797), or `None` to decline.
    ///
    /// A transcript's poly-A tail is post-transcriptional and has no genomic
    /// coordinate, so a derived exon→genome structure (#790) strips it and the
    /// CDS-aware mapper declines a `c.*` position landing there. mutalyzer maps
    /// such a position by a base-agnostic coordinate walk into downstream genome
    /// from the transcript's own 3'-terminal exon; this performs that walk via
    /// [`CdotTranscript::tx_to_genome_extending_polya`].
    ///
    /// Returns `None` (so the caller propagates the original decline) unless ALL
    /// hold, so normal projection is never perturbed and no wrong coordinate is
    /// emitted:
    /// - the position is a `c.*` 3'UTR position (`utr3`, no `special` sentinel),
    ///   optionally with a *positive* intronic-style `offset` (`c.*N+k`);
    /// - the base position `c.*N` itself lands in the poly-A region (past the
    ///   3'-terminal exon) — never a real intron between exons;
    /// - the transcript has a `cds_end` (needed to resolve `c.*N`);
    /// - the effective transcript position (`c.*N` plus any offset) is
    ///   `< true_tx_length` — i.e. it lies inside the transcript's real poly-A
    ///   tail, not past the transcript end. The true length comes from the
    ///   transcript FASTA, because the derived exon map's length excludes the
    ///   stripped tail.
    ///
    /// A `c.*N+k` whose base is in the poly-A region is reinterpreted linearly:
    /// the poly-A tail is the transcript's 3' terminus, not an intron, so the
    /// offset displaces `k` bases further along the contiguous downstream genome
    /// (this matches the normalizer's `c.*N+k → c.*(N+k)` folding and mutalyzer).
    /// Folding the offset into the effective tx coordinate is strand-correct
    /// because [`CdotTranscript::tx_to_genome_extending_polya`] already encodes
    /// the strand walk direction.
    fn try_extend_polya_to_genome(
        &self,
        cdot_tx: &crate::data::cdot::CdotTranscript,
        transcript_id: &str,
        p: &crate::hgvs::location::CdsPos,
    ) -> Option<u64> {
        if !p.utr3 || p.special.is_some() || p.base < 1 {
            return None;
        }
        // Only a non-negative downstream offset is meaningful in the poly-A
        // region (a negative offset would point back toward the genomic core,
        // which is mapped — leave it to the normal path / decline).
        let offset = match p.offset {
            None => 0,
            Some(k) if k >= 0 => k as u64,
            Some(_) => return None,
        };
        let cds_end = cdot_tx.cds_end?;
        // `c.*N` lives at tx position `cds_end + N - 1` (cds_end is 0-based
        // exclusive, so `c.*1` is at `cds_end`); mirrors `cds_to_tx_aware` in
        // `data::mapping`. Use checked arithmetic throughout.
        let base_tx_pos = cds_end.checked_add(p.base as u64)?.checked_sub(1)?;
        // The base position must itself be in the poly-A region (past the
        // 3'-terminal exon). If `c.*N` is still inside an exon, an offset would
        // be a genuine intron position — not ours to handle. (tx coords are
        // contiguous in real cdot, so "no exon" here means "past the terminal
        // exon", never a pre-terminal hole; `tx_to_genome_extending_polya`
        // independently rejects a pre-terminal hole as a second guard.)
        if cdot_tx.tx_to_genome(base_tx_pos).is_some() {
            return None;
        }
        // Fold the (downstream) offset into the effective tx coordinate.
        let tx_pos = base_tx_pos.checked_add(offset)?;
        // Gate to the transcript's TRUE length (FASTA) so the walk stays within
        // the real poly-A tail. Decline if the length is unavailable.
        let true_tx_length = self.provider.get_sequence_length(transcript_id).ok()?;
        if tx_pos >= true_tx_length {
            return None;
        }
        let g = cdot_tx.tx_to_genome_extending_polya(tx_pos)?;
        // Validate the walked coordinate against the contig bounds: the plus-strand
        // walk (`genome_end + delta`) could in principle run off the 3' end of the
        // contig, and a 1-based genome coordinate is never 0. The transcript-length
        // guard above bounds the walk to the real poly-A tail, but that is a
        // transcript-frame bound, not a contig-frame one — decline if the resulting
        // genome coordinate falls outside `1..=contig_len`.
        if g == 0 {
            return None;
        }
        if let Ok(contig_len) = self.provider.get_sequence_length(&cdot_tx.contig) {
            if g > contig_len {
                return None;
            }
        }
        Some(g)
    }

    fn project_to_genomic_nc(&self, variant: &HgvsVariant) -> Result<HgvsVariant, FerroError> {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::interval::GenomeInterval;
        use crate::hgvs::location::{CdsPos, GenomePos, TxPos};
        use crate::hgvs::variant::{GenomeVariant, HgvsVariant, LocEdit};
        use crate::reference::Strand as RefStrand;

        // The `UncertainBoundary<T>` → `T` resolver lives at module scope
        // (`resolve_uncertain_boundary`) so the direct bare-NM_ path shares the
        // exact `?`-vs-range classification used here.

        // 0. Reject variant kinds that have no genomic coordinate counterpart
        //    in this iteration (see #328 for allele support).
        match variant {
            HgvsVariant::Genome(_) => return Ok(variant.clone()),
            HgvsVariant::Protein(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support protein (p.) inputs".to_string(),
                })
            }
            HgvsVariant::Mt(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support mitochondrial (m.) inputs"
                        .to_string(),
                })
            }
            HgvsVariant::Circular(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support circular (o.) inputs".to_string(),
                })
            }
            HgvsVariant::RnaFusion(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support RNA fusion inputs".to_string(),
                })
            }
            HgvsVariant::GenomeRing(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support genome ring inputs".to_string(),
                })
            }
            HgvsVariant::Supernumerary(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support supernumerary (sup) inputs"
                        .to_string(),
                })
            }
            // Unreachable in practice: `project_to_genomic` handles `Allele`
            // before delegating here (#851). Kept for match exhaustiveness.
            HgvsVariant::Allele(_) => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic_nc does not handle allele inputs directly (#851)"
                        .to_string(),
                })
            }
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "project_to_genomic does not support null/unknown allele markers"
                        .to_string(),
                })
            }
            HgvsVariant::Cds(_) | HgvsVariant::Tx(_) | HgvsVariant::Rna(_) => {}
        }

        // 1. Pull the common pieces (accession, gene_symbol, edit, raw start/end
        //    as CdsPos) out of each transcript-coord variant kind. RNA and Tx
        //    positions have the same shape as CdsPos (base / offset / utr3-or-
        //    downstream), so we normalize to CdsPos here. Coding positions
        //    route through `data::mapping::CoordinateMapper::cds_to_genome`
        //    (which also handles intronic offsets and the `utr3` flag); the
        //    non-coding path uses `CdotTranscript::tx_to_genome` directly.
        // `input_is_transcript_coord` distinguishes `n.`/`r.` inputs (whose
        // `base` is a 1-based *transcript* position) from `c.` inputs (CDS
        // coordinates). For a coding transcript an `n.`/`r.` base must be
        // converted to a CDS coordinate before the CDS-aware genome mapping runs
        // (#693) — without it, `n.204` is mapped as if it were `c.204`.
        let (accession, gene_symbol, edit_mu, start_cds, end_cds, input_is_transcript_coord) =
            match variant {
                HgvsVariant::Cds(v) => {
                    let s = resolve_uncertain_boundary(&v.loc_edit.location.start, "c.", "start")?;
                    let e = resolve_uncertain_boundary(&v.loc_edit.location.end, "c.", "end")?;
                    (
                        v.accession.clone(),
                        v.gene_symbol.clone(),
                        v.loc_edit.edit.clone(),
                        s,
                        e,
                        false,
                    )
                }
                HgvsVariant::Tx(v) => {
                    let s = resolve_uncertain_boundary(&v.loc_edit.location.start, "n.", "start")?;
                    let e = resolve_uncertain_boundary(&v.loc_edit.location.end, "n.", "end")?;
                    let to_cds = |p: TxPos| CdsPos {
                        base: p.base,
                        offset: p.offset,
                        utr3: p.downstream,
                        special: None,
                    };
                    (
                        v.accession.clone(),
                        v.gene_symbol.clone(),
                        v.loc_edit.edit.clone(),
                        to_cds(s),
                        to_cds(e),
                        true,
                    )
                }
                HgvsVariant::Rna(v) => {
                    let s = resolve_uncertain_boundary(&v.loc_edit.location.start, "r.", "start")?;
                    let e = resolve_uncertain_boundary(&v.loc_edit.location.end, "r.", "end")?;
                    // RnaPos shape mirrors CdsPos exactly.
                    let to_cds = |p: crate::hgvs::location::RnaPos| CdsPos {
                        base: p.base,
                        offset: p.offset,
                        utr3: p.utr3,
                        special: None,
                    };
                    (
                        v.accession.clone(),
                        v.gene_symbol.clone(),
                        v.loc_edit.edit.clone(),
                        to_cds(s),
                        to_cds(e),
                        true,
                    )
                }
                _ => unreachable!("variant kind already filtered above"),
            };

        // 2. Reject `?` position sentinels explicitly — these have base == 0
        //    and no offset, which would otherwise propagate into the genomic
        //    mapper as a meaningless coordinate.
        if start_cds.is_unknown()
            || start_cds.is_special()
            || end_cds.is_unknown()
            || end_cds.is_special()
        {
            return Err(FerroError::UnsupportedProjection {
                reason: "project_to_genomic cannot resolve `?` or special (pter/qter/cen) \
                         position sentinels"
                    .to_string(),
            });
        }

        // 3. Resolve the genomic parent reference. Normally this is an explicit
        //    NG/NC parent on `genomic_context` — per #327 we do NOT synthesize a
        //    parent from cdot. The one exception is an LRG transcript/protein
        //    input (`LRG_<n>t<m>` / `LRG_<n>p<m>`): its genomic parent `LRG_<n>`
        //    is determined *structurally* by the accession itself, not inferred
        //    from cdot, so we derive it when no explicit parent is given (#480).
        let parent = match accession.genomic_context.as_deref().cloned() {
            Some(p) => p,
            None => Self::lrg_genomic_parent(&accession).ok_or_else(|| {
                FerroError::UnsupportedProjection {
                    reason: format!(
                        "input variant has no parent reference (genomic_context) on accession {}; \
                         project_to_genomic requires an explicit NG/NC parent (see #327)",
                        accession.full()
                    ),
                }
            })?,
        };

        // 4. Resolve the transcript via the cdot mapper (the same backing
        //    store used by the g. → c./n. path).  Working directly off cdot
        //    keeps us symmetric with `project_single_inner` and avoids a
        //    second `provider.get_transcript` load whose exon entries may
        //    lack genomic coordinates. When the NG/NC parent's version
        //    disambiguates a genome build (NC_*.10 → GRCh37, NC_*.11 →
        //    GRCh38) we ask cdot for that build's view so multi-build cdot
        //    loads project against the correct alignment (issue #389).
        let transcript_id = accession.transcript_accession();
        // Prefer the genome build encoded by the parent's chromosomal placement.
        // An `NG_`/`LRG_` parent accession carries no build, so without this the
        // cdot pivot below would fall back to cdot's primary build and could
        // compute coordinates on a different `NC_` build than the one the
        // placement — and the downstream re-anchor / de-anchor steps that stamp
        // `placement.nc` — assume (#646). Fall back to the parent's own build tag
        // (e.g. an explicit `NC_*` parent) when there is no placement.
        // Prefer the build the *input variant* carries (an explicit `NC_*.10`/`.11`
        // accession or `genomic_context`, or a `GRCh37(...)` assembly tag) — the
        // authoritative signal of intended build (#653/#713). Only when the input
        // is build-agnostic (e.g. a bare `NG_` parent) fall back to the build
        // implied by the GRCh38-preferred placement, then to the parent's own
        // version. Deriving the build from the variant first (rather than from the
        // chosen placement) lets the placement *selection* below honor it, instead
        // of being circularly fixed to GRCh38.
        let build_hint = self
            .build_hint_for_variant(variant)
            .or_else(|| {
                self.provider
                    .genomic_placement(&parent)
                    .and_then(|p| self.provider.infer_genome_build(&p.nc))
            })
            .or_else(|| self.provider.infer_genome_build(&parent));
        let cdot_mapper = self.projector.mapper();
        // Silent version-substitution gate (#785), c.→g. direction. Unlike
        // `project_variant`/`project_variant_all`, this path does not normalize
        // the original transcript-coordinate input first, so it would not inherit
        // the gate the normalizer applies in `normalize_core`. Apply the same
        // guarantee here: when the input names an EXPLICIT transcript version the
        // reference does not carry exactly, the lenient `get_transcript[_on_build]`
        // pivot below would silently fall back to a *sibling* version and project
        // onto the genome using that sibling's exon/CDS frame — a spec-invalid
        // result whose stated reference and coordinate frame disagree, with no
        // error. Decline with `TranscriptVersionNotExact` instead of fuzzing the
        // version. The predicate permits exact-version cdot-genome synthesis and
        // rejects only cross-version substitution; an LRG transcript carries no
        // version (`accession.version == None`), so its exact RefSeq alias path is
        // untouched. A transcript that is simply absent with no sibling to
        // substitute is left to the `ReferenceNotFound` miss below, unchanged.
        if accession.version.is_some() && is_transcript_accession(&transcript_id) {
            let exact = self.provider.has_transcript_version_exact(&transcript_id);
            let resolves = match build_hint {
                Some(b) => cdot_mapper
                    .cdot()
                    .get_transcript_on_build(&transcript_id, b)
                    .is_some(),
                None => cdot_mapper.cdot().get_transcript(&transcript_id).is_some(),
            };
            if !exact && resolves {
                return Err(FerroError::TranscriptVersionNotExact {
                    requested: transcript_id.to_string(),
                });
            }
        }
        let cdot_tx = match build_hint {
            Some(b) => cdot_mapper
                .cdot()
                .get_transcript_on_build(&transcript_id, b),
            None => cdot_mapper.cdot().get_transcript(&transcript_id),
        }
        .ok_or_else(|| FerroError::ReferenceNotFound {
            id: transcript_id.to_string(),
        })?;
        let strand = cdot_tx.strand;
        if strand == RefStrand::Unknown {
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "transcript {} has unknown strand; cannot project to g.",
                    transcript_id
                ),
            });
        }

        // 5. Pick the right primitive: for coding transcripts use the
        //    CDS-aware path from `data::mapping` (handles intronic offsets,
        //    UTR markers, and strand swap). For non-coding transcripts the
        //    input `base` is the 1-based tx_pos directly; convert it to the
        //    0-based form cdot exposes and use `CdotTranscript::tx_to_genome`.
        let is_coding = cdot_tx.cds_start.is_some() && cdot_tx.cds_end.is_some();
        // For an `n.`/`r.` input on a coding transcript, the carried `base` is a
        // 1-based *transcript* position, not a CDS coordinate; convert it to the
        // CDS frame before the CDS-aware genome mapping (#693). An exonic
        // position converts directly; an intronic-offset position converts its
        // exon-anchor `base` and keeps the offset/utr3 flags. A transcript
        // position that cannot be expressed in the CDS frame (`base < 1`, or a
        // `None` conversion for a position past the last exon) is rejected here:
        // leaving it unconverted would let the CDS-aware mapping interpret an
        // `n.`/`r.` coordinate as a `c.` coordinate and emit the wrong g.
        // position (#693).
        let to_cds_frame = |p: CdsPos| -> Result<CdsPos, FerroError> {
            // A position already flagged 3'UTR/downstream (`r.*N`, `n.*N`) carries
            // its `*N` distance in `base`, not an absolute 1-based transcript
            // position — feeding it to `cds_pos_from_tx_pos` would renumber it as a
            // transcript base. Leave such positions untouched so the downstream
            // CDS-aware mapping interprets the `*N` directly.
            if !(input_is_transcript_coord && is_coding) || p.utr3 {
                return Ok(p);
            }
            if p.base < 1 {
                return Err(FerroError::InvalidCoordinates {
                    msg: format!(
                        "transcript position {} is outside coding transcript {} and cannot be \
                         mapped as a CDS coordinate",
                        p.base, transcript_id
                    ),
                });
            }
            match cdot_tx.cds_pos_from_tx_pos((p.base - 1) as u64) {
                Some(cds) => Ok(CdsPos {
                    base: cds.base,
                    offset: p.offset,
                    // Keep the resolved 3'UTR flag even for intronic anchors
                    // (`c.*N±M`); dropping it when an offset is present would
                    // collapse a 3'UTR intron into a plain `c.N±M`.
                    utr3: cds.utr3,
                    special: p.special,
                }),
                None => Err(FerroError::InvalidCoordinates {
                    msg: format!(
                        "transcript position {} is outside coding transcript {} and cannot be \
                         mapped as a CDS coordinate",
                        p.base, transcript_id
                    ),
                }),
            }
        };
        // The poly-A 3'UTR fallback below is a `c.*`-only refinement (#797): it
        // reinterprets a `c.*` 3'UTR endpoint that lands in the transcript's
        // post-transcriptional poly-A tail. `n.`/`r.` inputs are reshaped into
        // `CdsPos` for this closure (with `utr3` carrying their own
        // downstream/3'UTR flag), so without this gate an off-exon `n.*`/`r.*`
        // endpoint would be reinterpreted as a `c.*` poly-A walk and emit a
        // genomic coordinate from a coordinate class that has no poly-A
        // semantics. Restrict the fallback to genuine `c.` (CDS) inputs.
        let is_cds_input = matches!(variant, HgvsVariant::Cds(_));
        let map_pos = |p: CdsPos| -> Result<u64, FerroError> {
            let p = to_cds_frame(p)?;
            if is_coding {
                let result =
                    match cdot_mapper.cds_to_genome_on_build(&transcript_id, &p, build_hint) {
                        Ok(r) => r,
                        Err(e) => {
                            // Poly-A-region 3'UTR fallback (#797). A `c.*` position
                            // can land in the transcript's post-transcriptional
                            // poly-A tail, which the (e.g. #790-derived) exon→genome
                            // map strips — so `cds_to_genome` declines. mutalyzer
                            // maps such positions by a contiguous coordinate walk into
                            // downstream genome from the transcript's own 3'-terminal
                            // exon. `try_extend_polya_to_genome` self-gates (a `c.*`
                            // position whose base lands in the poly-A region, an
                            // optional positive downstream offset, and an effective
                            // tx position within the transcript's TRUE length, i.e.
                            // inside the real poly-A tail) and returns None to
                            // propagate the original decline otherwise — so beyond
                            // the tail we never walk into non-transcript genome.
                            if is_cds_input {
                                if let Some(g) =
                                    self.try_extend_polya_to_genome(cdot_tx, &transcript_id, &p)
                                {
                                    return Ok(g);
                                }
                            }
                            return Err(e);
                        }
                    };
                // Sequence-aware correction (#644), mirroring the g.→c. path.
                // Only simple exonic positions are corrected; intronic / special
                // positions are derived by boundary arithmetic the realignment
                // doesn't model, so leave them to the naive result.
                if p.offset.is_some() || p.special.is_some() || p.utr3 {
                    return Ok(result.variant.base);
                }
                let tx_pos = match cdot_tx.cds_to_tx(p.base) {
                    Some(t) => t,
                    None => return Ok(result.variant.base),
                };
                self.correct_genome_for_exon_indel(
                    cdot_tx,
                    &transcript_id,
                    tx_pos,
                    result.variant.base,
                )
            } else {
                if p.offset.is_some() {
                    // Intronic n. offsets require coding-style exon-boundary
                    // arithmetic; we don't currently expose that primitive
                    // for non-coding transcripts. Surface a clear error.
                    return Err(FerroError::UnsupportedProjection {
                        reason: format!(
                            "project_to_genomic does not yet support intronic offsets on \
                             non-coding transcripts ({}); see #332",
                            transcript_id
                        ),
                    });
                }
                // `base` is 1-based on TxPos / RnaPos; cdot's tx coords are
                // 0-based, so subtract 1.
                if p.base <= 0 {
                    return Err(FerroError::InvalidCoordinates {
                        msg: format!(
                            "non-coding transcript {} position must be positive, got {}",
                            transcript_id, p.base
                        ),
                    });
                }
                let tx_pos_0based = (p.base - 1) as u64;
                let naive = cdot_tx.tx_to_genome(tx_pos_0based).ok_or_else(|| {
                    FerroError::InvalidCoordinates {
                        msg: format!(
                            "tx position {} not in any exon of {}",
                            p.base, transcript_id
                        ),
                    }
                })?;
                // Sequence-aware correction (#644): the inbound g.→n. path runs
                // `correct_cds_for_exon_indel` regardless of coding status, so a
                // non-coding `n.`/`r.` coordinate is corrected on the way in.
                // Apply the matching outbound correction here so it round-trips.
                self.correct_genome_for_exon_indel(cdot_tx, &transcript_id, tx_pos_0based, naive)
            }
        };

        let start_g = map_pos(start_cds)?;
        let end_g = map_pos(end_cds)?;

        // 6. Build the genomic interval. On minus strand the c./n./r. interval
        //    runs anti-parallel to the genome, so the smaller genomic coord
        //    becomes the start.
        let (g_start_pos, g_end_pos) = match strand {
            RefStrand::Plus => (GenomePos::new(start_g), GenomePos::new(end_g)),
            RefStrand::Minus => (GenomePos::new(end_g), GenomePos::new(start_g)),
            RefStrand::Unknown => {
                return Err(FerroError::UnsupportedProjection {
                    reason: format!(
                        "transcript {} has unknown strand; cannot project to g.",
                        transcript_id
                    ),
                });
            }
        };
        let g_interval = GenomeInterval::new(g_start_pos, g_end_pos);

        // 7. Transform the edit for the transcript strand. The c./n./r. edit
        //    reads on the transcript's sense strand; on minus strand we need
        //    to reverse-complement to put it back onto the genomic reference.
        //    We extract the concrete `NaEdit` (preserving `Mu` certainty) and
        //    apply `transform_edit_for_strand` from the project::edit module.
        let edit_inner =
            edit_mu
                .inner()
                .cloned()
                .ok_or_else(|| FerroError::UnsupportedProjection {
                    reason: "project_to_genomic requires a concrete edit (not `?`)".to_string(),
                })?;
        let g_edit_inner: NaEdit = transform_edit_for_strand(&edit_inner, strand);
        // Preserve the Mu certainty (Certain vs Uncertain).
        let g_edit_mu = edit_mu.map(|_| g_edit_inner);

        // This NC-frame result is the **pivot** used downstream (genome→CDS via
        // cdot needs chromosome coordinates). The public `project_to_genomic`
        // re-anchors it into an NG_/LRG_ parent's own frame for the genomic
        // *output*; the pivot must stay NC (#480).
        let g_variant = GenomeVariant {
            accession: parent,
            gene_symbol,
            loc_edit: LocEdit::with_uncertainty(g_interval, g_edit_mu),
        };
        Ok(HgvsVariant::Genome(g_variant))
    }

    /// Re-anchor an NC-frame genome variant into a genomic parent's own frame
    /// using its [`GenomicPlacement`] (#480): apply the affine NC→parent
    /// transform to the interval, and reverse-complement the edit when the
    /// parent runs antiparallel to the chromosome. If an endpoint falls outside
    /// the placed span the variant is returned unchanged (chromosome frame).
    fn reanchor_genome_to_parent(
        gv: crate::hgvs::variant::GenomeVariant,
        placement: &crate::reference::GenomicPlacement,
    ) -> Result<crate::hgvs::variant::GenomeVariant, FerroError> {
        use crate::hgvs::interval::GenomeInterval;
        use crate::hgvs::variant::GenomeVariant;
        let bounds = match (
            gv.loc_edit.location.start.inner(),
            gv.loc_edit.location.end.inner(),
        ) {
            (Some(s), Some(e)) => Some((s.base, e.base)),
            _ => None,
        };
        let Some((nc_lo, nc_hi)) = bounds else {
            // Uncertain/compound endpoint: no single position to re-anchor.
            // Decline rather than stamp the parent accession on an unresolved
            // chromosome coordinate, which would be invalid HGVS (#655).
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "cannot re-anchor {} into its NG_/LRG_ parent frame: an endpoint has no \
                     single resolved position (uncertain or compound boundary) (#480/#655)",
                    gv.accession.transcript_accession(),
                ),
            });
        };
        let (Some(a), Some(b)) = (placement.nc_to_parent(nc_lo), placement.nc_to_parent(nc_hi))
        else {
            // Endpoint falls outside the placed genomic span: the parent frame
            // cannot express this coordinate. Decline rather than emit a
            // chromosome coordinate under the parent accession (#655).
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "cannot re-anchor {} into its NG_/LRG_ parent frame: an endpoint falls \
                     outside the placed genomic span (#480/#655)",
                    gv.accession.transcript_accession(),
                ),
            });
        };
        // A minus-strand placement reverses coordinate order.
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        let interval = GenomeInterval::new(GenomePos::new(lo), GenomePos::new(hi));
        let edit = if placement.strand == crate::reference::Strand::Minus {
            gv.loc_edit
                .edit
                .map(|inner| transform_edit_for_strand(&inner, crate::reference::Strand::Minus))
        } else {
            gv.loc_edit.edit
        };
        Ok(GenomeVariant {
            accession: gv.accession,
            gene_symbol: gv.gene_symbol,
            loc_edit: LocEdit::with_uncertainty(interval, edit),
        })
    }

    /// If `variant` is a genome variant on an `NG_`/`LRG_` parent whose
    /// chromosomal [`GenomicPlacement`] is known, rewrite it into the chromosome
    /// (`NC_`) frame — the inverse of [`Self::reanchor_genome_to_parent`] (#480).
    ///
    /// cdot indexes and aligns transcripts on `NC_` contigs, so an `NG_`/`LRG_`
    /// genomic *input* (e.g. `NG_012337.1:g.4812_4813insTAC`) cannot be
    /// enumerated or projected in its own frame. Translating it into the `NC_`
    /// frame lets the normal fan-out path find and project onto the overlapping
    /// transcripts.
    ///
    /// Returns the rewritten variant and the original parent accession when a
    /// rewrite happened (so the per-transcript output can be re-framed under the
    /// parent); otherwise returns the input unchanged with `None`.
    fn deanchor_genomic_parent_input(
        &self,
        variant: &HgvsVariant,
    ) -> (HgvsVariant, Option<Accession>) {
        use crate::hgvs::interval::GenomeInterval;
        use crate::hgvs::variant::GenomeVariant;

        // The parent placement is selected by the build the input variant carries
        // (an explicit `NC_*.10`/`.11` or `GRCh37(...)` tag; `None` for a bare
        // `NG_`, which defaults to GRCh38). `None` keeps the prior GRCh38-preferred
        // behavior, so this is inert until an explicit build is present (#653/#713).
        let build = self.build_hint_for_variant(variant);

        // Transcript-coordinate inputs (c./n./r.) carrying an NG_/LRG_
        // genomic_context are taken into the chromosome (NC_) frame here too, so
        // the fan-out below enumerates the OTHER overlapping transcripts
        // (cross-isoform) and frames each result under the parent (#646). Without
        // this they project only onto the input transcript, unframed.
        let cnr_accession = match variant {
            HgvsVariant::Cds(c) => Some(&c.accession),
            HgvsVariant::Tx(t) => Some(&t.accession),
            HgvsVariant::Rna(r) => Some(&r.accession),
            _ => None,
        };
        if let Some(acc) = cnr_accession {
            if let Some(ctx) = acc.genomic_context.as_deref() {
                if let Some(placement) = (ctx.prefix.as_ref() == "NG" || ctx.is_lrg())
                    .then(|| self.provider.genomic_placement_on_build(ctx, build))
                    .flatten()
                {
                    // A poly-A `c.*` endpoint (#797) has no genomic alignment, so
                    // `project_to_genomic_nc` now derives it by a contiguous walk to
                    // an OFF-exon `NC_` coordinate. De-anchoring to that synthetic
                    // `Genome` would seed the fan-out stab off the transcript span
                    // (and then fail genome→CDS re-derivation). Keep the input's
                    // transcript-coordinate (`c.`) form instead and just signal the
                    // parent for re-framing: the `c.` fan-out seed path
                    // (`extract_contig_and_pos`) anchors poly-A endpoints to the
                    // 3'-terminal exon, so enumeration + the multi-axis short-circuit
                    // proceed correctly — mirroring how a bare-`NM_` poly-A `c.*`
                    // already fans out.
                    if self.cds_input_has_polya_endpoint(variant) {
                        return (variant.clone(), Some(ctx.clone()));
                    }
                    if let Ok(HgvsVariant::Genome(nc_gv)) = self.project_to_genomic_nc(variant) {
                        // `nc_gv` carries the parent (NG_/LRG_) accession but NC_
                        // coordinates. Stamp the placement's own `NC_` accession —
                        // the authoritative chromosome for this parent — so the
                        // projector enumerates transcripts at the locus. This
                        // mirrors the bare-genomic branch below and
                        // `reanchor_genome_to_parent`, which also key off
                        // `placement.nc`; deriving the contig from cdot instead can
                        // disagree when cdot stores a `chr`-style alias rather than
                        // the `NC_` accession.
                        let nc_variant = HgvsVariant::Genome(GenomeVariant {
                            accession: placement.nc.clone(),
                            gene_symbol: None,
                            loc_edit: nc_gv.loc_edit,
                        });
                        return (nc_variant, Some(ctx.clone()));
                    }
                }
            }
        }

        let HgvsVariant::Genome(gv) = variant else {
            return (variant.clone(), None);
        };
        let Some(placement) = self
            .provider
            .genomic_placement_on_build(&gv.accession, build)
        else {
            return (variant.clone(), None);
        };
        let bounds = match (
            gv.loc_edit.location.start.inner(),
            gv.loc_edit.location.end.inner(),
        ) {
            (Some(s), Some(e)) => Some((s.base, e.base)),
            _ => None,
        };
        let Some((p_lo, p_hi)) = bounds else {
            return (variant.clone(), None);
        };
        let (Some(a), Some(b)) = (placement.parent_to_nc(p_lo), placement.parent_to_nc(p_hi))
        else {
            return (variant.clone(), None);
        };
        // A minus-strand placement reverses coordinate order.
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        let interval = GenomeInterval::new(GenomePos::new(lo), GenomePos::new(hi));
        // Resolve any range-reference insertion (`ins<start>_<end>`, `delins…`)
        // against the *parent* (`NG_`/`LRG_`) sequence in the parent frame BEFORE
        // de-anchoring, so the de-anchored variant carries a frame-independent
        // literal. The range-ref endpoints are parent-frame coordinates; left
        // unresolved they would be fetched against the `NC_` accession stamped
        // below and come back as `insNNN…` (#652). Best-effort: if the parent
        // sequence can't be fetched, keep the original edit (no regression).
        let parent_resolved_edit = gv.loc_edit.edit.clone().map(|inner| {
            match crate::normalize::rules::canonicalize_insertion_expand(
                &inner,
                // The bare parent accession (strips any genomic_context), matching
                // the accessor `reanchor_genome_to_parent` uses — `full()` would
                // prepend a `CTX(…)` wrapper and fail the provider name lookup.
                &gv.accession.transcript_accession(),
                crate::normalize::rules::InsCoordKind::Direct,
                &self.provider,
            ) {
                // Range reference resolved to a literal — carry the frame-independent bases.
                Ok(Some(resolved)) => resolved,
                // Legitimate no-op: a plain literal / count / non-range payload returns
                // `Ok(None)` (src/normalize/rules.rs); keep the original edit.
                Ok(None) => inner,
                // A genuine resolution failure (out-of-bounds range, non-IUPAC base,
                // provider error). Do NOT silently re-emit the unresolved
                // `ins<start>_<end>` as `insNNN…` (#652) without a trace: surface it so
                // an unresolvable range reference is at least visible, even though the
                // `(HgvsVariant, Option<Accession>)` return type has no `Result` channel
                // to decline through.
                Err(e) => {
                    log::warn!(
                        "deanchor_genomic_parent_input: failed to resolve range-reference \
                         insertion against parent {}: {} — keeping the unresolved edit",
                        gv.accession.transcript_accession(),
                        e,
                    );
                    inner
                }
            }
        });
        let edit = if placement.strand == crate::reference::Strand::Minus {
            parent_resolved_edit
                .map(|inner| transform_edit_for_strand(&inner, crate::reference::Strand::Minus))
        } else {
            parent_resolved_edit
        };
        let parent = gv.accession.clone();
        let nc_variant = HgvsVariant::Genome(GenomeVariant {
            accession: placement.nc.clone(),
            gene_symbol: gv.gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(interval, edit),
        });
        (nc_variant, Some(parent))
    }

    /// Project an already-normalized g. variant onto the *curated* set of
    /// overlapping transcripts (#656), skipping re-normalization.
    ///
    /// The result is the curated enumerated set rather than every
    /// cdot-overlapping record: superseded transcript versions are collapsed
    /// (only the highest version per base accession is kept) and predicted
    /// `XM_`/`XR_` models are dropped when a curated `NM_`/`NR_` transcript
    /// covers the same locus; predicted models are kept only when they are the
    /// sole coverage. See [`select_enumerated_transcript_ids`] for the policy.
    ///
    /// Callers that pre-normalize once and then fan-out across transcripts
    /// should use this method.
    pub fn project_normalized_all(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Vec<VariantProjection>, FerroError> {
        // Surface a contradictory `--assembly` override like every other public
        // entry point (#715). The internal caller `project_variant_all` already
        // warns before delegating, so it routes through the non-warning
        // `project_normalized_all_inner` to avoid double-warning.
        self.warn_assembly_conflict(variant);
        self.project_normalized_all_inner(variant)
    }

    /// Fan-out projection without the `--assembly`-conflict warning, for callers
    /// (e.g. [`project_variant_all`]) that have already emitted it.
    fn project_normalized_all_inner(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Vec<VariantProjection>, FerroError> {
        // 2. Extract contig + first-base position from the normalized variant.
        //    For Allele variants we use the first inner variant's accession.
        //    c./n./r. inputs are resolved via cdot (contig from the
        //    transcript's exon table, position from the start-of-range
        //    transcript→genome projection) so PR #379's per-transcript
        //    widening reaches the fan-out path here (issue #389 item 3).
        let (contig, pos) = self.extract_contig_and_pos(variant)?;

        // 3. Find overlapping transcripts via the Projector (sorted by priority).
        let projection_result = self.projector.project(&contig, pos)?;

        // 3a. Apply the enumeration policy (#656): collapse superseded versions
        //     and prefer curated transcripts over predicted models, so the
        //     enumerated set matches mutalyzer's curated set rather than every
        //     cdot-overlapping record.
        let all_ids: Vec<&str> = projection_result
            .projections
            .iter()
            .map(|p| p.transcript_id.as_str())
            .collect();
        let keep: std::collections::HashSet<&str> = select_enumerated_transcript_ids(&all_ids)
            .into_iter()
            .collect();

        // 4. Project against each kept overlapping transcript (priority order preserved).
        let mut results = Vec::with_capacity(keep.len());
        for tx_proj in &projection_result.projections {
            if !keep.contains(tx_proj.transcript_id.as_str()) {
                continue;
            }
            match self.project_variant_inner(variant, &tx_proj.transcript_id) {
                Ok(vp) => results.push(vp),
                Err(e) => {
                    // Skip transcripts that fail — log at trace level only.
                    log::trace!(
                        "project_normalized_all: skipping {} for {}: {}",
                        tx_proj.transcript_id,
                        variant,
                        e
                    );
                }
            }
        }

        Ok(results)
    }

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    /// Core projection logic, operating on a variant that is assumed to be
    /// already normalized.  Does NOT re-normalize.
    fn project_variant_inner(
        &self,
        variant: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        // Dispatch on variant kind.
        match variant {
            HgvsVariant::Allele(allele) => {
                self.project_allele_inner(allele, variant, transcript_id)
            }
            _ => self.project_single_inner(variant, transcript_id),
        }
    }

    /// Project a compound [`AlleleVariant`] by recursively projecting each
    /// inner variant and combining the results.
    fn project_allele_inner(
        &self,
        allele: &AlleleVariant,
        original: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        // Empty allele: pass through with no coding/protein, but validate the
        // transcript ID — every other path through `project_variant_inner` errors
        // out for unknown transcripts (via `project_single_inner`), and this
        // branch must behave the same to avoid silently returning bogus
        // projections for typo'd or missing accessions.
        if allele.variants.is_empty() {
            // Honor the input's parent NG/NC (or assembly tag) when
            // resolving the gene_name from cdot, so an NG/NC-parented
            // allele on a multi-build cdot picks the correct build's
            // transcript record (issue #389).
            let build_hint = self.build_hint_for_variant(original);
            let gene_symbol = match build_hint {
                Some(b) => self
                    .projector
                    .mapper()
                    .cdot()
                    .get_transcript_on_build(transcript_id, b),
                None => self.projector.mapper().cdot().get_transcript(transcript_id),
            }
            .ok_or_else(|| FerroError::ReferenceNotFound {
                id: transcript_id.to_string(),
            })?
            .gene_name
            .clone();
            return Ok(VariantProjection {
                // No inner variant to derive a genomic form from; report `None`
                // rather than wrapping the (possibly non-genomic) input allele
                // into `.genomic`. Consistent with `coding`/`protein`, which are
                // already `None` for the empty allele (#508 review).
                genomic: None,
                coding: None,
                // Empty allele: no inner member to derive an n. form from.
                noncoding: None,
                protein: None,
                rna: None,
                transcript_id: transcript_id.to_string(),
                gene_symbol,
                is_frameshift: false,
                is_intronic: false,
                is_utr: false,
                // Empty allele: no inner edit, so no initiation-codon effect.
                affects_init: false,
            });
        }

        // Recursively project each inner variant.
        let mut inner_projections = Vec::with_capacity(allele.variants.len());
        for inner in &allele.variants {
            let proj = self.project_variant_inner(inner, transcript_id)?;
            inner_projections.push(proj);
        }

        // Aggregate flags.
        let is_frameshift = inner_projections.iter().any(|p| p.is_frameshift);
        let is_intronic = inner_projections.iter().any(|p| p.is_intronic);
        let is_utr = inner_projections.iter().any(|p| p.is_utr);

        // gene_symbol from any projection that has one.
        let gene_symbol = inner_projections.iter().find_map(|p| p.gene_symbol.clone());

        // Build the coding allele from inner c./n. variants.
        let coding_variants: Option<Vec<HgvsVariant>> =
            inner_projections.iter().map(|p| p.coding.clone()).collect();
        let coding = coding_variants
            .map(|variants| HgvsVariant::Allele(AlleleVariant::new(variants, allele.phase)));

        // Build the n. (transcript-relative) allele the same all-or-nothing way:
        // only when EVERY member carries an n. form (mirrors the coding rule).
        let noncoding_variants: Option<Vec<HgvsVariant>> = inner_projections
            .iter()
            .map(|p| p.noncoding.clone())
            .collect();
        let noncoding = noncoding_variants
            .map(|variants| HgvsVariant::Allele(AlleleVariant::new(variants, allele.phase)));

        // A cis allele whose member disrupts the translation initiation codon
        // has an unknown whole-protein consequence: once initiation is uncertain
        // the downstream members are moot (a protein product whose start codon
        // is unknown cannot be described). Collapse to that member's init-unknown
        // form — `p.(Met1?)` on a canonical `ATG` start, `p.?` on a
        // non-AUG-initiation transcript — mirroring the single-variant path
        // (#771) and the corpus row `NM_024426.4:c.1_4delinsATGA` →
        // `c.[1C>A;4C>A]` → `p.?`. This also rescues the common case where a
        // downstream non-initiation member declines under the #625 non-`ATG`
        // guard, which would otherwise drop the all-or-nothing allele protein to
        // `None`. Cis only: trans members are independent alleles, left to the
        // all-or-nothing rule below.
        let cis_init_unknown = if allele.phase == AllelePhase::Cis {
            inner_projections
                .iter()
                .find(|p| p.affects_init)
                .and_then(|p| p.protein.clone())
        } else {
            None
        };
        let protein = if let Some(init_unknown) = cis_init_unknown {
            Some(init_unknown)
        } else if allele.phase == AllelePhase::Cis
            && is_frameshift
            && inner_projections.iter().any(|p| p.protein.is_some())
        {
            // In-cis frameshift (#855): a downstream frameshift makes the combined
            // protein product uncertain — there is no spec form for "sense change +
            // frameshift" in one allele — so the whole-protein consequence collapses
            // to `p.?` ("an effect is expected but cannot be reliably predicted",
            // uncertain.md), reusing a member protein's accession/gene_symbol.
            inner_projections
                .iter()
                .find_map(|p| p.protein.as_ref())
                .map(build_cis_whole_protein_unknown)
        } else {
            // Build the protein allele only if ALL inner projections have a protein.
            let all_have_protein = inner_projections.iter().all(|p| p.protein.is_some());
            if all_have_protein {
                let protein_variants: Vec<HgvsVariant> = inner_projections
                    .iter()
                    .filter_map(|p| p.protein.clone())
                    .collect();
                // In-cis adjacent substitutions are described as a single combined
                // delins, not a bracketed list (#855; delins.md: `p.[Arg76Ser;
                // Cys77Trp]` is "not correct" for adjacent residues). Separated
                // substitutions and any non-substitution member keep the bracketed
                // allele (delins.md: "two variants separated by one or more amino
                // acids should be described individually"). Cis only — trans members
                // are independent alleles.
                if allele.phase == AllelePhase::Cis {
                    combine_cis_substitution_proteins(&protein_variants).or_else(|| {
                        Some(HgvsVariant::Allele(AlleleVariant::new(
                            protein_variants,
                            allele.phase,
                        )))
                    })
                } else {
                    Some(HgvsVariant::Allele(AlleleVariant::new(
                        protein_variants,
                        allele.phase,
                    )))
                }
            } else {
                None
            }
        };

        // Derive `.genomic` from the inner projections, not the raw `original`:
        // for a c./n./r. allele `original` is a non-genomic allele, and a
        // bare-NM_ inner projection has no genomic form at all. Mirror the
        // protein rule — build a g. allele only when EVERY member carries a
        // genomic form; otherwise report `None` so a genome-less projection
        // honestly says so rather than stuffing a transcript-coordinate allele
        // into `.genomic` (#508 review).
        let all_have_genomic = inner_projections.iter().all(|p| p.genomic.is_some());
        let genomic = if all_have_genomic {
            let genomic_variants: Vec<HgvsVariant> = inner_projections
                .iter()
                .filter_map(|p| p.genomic.clone())
                .collect();
            Some(HgvsVariant::Allele(AlleleVariant::new(
                genomic_variants,
                allele.phase,
            )))
        } else {
            None
        };

        // Predicted RNA allele (#693): predict from the assembled coding allele
        // — `predict_rna` handles the `Allele` shape directly, emitting one
        // outer predicted wrapper `r.([a;b])` (recommendations/RNA/alleles.md).
        // Re-frame under the input's parent so it matches the input framing.
        // Build-scope the transcript fetch by the input's NC/NG build (#843):
        // for a bare `g.` allele the member NC accession is build-bearing but
        // unparented, so a raw `cached_get_transcript_for_variant(original, …)`
        // would read the primary (GRCh38) transcript even for a GRCh37 input,
        // yielding a wrong-build `r.`. `predict_rna_transcript` stamps the
        // build-bearing accession as `genomic_context` first, mirroring the
        // single genome-pivot path's `cache_variant`.
        let rna = coding.as_ref().and_then(|coding_allele| {
            self.predict_rna_transcript(original, transcript_id)
                .ok()
                .and_then(|tx| crate::project::rna::predict_rna(coding_allele, &tx))
        });
        let rna = Self::reframe_rna_from_input(rna, original);

        Ok(VariantProjection {
            genomic,
            coding,
            noncoding,
            protein,
            rna,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift,
            is_intronic,
            is_utr,
            // Aggregate: true when any member's edit reaches the initiation codon.
            affects_init: inner_projections.iter().any(|p| p.affects_init),
        })
    }

    /// Project a single (non-allele) variant, assuming it has already been
    /// normalized.
    ///
    /// Accepts g./c./n./r. inputs. Non-genomic inputs are first
    /// projected back to their parent NG/NC reference via
    /// [`Self::project_to_genomic`], then re-enter the g.-anchored
    /// projection path. The parent reference must be present in the
    /// input's `Accession.genomic_context` (per #327); inputs without
    /// it surface a clear "no parent reference" diagnostic.
    ///
    /// **Transcript-id constraint.** For c./n./r. inputs the caller's
    /// `transcript_id` must equal the input's own
    /// `accession.transcript_accession()`. Mixing the two (e.g. passing
    /// `c.4C>A` on NM_FOO.1 but requesting projection against NM_BAR.1)
    /// is rejected with a clear error rather than silently projecting
    /// through one transcript's exons then re-projecting through
    /// another's. For g. inputs there is no such input-transcript, so
    /// the caller's `transcript_id` is the only target.
    ///
    /// **3'-shift asymmetry note (re #334).** This routine receives a
    /// pre-normalized variant. The c./n./r. normalizer respects the
    /// HGVS exon-junction exception (does not shift across exons),
    /// while the g. normalizer does not. Consequently, the same
    /// biological variant fed as c. vs. g. near an exon junction can
    /// project to different `(g., c., p.)` tuples. This is intentional:
    /// the input axis carries the spec-correct semantic for *that*
    /// axis. Callers who want g.-axis semantics should pass the g.
    /// variant. The c.-input projection preserves the c.-axis canonical
    /// form by construction.
    /// (#328)
    /// Predict the protein consequence of a coding variant from its 1-based
    /// CDS position(s) and edit. Shared by the genome-pivot path
    /// (`project_single_inner`) and the direct bare-NM_ path
    /// (`project_coding_direct`) so both compute protein identically.
    ///
    /// Returns `Ok(None)` when no protein consequence applies — most intronic
    /// edits, UTR, a non-coding transcript, or an edit shape with no
    /// in-frame/frameshift prediction. The intronic carve-out: a deletion whose
    /// endpoints are both intronic and bracket one or more complete coding exons
    /// (`whole_exon_deletion_span` returns `Some`) *is* predicted — it routes
    /// through the indel predictor on its clamped exonic CDS span (#498) — so
    /// only non-whole-exon intronic edits return `None`. `Err` is reserved for
    /// genuine failures (unknown reference, etc.), never for "no consequence"
    /// (mirrors the existing contract).
    ///
    /// `cds_start`/`cds_end` are 1-based CDS positions; `cache_variant` keys
    /// the per-(transcript, parent) sequence and ref-translation caches.
    #[allow(clippy::too_many_arguments)]
    fn predict_protein_consequence(
        &self,
        transcript_id: &str,
        cdot_protein: Option<String>,
        is_coding: bool,
        is_intronic: bool,
        is_utr: bool,
        c_edit: &NaEdit,
        cds_start: &CdsPos,
        cds_end: &CdsPos,
        cache_variant: &HgvsVariant,
    ) -> Result<Option<HgvsVariant>, FerroError> {
        // An edit whose span reaches the translation initiation codon (CDS
        // 1–3) has an unpredictable protein consequence, reported as
        // `p.(Met1?)` on a canonical-`ATG`-initiation transcript, or as `p.?`
        // on a non-AUG-initiation transcript where residue 1 is not `Met`
        // (#771; see the start-codon split at the `affects_init` branch below)
        // (HGVS recommendations/protein/{substitution.md:51,
        // deletion.md:62}; #512). This holds even when the edit's 5' end lies
        // in the 5'UTR (CDS base ≤ 0) — e.g. a boundary-spanning insertion or
        // duplication like `c.-1_2dup` — so the coarse `is_utr` gate must not
        // drop it to "no protein predicted" (#504). Intronic offsets never
        // reach the initiation codon, so an offset disqualifies the edit here.
        let affects_init = edit_reaches_initiation_codon(c_edit, cds_start, cds_end, is_intronic);
        // A deletion that removes one or more complete coding exons (both
        // endpoints intronic, bracketing the exon) has a predictable protein
        // consequence (HGVS protein/deletion.md "one or more exons"; #498).
        // Clamp it to the exonic CDS bases removed; if that yields a span, let
        // it through the intronic gate so the indel predictor runs on it below.
        let whole_exon = whole_exon_deletion_span(c_edit, cds_start, cds_end);
        // A deletion spanning the CDS→3'UTR boundary disrupts the termination
        // codon (stop-loss) and yields a C-terminal extension — the 3'-side
        // mirror of `affects_init`. Without this carve-out the coarse `is_utr`
        // gate (set because the 3' end is a `*N` position) would drop it to "no
        // protein predicted" before the extension predictor runs (#857).
        let affects_term = edit_spans_cds_into_3utr(cds_start, cds_end, is_intronic);
        if (is_intronic && whole_exon.is_none())
            || !is_coding
            || (is_utr && !affects_init && !affects_term)
        {
            return Ok(None);
        }
        // Resolve the protein accession: the authoritative cdot/reference value
        // if present, else the transcript id itself as the `p.` accession (the
        // HGVS p. grammar requires a sequence identifier; see #310).
        //
        // We deliberately do NOT infer `NP_*`/`XP_*` from `NM_*`/`XM_*` by
        // preserving the number: RefSeq does not guarantee the NM and NP
        // numbers match (e.g. `NM_000077.4` ↔ `NP_000068.1`), so that inference
        // is frequently wrong and would attach a real, translated prediction to
        // a non-existent protein accession (#808). Falling back to the
        // transcript id keeps the output honest — the identifier names the
        // transcript whose CDS we actually translated — and still satisfies the
        // grammar.
        // #860: when the input names an LRG transcript (LRG_<n>t<k>), echo its
        // LRG protein namespace (LRG_<n>p<k>) rather than the resolved NP_ — both
        // are spec-valid public references (general.md L16), and preserving the
        // input's namespace matches the reference the caller used. Non-LRG inputs
        // return None here and keep the cdot NP_ (or the honest transcript-id
        // fallback, #310/#808) — byte-identical to the prior behavior.
        let prot_acc = lrg_protein_accession(transcript_id)
            .or(cdot_protein)
            .unwrap_or_else(|| transcript_id.to_string());
        // Issue #505: protein prediction below reads this transcript's CDS
        // bases directly and translates them. If the reference only carries a
        // *different* version of the transcript (a silent version-strip
        // substitution), those bases may pair with a different CDS/exon
        // structure, so the prediction would be wrong yet attributed to the
        // requested version. Decline to predict rather than emit a wrong
        // protein — this gates *all* downstream prediction, including the
        // initiation-codon short-circuit below, whose CDS positions are
        // likewise mapped against the requested-version transcript. The
        // predicate still permits exact-version cdot-genome synthesis — it
        // rejects only cross-version substitution.
        if !self.provider.has_transcript_version_exact(transcript_id) {
            return Ok(None);
        }
        // Issue #625: CDS-start sanity. Every protein prediction below trusts
        // that `Transcript.cds_start` points at the first base of the
        // initiation codon and translates the CDS bases directly from there.
        // When the cdot CDS annotation is inconsistent with the transcript
        // FASTA — cds_start lands mid-frame, not on an `ATG` (e.g.
        // NM_000425.3 cds_start→"AGG", NM_152263.2 cds_start→"GAT") — that
        // assumption is false: translation reads the wrong frame and would
        // fabricate a garbage consequence (premature internal stop, Xaa,
        // wrong anchor). Decline to predict rather than emit a bogus protein.
        //
        // The discriminator is the START codon only. Selenoproteins (e.g.
        // SELENON / NM_020451.2) legitimately carry an in-frame internal `TGA`
        // recoded as selenocysteine yet begin with `ATG`; keying on internal
        // stops would wrongly decline them, so we never inspect them here.
        //
        // This is a deliberate over-approximation: a few RefSeq transcripts use
        // experimentally-supported non-AUG initiation (CUG/GUG with a
        // `transl_except`), whose CDS is consistent yet non-`ATG`; we decline
        // those too (for non-initiation edits) rather than risk wrong-frame
        // translation of the common inconsistent-annotation case. An
        // initiation-codon-affecting edit is the exception — it reports an
        // unknown form (`p.?`) without translating, so it is handled below
        // before this decline (#771).
        //
        // A failure to read the CDS (missing sequence, degenerate coords) is
        // left to the per-edit paths below, which already treat
        // `ProteinSequenceUnavailable` as a non-fatal "no protein" — so we
        // only act on a CDS we could read whose start is not `ATG`. Read just
        // the start codon (not the whole 1–10 kbp CDS) since that is all the
        // guard inspects.
        let tx_for_cds_check =
            self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
        // Read the start codon once; it drives two independent decisions below,
        // with deliberately *different* defaults when the CDS is unreadable
        // (missing sequence / degenerate coords): the strict ATG form-choice
        // defaults conservative (`false` → `p.?`), the downstream decline gate
        // defaults permissive (`true` → translate). Each default is documented at
        // its assignment.
        let start_codon = read_cds_start_codon(&tx_for_cds_check);
        // Strict `ATG`: drives the initiation-codon-unknown form choice (#771) —
        // an init-codon edit yields `p.(Met1?)` on an `ATG` start, `p.?` on a
        // non-AUG start (residue 1 is not `Met`). A CUG-start init-codon variant
        // must stay `p.?`, so the form choice keeps the strict `ATG` test.
        // Unreadable start ⇒ default `false`: no `ATG` was observed, so emitting
        // `p.(Met1?)` would assert a `Met1` start that was never seen. Fall to the
        // conservative whole-protein-unknown `p.?` instead.
        let start_is_atg = start_codon
            .as_deref()
            .map(cds_has_valid_start)
            .unwrap_or(false);
        // Recognized initiator (`ATG` or a near-cognate `CTG`/`GTG`/`TTG`): drives
        // the #625 decline gate for a *downstream* edit. A near-cognate start
        // translates in the correct frame with the initiator `Met` installed, so
        // a downstream variant on a legitimate non-AUG-initiation transcript (e.g.
        // WT1 `NM_024426.4`) is translated rather than declined (#780). Genuine
        // off-the-start drift (neither `ATG` nor near-cognate) still declines.
        let start_is_recognized_initiator = start_codon
            .as_deref()
            .map(cds_has_recognized_start)
            .unwrap_or(true);
        // An initiation-codon-affecting edit short-circuits to an unknown-protein
        // form here, before edit-type dispatch, so boundary-spanning edits whose
        // 5' end is in the 5'UTR (which the per-edit paths reject via their
        // `cds base > 0` guards) are still reported. The in-CDS start-codon cases
        // (e.g. `c.1A>G`) also funnel through here; the per-edit paths keep their
        // own `affects_initiation_codon` guard as defense in depth for direct
        // callers (#504, #512). The consequence is reported WITHOUT translating,
        // so the non-ATG #625 frame concern does not apply: a canonical ATG start
        // gives `p.(Met1?)`; a non-AUG-initiation transcript (#625) gives the
        // plain whole-protein unknown `p.?` (residue 1 is not `Met`, so the
        // `Met1?` form would be wrong) — #771.
        if affects_init {
            let tx_for_codon =
                self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
            let prot = if start_is_atg {
                build_initiator_unknown(&prot_acc, &tx_for_codon)
            } else {
                build_whole_protein_unknown(&prot_acc, &tx_for_codon)
            };
            return Ok(Some(prot));
        }
        // Issue #625 (non-initiation edits only — initiation-codon edits already
        // returned an unknown form above). Decline to translate when the CDS does
        // not begin with a recognized initiator (`ATG` or a near-cognate
        // `CTG`/`GTG`/`TTG`): a start codon outside that set means `cds_start` is
        // likely inconsistent with the FASTA (off-the-start drift), so translation
        // would read the wrong frame and fabricate a bogus consequence. A
        // legitimate non-AUG-initiation transcript (CUG/GUG) keeps a correct frame
        // with the initiator `Met` installed, so its downstream variants now
        // translate (#780) — only genuine drift is declined.
        if !start_is_recognized_initiator {
            log::trace!(
                "declining protein prediction for {transcript_id}: reference CDS does not begin \
                 with a recognized initiator (ATG/CTG/GTG/TTG); cds_start is inconsistent with \
                 the transcript FASTA",
            );
            return Ok(None);
        }
        // Issue #332: route per-codon lookups through the variant-aware path
        // so an NG/NC-parented `cache_variant` picks the build-correct
        // chromosome; falls back to the bare provider lookup on
        // ReferenceNotFound. For a bare-NM_ `cache_variant` (no parent) this
        // is just the plain provider lookup.
        let mut protein = None;
        match c_edit {
            NaEdit::Substitution { .. } => {
                let tx_for_codon =
                    self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
                protein = Some(predict_substitution_protein(
                    &tx_for_codon,
                    cds_start.base,
                    c_edit,
                    &prot_acc,
                )?);
            }
            // Whole-exon deletion (both endpoints intronic): predict from the
            // clamped exonic CDS span, reusing the indel predictor on a clean
            // exonic deletion of that span (#498).
            NaEdit::Deletion { .. } if whole_exon.is_some() => {
                let (lo, hi) = whole_exon.expect("checked is_some");
                let tx_for_codon =
                    self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
                let ref_bundle =
                    self.cached_ref_translation(cache_variant, transcript_id, &tx_for_codon)?;
                let exonic_del = NaEdit::Deletion { sequence: None, length: None };
                match predict_indel_protein(
                    &tx_for_codon,
                    &ref_bundle,
                    lo,
                    hi,
                    &exonic_del,
                    &prot_acc,
                ) {
                    Ok(pv) => protein = Some(pv),
                    Err(FerroError::UnsupportedProjection { .. })
                    | Err(FerroError::ProteinSequenceUnavailable { .. }) => {}
                    Err(other) => return Err(other),
                }
            }
            // Boundary-spanning deletion (CDS 5' end, 3'UTR `*N` 3' end):
            // stop-loss → C-terminal extension. Splices against the extended
            // (CDS + 3'UTR) reference; `predict_stop_region_extension` returns
            // `None` if it is not a clean stop-loss. Must precede the combined
            // arm below (whose `cds_end.base > 0` would otherwise match a `*1`
            // end and panic in the CDS-only splice) (#857).
            NaEdit::Deletion { .. }
                if cds_end.utr3
                    && !cds_start.utr3
                    && cds_start.base > 0
                    && cds_start.offset.is_none()
                    && cds_end.offset.is_none() =>
            {
                let tx_for_codon =
                    self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
                let ref_bundle =
                    self.cached_ref_translation(cache_variant, transcript_id, &tx_for_codon)?;
                match predict_stop_region_extension(
                    &tx_for_codon,
                    &ref_bundle,
                    cds_start.base,
                    cds_end.base,
                    c_edit,
                    &prot_acc,
                ) {
                    Ok(pv) => protein = pv,
                    Err(FerroError::UnsupportedProjection { .. })
                    | Err(FerroError::ProteinSequenceUnavailable { .. }) => {}
                    Err(other) => return Err(other),
                }
            }
            NaEdit::Deletion { .. }
            | NaEdit::Insertion { .. }
            | NaEdit::Duplication { .. }
            | NaEdit::Delins { .. }
            | NaEdit::Inversion { .. }
                // Only predict for concrete exonic CDS positions. This arm's
                // own `offset.is_none()` guards exclude intronic offsets — the
                // loosened intronic gate above no longer does, since whole-exon
                // deletions now pass it (they are handled by the arm above). The
                // `!utr3` guards keep a `*N` end (which has `base > 0`) out of
                // the CDS-only splice, where it would panic; such a 3'UTR-
                // spanning edit is handled by the boundary arm above or falls to
                // `_ => {}` (#857).
                if cds_start.offset.is_none()
                    && cds_end.offset.is_none()
                    && cds_start.base > 0
                    && cds_end.base > 0
                    && !cds_start.utr3
                    && !cds_end.utr3 =>
            {
                let tx_for_codon =
                    self.cached_get_transcript_for_variant(cache_variant, transcript_id)?;
                let ref_bundle =
                    self.cached_ref_translation(cache_variant, transcript_id, &tx_for_codon)?;
                match predict_indel_protein(
                    &tx_for_codon,
                    &ref_bundle,
                    cds_start.base,
                    cds_end.base,
                    c_edit,
                    &prot_acc,
                ) {
                    Ok(pv) => protein = Some(pv),
                    // Non-fatal: unsupported edits or missing sequence → None.
                    Err(FerroError::UnsupportedProjection { .. })
                    | Err(FerroError::ProteinSequenceUnavailable { .. }) => {}
                    Err(other) => return Err(other),
                }
            }
            _ => {}
        }
        Ok(protein)
    }

    /// Direct c.→p. projection for a bare transcript coding input (no
    /// `genomic_context` parent). The c.→g.→CDS roundtrip cannot run without a
    /// genome alignment, but protein prediction only needs the transcript's
    /// CDS sequence and the 1-based CDS position — which an exonic c. variant
    /// already provides. The resulting projection has `genomic = None` (no
    /// genomic representation is available for a bare-NM_ input) (#498).
    fn project_coding_direct(
        &self,
        cds: &CdsVariant,
        normalized: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        // Mirror the genome path's transcript_id-mismatch guard.
        let input_tx = cds.accession.transcript_accession();
        if input_tx != transcript_id {
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "transcript_id mismatch: input is on {} but projection requested \
                     against {}; transcript-coordinate inputs must be projected against \
                     their own transcript",
                    input_tx, transcript_id,
                ),
            });
        }

        // Coding/protein/gene metadata: prefer cdot, fall back to the
        // sequence provider's transcript record.
        let (is_coding, cdot_protein, gene_symbol) =
            match self.projector.mapper().cdot().get_transcript(transcript_id) {
                Some(t) => (
                    t.cds_start.is_some(),
                    t.protein.clone(),
                    t.gene_name.clone(),
                ),
                None => {
                    let tx = self.cached_get_transcript_for_variant(normalized, transcript_id)?;
                    (
                        tx.cds_start.is_some(),
                        tx.protein_id.clone(),
                        tx.gene_symbol.clone(),
                    )
                }
            };

        let edit = cds.loc_edit.edit.inner().cloned().ok_or_else(|| {
            FerroError::UnsupportedProjection {
                reason: "coding variant has no concrete edit".to_string(),
            }
        })?;
        // Resolve the c. boundaries through the shared resolver so a bare `c.?`
        // (`Single(Mu::Unknown)`) surfaces as `UnsupportedProjection` — matching
        // the genome-pivot `project_to_genomic` contract — while a genuine range
        // boundary still surfaces as `InvalidCoordinates` (#508 review).
        let cds_start = resolve_uncertain_boundary(&cds.loc_edit.location.start, "c.", "start")?;
        let cds_end = resolve_uncertain_boundary(&cds.loc_edit.location.end, "c.", "end")?;

        // `resolve_uncertain_boundary` only rejects the parsed `?` sentinel
        // (`Single(Mu::Unknown)`). A programmatically-built `CdsPos::unknown`
        // carries its `base == 0` sentinel inside a `Mu::Certain`, so it slips
        // through above; without this guard `base <= 0` would misclassify it as
        // 5'UTR and return `Ok`. Reject it here to match the
        // `project_to_genomic` contract (#508 review).
        if cds_start.is_unknown()
            || cds_start.is_special()
            || cds_end.is_unknown()
            || cds_end.is_special()
        {
            return Err(FerroError::UnsupportedProjection {
                reason: "cannot resolve `?` or special (pter/qter/cen) position sentinels on \
                         direct c.→p. projection"
                    .to_string(),
            });
        }

        // Flags derived straight from the input c. position (no genome):
        // intronic = any offset; UTR = 5'UTR (base ≤ 0) or 3'UTR (`*`).
        let is_intronic = cds_start.offset.is_some() || cds_end.offset.is_some();
        let is_utr = !is_intronic
            && (cds_start.base <= 0 || cds_end.base <= 0 || cds_start.utr3 || cds_end.utr3);

        let protein = self.predict_protein_consequence(
            transcript_id,
            cdot_protein,
            is_coding,
            is_intronic,
            is_utr,
            &edit,
            &cds_start,
            &cds_end,
            normalized,
        )?;

        let frameshift = is_frameshift(normalized);

        // c↔n axis: derive the n. (transcript-relative) form from the resolved
        // CDS positions via the exon/CIGAR-aware mapper. Carry the *input's*
        // gene symbol (`cds.gene_symbol`) so the n. form matches the `coding`
        // form (which is the raw input here) rather than disagreeing on the
        // gene tag. Best-effort — a derivation edge (e.g. legacy c.0) is logged
        // and leaves the axis None rather than failing the whole projection.
        let noncoding = match self.cached_get_transcript_for_variant(normalized, transcript_id) {
            Ok(tx) => match crate::project::transcript_axis::noncoding_from_coding(
                &cds_start,
                &cds_end,
                &tx,
                &edit,
                transcript_id,
                cds.gene_symbol.clone(),
            ) {
                Ok(v) => Some(v),
                Err(e) => {
                    log::trace!("c→n derivation failed for {transcript_id}: {e}");
                    None
                }
            },
            Err(_) => None,
        };

        // Build-scope the `predict_rna` transcript fetch via the shared helper
        // for consistency with the allele site and the genome-pivot path (#843).
        // NOTE: this site is reached ONLY for a bare-`NM_` c. input
        // (`genomic_context.is_none()` gates the dispatch above), which carries no
        // genome build at all — so `rna_build_context` returns `None` and the
        // fetch is byte-identical to the prior build-agnostic lookup. The change
        // is therefore a deliberate latent no-op here: it makes every
        // `predict_rna` fetch flow through one build-aware helper, so a future
        // build-bearing input reaching this path would be scoped automatically
        // rather than silently reading the primary build. (The reachable
        // build-scoping fix is at the allele site.)
        let rna = self
            .predict_rna_transcript(normalized, transcript_id)
            .ok()
            .and_then(|tx| crate::project::rna::predict_rna(normalized, &tx));
        let rna = Self::reframe_rna_from_input(rna, normalized);

        Ok(VariantProjection {
            genomic: None,
            coding: Some(normalized.clone()),
            noncoding,
            protein,
            rna,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift: frameshift,
            is_intronic,
            is_utr,
            affects_init: edit_reaches_initiation_codon(&edit, &cds_start, &cds_end, is_intronic),
        })
    }

    /// Direct n./r.→CDS→p. projection for a bare transcript non-coding (`n.`)
    /// or RNA (`r.`) input (no `genomic_context` parent). The c.→g.→CDS
    /// roundtrip cannot run without a genome alignment, but an `n.` base is a
    /// 1-based transcript position: convert it to a CDS position via the
    /// transcript's `cds_start` (the exon/CIGAR-aware
    /// [`CoordinateMapper::tx_to_cds`]), then run the same protein-consequence
    /// prediction the coding path uses. The resulting projection has
    /// `genomic = None` (no genomic representation is available for a bare-NM_
    /// input) (#506).
    ///
    /// `normalized` must be an [`HgvsVariant::Tx`] or [`HgvsVariant::Rna`].
    fn project_noncoding_direct(
        &self,
        normalized: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        use crate::reference::Strand as RefStrand;

        // Pull the common pieces (accession, gene_symbol, edit, raw start/end
        // as a transcript TxPos) out of the n. / r. input. RNA positions share
        // the same `base`/`offset`/`utr3` shape as transcript positions, so we
        // normalize both to `TxPos`. RNA edits carry `Base::U`; the protein
        // machinery reads the transcript's DNA CDS, so we translate U→T below.
        let (accession, gene_symbol, edit_mu, axis_label, start_tx, end_tx) = match normalized {
            HgvsVariant::Tx(v) => {
                let s = resolve_uncertain_boundary(&v.loc_edit.location.start, "n.", "start")?;
                let e = resolve_uncertain_boundary(&v.loc_edit.location.end, "n.", "end")?;
                (
                    v.accession.clone(),
                    v.gene_symbol.clone(),
                    v.loc_edit.edit.clone(),
                    "n.",
                    s,
                    e,
                )
            }
            HgvsVariant::Rna(v) => {
                let s = resolve_uncertain_boundary(&v.loc_edit.location.start, "r.", "start")?;
                let e = resolve_uncertain_boundary(&v.loc_edit.location.end, "r.", "end")?;
                // RnaPos shape mirrors TxPos: base / offset / utr3↔downstream.
                let to_tx = |p: crate::hgvs::location::RnaPos| TxPos {
                    base: p.base,
                    offset: p.offset,
                    downstream: p.utr3,
                };
                (
                    v.accession.clone(),
                    v.gene_symbol.clone(),
                    v.loc_edit.edit.clone(),
                    "r.",
                    to_tx(s),
                    to_tx(e),
                )
            }
            _ => {
                return Err(FerroError::UnsupportedProjection {
                    reason: format!(
                        "project_noncoding_direct only accepts n./r. inputs, got {}",
                        normalized.variant_type()
                    ),
                })
            }
        };

        // Mirror the genome path's transcript_id-mismatch guard: an n./r. input
        // must be projected against its own transcript.
        let input_tx = accession.transcript_accession();
        if input_tx != transcript_id {
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "transcript_id mismatch: input is on {} but projection requested \
                     against {}; transcript-coordinate inputs must be projected against \
                     their own transcript",
                    input_tx, transcript_id,
                ),
            });
        }

        // `?` sentinel (base == 0) has no concrete transcript position; reject
        // it up front to match the coding-path contract. (Transcript positions
        // cannot carry chromosome-arm markers — only `CdsPos` has a `special`
        // field — so there is no pter/qter/cen case to handle here.)
        if start_tx.base == 0 || end_tx.base == 0 {
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "cannot resolve `?` position sentinel on direct {axis_label}→p. projection"
                ),
            });
        }

        // 3'-downstream (`*N`) positions carry the `*N` *relative* 3'UTR offset
        // in `base` with `downstream = true`. `CoordinateMapper::tx_to_cds`
        // (used below) classifies purely by comparing `base` against the CDS
        // bounds and never reads this flag, so it would misread `n.*5` / `r.*5`
        // as in-transcript position 5 — yielding a bogus `c.` form and a
        // spurious protein prediction. The genome-pivot path reshapes
        // `downstream` into `CdsPos.utr3` (which `cds_to_genome` honors), but the
        // direct path has no exon-aware 3'UTR translation, so reject these up
        // front rather than mis-project them.
        if start_tx.downstream || end_tx.downstream {
            return Err(FerroError::UnsupportedProjection {
                reason: format!(
                    "cannot project 3'-downstream (`*N`) positions on direct {axis_label}→p. \
                     projection without a genome alignment"
                ),
            });
        }

        // The transcript record is needed both for its CDS structure (to
        // convert tx→CDS) and for protein metadata. Prefer the cdot view for
        // gene/protein metadata, but the full `Transcript` is what carries the
        // exon table and `cds_start`/`cds_end` used by `tx_to_cds`.
        let tx = self.cached_get_transcript_for_variant(normalized, transcript_id)?;
        // Coding/protein/gene metadata: prefer cdot, fall back to the sequence
        // provider's transcript record. Mirror the coding-path sibling
        // (`project_coding_direct`) by resolving `is_coding` cdot-first too — the
        // two transcript stores can disagree (cdot lists a CDS, the provider
        // record doesn't), and reading `is_coding` from a different source than
        // the sibling would silently drop protein where the sibling predicts it.
        let (is_coding, cdot_protein, gene_symbol_meta) =
            match self.projector.mapper().cdot().get_transcript(transcript_id) {
                Some(t) => (
                    t.cds_start.is_some(),
                    t.protein.clone(),
                    t.gene_name.clone(),
                ),
                None => (
                    tx.cds_start.is_some() && tx.cds_end.is_some(),
                    tx.protein_id.clone(),
                    tx.gene_symbol.clone(),
                ),
            };
        // Carry the input's gene symbol through to the derived axes when it has
        // one, else the metadata gene symbol; mirrors the coding path.
        let gene_symbol = gene_symbol.or(gene_symbol_meta);

        let edit = edit_mu
            .inner()
            .cloned()
            .ok_or_else(|| FerroError::UnsupportedProjection {
                reason: format!("{axis_label} variant has no concrete edit"),
            })?;
        // Translate any RNA `U` bases to DNA `T` so the edit reads against the
        // transcript's DNA CDS. On the transcript sense strand no reverse-
        // complement is needed (the c./n. form already reads in transcript
        // orientation), so we pass `Strand::Plus`, whose only effect is the
        // U→T mapping. This is a no-op for n. (DNA) inputs.
        let c_edit = transform_edit_for_strand(&edit, RefStrand::Plus);

        // Non-coding transcript: there is no CDS to convert into, so there is
        // no protein consequence and no c. form. The n. form is the input
        // itself; report it on the non-coding axis and stop.
        if !is_coding {
            return Ok(VariantProjection {
                genomic: None,
                coding: None,
                noncoding: Some(normalized.clone()),
                protein: None,
                // No coding form on a non-coding transcript; the RNA prediction
                // surface (#485) predicts from the c. form, so it's unavailable here.
                rna: None,
                transcript_id: transcript_id.to_string(),
                gene_symbol,
                is_frameshift: false,
                is_intronic: start_tx.offset.is_some() || end_tx.offset.is_some(),
                is_utr: false,
                // Non-coding transcript: no CDS, so no initiation codon to reach.
                affects_init: false,
            });
        }

        // Convert the 1-based transcript positions to CDS positions via the
        // exon/CIGAR-aware mapper (the inverse of `noncoding_from_coding`'s
        // `cds_to_tx`). The intronic offset is carried through.
        let mapper = crate::convert::mapper::CoordinateMapper::new(&tx);
        let cds_start = mapper.tx_to_cds(&start_tx)?;
        let cds_end = mapper.tx_to_cds(&end_tx)?;

        // Flags derived from the resolved CDS positions (no genome): intronic =
        // any offset; UTR = 5'UTR (base ≤ 0) or 3'UTR (`*`). Mirrors the coding
        // path so the protein gate behaves identically.
        let is_intronic = cds_start.offset.is_some() || cds_end.offset.is_some();
        let is_utr = !is_intronic
            && (cds_start.base <= 0 || cds_end.base <= 0 || cds_start.utr3 || cds_end.utr3);

        let protein = self.predict_protein_consequence(
            transcript_id,
            cdot_protein,
            is_coding,
            is_intronic,
            is_utr,
            &c_edit,
            &cds_start,
            &cds_end,
            normalized,
        )?;

        // Build the c. (coding) form from the resolved CDS positions and the
        // (U→T-translated) edit. Best-effort: keep the bare transcript
        // accession so it renders without an NC_* wrapper, matching the
        // genome-pivot `coding` field.
        let coding = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession(transcript_id),
            gene_symbol: gene_symbol.clone(),
            loc_edit: LocEdit::new(CdsInterval::new(cds_start, cds_end), c_edit.clone()),
        });

        let frameshift = is_frameshift(normalized);

        let rna = crate::project::rna::predict_rna(&coding, &tx);
        // The synthesized `coding` form carries a cdot-looked-up gene symbol the
        // input did not have; re-frame the predicted `r.` from the input so it
        // renders bare (or under the input's parent for an `NG_(NM_)` input).
        let rna = Self::reframe_rna_from_input(rna, normalized);
        Ok(VariantProjection {
            genomic: None,
            coding: Some(coding),
            // The non-coding axis is the input itself.
            noncoding: Some(normalized.clone()),
            protein,
            rna,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift: frameshift,
            is_intronic,
            is_utr,
            affects_init: edit_reaches_initiation_codon(&c_edit, &cds_start, &cds_end, is_intronic),
        })
    }

    /// Fetch the exon's genome and transcript sequence slices in transcript
    /// 5'→3' orientation, ready to feed [`crate::data::exon_realign`].
    ///
    /// Returns `None` (decline) on any condition where a sequence-aware
    /// correction should not be attempted: the exon already carries a cdot
    /// CIGAR (its alignment is modelled — leave the existing E3007 / arithmetic
    /// behavior untouched), genomic data is unavailable, or a slice cannot be
    /// fetched. The two slices are returned `(genome, tx)` with the genome slice
    /// reverse-complemented on a minus-strand transcript so both read in
    /// transcript orientation.
    fn exon_slices_for_realign(
        &self,
        cdot_tx: &CdotTranscript,
        transcript_id: &str,
        exon_idx: usize,
    ) -> Option<(Vec<u8>, Vec<u8>)> {
        // Only the ungapped-but-divergent class is in scope (#644). Exons cdot
        // models with a CIGAR keep their existing semantics.
        match cdot_tx.exon_cigars.get(exon_idx) {
            Some(None) | None => {}
            Some(Some(_)) => return None,
        }
        if !self.provider.has_genomic_data() {
            return None;
        }
        let exon = cdot_tx.exons.get(exon_idx)?;
        let (genome_start, genome_end, tx_start, tx_end) = (exon[0], exon[1], exon[2], exon[3]);
        if genome_end <= genome_start || tx_end <= tx_start {
            return None;
        }

        // Transcript slice: cdot tx coords are 0-based half-open; the provider's
        // `get_sequence` is also 0-based half-open.
        let tx_seq = self
            .provider
            .get_sequence(transcript_id, tx_start, tx_end)
            .ok()?;
        // Genome slice: cdot genome coords are HGVS-value-based, so HGVS g.X is
        // 0-based X-1. The exon spans HGVS [genome_start, genome_end), i.e.
        // 0-based [genome_start - 1, genome_end - 1).
        let genome_seq = self
            .provider
            .get_genomic_sequence(&cdot_tx.contig, genome_start - 1, genome_end - 1)
            .ok()?;

        let genome_bytes = if cdot_tx.strand == crate::reference::Strand::Minus {
            crate::sequence::reverse_complement(&genome_seq).into_bytes()
        } else {
            genome_seq.into_bytes()
        };
        Some((genome_bytes, tx_seq.into_bytes()))
    }

    /// Apply the sequence-aware exon-indel correction to a genome→CDS mapping
    /// (issue #644).
    ///
    /// `naive` is the `CdsPos` the cdot arithmetic produced for `genome_pos`.
    /// When the containing exon is reported ungapped by cdot but its transcript
    /// and genome sequences differ by an indel, the naive position is off by the
    /// indel size; this re-derives the position from a local realignment. Only
    /// simple exonic positions are corrected — intronic / special / UTR-range
    /// positions (which carry an `offset` or a `special` marker) are returned
    /// unchanged, as is any position the realignment cannot confidently place.
    ///
    /// Returns `AlignmentGap` when the corrected position would fall strictly
    /// inside a genome-only gap discovered by the realignment (a genome base
    /// with no transcript counterpart) — consistent with the #603 E3007
    /// semantics for modelled CIGAR deletions.
    fn correct_cds_for_exon_indel(
        &self,
        cdot_tx: &CdotTranscript,
        transcript_id: &str,
        info: &MappingInfo,
        genome_pos: &GenomePos,
        naive: &CdsPos,
    ) -> Result<CdsPos, FerroError> {
        use crate::data::exon_realign::{correct_genome_offset, ExonOffsetCorrection};

        // Only correct simple exonic positions. Intronic offsets and special
        // sentinels are derived by boundary arithmetic the realignment doesn't
        // model; leave them alone.
        if info.is_intronic || naive.offset.is_some() || naive.special.is_some() {
            return Ok(*naive);
        }
        // The mapping records the 1-based exon number; the parallel exon /
        // exon_cigars index is `number - 1`.
        let exon_number = match info.exon_numbers.first() {
            Some(n) => *n,
            None => return Ok(*naive),
        };
        let exon_idx = (exon_number as usize).saturating_sub(1);
        let exon = match cdot_tx.exons.get(exon_idx) {
            Some(e) => e,
            None => return Ok(*naive),
        };
        let (genome_start, genome_end) = (exon[0], exon[1]);

        let (genome_seq, tx_seq) =
            match self.exon_slices_for_realign(cdot_tx, transcript_id, exon_idx) {
                Some(slices) => slices,
                None => return Ok(*naive),
            };

        // Offset of the queried genome position from the exon's transcript-5'
        // end, measured on the genome axis — the same orientation the slices
        // are in (see `cigar_deletion_gap_at_genome_pos`).
        let genome_offset = match cdot_tx.strand {
            crate::reference::Strand::Plus => {
                if genome_pos.base < genome_start || genome_pos.base >= genome_end {
                    return Ok(*naive);
                }
                (genome_pos.base - genome_start) as usize
            }
            crate::reference::Strand::Minus => {
                if genome_pos.base < genome_start || genome_pos.base >= genome_end {
                    return Ok(*naive);
                }
                (genome_end - 1 - genome_pos.base) as usize
            }
            crate::reference::Strand::Unknown => return Ok(*naive),
        };

        match correct_genome_offset(&genome_seq, &tx_seq, genome_offset) {
            None => Ok(*naive),
            Some(ExonOffsetCorrection::InsideGenomeOnlyGap) => Err(FerroError::AlignmentGap {
                msg: format!(
                    "genomic position {} falls in a transcript-genome alignment gap \
                     (sequence-detected genome-only indel) in exon {}",
                    genome_pos.base, exon_number
                ),
            }),
            Some(ExonOffsetCorrection::Delta(0)) => Ok(*naive),
            Some(ExonOffsetCorrection::Delta(delta)) => {
                // Apply the correction in transcript space and re-derive the
                // CDS position so a correction crossing a CDS/UTR boundary is
                // resolved correctly.
                let naive_tx = (exon[2] as i64) + (genome_offset as i64);
                let corrected_tx = naive_tx + delta;
                if corrected_tx < 0 {
                    return Ok(*naive);
                }
                cdot_tx
                    .cds_pos_from_tx_pos(corrected_tx as u64)
                    .map(Ok)
                    .unwrap_or(Ok(*naive))
            }
        }
    }

    /// Apply the sequence-aware exon-indel correction to a transcript→genome
    /// mapping — the c.→g. mirror of [`Self::correct_cds_for_exon_indel`]
    /// (issue #644).
    ///
    /// `tx_pos` is the 0-based transcript position being mapped and `naive` is
    /// the genome coordinate the cdot arithmetic produced for it. When the
    /// containing exon is reported ungapped by cdot but its sequences differ by
    /// an indel, the naive genome coordinate is off by the indel size; this
    /// re-derives it from the local realignment. Returns `naive` unchanged when
    /// no confident correction applies, and `AlignmentGap` when the transcript
    /// position has no genome counterpart (a transcript-only indel base).
    fn correct_genome_for_exon_indel(
        &self,
        cdot_tx: &CdotTranscript,
        transcript_id: &str,
        tx_pos: u64,
        naive: u64,
    ) -> Result<u64, FerroError> {
        use crate::data::exon_realign::{correct_tx_offset, TxOffsetCorrection};

        let exon = match cdot_tx.exon_for_tx_pos(tx_pos) {
            Some(e) => e,
            None => return Ok(naive),
        };
        let exon_idx = (exon.number as usize).saturating_sub(1);
        let (genome_seq, tx_seq) =
            match self.exon_slices_for_realign(cdot_tx, transcript_id, exon_idx) {
                Some(slices) => slices,
                None => return Ok(naive),
            };

        let tx_offset = (tx_pos - exon.tx_start) as usize;
        match correct_tx_offset(&genome_seq, &tx_seq, tx_offset) {
            None | Some(TxOffsetCorrection::Delta(0)) => Ok(naive),
            Some(TxOffsetCorrection::InsideTxOnlyGap) => Err(FerroError::AlignmentGap {
                msg: format!(
                    "transcript position {} falls in a transcript-genome alignment gap \
                     (sequence-detected transcript-only indel) in exon {}",
                    tx_pos, exon.number
                ),
            }),
            Some(TxOffsetCorrection::Delta(delta)) => {
                // The delta is in the genome axis, measured from the exon's
                // transcript-5' end. Apply it to the naive genome coordinate in
                // the strand's genomic direction (the naive coord already came
                // from `genome_start + tx_offset` on plus, `genome_end - 1 -
                // tx_offset` on minus, so on minus the axis runs the other way).
                let corrected = match cdot_tx.strand {
                    crate::reference::Strand::Plus => naive as i64 + delta,
                    crate::reference::Strand::Minus => naive as i64 - delta,
                    crate::reference::Strand::Unknown => return Ok(naive),
                };
                if corrected < 0 {
                    return Ok(naive);
                }
                Ok(corrected as u64)
            }
        }
    }

    fn project_single_inner(
        &self,
        normalized: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        // #887: resolve a whole-arm terminus (pter/qter) to the parent genome
        // frame via the same handler the raw project_to_genomic pivot uses
        // (1310), BEFORE the direct-path gate and the nc pivot — both reject the
        // sentinel. Fires only for a CDS special-terminus with a derivable
        // parent; otherwise falls through to the existing (correct) rejection.
        // Enforce the transcript_id match first (mirrors the pivot arm at
        // 3572-3584) so a mismatched target still errors.
        if let HgvsVariant::Cds(c) = normalized {
            // `.is_special()` lives on `CdsPos`, not on `UncertainBoundary<CdsPos>`,
            // so reach through `.inner()` — the exact idiom normalize uses
            // (normalize/mod.rs:1993).
            let start_special = c
                .loc_edit
                .location
                .start
                .inner()
                .is_some_and(|p| p.is_special());
            let end_special = c
                .loc_edit
                .location
                .end
                .inner()
                .is_some_and(|p| p.is_special());
            if start_special || end_special {
                if let Some(acc) = normalized.accession() {
                    let input_tx = acc.transcript_accession();
                    if input_tx != transcript_id {
                        return Err(FerroError::UnsupportedProjection {
                            reason: format!(
                                "transcript_id mismatch: input is on {} but projection \
                                 requested against {}; transcript-coordinate inputs must \
                                 be projected against their own transcript",
                                input_tx, transcript_id,
                            ),
                        });
                    }
                }
                if let Some(result) = self.project_cds_terminus_to_parent(normalized) {
                    // `result?` propagates a genuine terminus-resolution error
                    // (parent known but its length lookup failed). Normalize the
                    // resolved terminus to match the raw pivot + harness normalize;
                    // on a normalize failure fall back to the un-normalized parent
                    // terminus (still a valid coordinate) rather than aborting the
                    // whole projection — the genomic axis degrades gracefully,
                    // mirroring reanchor_and_normalize_genomic.
                    let raw = result?;
                    let genomic = self.normalizer.normalize(&raw).unwrap_or(raw);
                    return Ok(self.terminus_multiaxis_projection(
                        normalized,
                        transcript_id,
                        genomic,
                    ));
                }
                // No derivable parent (bare NM_/NR_) or cen/?/reversed -> fall
                // through; the direct-path gate / nc pivot rejects it as today.
            }
        }
        // Direct c.→p. path: a bare coding input (no genomic_context parent)
        // has no genome alignment to roundtrip through, but protein can be
        // predicted straight from the CDS (#498).
        if let HgvsVariant::Cds(c) = normalized {
            if c.accession.genomic_context.is_none() {
                return self.project_coding_direct(c, normalized, transcript_id);
            }
        }
        // Direct n./r.→CDS→p. path: a bare non-coding (`n.`) or RNA (`r.`)
        // input likewise has no genome alignment. An `n.` base is a 1-based
        // transcript position; convert it to a CDS position via the
        // transcript's `cds_start` and run the same protein prediction the
        // coding path uses (#506).
        match normalized {
            HgvsVariant::Tx(t) if t.accession.genomic_context.is_none() => {
                return self.project_noncoding_direct(normalized, transcript_id);
            }
            HgvsVariant::Rna(r) if r.accession.genomic_context.is_none() => {
                return self.project_noncoding_direct(normalized, transcript_id);
            }
            _ => {}
        }
        // Track which g. variant feeds the downstream pipeline AND will
        // be reported in `VariantProjection.genomic`. For g. input the
        // two are the same; for c./n./r. input we project once and use
        // the projected form for both — the `.genomic` field is the
        // canonical g. representation of the variant, not the input
        // axis.
        let projected_genome: HgvsVariant = match normalized {
            HgvsVariant::Genome(_) => normalized.clone(),
            HgvsVariant::Cds(_) | HgvsVariant::Tx(_) | HgvsVariant::Rna(_) => {
                // Reject transcript_id mismatch up front: projecting a
                // c.-on-NM_FOO through NM_BAR's exons and back through
                // NM_FOO's would silently produce nonsensical results.
                if let Some(acc) = normalized.accession() {
                    let input_tx = acc.transcript_accession();
                    if input_tx != transcript_id {
                        return Err(FerroError::UnsupportedProjection {
                            reason: format!(
                                "transcript_id mismatch: input is on {} but projection \
                                 requested against {}; transcript-coordinate inputs must \
                                 be projected against their own transcript",
                                input_tx, transcript_id,
                            ),
                        });
                    }
                }
                // Pivot in the chromosome (NC_) frame: genome→CDS/protein
                // derivation below runs through cdot, which is NC-based. The
                // re-anchored (NG_/LRG_) form is the genomic *output* only and
                // would mis-map back through cdot (#480).
                let g = self.project_to_genomic_nc(normalized)?;
                match &g {
                    HgvsVariant::Genome(_) => g,
                    _ => {
                        // Defensive: project_to_genomic is documented to
                        // return Genome for these axes. If a future
                        // refactor changes that contract, surface a
                        // clear error rather than silently misprojecting.
                        return Err(FerroError::UnsupportedProjection {
                            reason: format!(
                                "project_to_genomic returned a non-Genome variant for {} input",
                                normalized.variant_type()
                            ),
                        });
                    }
                }
            }
            _ => {
                return Err(FerroError::UnsupportedProjection {
                    reason: format!(
                        "VariantProjector does not accept {} inputs",
                        normalized.variant_type()
                    ),
                });
            }
        };
        let genome_variant = match &projected_genome {
            HgvsVariant::Genome(g) => g.clone(),
            _ => unreachable!("projected_genome is always Genome by construction above"),
        };

        let edit = genome_variant
            .loc_edit
            .edit
            .inner()
            .cloned()
            .ok_or_else(|| FerroError::UnsupportedProjection {
                reason: "g. variant has no concrete edit".to_string(),
            })?;

        // Look up the transcript in the cdot mapper. Prefer the build the
        // *input* carries — an explicit `NC_*.10`/`.11`, a `genomic_context`
        // parent, or the `--assembly` override that `project_to_genomic_nc`
        // already used to pick the placement (#716). Re-deriving the build
        // solely from the projected NC_ accession loses that signal when the
        // pivot anchored onto a contig whose version the heuristic does not
        // recognize (it would fall back to cdot's primary view and reject or
        // misproject the NC-frame pivot). Fall back to inferring from the
        // projected accession only when the input is build-agnostic (issue
        // #389: NC_*.10 → GRCh37, NC_*.11 → GRCh38).
        let build_hint = self
            .build_hint_for_variant(normalized)
            .or_else(|| self.provider.infer_genome_build(&genome_variant.accession));
        let cdot_tx = match build_hint {
            Some(b) => self
                .projector
                .mapper()
                .cdot()
                .get_transcript_on_build(transcript_id, b),
            None => self.projector.mapper().cdot().get_transcript(transcript_id),
        }
        .ok_or_else(|| FerroError::ReferenceNotFound {
            id: transcript_id.to_string(),
        })?;
        let strand = cdot_tx.strand;
        let gene_symbol = cdot_tx.gene_name.clone();
        let cdot_protein = cdot_tx.protein.clone();
        let is_coding = cdot_tx.cds_start.is_some();

        // Extract start and end genomic positions from the variant interval.
        let g_start = genome_variant
            .loc_edit
            .location
            .start
            .inner()
            .cloned()
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: "genomic interval start is unknown".to_string(),
            })?;
        let g_end = genome_variant
            .loc_edit
            .location
            .end
            .inner()
            .cloned()
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: "genomic interval end is unknown".to_string(),
            })?;

        let mapper = self.projector.mapper();

        // Genomic extent of the transcript — used to short-circuit the
        // genome-to-CDS mapping when the variant is outside the transcript's
        // span. The cached primary-build side-table is consulted via
        // `transcript_genome_span_on_build`; for alt builds the span is
        // computed on demand against the alt-build exon table (#389).
        let (tx_genome_start, tx_genome_end) = match build_hint {
            Some(b) => mapper
                .cdot()
                .transcript_genome_span_on_build(transcript_id, b),
            None => mapper.cdot().transcript_genome_span(transcript_id),
        }
        .ok_or_else(|| FerroError::ReferenceNotFound {
            id: transcript_id.to_string(),
        })?;

        // #797 (multi-axis): a `c.*` poly-A input maps to a genome coordinate that
        // is outside the cdot exon span by construction (the post-transcriptional
        // poly-A tail has no genomic alignment; `project_to_genomic_nc` derives it
        // by a contiguous downstream walk). Re-deriving the coding axis from that
        // synthetic coordinate via genome→CDS would fail the span check below
        // (`TranscriptNotOverlapping`), so short-circuit: keep the input's own
        // coding form and report the walked genome coordinate as the genomic axis.
        // Only a `c.*` (CDS) input whose pivot lands past the span qualifies: the
        // poly-A tail fallback is keyed on a `c.*` 3'UTR endpoint
        // (`try_extend_polya_to_genome` self-gates on `p.utr3`), so it has no
        // meaning for `n.`/`r.` inputs. A genuine off-transcript `g.` input — or a
        // non-overlapping `n.`/`r.` input — must still raise
        // `TranscriptNotOverlapping` via the normal path below.
        //
        // Gate on ANY outside endpoint, not both: a boundary-spanning interval
        // like `c.*3_*4del` has one endpoint inside the stripped terminal exon
        // (the 3'UTR genomic core) and the other in the poly-A tail past the span.
        // Requiring BOTH outside would miss that case — `map_position` would then
        // reject the single outside (poly-A) endpoint as `TranscriptNotOverlapping`
        // even though `project_to_genomic_nc` already walked it. This is safe for
        // exactly this `c.` branch: `projected_genome` only ever holds an off-span
        // coordinate that `project_to_genomic_nc` produced via the `c.`-gated
        // poly-A walk (`cds_to_genome` otherwise yields an in-exon coordinate), so
        // an outside endpoint here is necessarily a validated poly-A endpoint —
        // `polya_multiaxis_projection` re-anchors that already-walked genome
        // coordinate rather than re-deriving anything from it.
        if matches!(normalized, HgvsVariant::Cds(_)) {
            let outside = |p: &GenomePos| p.base < tx_genome_start || p.base >= tx_genome_end;
            if outside(&g_start) || outside(&g_end) {
                if let Some(projection) = self.polya_multiaxis_projection(
                    normalized,
                    transcript_id,
                    &projected_genome,
                    gene_symbol.clone(),
                ) {
                    return Ok(projection);
                }
            }
        }

        // Helper: map one GenomePos → CdsPos, converting out-of-range errors.
        // `normalized.to_string()` is only consumed by the error message and
        // the happy path doesn't touch it; building it eagerly was ~2% of
        // SNP CPU per `fmt::write` self-time, so defer it to the moment we
        // actually construct a `TranscriptNotOverlapping` error.
        let map_position = |gp: &GenomePos| -> Result<(CdsPos, MappingInfo), FerroError> {
            if gp.base < tx_genome_start || gp.base >= tx_genome_end {
                return Err(FerroError::TranscriptNotOverlapping {
                    variant: normalized.to_string(),
                    transcript_id: transcript_id.to_string(),
                });
            }
            match mapper.genome_to_cds_on_build(transcript_id, gp, build_hint) {
                Ok(res) => {
                    // Sequence-aware correction (#644): when cdot reports the
                    // containing exon ungapped but the transcript and genome
                    // FASTAs genuinely differ by an indel, the naive arithmetic
                    // `res.variant` is off by the indel's size. Re-derive the
                    // CDS position from a local realignment of the exon's two
                    // sequences. A no-op (and unchanged behavior) when the
                    // sequences agree, the exon already carries a CIGAR, or the
                    // alignment is not confident enough to correct.
                    let corrected = self.correct_cds_for_exon_indel(
                        cdot_tx,
                        transcript_id,
                        &res.info,
                        gp,
                        &res.variant,
                    )?;
                    Ok((corrected, res.info))
                }
                Err(FerroError::InvalidCoordinates { .. })
                | Err(FerroError::ConversionError { .. }) => {
                    Err(FerroError::TranscriptNotOverlapping {
                        variant: normalized.to_string(),
                        transcript_id: transcript_id.to_string(),
                    })
                }
                Err(other) => Err(other),
            }
        };

        let (cds_start_raw, info_start) = map_position(&g_start)?;
        let (cds_end_raw, info_end) = map_position(&g_end)?;

        // On minus strand the start and end of the c. interval are swapped.
        let (cds_start, cds_end) = match strand {
            crate::reference::Strand::Plus => (cds_start_raw, cds_end_raw),
            crate::reference::Strand::Minus => (cds_end_raw, cds_start_raw),
            crate::reference::Strand::Unknown => {
                return Err(FerroError::UnsupportedProjection {
                    reason: format!(
                        "transcript {} has unknown strand; cannot project g. → c./n.",
                        transcript_id
                    ),
                });
            }
        };

        // Transform the edit for the transcript strand.
        let c_edit = transform_edit_for_strand(&edit, strand);

        // Build the c./n. HGVS variant returned in `VariantProjection.coding`.
        // Kept bare-accession (no NC_* wrapper) for caller-visible rendering;
        // build-awareness for the parent-aware caches goes through
        // `cache_variant` below.
        let coding = if is_coding {
            let interval = CdsInterval::new(cds_start, cds_end);
            HgvsVariant::Cds(CdsVariant {
                accession: parse_accession(transcript_id),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, c_edit.clone()),
            })
        } else {
            let tx_start = TxPos::new(cds_start.base);
            let tx_end = TxPos::new(cds_end.base);
            HgvsVariant::Tx(TxVariant {
                accession: parse_accession(transcript_id),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(TxInterval::new(tx_start, tx_end), c_edit.clone()),
            })
        };

        // Parallel "cache key" variant that stamps the projected g.
        // accession (NC_*.10 / NC_*.11) as `genomic_context`. Used only
        // for `cached_get_transcript_for_variant` and
        // `cached_ref_translation` so that:
        //   1. `parent_key_for` reads the build-bearing NC accession,
        //      partitioning the (transcript_id, parent) caches by
        //      build — NC_*.10 vs NC_*.11 inputs no longer collide
        //      under the same key.
        //   2. `MultiFastaProvider::get_transcript_for_variant` sees the
        //      NC parent and probes the matching cdot build first via
        //      `get_transcript_on_build`, so the cached `Transcript` and
        //      `RefProteinBundle` are built against the correct
        //      alignment.
        // The returned `coding` field is unaffected (#389 follow-up).
        let cache_variant = if is_coding {
            HgvsVariant::Cds(CdsVariant {
                accession: parse_accession(transcript_id)
                    .with_genomic_context(genome_variant.accession.clone()),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(CdsInterval::new(cds_start, cds_end), c_edit.clone()),
            })
        } else {
            HgvsVariant::Tx(TxVariant {
                accession: parse_accession(transcript_id)
                    .with_genomic_context(genome_variant.accession.clone()),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(
                    TxInterval::new(TxPos::new(cds_start.base), TxPos::new(cds_end.base)),
                    c_edit.clone(),
                ),
            })
        };

        // Derive per-position flags from MappingInfo.
        let is_intronic = info_start.is_intronic || info_end.is_intronic;
        let is_utr = !is_intronic
            && (info_start.in_5utr || info_start.in_3utr || info_end.in_5utr || info_end.in_3utr);

        // Predict protein consequence for CDS variants. Shared with the
        // direct bare-NM_ c.→p. path (`project_coding_direct`) so both entry
        // points use identical CDS→protein logic (#498).
        let protein = self.predict_protein_consequence(
            transcript_id,
            cdot_protein,
            is_coding,
            is_intronic,
            is_utr,
            &c_edit,
            &cds_start,
            &cds_end,
            &cache_variant,
        )?;

        let frameshift = is_frameshift(&coding);

        // c↔n axis. For a coding transcript derive the n. form genome-free from
        // the resolved CDS positions (same edit, reframed coordinates);
        // best-effort — never fail the whole projection on an n.-derivation
        // edge. For a non-coding transcript `coding` already IS the n. (Tx)
        // form, so reuse it (forward-compatible: the non-coding full-projection
        // path is currently pinned unsupported — see
        // `tests/issue_328_projector_accepts_cnr.rs` — so this arm is not yet
        // reached, but carries the correct semantics for when it is).
        let noncoding = if is_coding {
            match self.cached_get_transcript_for_variant(&cache_variant, transcript_id) {
                Ok(tx) => match crate::project::transcript_axis::noncoding_from_coding(
                    &cds_start,
                    &cds_end,
                    &tx,
                    &c_edit,
                    transcript_id,
                    gene_symbol.clone(),
                ) {
                    Ok(v) => Some(v),
                    Err(e) => {
                        log::trace!("c→n derivation failed for {transcript_id}: {e}");
                        None
                    }
                },
                Err(_) => None,
            }
        } else {
            Some(coding.clone())
        };

        // `.genomic` is the canonical g. *output*. For a g. input it is the
        // (normalized) input itself, emitted as-is. For a c./n./r. input it is
        // the chromosome-frame pivot re-anchored into an `NG_`/`LRG_` parent's
        // own frame (#480) — the derivation above already ran on the `NC_` pivot
        // (`genome_variant`), which must stay in the chromosome frame for cdot.
        // When the parent frame cannot be resolved (no placement, an endpoint
        // outside the placed span, or an uncertain endpoint) the genomic axis is
        // dropped to `None` rather than stamping the parent accession on a
        // chromosome coordinate (invalid HGVS): unlike the g.-only
        // `project_to_genomic`, this path returns a multi-axis projection, so
        // dropping the unframable genomic axis preserves the still-valid
        // coding/noncoding/protein axes (#655/#702, mirroring `frame_projection_owned`).
        let genomic = match normalized {
            HgvsVariant::Genome(_) => Some(projected_genome),
            _ => match projected_genome {
                // Re-anchor + genomic-frame renormalize (#737); decline degrades the
                // genomic axis to `None` so the other axes survive (#655/#702).
                HgvsVariant::Genome(gv) => self
                    .reanchor_and_normalize_genomic(gv, self.build_hint_for_variant(normalized))
                    .ok(),
                other => Some(other),
            },
        };

        // Predicted RNA consequence (#485): derived best-effort from the c. form
        // and the transcript sequence; `None` when not representable.
        let rna = self
            .cached_get_transcript_for_variant(&cache_variant, transcript_id)
            .ok()
            .and_then(|tx| crate::project::rna::predict_rna(&coding, &tx));
        // `coding` carries the bare transcript accession plus a synthesized gene
        // symbol; re-frame the predicted `r.` under the input's parent (the
        // `NG_(NM_)` genomic context) so the RNA axis matches the input framing
        // rather than dropping the parent and substituting a gene symbol (#693).
        let rna = Self::reframe_rna_from_input(rna, normalized);

        Ok(VariantProjection {
            genomic,
            coding: Some(coding),
            noncoding,
            protein,
            rna,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift: frameshift,
            is_intronic,
            is_utr,
            affects_init: edit_reaches_initiation_codon(&c_edit, &cds_start, &cds_end, is_intronic),
        })
    }

    /// Build the multi-axis projection for a `c.*` poly-A input whose genome
    /// pivot lands outside the cdot exon span (#797).
    ///
    /// The post-transcriptional poly-A tail has no genomic alignment, so the
    /// synthetic genome coordinate produced by `project_to_genomic_nc`'s
    /// contiguous downstream walk cannot be re-derived back to a CDS position via
    /// the normal genome→CDS path (it would be rejected as
    /// `TranscriptNotOverlapping`). The coding axis is therefore the input's own
    /// (bare-accession) form, and the genomic axis is the re-anchored walked
    /// coordinate. A `c.*` 3'UTR/poly-A position has no protein consequence, so
    /// `protein`/`rna`/`noncoding` are `None` (best-effort: nothing to derive
    /// from a position with no genomic exon context).
    ///
    /// Returns `None` (so the caller falls through to the normal path, which
    /// raises the canonical decline) when the input is not a `c.*` (CDS)
    /// variant or when the genomic axis cannot be re-anchored into the
    /// parent frame — never stamping a chromosome coordinate under a parent
    /// accession (invalid HGVS).
    fn polya_multiaxis_projection(
        &self,
        normalized: &HgvsVariant,
        transcript_id: &str,
        projected_genome: &HgvsVariant,
        gene_symbol: Option<String>,
    ) -> Option<VariantProjection> {
        // The coding axis is the input's own form, rendered bare (no NC_ wrapper),
        // matching how `coding` is reported on the normal path. Only `c.*` (CDS)
        // inputs qualify for the poly-A fallback (the caller gates on
        // `HgvsVariant::Cds`); `n.`/`r.` inputs have no poly-A 3'UTR endpoint and
        // must fall through to the canonical `TranscriptNotOverlapping` decline.
        let coding = match normalized {
            HgvsVariant::Cds(c) => {
                let mut c = c.clone();
                c.accession.genomic_context = None;
                HgvsVariant::Cds(c)
            }
            _ => return None,
        };

        // Re-anchor the walked genome coordinate into the parent's own frame, or
        // drop the genomic axis if it cannot be framed (rather than emit a
        // chromosome coordinate under the parent accession; mirrors the normal
        // path's #655/#702 handling).
        let genomic = match projected_genome {
            HgvsVariant::Genome(gv) => self
                .reanchor_and_normalize_genomic(gv.clone(), self.build_hint_for_variant(normalized))
                .ok(),
            _ => return None,
        };

        Some(VariantProjection {
            genomic,
            coding: Some(coding),
            noncoding: None,
            protein: None,
            rna: None,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift: false,
            is_intronic: false,
            // A `c.*` poly-A position is in the 3'UTR (its post-transcriptional
            // extension), so flag the UTR axis.
            is_utr: true,
            affects_init: false,
        })
    }

    /// Genomic-only projection for a whole-arm terminus (pter/qter) `c.` input
    /// whose parent genome frame IS resolvable (#887). Mirrors
    /// `polya_multiaxis_projection`'s shape: echo the input's bare coding form,
    /// no protein/noncoding/rna axis (a whole-arm marker does not round-trip
    /// through cdot). `genomic` is the already-normalized parent-frame
    /// (NG_/LRG_) coordinate.
    fn terminus_multiaxis_projection(
        &self,
        normalized: &HgvsVariant,
        transcript_id: &str,
        genomic: HgvsVariant,
    ) -> VariantProjection {
        let coding = match normalized {
            HgvsVariant::Cds(c) => {
                let mut c = c.clone();
                c.accession.genomic_context = None;
                (Some(HgvsVariant::Cds(c.clone())), c.gene_symbol.clone())
            }
            _ => (None, None),
        };
        VariantProjection {
            genomic: Some(genomic),
            coding: coding.0,
            noncoding: None,
            protein: None,
            rna: None,
            transcript_id: transcript_id.to_string(),
            gene_symbol: coding.1,
            is_frameshift: false,
            is_intronic: false,
            // A whole-arm terminus is not a canonical UTR position.
            is_utr: false,
            affects_init: false,
        }
    }
}

/// Compute a representative 1-based genomic position for a non-coding
/// (n.) or RNA (r.) transcript coordinate, suitable for seeding
/// `Projector::project`'s stab query.
///
/// - For "simple" exonic positions (`offset.is_none() && !utr3 &&
///   base >= 1`), the exact `tx_to_genome` mapping is returned.
/// - For intronic (`offset.is_some()`), downstream/3'UTR (`utr3`), or
///   non-positive `base`, the conversion is best-effort: the
///   transcript's first-exon genome_start is returned. The stab query
///   at that position still lands inside the transcript, and the
///   per-transcript projection downstream re-derives the precise
///   coordinate. Falling back here keeps the fan-out path lenient for
///   inputs the precise tx-to-genome path can't resolve (which is also
///   what `project_to_genomic` itself rejects elsewhere for intronic
///   non-coding offsets).
fn nr_representative_genome_pos(
    cdot_tx: &CdotTranscript,
    base: i64,
    offset: Option<i64>,
    utr3: bool,
) -> Result<u64, FerroError> {
    if offset.is_some() || utr3 || base < 1 {
        return cdot_tx
            .exons
            .first()
            .map(|e| e[0])
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: "transcript has no exons; cannot derive stab position".to_string(),
            });
    }
    let tx_pos_0based = (base - 1) as u64;
    cdot_tx
        .tx_to_genome(tx_pos_0based)
        .ok_or_else(|| FerroError::InvalidCoordinates {
            msg: format!(
                "tx position {} not in any exon of {} (cannot derive stab position)",
                base, cdot_tx.contig
            ),
        })
}

/// Resolve an [`UncertainBoundary<T>`] to an owned `T`, distinguishing the two
/// `None` returns of [`UncertainBoundary::inner`]:
///   - `Single(Mu::Unknown)` (parsed `?` sentinel) → [`FerroError::UnsupportedProjection`]
///   - `Range { .. }`        (parsed `(a_b)` range) → [`FerroError::InvalidCoordinates`]
///
/// Without this split, a parsed `c.?` slips through as `InvalidCoordinates` and
/// contradicts the documented projection contract. Shared by both projection
/// entry points — `project_to_genomic` (genome-pivot path) and
/// `project_coding_direct` (direct bare-NM_ path) — so they classify `?` vs.
/// ranges identically.
fn resolve_uncertain_boundary<T: Copy>(
    boundary: &crate::hgvs::interval::UncertainBoundary<T>,
    coord_label: &str,
    end_label: &str,
) -> Result<T, FerroError> {
    use crate::hgvs::interval::UncertainBoundary;
    use crate::hgvs::uncertainty::Mu;
    match boundary {
        UncertainBoundary::Single(Mu::Certain(t) | Mu::Uncertain(t)) => Ok(*t),
        UncertainBoundary::Single(Mu::Unknown) => Err(FerroError::UnsupportedProjection {
            reason: format!("cannot resolve unknown position (`?`) at {coord_label} {end_label}"),
        }),
        UncertainBoundary::Range { .. } => Err(FerroError::InvalidCoordinates {
            msg: format!("{coord_label} interval {end_label} is a range"),
        }),
    }
}

/// Split a transcript accession into its version-stripped base and version:
/// `"NM_000532.5" → ("NM_000532", Some(5))`, `"NM_000532" → ("NM_000532", None)`.
fn split_accession_version(id: &str) -> (&str, Option<u32>) {
    match id.rsplit_once('.') {
        Some((base, ver)) if !ver.is_empty() && ver.bytes().all(|b| b.is_ascii_digit()) => {
            (base, ver.parse::<u32>().ok())
        }
        _ => (id, None),
    }
}

/// Whether a transcript id is a *predicted* RefSeq model (`XM_`/`XR_`), as
/// opposed to a curated `NM_`/`NR_` transcript.
fn is_predicted_model(id: &str) -> bool {
    id.starts_with("XM_") || id.starts_with("XR_")
}

/// Whether `id` names a transcript reading frame (RefSeq `NM_`/`NR_`/`XM_`/`XR_`
/// or Ensembl `ENST`), gating the #785 silent version-substitution refusal to
/// genuine transcript versions. Mirrors the normalizer's predicate of the same
/// name: a genomic/gene reference (`NG_`/`NC_`/`NW_`, `ENSG`) is *not* a
/// transcript version and must not be gated as one (e.g. an unrewritten
/// `NG_(GENE):c.…` selector).
fn is_transcript_accession(id: &str) -> bool {
    id.starts_with("NM_")
        || id.starts_with("NR_")
        || id.starts_with("XM_")
        || id.starts_with("XR_")
        || id.starts_with("ENST")
}

/// Apply the `project_*_all` enumeration policy (#656) to a priority-ordered
/// list of overlapping transcript ids:
///
/// 1. **Collapse superseded versions** — keep only the highest version per base
///    accession (drop `NM_000532.4` when `NM_000532.5` overlaps the same locus).
/// 2. **Prefer curated transcripts** — drop predicted `XM_`/`XR_` models when a
///    curated `NM_`/`NR_` transcript also covers the locus (matching mutalyzer's
///    curated set), but keep predicted models when they are the *sole* coverage,
///    so a locus a predicted model would have covered never returns empty.
///
/// Input order (clinical priority) is preserved among the kept ids.
pub(crate) fn select_enumerated_transcript_ids<'a>(ids: &[&'a str]) -> Vec<&'a str> {
    use std::collections::HashMap;
    // 1. Highest version seen per base accession.
    let mut max_ver: HashMap<&str, u32> = HashMap::new();
    for id in ids {
        let (base, ver) = split_accession_version(id);
        if let Some(v) = ver {
            max_ver
                .entry(base)
                .and_modify(|m| *m = (*m).max(v))
                .or_insert(v);
        }
    }
    let highest: Vec<&str> = ids
        .iter()
        .copied()
        .filter(|id| {
            let (base, ver) = split_accession_version(id);
            match ver {
                Some(v) => max_ver.get(base) == Some(&v),
                // A bare (unversioned) accession is superseded by any versioned
                // form of the same base; keep it only when no versioned id for
                // that base was seen.
                None => !max_ver.contains_key(base),
            }
        })
        .collect();
    // 2. Curated-preferred, never-empty.
    if highest.iter().any(|id| !is_predicted_model(id)) {
        highest
            .into_iter()
            .filter(|id| !is_predicted_model(id))
            .collect()
    } else {
        highest
    }
}

#[cfg(test)]
mod enumeration_policy_tests {
    use super::select_enumerated_transcript_ids;

    #[test]
    fn collapses_superseded_versions() {
        // The #656 example: both versions of two transcripts → keep the highest.
        let ids = [
            "NM_000532.4",
            "NM_000532.5",
            "NM_001178014.1",
            "NM_001178014.2",
        ];
        assert_eq!(
            select_enumerated_transcript_ids(&ids),
            vec!["NM_000532.5", "NM_001178014.2"]
        );
    }

    #[test]
    fn drops_predicted_models_when_curated_present() {
        let ids = ["NM_000532.5", "XM_011512873.2", "XM_005247508.1"];
        assert_eq!(select_enumerated_transcript_ids(&ids), vec!["NM_000532.5"]);
    }

    #[test]
    fn keeps_predicted_models_when_sole_coverage() {
        // No curated transcript overlaps → keep predicted models (never empty).
        let ids = ["XM_005247508.1", "XM_011512873.2"];
        assert_eq!(
            select_enumerated_transcript_ids(&ids),
            vec!["XM_005247508.1", "XM_011512873.2"]
        );
    }

    #[test]
    fn collapses_predicted_versions_then_drops_when_curated_present() {
        let ids = ["NM_000532.5", "XM_011512873.1", "XM_011512873.2"];
        assert_eq!(select_enumerated_transcript_ids(&ids), vec!["NM_000532.5"]);
    }

    #[test]
    fn preserves_priority_order_among_kept() {
        let ids = ["NM_B.1", "NM_A.2", "NM_A.1"];
        assert_eq!(
            select_enumerated_transcript_ids(&ids),
            vec!["NM_B.1", "NM_A.2"]
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::cdot::{CdotMapper, CdotTranscript};
    use crate::data::projection::Projector;
    use crate::hgvs::variant::{Accession, AllelePhase, AlleleVariant};
    use crate::reference::mock::MockProvider;
    use crate::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
    use crate::reference::Strand as ProvStrand;
    use std::sync::OnceLock;

    #[test]
    fn lrg_protein_accession_swaps_t_to_p() {
        assert_eq!(lrg_protein_accession("LRG_1t1").as_deref(), Some("LRG_1p1"));
        assert_eq!(
            lrg_protein_accession("LRG_24t2").as_deref(),
            Some("LRG_24p2")
        );
        assert_eq!(
            lrg_protein_accession("LRG_1t10").as_deref(),
            Some("LRG_1p10")
        );
        // Not an LRG transcript accession:
        assert_eq!(lrg_protein_accession("NM_000077.4"), None);
        assert_eq!(lrg_protein_accession("LRG_24"), None); // bare genomic LRG
        assert_eq!(lrg_protein_accession("LRG_1p1"), None); // already a protein
        assert_eq!(lrg_protein_accession("NP_000068.1"), None);
    }

    fn make_test_provider_and_projector() -> (Projector, MockProvider) {
        // A single coding transcript on chr1, plus strand.
        //
        // Genomic layout (cdot 0-based coords):
        //   Exon: genome [1000, 1009), tx [0, 9), so 9 bases total.
        //   cds_start = 0 (0-based, no 5'UTR), cds_end = 9.
        //
        // Sequence: "ATGCGCTAA" = Met-Arg-Stop (3 codons).
        //
        // Coordinate mapping (cdot genome 0-based → HGVS c. 1-based):
        //   g.1000 → tx_pos=0 → c.1 (A = first base of Met)
        //   g.1001 → tx_pos=1 → c.2 (T)
        //   g.1002 → tx_pos=2 → c.3 (G)
        //   g.1003 → tx_pos=3 → c.4 (C ← ref base for test substitution)
        //   g.1004 → tx_pos=4 → c.5 (G)
        //   ...
        //
        // c.4C>A: codon 2 CGC (Arg) → AGC (Ser), missense.
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_TEST.1".to_string(),
            CdotTranscript {
                gene_name: Some("TESTGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                // [genome_start(0-based), genome_end(0-based excl), tx_start(1-based), tx_end(1-based)]
                exons: vec![[1000, 1009, 0, 9]],
                cds_start: Some(0), // 0-based: CDS starts at tx_pos 0 (no 5'UTR)
                cds_end: Some(9),   // 0-based exclusive: CDS ends at tx_pos 9
                gene_id: None,
                protein: Some("NP_TEST.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        // A second transcript with a NON-ATG start codon ("CTG…"), for the #771
        // initiation-codon-on-a-non-AUG-transcript path. Bare-NM_ direct
        // projection needs only the provider transcript + cdot version check.
        cdot.add_transcript(
            "NM_NOATG.1".to_string(),
            CdotTranscript {
                gene_name: Some("NOATGGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[2000, 2009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: Some("NP_NOATG.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        // A third transcript with a genuine off-the-start drift codon ("AGT…",
        // not a recognized near-cognate initiator), for the #625/#780 boundary:
        // a downstream variant on this must STILL decline.
        cdot.add_transcript(
            "NM_DRIFT.1".to_string(),
            CdotTranscript {
                gene_name: Some("DRIFTGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[3000, 3009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: Some("NP_DRIFT.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        // A fourth transcript whose CDS start codon is UNREADABLE: the provider
        // copy carries no sequence (`sequence: None`), so `read_cds_start_codon`
        // returns an error rather than three bases. Used to pin the #780/#771
        // strict-form-choice default: an init-codon edit on an unreadable start
        // must fall to the conservative `p.?`, not `p.(Met1?)`.
        cdot.add_transcript(
            "NM_NOSEQ.1".to_string(),
            CdotTranscript {
                gene_name: Some("NOSEQGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[4000, 4009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: Some("NP_NOSEQ.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        // Unreadable start codon: a transcript whose sequence is absent
        // (`sequence: None`), so `read_cds_start_codon` returns
        // `ProteinSequenceUnavailable` and the start codon cannot be observed.
        // The CDS coords are otherwise well-formed (matching the cdot copy).
        provider.add_transcript(Transcript {
            id: "NM_NOSEQ.1".to_string(),
            gene_symbol: Some("NOSEQGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(4000),
            genomic_end: Some(4008),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        // Non-ATG transcript: "CTGCGCTAA" (CTG start), cds_start=1.
        provider.add_transcript(Transcript {
            id: "NM_NOATG.1".to_string(),
            gene_symbol: Some("NOATGGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("CTGCGCTAA".to_string()),
            cds_start: Some(1),
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(2000),
            genomic_end: Some(2008),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        // Genuine off-the-start drift: "AGTCGCTAA" (AGT start — not ATG and not a
        // recognized near-cognate initiator), cds_start=1. A downstream variant
        // on this must still decline under #780 (the #625 protection).
        provider.add_transcript(Transcript {
            id: "NM_DRIFT.1".to_string(),
            gene_symbol: Some("DRIFTGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("AGTCGCTAA".to_string()),
            cds_start: Some(1),
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(3000),
            genomic_end: Some(3008),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        // Transcript: sequence "ATGCGCTAA", cds_start=1 (1-based, first base).
        provider.add_transcript(Transcript {
            id: "NM_TEST.1".to_string(),
            gene_symbol: Some("TESTGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("ATGCGCTAA".to_string()),
            cds_start: Some(1), // 1-based inclusive per Transcript convention
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1008),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        // Genomic sequence: 999 N's + "ATGCGCTAA" + 100 N's.
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ATGCGCTAA", suffix));
        (projector, provider)
    }

    /// Provider with two longer plus-strand CDS transcripts for the #855 in-cis
    /// protein-combination tests:
    /// - `NM_DELINS.1` = "ATGGATTATTAA" (Met-Asp-Tyr-Stop) — adjacent-codon subs
    ///   and an in-cis frameshift.
    /// - `NM_SEP.1` = "ATGGATTATTGCTAA" (Met-Asp-Tyr-Cys-Stop) — subs separated by
    ///   an unchanged residue.
    fn make_multicodon_provider_and_projector() -> (Projector, MockProvider) {
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_DELINS.1".to_string(),
            CdotTranscript {
                gene_name: Some("DELINSGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1012, 0, 12]],
                cds_start: Some(0),
                cds_end: Some(12),
                gene_id: None,
                protein: Some("NP_DELINS.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        cdot.add_transcript(
            "NM_SEP.1".to_string(),
            CdotTranscript {
                gene_name: Some("SEPGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[2000, 2015, 0, 15]],
                cds_start: Some(0),
                cds_end: Some(15),
                gene_id: None,
                protein: Some("NP_SEP.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_DELINS.1".to_string(),
            gene_symbol: Some("DELINSGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("ATGGATTATTAA".to_string()),
            cds_start: Some(1),
            cds_end: Some(12),
            exons: vec![Exon::new(1, 1, 12)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1011),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        provider.add_transcript(Transcript {
            id: "NM_SEP.1".to_string(),
            gene_symbol: Some("SEPGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("ATGGATTATTGCTAA".to_string()),
            cds_start: Some(1),
            cds_end: Some(15),
            exons: vec![Exon::new(1, 1, 15)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(2000),
            genomic_end: Some(2014),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        // Genomic sequence so the genomic-projection step of each member resolves:
        // N*1000 + DELINS CDS + filler + SEP CDS + tail (plus strand → genomic
        // bases equal the transcript bases at each exon offset).
        let mut genome = String::new();
        genome.push_str(&"N".repeat(1000));
        genome.push_str("ATGGATTATTAA"); // [1000, 1012)
        genome.push_str(&"N".repeat(2000 - 1012));
        genome.push_str("ATGGATTATTGCTAA"); // [2000, 2015)
        genome.push_str(&"N".repeat(100));
        provider.add_genomic_sequence("chr1", genome);
        (projector, provider)
    }

    fn make_minus_strand_provider_and_projector() -> (Projector, MockProvider) {
        // Same 9bp CDS on chr1, but transcript is on the minus strand.
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_TEST_MINUS.1".to_string(),
            CdotTranscript {
                gene_name: Some("TESTGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Minus,
                exons: vec![[1000, 1009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: Some("NP_TEST_MINUS.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_TEST_MINUS.1".to_string(),
            gene_symbol: Some("TESTGENE".to_string()),
            strand: TxStrand::Minus,
            sequence: Some("ATGCGCTAA".to_string()), // CDS as the transcript reads it
            cds_start: Some(1),                      // 1-based inclusive per Transcript convention
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1008),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}TTAGCGCAT{}", prefix, suffix));
        (projector, provider)
    }

    /// Build a two-transcript setup for project_all tests.
    ///
    /// NM_TX1.1: chr1 [1000,1009), plus strand, 9bp CDS "ATGCGCTAA"
    /// NM_TX2.1: chr1 [1000,1009), plus strand, 9bp CDS "ATGCGCTAA" (same region)
    ///           NM_TX2.1 is registered as MANE Select so it sorts first.
    fn make_two_transcript_setup() -> (Projector, MockProvider) {
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_TX1.1".to_string(),
            CdotTranscript {
                gene_name: Some("GENE1".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: Some("NP_TX1.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        cdot.add_transcript(
            "NM_TX2.1".to_string(),
            CdotTranscript {
                gene_name: Some("GENE1".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: Some("NP_TX2.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot)
            // Make NM_TX2.1 the MANE Select so it sorts first.
            .with_mane(vec!["NM_TX2.1".to_string()], vec![]);

        let mut provider = MockProvider::new();
        for id in ["NM_TX1.1", "NM_TX2.1"] {
            provider.add_transcript(Transcript {
                id: id.to_string(),
                gene_symbol: Some("GENE1".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("ATGCGCTAA".to_string()),
                cds_start: Some(1),
                cds_end: Some(9),
                exons: vec![Exon::new(1, 1, 9)],
                chromosome: Some("chr1".to_string()),
                genomic_start: Some(1000),
                genomic_end: Some(1008),
                genome_build: Default::default(),
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                protein_id: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
            });
        }
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ATGCGCTAA", suffix));
        (projector, provider)
    }

    /// Like [`make_two_transcript_setup`] but seeds three transcripts on the
    /// *same* locus to exercise the #656 enumeration policy end-to-end through
    /// the `project_*_all` fan-out:
    /// - `NM_TX1.1` and `NM_TX1.2` — superseded/current versions of one curated
    ///   base accession (only `.2` should survive the version collapse);
    /// - `XM_TX9.1` — a predicted model that should be dropped because a curated
    ///   transcript covers the same locus.
    ///
    /// So a curated set of exactly `{NM_TX1.2}` is expected, versus the three
    /// overlapping records cdot reports.
    #[cfg(test)]
    fn make_curated_enumeration_setup() -> (Projector, MockProvider) {
        let mut cdot = CdotMapper::new();
        for id in ["NM_TX1.1", "NM_TX1.2", "XM_TX9.1"] {
            cdot.add_transcript(
                id.to_string(),
                CdotTranscript {
                    gene_name: Some("GENE1".to_string()),
                    contig: "chr1".to_string(),
                    strand: ProvStrand::Plus,
                    exons: vec![[1000, 1009, 0, 9]],
                    cds_start: Some(0),
                    cds_end: Some(9),
                    gene_id: None,
                    protein: None,
                    exon_cigars: Vec::new(),
                },
            );
        }
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        for id in ["NM_TX1.1", "NM_TX1.2", "XM_TX9.1"] {
            provider.add_transcript(Transcript {
                id: id.to_string(),
                gene_symbol: Some("GENE1".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("ATGCGCTAA".to_string()),
                cds_start: Some(1),
                cds_end: Some(9),
                exons: vec![Exon::new(1, 1, 9)],
                chromosome: Some("chr1".to_string()),
                genomic_start: Some(1000),
                genomic_end: Some(1008),
                genome_build: Default::default(),
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                protein_id: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
            });
        }
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ATGCGCTAA", suffix));
        (projector, provider)
    }

    /// Like [`make_two_transcript_setup`] but with the transcripts placed on an
    /// `NC_` accession contig (`NC_000001.11`) rather than the `chr1` alias, so
    /// an `NG_` `GenomicPlacement` (whose `nc` field is an `NC_` accession) maps
    /// onto the same contig the projector indexes. Used to exercise the
    /// NG_-genomic-input fan-out path (#480 inverse).
    #[cfg(test)]
    fn make_nc_two_transcript_setup() -> (Projector, MockProvider) {
        let mut cdot = CdotMapper::new();
        for (id, np) in [("NM_TX1.1", "NP_TX1.1"), ("NM_TX2.1", "NP_TX2.1")] {
            cdot.add_transcript(
                id.to_string(),
                CdotTranscript {
                    gene_name: Some("GENE1".to_string()),
                    contig: "NC_000001.11".to_string(),
                    strand: ProvStrand::Plus,
                    exons: vec![[1000, 1009, 0, 9]],
                    cds_start: Some(0),
                    cds_end: Some(9),
                    gene_id: None,
                    protein: Some(np.to_string()),
                    exon_cigars: Vec::new(),
                },
            );
        }
        let projector = Projector::new(cdot).with_mane(vec!["NM_TX2.1".to_string()], vec![]);

        let mut provider = MockProvider::new();
        for id in ["NM_TX1.1", "NM_TX2.1"] {
            provider.add_transcript(Transcript {
                id: id.to_string(),
                gene_symbol: Some("GENE1".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("ATGCGCTAA".to_string()),
                cds_start: Some(1),
                cds_end: Some(9),
                exons: vec![Exon::new(1, 1, 9)],
                chromosome: Some("NC_000001.11".to_string()),
                genomic_start: Some(1000),
                genomic_end: Some(1008),
                genome_build: Default::default(),
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                protein_id: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
            });
        }
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence(
            "NC_000001.11",
            format!("{}{}{}", prefix, "ATGCGCTAA", suffix),
        );
        // The NG_ parent's own sequence (NG_900.1 base 1..9 == NC 1000..1008),
        // so normalization of an NG_-relative input resolves its reference.
        provider.add_genomic_sequence("NG_900.1", "ATGCGCTAA".to_string());
        (projector, provider)
    }

    /// A bare `NG_<n>:g.` input must enumerate the transcripts overlapping that
    /// NG_ position by translating the NG_-relative coordinate into the
    /// chromosome frame via the parent's [`GenomicPlacement`] (#480 inverse),
    /// then projecting onto each overlapping transcript. Pre-fix
    /// `extract_contig_and_pos` used the bare NG_ accession as the contig, which
    /// the projector cannot map, so `project_variant_all` returned `[]`.
    #[test]
    fn project_variant_all_enumerates_transcripts_for_ng_genomic_input() {
        let (projector, mut provider) = make_nc_two_transcript_setup();
        // NG_900.1 placed on NC_000001.11, plus strand, NG_ base 1 == nc 1000.
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1008,
                strand: crate::reference::Strand::Plus,
            },
        );
        let vp = VariantProjector::new(projector, provider);

        // NG_900.1:g.4C>A — NG_ base 4 == nc 1003 == c.4 of the CDS.
        let results = vp
            .project_all("NG_900.1:g.4C>A")
            .expect("project_all should enumerate transcripts for an NG_ genomic input");

        assert_eq!(
            results.len(),
            2,
            "expected both transcripts overlapping the NG_ position, got {:?}",
            results.iter().map(|p| &p.transcript_id).collect::<Vec<_>>()
        );
    }

    /// #637: a legacy gene-model selector on a genomic reference
    /// (`NG_900.1(GENE1_v001):c.…`) resolves to the gene's transcript via the
    /// normalize-first step in projection, so it no longer fails with
    /// "Reference not found: NG_900.1" — it projects onto the resolved `NM_`.
    #[test]
    fn legacy_gene_model_selector_projects_via_normalize_rewrite() {
        let (projector, mut provider) = make_nc_two_transcript_setup();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1008,
                strand: crate::reference::Strand::Plus,
            },
        );
        provider.add_legacy_gene_model("GENE1", "NM_TX2.1");
        let vp = VariantProjector::new(projector, provider);

        let proj = vp
            .project("NG_900.1(GENE1_v001):c.4C>A", "NM_TX2.1")
            .expect("legacy gene-model selector must resolve and project, not error on NG_");
        let coding = proj
            .coding
            .as_ref()
            .expect("coding axis present")
            .to_string();
        assert!(
            coding.contains("NM_TX2.1"),
            "the GENE1_v001 selector must resolve to the transcript: {coding}"
        );
    }

    #[test]
    fn legacy_gene_model_selector_projects_against_self_derived_target() {
        // Realistic CLI / conformance-harness path (#783): the projection target
        // is derived from the *raw* input accession via `transcript_accession()`.
        // For an `NG_(GENE_v001)` legacy selector the transcript slot is the
        // genomic ref itself (`NG_900.1`), not a transcript — unlike `NG_(NM_)`,
        // whose `transcript_accession()` already returns the inner `NM_`.
        // Normalization resolves the selector to `NM_TX2.1` (#637); without
        // reconciling the target against the resolved transcript, this
        // self-projection tripped the transcript_id-mismatch guard.
        let (projector, mut provider) = make_nc_two_transcript_setup();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1008,
                strand: crate::reference::Strand::Plus,
            },
        );
        provider.add_legacy_gene_model("GENE1", "NM_TX2.1");
        let vp = VariantProjector::new(projector, provider);

        let variant = crate::parse_hgvs("NG_900.1(GENE1_v001):c.4C>A").unwrap();
        // Derive the target exactly as the CLI and conformance harness do.
        let target = variant.accession().unwrap().transcript_accession();
        assert_eq!(
            target, "NG_900.1",
            "the raw LOVD selector form's transcript slot is the genomic ref"
        );

        let proj = vp
            .project_variant(&variant, &target)
            .expect("self-projection of a legacy selector must resolve, not mismatch");
        let coding = proj
            .coding
            .as_ref()
            .expect("coding axis present")
            .to_string();
        assert!(
            coding.contains("NM_TX2.1"),
            "the GENE1_v001 selector must resolve to the transcript: {coding}"
        );
    }

    /// #652: a range-reference insertion (`ins<start>_<end>`) on an `NG_`
    /// genomic input must resolve to the literal parent bases — the range-ref
    /// endpoints are in the `NG_` parent frame and, if left unresolved through
    /// the de-anchor into the `NC_` frame, would be fetched against the `NC_`
    /// accession and come back as `insNNN…`.
    #[test]
    fn deanchor_resolves_ng_range_reference_insertion_to_literal() {
        let mut provider = MockProvider::new();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 30000,
                strand: crate::reference::Strand::Plus,
            },
        );
        // NG_900.1 parent sequence: bases 4300..=4320 (1-based) are the insert.
        let insert = "GTCCTGTGCTCATTATCTGGC"; // 21 bases
        let mut seq = vec![b'A'; 4299];
        seq.extend_from_slice(insert.as_bytes());
        provider.add_genomic_sequence("NG_900.1", String::from_utf8(seq).unwrap());

        let vp = VariantProjector::new(Projector::new(CdotMapper::new()), provider);
        let variant = crate::parse_hgvs("NG_900.1:g.5207_5208ins4300_4320").unwrap();
        let (deanchored, parent) = vp.deanchor_genomic_parent_input(&variant);

        assert_eq!(
            parent.as_ref().map(|a| a.full()),
            Some("NG_900.1".to_string())
        );
        let rendered = deanchored.to_string();
        assert!(
            rendered.contains(&format!("ins{insert}")),
            "expected the literal parent bases, got: {rendered}"
        );
        assert!(
            !rendered.contains("insN"),
            "must not emit insNNN (unresolved range reference): {rendered}"
        );
    }

    /// #652, minus-strand placement: the range reference is resolved to the
    /// literal parent bases in the parent frame, then the existing strand
    /// transform reverse-complements that literal for the `NC_` frame (the same
    /// path a literal `insACGT` already takes), never `insNNN`.
    #[test]
    fn deanchor_resolves_ng_range_reference_insertion_minus_strand() {
        let mut provider = MockProvider::new();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 30000,
                strand: crate::reference::Strand::Minus,
            },
        );
        let insert = "GTCCTGTGCTCATTATCTGGC"; // parent bases 4300..=4320
        let revcomp = "GCCAGATAATGAGCACAGGAC"; // reverse-complement of `insert`
        let mut seq = vec![b'A'; 4299];
        seq.extend_from_slice(insert.as_bytes());
        provider.add_genomic_sequence("NG_900.1", String::from_utf8(seq).unwrap());

        let vp = VariantProjector::new(Projector::new(CdotMapper::new()), provider);
        let variant = crate::parse_hgvs("NG_900.1:g.5207_5208ins4300_4320").unwrap();
        let (deanchored, _) = vp.deanchor_genomic_parent_input(&variant);
        let rendered = deanchored.to_string();
        assert!(
            rendered.contains(&format!("ins{revcomp}")),
            "expected the reverse-complemented literal, got: {rendered}"
        );
        assert!(
            !rendered.contains("insN"),
            "must not emit insNNN: {rendered}"
        );
    }

    /// #652, `delins` range reference: `canonicalize_insertion_expand` and the
    /// production doc-comment explicitly handle `delins<start>_<end>`, so a
    /// `delins`-shaped range-ref input must also resolve to the literal parent
    /// bases (never `insNNN…`), exercising the named-but-otherwise-untested edit
    /// kind alongside the plain `Insertion` cases above.
    #[test]
    fn deanchor_resolves_ng_range_reference_delins_to_literal() {
        let mut provider = MockProvider::new();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 30000,
                strand: crate::reference::Strand::Plus,
            },
        );
        // NG_900.1 parent sequence: bases 4300..=4320 (1-based) are the insert.
        let insert = "GTCCTGTGCTCATTATCTGGC"; // 21 bases
        let mut seq = vec![b'A'; 4299];
        seq.extend_from_slice(insert.as_bytes());
        provider.add_genomic_sequence("NG_900.1", String::from_utf8(seq).unwrap());

        let vp = VariantProjector::new(Projector::new(CdotMapper::new()), provider);
        let variant = crate::parse_hgvs("NG_900.1:g.5207_5208delins4300_4320").unwrap();
        let (deanchored, parent) = vp.deanchor_genomic_parent_input(&variant);

        assert_eq!(
            parent.as_ref().map(|a| a.full()),
            Some("NG_900.1".to_string())
        );
        let rendered = deanchored.to_string();
        assert!(
            rendered.contains(&format!("delins{insert}")),
            "expected the literal parent bases, got: {rendered}"
        );
        assert!(
            !rendered.contains("insN"),
            "must not emit insNNN (unresolved range reference): {rendered}"
        );
    }

    /// For an `NG_`-genomic input, each per-transcript coding and protein
    /// description must be framed under the original `NG_` parent
    /// (`NG_900.1(NM_TX1.1):c.…` / `NG_900.1(NP_TX1.1):p.…`), matching the
    /// mutalyzer corpus, rather than the bare transcript or the de-anchored
    /// `NC_` accession used internally for projection.
    #[test]
    fn project_variant_all_frames_coding_protein_under_ng_parent() {
        let (projector, mut provider) = make_nc_two_transcript_setup();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1008,
                strand: crate::reference::Strand::Plus,
            },
        );
        let vp = VariantProjector::new(projector, provider);

        let results = vp
            .project_all("NG_900.1:g.4C>A")
            .expect("project_all should succeed for an NG_ genomic input");

        let coding_strings: Vec<String> = results
            .iter()
            .filter_map(|r| r.coding.as_ref().map(|c| c.to_string()))
            .collect();
        assert!(
            !coding_strings.is_empty(),
            "expected at least one coding description"
        );
        for c in &coding_strings {
            assert!(
                c.starts_with("NG_900.1("),
                "coding description should be framed under the NG_ parent, got {c:?}"
            );
        }
        for r in &results {
            if let Some(p) = r.protein.as_ref().map(|p| p.to_string()) {
                assert!(
                    p.starts_with("NG_900.1("),
                    "protein description should be framed under the NG_ parent, got {p:?}"
                );
            }
            // The genomic description is re-anchored back into the NG_ frame
            // (NG_ base 4 == nc 1003), so it reads `NG_900.1:g.4…`.
            if let Some(g) = r.genomic.as_ref().map(|g| g.to_string()) {
                assert!(
                    g.starts_with("NG_900.1:g.4"),
                    "genomic description should be in the NG_ frame, got {g:?}"
                );
            }
        }
    }

    /// A transcript-coordinate input carrying an `NG_` `genomic_context`
    /// (`NG_900.1(NM_TX1.1):c.4C>A`) must, like a bare `NG_:g.` input, enumerate
    /// the overlapping transcripts and frame each description under the `NG_`
    /// parent — not return only the input transcript, bare (#646 / Mode A).
    #[test]
    fn project_variant_all_frames_coding_input_with_ng_context() {
        let (projector, mut provider) = make_nc_two_transcript_setup();
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1008,
                strand: crate::reference::Strand::Plus,
            },
        );
        let vp = VariantProjector::new(projector, provider);

        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::CdsVariant;
        // NG_900.1(NM_TX1.1):c.4C>A — c. input carrying an NG_ genomic_context.
        let cds = CdsVariant {
            accession: parse_accession("NM_TX1.1").with_genomic_context(Accession::new(
                "NG",
                "900",
                Some(1),
            )),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let v = HgvsVariant::Cds(cds);
        let results = vp
            .project_variant_all(&v)
            .expect("project_all should succeed for a c. input with NG_ context");

        // Both overlapping transcripts are enumerated.
        let txs: Vec<&str> = results.iter().map(|r| r.transcript_id.as_str()).collect();
        assert!(
            txs.contains(&"NM_TX1.1") && txs.contains(&"NM_TX2.1"),
            "expected both transcripts enumerated, got {txs:?}"
        );
        // Every coding description is framed under the NG_ parent.
        for r in &results {
            if let Some(c) = r.coding.as_ref().map(|c| c.to_string()) {
                assert!(
                    c.starts_with("NG_900.1("),
                    "coding should be framed under the NG_ parent, got {c:?}"
                );
            }
        }
    }

    /// Coexistence of the #646 cross-isoform enumeration with the #655 graceful
    /// decline. A `c.` input carrying an `NG_` `genomic_context` is de-anchored
    /// into the chromosome frame and enumerated against the overlapping
    /// transcripts, each re-framed under the parent (#646). Here the registered
    /// placement spans a region disjoint from where the transcripts actually map
    /// (chr 1003), so the genomic axis cannot be re-anchored into the parent
    /// frame. Rather than failing the whole projection (the #646 enumerate path)
    /// or stamping the parent accession on an out-of-frame chromosome coordinate
    /// (invalid HGVS, #655), `project_variant_all` returns the transcript-relative
    /// coding/protein forms framed under the parent and drops only the genomic
    /// axis (`None`). This locks in that #646 (enumerate) and #655 (graceful
    /// decline) coexist: `frame_projection_owned` degrades the unframable axis
    /// via `reanchor_genome_to_parent` instead of hard-erroring.
    #[test]
    fn project_variant_all_drops_genomic_axis_but_still_enumerates_under_ng_parent() {
        let (projector, mut provider) = make_nc_two_transcript_setup();
        // Placement exists (so de-anchor enumerates) but its span [2000, 2010]
        // is disjoint from chr 1003 where the transcripts map, so the genomic
        // re-anchor `nc_to_parent(1003)` declines and the axis is dropped.
        provider.add_genomic_placement(
            "NG_900.1",
            crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 2000,
                nc_end: 2010,
                strand: crate::reference::Strand::Plus,
            },
        );
        let vp = VariantProjector::new(projector, provider);

        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::CdsVariant;
        let cds = CdsVariant {
            accession: parse_accession("NM_TX1.1").with_genomic_context(Accession::new(
                "NG",
                "900",
                Some(1),
            )),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let v = HgvsVariant::Cds(cds);
        let results = vp
            .project_variant_all(&v)
            .expect("enumerate path must degrade the unframable genomic axis, not hard-error");

        // Enumeration still finds the overlapping transcripts (#646 path alive).
        let txs: Vec<&str> = results.iter().map(|r| r.transcript_id.as_str()).collect();
        assert!(
            txs.contains(&"NM_TX1.1") && txs.contains(&"NM_TX2.1"),
            "expected both transcripts enumerated, got {txs:?}"
        );
        for r in &results {
            // The genomic axis is dropped — it cannot be expressed in the parent
            // frame — rather than emitting a chromosome coordinate under NG_900.1.
            assert!(
                r.genomic.is_none(),
                "genomic axis should be dropped when it cannot be re-anchored, got: {:?}",
                r.genomic
            );
            // The coding axis survives, re-framed under the NG_ parent.
            if let Some(c) = r.coding.as_ref().map(|c| c.to_string()) {
                assert!(
                    c.starts_with("NG_900.1("),
                    "coding should be framed under the NG_ parent, got {c:?}"
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    // Cache identity tests
    // -------------------------------------------------------------------------

    /// Build a minimal coding HgvsVariant pointing at `transcript_id` for use
    /// in `cached_ref_translation` tests. The variant's `accession` carries
    /// `genomic_context = ctx`, so `parent_key_for` returns the same string
    /// the production code would derive at projection time.
    #[cfg(test)]
    fn make_coding_variant(transcript_id: &str, ctx: Option<Accession>) -> HgvsVariant {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};

        let mut accession = parse_accession(transcript_id);
        if let Some(parent) = ctx {
            accession = accession.with_genomic_context(parent);
        }
        HgvsVariant::Cds(CdsVariant {
            accession,
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        })
    }

    /// Two consecutive `cached_ref_translation` calls for the same
    /// `(transcript_id, parent_accession)` pair must return the same
    /// `Arc<RefProteinBundle>` (pointer-equal). Catches any future
    /// regression that accidentally rebuilds the bundle per call.
    #[test]
    fn cached_ref_translation_returns_same_arc() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let tx = vp
            .cached_get_transcript("NM_TEST.1")
            .expect("transcript fetch");
        let variant = make_coding_variant("NM_TEST.1", None);
        let a = vp
            .cached_ref_translation(&variant, "NM_TEST.1", &tx)
            .expect("first translation");
        let b = vp
            .cached_ref_translation(&variant, "NM_TEST.1", &tx)
            .expect("second translation");
        assert!(
            Arc::ptr_eq(&a, &b),
            "ref-protein cache must return the same Arc for repeated lookups"
        );
        // ATGCGCTAA = Met-Arg-Ter; translate_full_cds stops before the Ter.
        assert_eq!(a.ref_protein.len(), 2);
        assert_eq!(a.ref_protein_with_stop.len(), 3);
    }

    /// Two `cached_ref_translation` lookups for the same `transcript_id`
    /// but different parent `genomic_context` accessions must produce
    /// distinct cache entries (different `Arc`s). Without parent-aware
    /// keying the second call would reuse the bundle built for the first
    /// parent, which can poison indel protein predictions once the
    /// underlying transcript resolves to a different cdot build under
    /// each parent (issue #332).
    #[test]
    fn cached_ref_translation_is_parent_aware() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let tx = vp
            .cached_get_transcript("NM_TEST.1")
            .expect("transcript fetch");

        let v_none = make_coding_variant("NM_TEST.1", None);
        let v_ng1 = make_coding_variant("NM_TEST.1", Some(Accession::new("NG", "TEST", Some(1))));
        let v_ng2 = make_coding_variant("NM_TEST.1", Some(Accession::new("NG", "TEST", Some(2))));

        let b_none = vp
            .cached_ref_translation(&v_none, "NM_TEST.1", &tx)
            .expect("no-parent bundle");
        let b_ng1 = vp
            .cached_ref_translation(&v_ng1, "NM_TEST.1", &tx)
            .expect("NG_TEST.1 bundle");
        let b_ng2 = vp
            .cached_ref_translation(&v_ng2, "NM_TEST.1", &tx)
            .expect("NG_TEST.2 bundle");

        // Each parent identity must occupy its own cache slot.
        assert!(
            !Arc::ptr_eq(&b_none, &b_ng1),
            "no-parent and NG-parented entries must not collide"
        );
        assert!(
            !Arc::ptr_eq(&b_ng1, &b_ng2),
            "two distinct NG parents must produce distinct cache entries"
        );

        // Re-querying with the same parent must still hit the cache.
        let b_ng1_again = vp
            .cached_ref_translation(&v_ng1, "NM_TEST.1", &tx)
            .expect("NG_TEST.1 bundle (repeat)");
        assert!(
            Arc::ptr_eq(&b_ng1, &b_ng1_again),
            "repeated lookup with the same parent must hit the cache"
        );
    }

    // -------------------------------------------------------------------------
    // Existing tests (preserved)
    // -------------------------------------------------------------------------

    #[test]
    fn project_substitution_minus_strand_revcomps_ref_alt() {
        let (projector, provider) = make_minus_strand_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1005G>A", "NM_TEST_MINUS.1")
            .expect("minus-strand projection should succeed");
        let c = result
            .coding
            .as_ref()
            .expect("c. should be present")
            .to_string();
        assert!(
            c.contains(":c.4C>T"),
            "expected revcomp c.4C>T for minus strand, got: {}",
            c
        );
        let p = result
            .protein
            .as_ref()
            .expect("p. should be present")
            .to_string();
        assert_eq!(p, "NP_TEST_MINUS.1:p.(Arg2Cys)");
        assert!(!result.is_frameshift);
        assert!(!result.is_intronic);
        assert!(!result.is_utr);
    }

    fn make_intronic_test_data() -> (Projector, MockProvider) {
        // Two-exon coding transcript on chr1, plus strand.
        //
        //   Exon 1: genome [1000, 1010), tx [0, 10)  — 10 bp
        //   Intron: genome [1010, 2000)               — 990 bp
        //   Exon 2: genome [2000, 2010), tx [10, 20) — 10 bp (last 2 are pad)
        //
        // CDS: tx [0, 18) → 18 bp → 6 codons: ATG-CGC-AAA-GGG-TAA-CCC
        //   (Met-Arg-Lys-Gly-Stop-Pro)
        //
        // Test position g.1015 is in the intron, 6 bases after the end of exon 1
        // (exon 1 ends at genome 1010, exclusive; last exonic base is g.1009;
        //  distance = 1015-1010+1 = 6).  Maps to c.9+6 or c.10+6 depending on
        // whether tx_end is inclusive or exclusive in the boundary calculation.
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_INTR.1".to_string(),
            CdotTranscript {
                gene_name: Some("INTRGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1010, 0, 10], [2000, 2010, 10, 20]],
                cds_start: Some(0),
                cds_end: Some(18),
                gene_id: None,
                protein: Some("NP_INTR.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_INTR.1".to_string(),
            gene_symbol: Some("INTRGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("ATGCGCAAAGGGTAACCC".to_string()), // 18 bp
            cds_start: Some(1),
            cds_end: Some(18),
            exons: vec![Exon::new(1, 1, 10), Exon::new(2, 11, 20)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2009),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let exon1 = "ATGCGCAAAG"; // 10 bp (genomic positions 1000..1009)
        let intron = "N".repeat(990); // genomic positions 1010..1999 (990 bp)
        let exon2 = "GGTAACCCNN"; // 10 bp pad (genomic positions 2000..2009)
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{prefix}{exon1}{intron}{exon2}{suffix}"));
        (projector, provider)
    }

    #[test]
    fn project_intronic_substitution_no_protein() {
        let (projector, provider) = make_intronic_test_data();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1015A>G", "NM_INTR.1")
            .expect("intronic substitution should project to c. with offset");
        assert!(result.is_intronic, "expected is_intronic=true");
        assert!(result.protein.is_none(), "no p. for intronic substitutions");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(
            c.contains('+'),
            "expected intronic offset notation (e.g. c.9+6 or c.10+6), got: {}",
            c
        );
    }

    #[test]
    fn project_no_overlap_returns_transcript_not_overlapping() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let err = vp
            .project("NC_000001.11:g.5000A>G", "NM_TEST.1")
            .expect_err("should fail to project outside the transcript");
        assert!(
            matches!(err, FerroError::TranscriptNotOverlapping { .. }),
            "expected TranscriptNotOverlapping, got: {:?}",
            err
        );
    }

    // ------------------------------------------------------------------
    // Sequence-aware exon-indel correction (issue #644)
    // ------------------------------------------------------------------

    /// Build a projector + provider for a plus-strand transcript whose cdot
    /// record is **ungapped** (`exon_cigars` empty) but whose transcript and
    /// genome sequences genuinely differ by an indel — the cdot 0.2.32 /
    /// `NM_000532.5` (PCCB) failure shape, reduced to a CI-runnable fixture.
    ///
    /// Single exon, genome HGVS `[1000, 1012)` (12 bp) ↔ tx `[0, 11)` (11 bp),
    /// `cds_start = 0`. The true alignment (recoverable from the sequences) is:
    /// 2 genome-only bases at the exon start, then 10 matched bases, then 1
    /// transcript-only base at the end (net genome − tx = +1). So a genomic
    /// position in the matched region maps to a transcript position 2 lower than
    /// the naive cdot arithmetic would give.
    fn make_ungapped_indel_provider_and_projector() -> (Projector, MockProvider) {
        let common = "ATGCGCTAAC"; // 10 bases, matched region
        let genome_exon = format!("CA{common}"); // 12 bp: 2 genome-only + 10 matched
        let tx_exon = format!("{common}T"); // 11 bp: 10 matched + 1 tx-only
        debug_assert_eq!(genome_exon.len(), 12);
        debug_assert_eq!(tx_exon.len(), 11);

        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_GAP644.1".to_string(),
            CdotTranscript {
                gene_name: Some("GAPGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                // genome [1000, 1012) (HGVS-value), tx [0, 11) — lengths match
                // the FASTAs, but cdot records NO CIGAR (the bug shape).
                exons: vec![[1000, 1012, 0, 11]],
                cds_start: Some(0),
                cds_end: Some(11),
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_GAP644.1".to_string(),
            gene_symbol: Some("GAPGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some(tx_exon.clone()),
            cds_start: Some(1), // 1-based per Transcript convention
            cds_end: Some(11),
            exons: vec![Exon::new(1, 1, 11)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1011),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        // Genomic contig: 999 N's (so HGVS g.1000 == 0-based index 999) + the
        // exon sequence + trailing N's.
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(50);
        provider.add_genomic_sequence("chr1", format!("{prefix}{genome_exon}{suffix}"));
        (projector, provider)
    }

    #[test]
    fn project_corrects_genome_to_cds_across_unmodelled_indel() {
        // The reproducer class: a g. SNV in the matched region must project to a
        // c. coordinate corrected by the 2-bp genome-only prefix.
        //
        // First matched genome base is HGVS g.1002 (0-based exon offset 2). It
        // maps to tx offset 0 → c.1. Without the correction the naive arithmetic
        // would give tx offset 2 → c.3.
        let (projector, provider) = make_ungapped_indel_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        // g.1006 is exon offset 6 → naive tx 6 (c.7); corrected tx 4 (c.5).
        let result = vp
            .project("chr1:g.1006A>T", "NM_GAP644.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. present").to_string();
        assert!(
            c.contains(":c.5"),
            "expected corrected ':c.5' (not naive ':c.7') in '{c}'"
        );
    }

    #[test]
    fn project_first_matched_base_maps_to_c1() {
        let (projector, provider) = make_ungapped_indel_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // g.1002 is the first matched base → c.1.
        let result = vp
            .project("chr1:g.1002A>T", "NM_GAP644.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. present").to_string();
        assert!(c.contains(":c.1A>T"), "expected ':c.1A>T' in '{c}'");
    }

    /// Non-coding (`cds_start`/`cds_end == None`) sibling of
    /// `make_ungapped_indel_provider_and_projector`: same 2-bp genome-only
    /// prefix + 10-bp matched + 1-bp tx-only exon shape, but the transcript
    /// carries no CDS so it projects on the `n.` axis.
    fn make_noncoding_ungapped_indel_provider_and_projector() -> (Projector, MockProvider) {
        let common = "ATGCGCTAAC"; // 10 bases, matched region
        let genome_exon = format!("CA{common}"); // 12 bp: 2 genome-only + 10 matched
        let tx_exon = format!("{common}T"); // 11 bp: 10 matched + 1 tx-only

        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NR_GAP644.1".to_string(),
            CdotTranscript {
                gene_name: Some("GAPGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1012, 0, 11]],
                // Non-coding: no CDS bounds.
                cds_start: None,
                cds_end: None,
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NR_GAP644.1".to_string(),
            gene_symbol: Some("GAPGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some(tx_exon.clone()),
            cds_start: None,
            cds_end: None,
            exons: vec![Exon::new(1, 1, 11)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1011),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(50);
        provider.add_genomic_sequence("chr1", format!("{prefix}{genome_exon}{suffix}"));
        (projector, provider)
    }

    #[test]
    fn project_to_genomic_noncoding_corrects_across_unmodelled_indel() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::TxInterval;
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{Accession, TxVariant};

        // The outbound n.→g. path must apply the same #644 sequence-aware
        // correction the inbound g.→n. path (`correct_cds_for_exon_indel`, which
        // is coding-agnostic) applies, so the two round-trip. Without the
        // correction in the non-coding `else` branch of `map_pos`, n.5 would map
        // to the naive g.1004; with it, the 2-bp genome-only prefix shifts it to
        // g.1006 (tx offset 4 → tx-oriented genome offset 6 → HGVS g.1006).
        let (projector, provider) = make_noncoding_ungapped_indel_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        let n = TxVariant {
            accession: parse_accession("NR_GAP644.1").with_genomic_context(Accession::new(
                "NC",
                "000001",
                Some(11),
            )),
            gene_symbol: Some("GAPGENE".to_string()),
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(5)),
                NaEdit::Substitution {
                    reference: Base::G,
                    alternative: Base::T,
                },
            ),
        };
        let g = vp
            .project_to_genomic(&HgvsVariant::Tx(n))
            .expect("n. → g. projection should succeed");
        assert!(
            g.to_string().contains("g.1006"),
            "expected corrected 'g.1006' (not naive 'g.1004'), got '{g}'"
        );
    }

    #[test]
    fn project_to_genomic_rejects_out_of_transcript_noncoding_on_coding_tx() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::TxInterval;
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{Accession, TxVariant};

        // An `n.` (transcript-coordinate) position on a CODING transcript must be
        // converted into the CDS frame via `cds_pos_from_tx_pos` before the
        // CDS-aware genome mapping runs (#693). When that conversion has no answer
        // — here `n.100` is far past the single 9-base exon of `NM_TEST.1`, so
        // `cds_pos_from_tx_pos` returns `None` — the position must be REJECTED with
        // an `InvalidCoordinates` error rather than falling through unconverted and
        // being silently re-interpreted as a `c.` coordinate (which would emit a
        // wrong g. position or a less-specific downstream error). The `NM_TEST.1`
        // fixture is a coding transcript (`cds_start`/`cds_end` set) on chr1.
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        let n = TxVariant {
            // `project_to_genomic` requires an explicit NC/NG parent (#327); the
            // reject under test fires before any re-anchoring against that parent.
            accession: parse_accession("NM_TEST.1").with_genomic_context(Accession::new(
                "NC",
                "000001",
                Some(11),
            )),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };
        let err = vp
            .project_to_genomic(&HgvsVariant::Tx(n))
            .expect_err("an out-of-transcript n. position on a coding transcript must error");
        // Pin the *early-reject* contract specifically: the error must be the
        // `to_cds_frame` rejection, identified by its message. The pre-#693
        // fallback (`None => p`) also eventually errored — but with a generic
        // downstream message — so asserting only `InvalidCoordinates` would not
        // distinguish the early-reject from the fall-through. The dedicated
        // phrase below is emitted only by `to_cds_frame`, so this assertion fails
        // if the fall-through is ever reintroduced.
        match err {
            FerroError::InvalidCoordinates { msg } => assert!(
                msg.contains("cannot be mapped as a CDS coordinate"),
                "expected the to_cds_frame early-reject message, got {msg:?}"
            ),
            other => panic!("expected InvalidCoordinates from to_cds_frame, got {other:?}"),
        }
    }

    /// Minus-strand sibling of `make_ungapped_indel_provider_and_projector`.
    /// The transcript-oriented exon is the same shape (2 genome-only prefix +
    /// 10 matched + 1 tx-only), but the transcript reads off the minus strand,
    /// so the genome FASTA stores the reverse complement of the
    /// transcript-oriented exon. This is the only fixture that drives the
    /// minus-strand axis-flip arithmetic in `correct_cds_for_exon_indel`
    /// (`genome_end - 1 - genome_pos.base`) and `correct_genome_for_exon_indel`
    /// (`naive - delta`).
    fn make_minus_ungapped_indel_provider_and_projector() -> (Projector, MockProvider) {
        // Transcript-oriented exon (5'→3' as the transcript reads it):
        let common = "ATGCGCTAAC"; // 10 matched bases
        let tx_oriented_genome_exon = format!("CA{common}"); // 12 bp: 2 genome-only + 10 matched
        let tx_exon = format!("{common}T"); // 11 bp: 10 matched + 1 tx-only
                                            // The genome FASTA stores the minus-strand (reference) orientation:
        let genome_fasta_exon = crate::sequence::reverse_complement(&tx_oriented_genome_exon);

        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_GAP644M.1".to_string(),
            CdotTranscript {
                gene_name: Some("GAPGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Minus,
                // genome HGVS [1000, 1012); tx [0, 11). No CIGAR (the bug shape).
                exons: vec![[1000, 1012, 0, 11]],
                cds_start: Some(0),
                cds_end: Some(11),
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_GAP644M.1".to_string(),
            gene_symbol: Some("GAPGENE".to_string()),
            strand: TxStrand::Minus,
            sequence: Some(tx_exon.clone()),
            cds_start: Some(1),
            cds_end: Some(11),
            exons: vec![Exon::new(1, 1, 11)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1011),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(50);
        provider.add_genomic_sequence("chr1", format!("{prefix}{genome_fasta_exon}{suffix}"));
        (projector, provider)
    }

    #[test]
    fn project_minus_strand_corrects_and_roundtrips_across_unmodelled_indel() {
        use crate::hgvs::variant::Accession;

        // On minus strand the transcript-5' end is the exon's genome-3' end
        // (HGVS g.1011). The first matched base is tx-oriented genome offset 2
        // → HGVS g.1009 → c.1. A SNV at HGVS g.1005 is tx-oriented genome
        // offset 6 → naive tx 6 (c.7); the 2-bp genome-only prefix corrects it
        // to tx 4 (c.5).
        let (projector, provider) = make_minus_ungapped_indel_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        // g.→c. correction (drives the minus-strand `correct_cds_for_exon_indel`
        // axis flip). Base at g.1005 is 'C' in reference orientation.
        let result = vp
            .project("chr1:g.1005C>A", "NM_GAP644M.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. present").clone();
        assert!(
            c.to_string().contains(":c.5"),
            "expected corrected ':c.5' (not naive ':c.7') in '{c}'"
        );

        // c.→g. round-trip (drives the minus-strand
        // `correct_genome_for_exon_indel` `naive - delta` branch).
        let c_with_parent = match c {
            HgvsVariant::Cds(mut cds) => {
                cds.accession =
                    cds.accession
                        .with_genomic_context(Accession::new("NC", "000001", Some(11)));
                HgvsVariant::Cds(cds)
            }
            other => panic!("expected a Cds (c.) variant, got {other:?}"),
        };
        let g = vp
            .project_to_genomic(&c_with_parent)
            .expect("c. → g. round-trip should succeed");
        assert!(
            g.to_string().contains("g.1005"),
            "expected round-trip back to 'g.1005', got '{g}'"
        );
    }

    #[test]
    fn project_inside_genome_only_gap_errors_alignment_gap() {
        // A variant on a genome-only base (the 2-bp prefix) has no transcript
        // counterpart and must surface as AlignmentGap, not a wrong coordinate.
        let (projector, provider) = make_ungapped_indel_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // g.1000 / g.1001 are the genome-only prefix bases.
        let err = vp
            .project("chr1:g.1000C>A", "NM_GAP644.1")
            .expect_err("expected AlignmentGap for a genome-only base");
        assert!(
            matches!(err, FerroError::AlignmentGap { .. }),
            "expected AlignmentGap, got {err:?}"
        );
    }

    #[test]
    fn project_identical_sequences_unaffected_by_correction() {
        // Control: the normal fixture (transcript == genome sequence) must be
        // byte-identical to its pre-#644 behavior — no spurious correction.
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003C>A", "NM_TEST.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. present").to_string();
        assert!(
            c.contains(":c.4C>A"),
            "expected uncorrected ':c.4C>A' in '{c}'"
        );
    }

    #[test]
    fn project_substitution_plus_strand_missense() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        let result = vp
            .project("chr1:g.1003C>A", "NM_TEST.1")
            .expect("projection should succeed");

        assert_eq!(result.transcript_id, "NM_TEST.1");
        assert_eq!(result.gene_symbol.as_deref(), Some("TESTGENE"));

        let c = result
            .coding
            .as_ref()
            .expect("c. should be present")
            .to_string();
        assert!(c.contains(":c.4C>A"), "expected ':c.4C>A' in '{}' ", c);

        let p = result
            .protein
            .as_ref()
            .expect("p. should be present")
            .to_string();
        assert_eq!(p, "NP_TEST.1:p.(Arg2Ser)");

        assert!(!result.is_frameshift);
        assert!(!result.is_intronic);
        assert!(!result.is_utr);
    }

    // -------------------------------------------------------------------------
    // CDS-start sanity (issue #625)
    //
    // When a transcript's cdot CDS coordinates are inconsistent with its FASTA
    // (cds_start does not land on an ATG start codon), the reference CDS reads
    // the wrong frame. ferro must DECLINE protein prediction (Ok(None)) rather
    // than emit a fabricated consequence translated from the wrong frame.
    // -------------------------------------------------------------------------

    /// Build a single-transcript projector whose CDS sequence and length are
    /// chosen by the caller. The whole transcript is one plus-strand exon at
    /// genome `[1000, 1000+len)`, with `cds_start = 1` and `cds_end = len`
    /// (1-based, full-length CDS, no UTR). The genomic sequence mirrors the
    /// transcript so a `g.(1000+k)` SNV maps to `c.(k+1)`.
    ///
    /// `seq` is the transcript/CDS sequence; its first codon is what the
    /// CDS-start sanity guard inspects.
    #[cfg(test)]
    fn make_custom_cds_provider_and_projector(seq: &str) -> (Projector, MockProvider) {
        let len_u = seq.len() as u64;
        let genome_end_excl_u = 1000_u64 + len_u; // 0-based exclusive

        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_CUSTOM.1".to_string(),
            CdotTranscript {
                gene_name: Some("CUSTOMGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                // [genome_start(0-based), genome_end(0-based excl), tx_start, tx_end]
                exons: vec![[1000, genome_end_excl_u, 0, len_u]],
                cds_start: Some(0),   // 0-based: CDS starts at tx_pos 0 (no 5'UTR)
                cds_end: Some(len_u), // 0-based exclusive
                gene_id: None,
                protein: Some("NP_CUSTOM.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_CUSTOM.1".to_string(),
            gene_symbol: Some("CUSTOMGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some(seq.to_string()),
            cds_start: Some(1), // 1-based inclusive per Transcript convention
            cds_end: Some(len_u),
            exons: vec![Exon::new(1, 1, len_u)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(genome_end_excl_u - 1),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        // cdot genome coords are HGVS-value-based (HGVS g.X == 0-based index
        // X-1), matching every other fixture here, so the exon's genome_start
        // 1000 lands the CDS at 0-based index 999 (999 N's of prefix). The
        // sequence-aware exon-indel correction (#644) reads the genome FASTA
        // through these cdot coords, so the FASTA must be placed consistently or
        // the realign would fabricate a phantom 1-bp indel against the identical
        // transcript sequence.
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{prefix}{seq}{suffix}"));
        (projector, provider)
    }

    /// (a) RED-FIRST: a transcript whose `cds_start` lands on a non-ATG codon
    /// ("AGG", not a start) is FASTA/CDS-inconsistent. A coding SNV must
    /// DECLINE protein prediction (no fabricated consequence from the wrong
    /// frame), even though c. and n. axes are still emitted.
    #[test]
    fn project_non_atg_start_declines_protein() {
        // CDS "AGGGGCGCTAA": first codon AGG (Arg) is NOT a start codon.
        // g.1003 → c.4 (4th base 'G'), a plainly coding position.
        let (projector, provider) = make_custom_cds_provider_and_projector("AGGGGCGCTAA");
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003G>A", "NM_CUSTOM.1")
            .expect("projection should still succeed at the c./n. level");

        // c. axis must still be present (the inconsistency is CDS-only).
        assert!(
            result.coding.is_some(),
            "c. axis should still be emitted for a coding SNV"
        );
        // Protein must be DECLINED: the CDS does not start with ATG, so any
        // translation would be from the wrong frame.
        assert!(
            result.protein.is_none(),
            "expected NO protein for a non-ATG-start (inconsistent) CDS, got: {:?}",
            result.protein
        );
    }

    /// (b) CONTROL: a normal ATG-start CDS predicts protein as before. This
    /// must stay green after the guard is added.
    #[test]
    fn project_atg_start_predicts_protein_control() {
        // CDS "ATGCGCTAA": Met-Arg-Ter, a valid ATG start.
        let (projector, provider) = make_custom_cds_provider_and_projector("ATGCGCTAA");
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003C>A", "NM_CUSTOM.1")
            .expect("projection should succeed");
        let p = result
            .protein
            .as_ref()
            .expect("p. should be present for a valid ATG-start CDS")
            .to_string();
        assert_eq!(p, "NP_CUSTOM.1:p.(Arg2Ser)");
    }

    /// (c) SELENOPROTEIN CONTROL: a valid ATG-start CDS that contains an
    /// in-frame internal TGA (recoded as selenocysteine in real selenoproteins
    /// such as SELENON / NM_020451.2) before the true stop. The guard MUST key
    /// on the START codon only — an internal stop is legitimate, so protein
    /// prediction must still succeed (NOT be declined).
    #[test]
    fn project_selenoprotein_internal_tga_predicts_protein() {
        // CDS "ATGAAATGAAAATAA":
        //   codon0 ATG (Met, start)  codon1 AAA (Lys)  codon2 TGA (internal,
        //   Sec in vivo / Ter for the naive translator)  codon3 AAA (Lys)
        //   codon4 TAA (true Ter).
        // A naive translator stops at the internal TGA, but the START is a
        // valid ATG, so the CDS-start guard must NOT decline.
        let (projector, provider) = make_custom_cds_provider_and_projector("ATGAAATGAAAATAA");
        let vp = VariantProjector::new(projector, provider);
        // g.1004 → c.5: 2nd base of codon1 (AAA, Lys) — squarely coding and
        // upstream of the internal TGA so prediction is well-defined.
        let result = vp
            .project("chr1:g.1004A>G", "NM_CUSTOM.1")
            .expect("projection should succeed");
        // Assert the concrete consequence, not just presence: codon1 AAA (Lys2)
        // → AGA (Arg2). A regression that declined selenoproteins — or emitted a
        // wrong-but-present protein — would slip past a bare is_some() check.
        let p = result
            .protein
            .as_ref()
            .expect(
                "selenoprotein-style internal TGA must NOT decline protein \
                 (guard keys on START codon, not internal stops)",
            )
            .to_string();
        assert_eq!(p, "NP_CUSTOM.1:p.(Lys2Arg)");
    }

    #[test]
    fn project_genome_pivot_coding_populates_noncoding_axis() {
        // c↔n axis: a g. SNV projected onto a coding transcript carries BOTH
        // the c. form and the n. form, in distinct coordinate frames.
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003C>A", "NM_TEST.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. present");
        assert!(
            matches!(c, HgvsVariant::Cds(_)),
            "coding should be a c. (Cds) variant, got {c:?}"
        );
        let n = result
            .noncoding
            .as_ref()
            .expect("n. axis populated for coding tx");
        assert!(
            matches!(n, HgvsVariant::Tx(_)),
            "noncoding should be an n. (Tx) variant, got {n:?}"
        );
        // NM_TEST.1's CDS begins at transcript position 1, so c.4 → n.4.
        assert!(
            n.to_string().contains(":n.4C>A"),
            "expected n.4C>A, got {n}"
        );
    }

    #[test]
    fn project_single_base_deletion_is_frameshift() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1004del", "NM_TEST.1")
            .expect("deletion should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains(":c."), "expected c. variant, got: {}", c);
        assert!(c.contains("del"), "expected del notation in c., got: {}", c);
        let p = result
            .protein
            .as_ref()
            .expect("p. expected for CDS del")
            .to_string();
        assert!(p.contains("Arg2"), "expected Arg2 in p.: {}", p);
        assert!(p.contains("fs"), "expected fs in p.: {}", p);
        assert!(result.is_frameshift, "1-base del should be frameshift");
        assert!(!result.is_intronic);
    }

    #[test]
    fn project_three_base_deletion_in_frame() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1003_1005del", "NM_TEST.1")
            .expect("3-base del should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains("del"), "expected del notation, got: {}", c);
        let p = result
            .protein
            .as_ref()
            .expect("p. expected for in-frame del")
            .to_string();
        assert!(p.contains("Arg2del"), "expected Arg2del in p.: {}", p);
        assert!(!result.is_frameshift, "3-base deletion is in-frame");
    }

    #[test]
    fn project_single_base_insertion_is_frameshift() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1003_1004insA", "NM_TEST.1")
            .expect("insertion should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains("ins"), "expected ins notation, got: {}", c);
        assert!(result.protein.is_some(), "p. expected for CDS insertion");
        assert!(result.is_frameshift, "1-base insertion is frameshift");
    }

    #[test]
    fn project_duplication() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1004dup", "NM_TEST.1")
            .expect("dup should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains("dup"), "expected dup notation, got: {}", c);
        assert!(result.protein.is_some(), "p. expected for CDS dup");
        assert!(result.is_frameshift, "1-base dup is frameshift");
    }

    #[test]
    fn project_delins() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1003_1004delinsAT", "NM_TEST.1")
            .expect("delins should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains("delins"), "expected delins notation, got: {}", c);
        let p = result
            .protein
            .as_ref()
            .expect("p. expected for delins")
            .to_string();
        // c.4_5delinsAT changes codon 2 CGC→ATC, i.e. Arg2→Ile (1 AA → 1 AA).
        // Per delins.md that is a SUBSTITUTION on the protein level, not a delins,
        // even though the DNA-level edit stays a delins (#855).
        assert_eq!(p, "NP_TEST.1:p.(Arg2Ile)");
        assert!(
            !result.is_frameshift,
            "delins of equal length is not a frameshift"
        );
    }

    #[test]
    fn project_inversion() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("NC_000001.11:g.1003_1005inv", "NM_TEST.1")
            .expect("inv should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains("inv"), "expected inv notation, got: {}", c);
        let p = result
            .protein
            .as_ref()
            .expect("p. expected for inv")
            .to_string();
        // c.4_6inv inverts codon 2 CGC→GCG, i.e. Arg2→Ala (1 AA → 1 AA). A
        // single-residue whole-codon inversion is a SUBSTITUTION on the protein
        // level per delins.md (1↔1), not a delins (#855).
        assert_eq!(p, "NP_TEST.1:p.(Arg2Ala)");
        assert!(!result.is_frameshift, "inversion is not a frameshift");
    }

    fn make_ensembl_provider_and_projector() -> (Projector, MockProvider) {
        // Same CDS as make_test_provider_and_projector but the transcript ID is
        // an Ensembl accession (no NM_/XM_ prefix) and cdot.protein is absent.
        // The projector must not fabricate a bogus accession via substring
        // substitution; per #310 it falls back to using the transcript id as
        // the `p.` accession prefix.
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "ENST00000000001.1".to_string(),
            CdotTranscript {
                gene_name: Some("TESTGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "ENST00000000001.1".to_string(),
            gene_symbol: Some("TESTGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("ATGCGCTAA".to_string()),
            cds_start: Some(1),
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1008),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ATGCGCTAA", suffix));
        (projector, provider)
    }

    #[test]
    fn project_substitution_falls_back_to_transcript_id_as_protein_accession() {
        // Per #310: when no cdot.protein and no NM_/XM_ prefix is available
        // to infer NP_/XP_ from, the projector falls back to using the
        // transcript id directly as the `p.` accession prefix rather than
        // silently dropping the protein prediction. The HGVS protein
        // grammar requires a sequence_identifier, so the alternative would
        // be a non-conformant accession-less `p.` form.
        let (projector, provider) = make_ensembl_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003C>A", "ENST00000000001.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains(":c.4C>A"), "expected c.4C>A, got: {}", c);
        let p = result
            .protein
            .as_ref()
            .expect("p. should be present via transcript-id fallback")
            .to_string();
        assert_eq!(p, "ENST00000000001.1:p.(Arg2Ser)");
    }

    /// Build a coding transcript on chr1 (plus strand, "ATGCGCTAA" = Met-Arg-Stop)
    /// keyed by an arbitrary accession with **no** authoritative protein
    /// accession (cdot `protein` and `Transcript.protein_id` both `None`), so
    /// the projector's protein-accession resolution must fall through to the
    /// transcript-id fallback rather than fabricate one. Mirrors
    /// `make_test_provider_and_projector` otherwise.
    #[cfg(test)]
    fn make_protein_less_provider_and_projector(tx_id: &str) -> (Projector, MockProvider) {
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            tx_id.to_string(),
            CdotTranscript {
                gene_name: Some("TESTGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1009, 0, 9]],
                cds_start: Some(0),
                cds_end: Some(9),
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: tx_id.to_string(),
            gene_symbol: Some("TESTGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("ATGCGCTAA".to_string()),
            cds_start: Some(1),
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1008),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{prefix}ATGCGCTAA{suffix}"));
        (projector, provider)
    }

    #[test]
    fn project_substitution_nm_without_protein_uses_transcript_id_not_inferred_np() {
        // #808: a coding `NM_` transcript with no authoritative protein
        // accession must NOT have one fabricated by number-preserving
        // substitution (`NM_000077`→`NP_000077`), which is frequently wrong
        // (RefSeq does not guarantee the NM and NP numbers match). Instead the
        // projector falls back to the transcript id itself (#310) — honest, and
        // still grammar-valid. The bug here would emit `NP_NOPROT.1:...`.
        let (projector, provider) = make_protein_less_provider_and_projector("NM_NOPROT.1");
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003C>A", "NM_NOPROT.1")
            .expect("projection should succeed");
        let p = result
            .protein
            .as_ref()
            .expect("p. should be present via transcript-id fallback")
            .to_string();
        assert_eq!(
            p, "NM_NOPROT.1:p.(Arg2Ser)",
            "must use the transcript id, never a fabricated NP_"
        );
        assert!(
            !p.starts_with("NP_"),
            "must not fabricate a number-preserving NP_ accession: {p}"
        );
    }

    #[test]
    fn project_substitution_xm_without_protein_uses_transcript_id_not_inferred_xp() {
        // #808, XM_/XP_ arm: same rule for predicted-model RefSeq transcripts.
        let (projector, provider) = make_protein_less_provider_and_projector("XM_NOPROT.1");
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003C>A", "XM_NOPROT.1")
            .expect("projection should succeed");
        let p = result
            .protein
            .as_ref()
            .expect("p. should be present via transcript-id fallback")
            .to_string();
        assert_eq!(p, "XM_NOPROT.1:p.(Arg2Ser)");
        assert!(
            !p.starts_with("XP_"),
            "must not fabricate a number-preserving XP_ accession: {p}"
        );
    }

    // -------------------------------------------------------------------------
    // project_normalized tests
    // -------------------------------------------------------------------------

    #[test]
    fn project_normalized_same_result_as_project_variant() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        let variant = crate::parse_hgvs("chr1:g.1003C>A").expect("parse should succeed");
        let via_project = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("project_variant should succeed");
        // Pre-normalize via the cached normalizer, then call project_normalized.
        let normalized = vp.normalizer.normalize(&variant).expect("normalize failed");
        let via_normalized = vp
            .project_normalized(&normalized, "NM_TEST.1")
            .expect("project_normalized should succeed");

        assert_eq!(
            via_project.coding.as_ref().map(|v| v.to_string()),
            via_normalized.coding.as_ref().map(|v| v.to_string()),
            "project_normalized should produce same c. as project_variant"
        );
        assert_eq!(
            via_project.protein.as_ref().map(|v| v.to_string()),
            via_normalized.protein.as_ref().map(|v| v.to_string()),
            "project_normalized should produce same p. as project_variant"
        );
    }

    // -------------------------------------------------------------------------
    // project_all / project_variant_all tests
    // -------------------------------------------------------------------------

    #[test]
    fn project_variant_all_returns_both_transcripts() {
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        let results = vp
            .project_all("chr1:g.1003C>A")
            .expect("project_all should succeed");

        assert_eq!(
            results.len(),
            2,
            "expected projections onto both transcripts, got: {}",
            results.len()
        );
    }

    #[test]
    fn project_all_mane_select_sorts_first() {
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        let results = vp
            .project_all("chr1:g.1003C>A")
            .expect("project_all should succeed");

        // NM_TX2.1 is MANE Select → must be first.
        assert_eq!(
            results[0].transcript_id, "NM_TX2.1",
            "MANE Select should be first"
        );
    }

    #[test]
    fn project_all_applies_curated_enumeration_policy() {
        // End-to-end: the #656 curation must be wired into the fan-out loop, not
        // just unit-tested on `select_enumerated_transcript_ids`. cdot reports
        // three overlapping records (NM_TX1.1, NM_TX1.2, XM_TX9.1) but the
        // curated set is exactly {NM_TX1.2}: the superseded .1 is collapsed and
        // the predicted XM_ model is dropped in favor of the curated transcript.
        let (projector, provider) = make_curated_enumeration_setup();
        let vp = VariantProjector::new(projector, provider);

        let results = vp
            .project_all("chr1:g.1003C>A")
            .expect("project_all should succeed");
        let ids: Vec<&str> = results.iter().map(|p| p.transcript_id.as_str()).collect();
        assert_eq!(
            ids,
            vec!["NM_TX1.2"],
            "expected only the curated, current-version transcript, got {:?}",
            ids
        );
    }

    #[test]
    fn project_normalized_all_applies_curated_enumeration_policy() {
        // Same curation, driven through the pre-normalized fan-out entrypoint.
        let (projector, provider) = make_curated_enumeration_setup();
        let vp = VariantProjector::new(projector, provider);

        let variant = crate::parse_hgvs("chr1:g.1003C>A").expect("parse should succeed");
        let normalized = vp.normalizer.normalize(&variant).expect("normalize failed");
        let results = vp
            .project_normalized_all(&normalized)
            .expect("project_normalized_all should succeed");
        let ids: Vec<&str> = results.iter().map(|p| p.transcript_id.as_str()).collect();
        assert_eq!(
            ids,
            vec!["NM_TX1.2"],
            "expected only the curated, current-version transcript, got {:?}",
            ids
        );
    }

    // -------------------------------------------------------------------------
    // c./n./r. fan-out via project_variant_all (closes #389 item 3)
    //
    // Pre-fix, `extract_contig_and_pos` matched only `HgvsVariant::Genome`
    // and returned `UnsupportedProjection { reason: "project_all
    // currently only accepts g. variants" }` for every other axis, which
    // defeated PR #379's widening of `project()` to c./n./r. inputs in
    // the fan-out entrypoints.
    // -------------------------------------------------------------------------

    /// `project_variant_all` on a c. input with an NG/NC parent must
    /// reach the fan-out path and produce a per-transcript projection
    /// for the matching transcript. Pre-#389 this errored out at the
    /// `extract_contig_and_pos` g.-only gate.
    #[test]
    fn project_variant_all_accepts_coding_input_with_parent() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{Accession, CdsVariant, HgvsVariant, LocEdit};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        // NM_TX1.1 c.4C>A with an NC_000001.11 parent; per the
        // make_two_transcript_setup fixture this should project to
        // g.1003 (which is c.4 on the plus-strand 9-base CDS).
        let cds = CdsVariant {
            accession: parse_accession("NM_TX1.1").with_genomic_context(Accession::new(
                "NC",
                "000001",
                Some(11),
            )),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let results = vp
            .project_variant_all(&HgvsVariant::Cds(cds))
            .expect("project_variant_all should accept c. input post-#389");
        assert!(
            results.iter().any(|p| p.transcript_id == "NM_TX1.1"),
            "expected at least the input transcript's projection, got {:?}",
            results.iter().map(|p| &p.transcript_id).collect::<Vec<_>>()
        );
    }

    /// Same as above but for the n. axis: `project_variant_all` on a
    /// non-coding (or n.-axis) variant with an NG/NC parent must also
    /// reach the fan-out path.
    #[test]
    fn project_variant_all_accepts_noncoding_n_input_with_parent() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::TxInterval;
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{Accession, HgvsVariant, LocEdit, TxVariant};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        // n.4 maps to genome 1003 on NM_TX1.1 in this fixture (tx_pos 3
        // → exon.genome_start + 3 = 1003).
        let tx = TxVariant {
            accession: parse_accession("NM_TX1.1").with_genomic_context(Accession::new(
                "NC",
                "000001",
                Some(11),
            )),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let results = vp
            .project_variant_all(&HgvsVariant::Tx(tx))
            .expect("project_variant_all should accept n. input post-#389");
        assert!(
            results.iter().any(|p| p.transcript_id == "NM_TX1.1"),
            "expected at least the input transcript's projection, got {:?}",
            results.iter().map(|p| &p.transcript_id).collect::<Vec<_>>()
        );
    }

    /// `project_variant_all` on an r. (RNA) input — parallel to n. but
    /// on the r. axis. Pre-#389 the gate rejected this too.
    #[test]
    fn project_variant_all_accepts_rna_input_with_parent() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::location::RnaPos;
        use crate::hgvs::variant::{Accession, HgvsVariant, LocEdit, RnaVariant};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        let rna = RnaVariant {
            accession: parse_accession("NM_TX1.1").with_genomic_context(Accession::new(
                "NC",
                "000001",
                Some(11),
            )),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let results = vp
            .project_variant_all(&HgvsVariant::Rna(rna))
            .expect("project_variant_all should accept r. input post-#389");
        assert!(
            results.iter().any(|p| p.transcript_id == "NM_TX1.1"),
            "expected at least the input transcript's projection, got {:?}",
            results.iter().map(|p| &p.transcript_id).collect::<Vec<_>>()
        );
    }

    /// c. input *without* a `genomic_context` parent: the fan-out
    /// rejects the input up front with `UnsupportedProjection`
    /// (#389 follow-up).
    ///
    /// Pre-fix the stab query still seeded and each per-transcript
    /// projection silently failed inside `project_to_genomic`'s parent
    /// requirement; the fan-out loop's `log::trace!`-and-drop converted
    /// the user error into `Ok([])`, turning "you forgot the parent"
    /// into "no overlaps." Post-fix the error is surfaced.
    #[test]
    fn project_variant_all_rejects_coding_input_without_parent() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{CdsVariant, HgvsVariant, LocEdit};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        // No `.with_genomic_context(...)` — bare NM accession on a c. input.
        let cds = CdsVariant {
            accession: parse_accession("NM_TX1.1"),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let err = vp
            .project_variant_all(&HgvsVariant::Cds(cds))
            .expect_err("c. input without parent must error, not silently return empty");
        match err {
            FerroError::UnsupportedProjection { reason } => {
                assert!(
                    reason.contains("project_all requires a parent reference")
                        && reason.contains("c."),
                    "unexpected error message: {reason}"
                );
            }
            other => panic!("expected UnsupportedProjection, got {other:?}"),
        }
    }

    /// Same as above but for n. — the Tx-axis arm of the
    /// `require_parent_for_fanout` gate (#389 follow-up).
    #[test]
    fn project_variant_all_rejects_noncoding_n_input_without_parent() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::TxInterval;
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{HgvsVariant, LocEdit, TxVariant};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        let tx = TxVariant {
            accession: parse_accession("NM_TX1.1"),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let err = vp
            .project_variant_all(&HgvsVariant::Tx(tx))
            .expect_err("n. input without parent must error");
        match err {
            FerroError::UnsupportedProjection { reason } => {
                assert!(reason.contains("n."), "unexpected error message: {reason}");
            }
            other => panic!("expected UnsupportedProjection, got {other:?}"),
        }
    }

    /// Same as above but for r. — the Rna-axis arm of the gate.
    #[test]
    fn project_variant_all_rejects_rna_input_without_parent() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::location::RnaPos;
        use crate::hgvs::variant::{HgvsVariant, LocEdit, RnaVariant};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        let rna = RnaVariant {
            accession: parse_accession("NM_TX1.1"),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let err = vp
            .project_variant_all(&HgvsVariant::Rna(rna))
            .expect_err("r. input without parent must error");
        match err {
            FerroError::UnsupportedProjection { reason } => {
                assert!(reason.contains("r."), "unexpected error message: {reason}");
            }
            other => panic!("expected UnsupportedProjection, got {other:?}"),
        }
    }

    /// `HgvsVariant::Allele` whose first inner variant is c. — the
    /// fan-out path should unwrap and process the inner correctly,
    /// matching the behavior of g.-inner alleles.
    #[test]
    fn project_variant_all_accepts_allele_of_coding_input() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{Accession, CdsVariant, HgvsVariant, LocEdit};
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        let cds = CdsVariant {
            accession: parse_accession("NM_TX1.1").with_genomic_context(Accession::new(
                "NC",
                "000001",
                Some(11),
            )),
            gene_symbol: Some("GENE1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        };
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![HgvsVariant::Cds(cds)]));
        let results = vp
            .project_variant_all(&allele)
            .expect("Allele wrapping a c. inner must reach the fan-out path");
        assert!(
            results.iter().any(|p| p.transcript_id == "NM_TX1.1"),
            "expected the input transcript's projection from the allele's inner c., got {:?}",
            results.iter().map(|p| &p.transcript_id).collect::<Vec<_>>()
        );
    }

    #[test]
    fn project_all_no_overlap_returns_empty() {
        let (projector, provider) = make_two_transcript_setup();
        let vp = VariantProjector::new(projector, provider);

        // g.5000 is far outside all transcripts.
        let results = vp
            .project_all("NC_000001.11:g.5000A>G")
            .expect("project_all should return Ok for no overlaps");

        assert!(
            results.is_empty(),
            "expected empty result for non-overlapping variant"
        );
    }

    // -------------------------------------------------------------------------
    // Allele compound projection tests
    // -------------------------------------------------------------------------

    /// Helper: build a cis allele `[chr1:g.1003C>A;chr1:g.1006T>A]` and
    /// project it onto NM_TEST.1.
    fn project_cis_allele() -> VariantProjection {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        // Use parse → project_variant to build the allele naturally.
        // g.1003 is C (c.4) and g.1006 is T (c.7) in the test fixture
        // "ATGCGCTAA"; using the correct ref bases keeps the test
        // valid HGVS rather than just exercising control flow.
        let v1 = crate::parse_hgvs("chr1:g.1003C>A").expect("v1 parse");
        let v2 = crate::parse_hgvs("chr1:g.1006T>A").expect("v2 parse");
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![v1, v2]));
        vp.project_variant(&allele, "NM_TEST.1")
            .expect("cis allele projection should succeed")
    }

    #[test]
    fn project_cis_allele_produces_cis_coding_allele() {
        let result = project_cis_allele();
        let coding = result.coding.as_ref().expect("c. allele should be present");
        // The coding variant should itself be an Allele.
        assert!(
            matches!(coding, HgvsVariant::Allele(av) if av.phase == AllelePhase::Cis),
            "expected Cis coding allele, got: {}",
            coding
        );
        // Two inner c. variants.
        if let HgvsVariant::Allele(av) = coding {
            assert_eq!(
                av.variants.len(),
                2,
                "expected 2 inner c. variants in cis allele"
            );
        }
    }

    #[test]
    fn project_cis_allele_aggregates_noncoding_axis() {
        // c↔n axis: every inner member carries an n. form, so the allele
        // projection exposes a cis n. allele mirroring the c. allele.
        let result = project_cis_allele();
        let noncoding = result
            .noncoding
            .as_ref()
            .expect("n. allele should be aggregated");
        assert!(
            matches!(noncoding, HgvsVariant::Allele(av) if av.phase == AllelePhase::Cis),
            "expected Cis noncoding allele, got: {noncoding}"
        );
        if let HgvsVariant::Allele(av) = noncoding {
            assert_eq!(
                av.variants.len(),
                2,
                "expected 2 inner n. variants in cis allele"
            );
            assert!(
                av.variants.iter().all(|v| matches!(v, HgvsVariant::Tx(_))),
                "every inner member of the n. allele should be an n. (Tx) variant"
            );
            // The fixture's inner members are g.1003→c.4 and g.1006→c.7 on
            // NM_TEST.1 (CDS at tx pos 1), so the aggregated n. coordinates are
            // n.4 and n.7 — pin them so an off-by-one in the derivation fails.
            let rendered: Vec<String> = av.variants.iter().map(|v| v.to_string()).collect();
            assert!(
                rendered.iter().any(|s| s.contains(":n.4")),
                "expected an n.4 member, got {rendered:?}"
            );
            assert!(
                rendered.iter().any(|s| s.contains(":n.7")),
                "expected an n.7 member, got {rendered:?}"
            );
        }
    }

    #[test]
    fn project_trans_allele_preserves_trans_phase() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        // Ref bases match the fixture: g.1003 = C, g.1006 = T.
        let v1 = crate::parse_hgvs("chr1:g.1003C>A").expect("v1 parse");
        let v2 = crate::parse_hgvs("chr1:g.1006T>A").expect("v2 parse");
        let allele = HgvsVariant::Allele(AlleleVariant::trans(vec![v1, v2]));
        let result = vp
            .project_variant(&allele, "NM_TEST.1")
            .expect("trans allele projection should succeed");

        let coding = result.coding.as_ref().expect("c. expected");
        assert!(
            matches!(coding, HgvsVariant::Allele(av) if av.phase == AllelePhase::Trans),
            "expected Trans coding allele, got: {}",
            coding
        );
    }

    #[test]
    fn project_allele_with_frameshift_inner_sets_is_frameshift() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        // g.1004del is a 1-base deletion → frameshift.
        let v1 = crate::parse_hgvs("chr1:g.1003C>A").expect("v1 parse");
        let v_fs = crate::parse_hgvs("NC_000001.11:g.1004del").expect("fs parse");
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![v1, v_fs]));
        let result = vp
            .project_variant(&allele, "NM_TEST.1")
            .expect("allele with frameshift should project");

        assert!(
            result.is_frameshift,
            "allele containing a frameshift inner variant should set is_frameshift=true"
        );
    }

    #[test]
    fn project_allele_with_non_protein_inner_has_no_protein() {
        // An allele where one inner variant has no protein (intronic) →
        // the whole allele protein should be None.
        let (projector, provider) = make_intronic_test_data();
        let vp = VariantProjector::new(projector, provider);

        // g.1003 (exonic, ref C in NM_INTR.1 fixture) + g.1015 (intronic, ref N
        // since g.1015 is in the intron gap — the projector doesn't validate
        // intronic ref bases against the genomic reference).
        let v_exon = crate::parse_hgvs("NC_000001.11:g.1003C>G").expect("exon parse");
        let v_intron = crate::parse_hgvs("NC_000001.11:g.1015N>G").expect("intron parse");
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![v_exon, v_intron]));
        let result = vp
            .project_variant(&allele, "NM_INTR.1")
            .expect("allele with intronic variant should project");

        assert!(
            result.protein.is_none(),
            "allele with intronic inner variant should have no protein"
        );
        assert!(result.is_intronic, "should be marked intronic");
    }

    // -------------------------------------------------------------------------
    // project_to_genomic tests (#327)
    // -------------------------------------------------------------------------

    mod project_to_genomic_tests {
        use super::*;
        use crate::hgvs::edit::{Base, InsertedSequence, NaEdit};
        use crate::hgvs::interval::{CdsInterval, GenomeInterval, TxInterval};
        use crate::hgvs::location::{CdsPos, GenomePos, TxPos};
        use crate::hgvs::variant::{
            Accession, CdsVariant, GenomeVariant, HgvsVariant, LocEdit, TxVariant,
        };

        /// Build an `NG_*` parent `Accession`.
        fn ng_parent(num: &str, version: u32) -> Accession {
            Accession::new("NG", num, Some(version))
        }

        /// An `NC_` chromosome parent. Unlike an `NG_`/`LRG_` parent it shares
        /// cdot's coordinate frame, so `project_to_genomic` emits it unchanged
        /// (no re-anchoring, no decline) — the right vehicle for exercising the
        /// c./n./r. → g. coordinate math without tripping the NG_/LRG_ decline
        /// (#655). GRCh38 primary (`.11`) to match the test cdot fixture.
        fn nc_parent() -> Accession {
            Accession::new("NC", "000001", Some(11))
        }

        /// Attach a genomic_context to an existing `CdsVariant`.
        fn attach_genomic_context_cds(mut v: CdsVariant, ctx: Accession) -> CdsVariant {
            v.accession = v.accession.with_genomic_context(ctx);
            v
        }

        /// Attach a genomic_context to an existing `TxVariant` (n.).
        fn attach_genomic_context_tx(mut v: TxVariant, ctx: Accession) -> TxVariant {
            v.accession = v.accession.with_genomic_context(ctx);
            v
        }

        /// Build a non-coding transcript test fixture for `n.` projection.
        ///
        /// `NR_TEST.1` is a 9-base non-coding RNA on chr1, plus strand,
        /// genome [1000, 1009).  Without a CDS the transcript uses pure n.
        /// coordinates (1-based: n.1..n.9). n.5 maps to genome 1004.
        fn make_noncoding_test_data() -> (Projector, MockProvider) {
            let mut cdot = CdotMapper::new();
            cdot.add_transcript(
                "NR_TEST.1".to_string(),
                CdotTranscript {
                    gene_name: Some("TESTGENE".to_string()),
                    contig: "chr1".to_string(),
                    strand: ProvStrand::Plus,
                    exons: vec![[1000, 1009, 0, 9]],
                    // No CDS for n. transcripts.
                    cds_start: None,
                    cds_end: None,
                    gene_id: None,
                    protein: None,
                    exon_cigars: Vec::new(),
                },
            );
            let projector = Projector::new(cdot);

            let mut provider = MockProvider::new();
            provider.add_transcript(Transcript {
                id: "NR_TEST.1".to_string(),
                gene_symbol: Some("TESTGENE".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("ATGCGCTAA".to_string()),
                cds_start: None,
                cds_end: None,
                exons: vec![Exon::with_genomic(1, 1, 9, 1000, 1008)],
                chromosome: Some("chr1".to_string()),
                genomic_start: Some(1000),
                genomic_end: Some(1008),
                genome_build: Default::default(),
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
                protein_id: None,
            });
            let prefix = "N".repeat(999);
            let suffix = "N".repeat(100);
            provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ATGCGCTAA", suffix));
            (projector, provider)
        }

        // Note: the non-coding (`else`) branch of `project_single_inner` — where
        // `noncoding` mirrors the `coding` `Tx` form — is currently unreachable
        // via full projection: the downstream g.→c. path errors with
        // `TranscriptNotOverlapping` for non-coding transcripts, pinned in
        // `tests/issue_328_projector_accepts_cnr.rs`. When that path is expanded
        // to support non-coding tx, add a `noncoding == coding` regression test
        // here. The branch stays for forward-compatibility (correct semantics
        // the moment the path is reachable).

        // -- 1: plus-strand substitution ---------------------------------------

        #[test]
        fn project_to_genomic_plus_strand_substitution() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // Build NG_TEST.1(NM_TEST.1):c.4C>A.
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("plus-strand c. → g. projection should succeed");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant, got: {}", out),
            };
            // Parent NC_000001.11, plus-strand → g.1003C>A.
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let s = g.to_string();
            assert!(
                s.starts_with("NC_000001.11") && s.contains(":g.1003C>A"),
                "expected NC_000001.11...:g.1003C>A, got: {}",
                s
            );
        }

        // -- 1b: versioned silent-substitution refusal in the c.→g. direction --

        #[test]
        fn project_to_genomic_declines_silent_version_substitution() {
            // Issue #785, c.→g. direction. `project_to_genomic` does not normalize
            // the original transcript-coordinate input first, so without an
            // explicit gate it would resolve a versioned request through the
            // *lenient* cdot pivot and silently project onto the genome using a
            // *sibling* version's exon/CDS frame. Here the reference carries only
            // `NM_TEST.1`; a request for the absent `NM_TEST.2` must DECLINE with
            // `TranscriptVersionNotExact` rather than silently use `.1`'s frame —
            // the same guarantee `project_variant`/`project_variant_all` inherit
            // by normalizing first.
            let (projector, mut provider) = make_test_provider_and_projector();
            // NM_TEST.2 is absent; a lenient lookup serves the .1 bases/frame.
            provider.mark_version_substitution("NM_TEST.2", "NM_TEST.1");
            let vp = VariantProjector::new(projector, provider);

            // NC_000001.11(NM_TEST.2):c.4C>A — explicitly versioned, sibling-only.
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.2"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);
            let err = vp
                .project_to_genomic(&input)
                .expect_err("a versioned silent substitution must decline in c.→g.");
            assert!(
                matches!(err, FerroError::TranscriptVersionNotExact { ref requested } if requested == "NM_TEST.2"),
                "expected TranscriptVersionNotExact for NM_TEST.2, got {err:?}"
            );
        }

        // -- 1c: exact-version request still projects (no over-decline) ---------

        #[test]
        fn project_to_genomic_allows_exact_version() {
            // Guard against over-declining: when the EXACT requested version is
            // served, the #785 gate must NOT fire. `NM_TEST.1` is present, so
            // `NC_000001.11(NM_TEST.1):c.4C>A` projects to g. unchanged.
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("exact-version c.→g. projection must not be declined");
            assert!(
                out.to_string().contains(":g.1003C>A"),
                "expected NC_000001.11...:g.1003C>A, got: {out}"
            );
        }

        // -- 2: minus-strand substitution --------------------------------------

        #[test]
        fn project_to_genomic_minus_strand_substitution_revcomps() {
            let (projector, provider) = make_minus_strand_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // c.4C>T on minus-strand NM_TEST_MINUS.1 → tx_pos 4 on minus,
            // c.1 is the last base of the exon [1000, 1009) at genome 1008,
            // so c.4 = genome 1005.  Ref C on transcript = G on genome;
            // alt T on transcript = A on genome → g.1005G>A.
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST_MINUS.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::T,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("minus-strand c. → g. should succeed");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let s = g.to_string();
            assert!(
                s.starts_with("NC_000001.11") && s.contains(":g.1005G>A"),
                "minus-strand c.4C>T should revcomp to NC_000001.11...:g.1005G>A, got: {}",
                s
            );
        }

        // -- 3: minus-strand insertion (revcomp inserted seq) ------------------

        #[test]
        fn project_to_genomic_minus_strand_insertion_revcomps_inserted_seq() {
            let (projector, provider) = make_minus_strand_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // c.4_5insATC on minus strand.  c.4 = g.1005, c.5 = g.1004 (minus).
            // Inserted "ATC" on the transcript → revcomp "GAT" on the genome.
            // g. start/end are min/max of the two genome coords, so 1004_1005.
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST_MINUS.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::new(CdsPos::new(4), CdsPos::new(5)),
                    NaEdit::Insertion {
                        sequence: InsertedSequence::Literal("ATC".parse().unwrap()),
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("minus-strand insertion should project");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            let s = g.to_string();
            // Pin the exact interval — a start/end swap regression would
            // otherwise still pass via `contains("ins")` + `contains("GAT")`.
            // c.4 = g.1005, c.5 = g.1004; canonical g. interval is low→high.
            assert!(
                s.contains(":g.1004_1005ins"),
                "expected interval :g.1004_1005ins, got: {}",
                s
            );
            assert!(
                s.contains("GAT"),
                "expected revcomp(ATC) = GAT in g., got: {}",
                s
            );
        }

        // -- 3b: minus-strand deletion is re-normalized in genomic space (#737) -

        /// A minus-strand transcript on `NC_000001.11` whose CDS reads
        /// `GTAAAAAAA` (cdna 1..9), so the genome (its reverse complement) carries
        /// a `T` homopolymer. cdot exon `[1000, 1009)` minus-strand maps c.N →
        /// genome `1009 - N`, so the transcript's poly-A run (cdna 3..9) is the
        /// genome poly-T run at g.1000..1006. A single-base deletion anywhere in
        /// that run is the *same* variant; its spec-canonical genomic form is the
        /// 3'-most (highest-coordinate) position, g.1006del.
        ///
        /// Keyed under `NC_000001.11` (not `chr1`) so `normalize` can fetch the
        /// genomic bases back through the same accession the re-anchored output
        /// carries — exercising the #737 genomic-frame renormalization end to end.
        fn make_minus_homopolymer_provider_and_projector() -> (Projector, MockProvider) {
            let mut cdot = CdotMapper::new();
            cdot.add_transcript(
                "NM_HOMO_MINUS.1".to_string(),
                CdotTranscript {
                    gene_name: Some("HOMOGENE".to_string()),
                    contig: "NC_000001.11".to_string(),
                    strand: ProvStrand::Minus,
                    exons: vec![[1000, 1009, 0, 9]],
                    cds_start: Some(0),
                    cds_end: Some(9),
                    gene_id: None,
                    protein: Some("NP_HOMO_MINUS.1".to_string()),
                    exon_cigars: Vec::new(),
                },
            );
            let projector = Projector::new(cdot);

            let mut provider = MockProvider::new();
            provider.add_transcript(Transcript {
                id: "NM_HOMO_MINUS.1".to_string(),
                gene_symbol: Some("HOMOGENE".to_string()),
                strand: TxStrand::Minus,
                sequence: Some("GTAAAAAAA".to_string()),
                cds_start: Some(1),
                cds_end: Some(9),
                exons: vec![Exon::new(1, 1, 9)],
                chromosome: Some("NC_000001.11".to_string()),
                genomic_start: Some(1000),
                genomic_end: Some(1008),
                genome_build: Default::default(),
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                protein_id: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
            });
            // Genome forward: prefix + "TTTTTTTAC" (g.1000..1008) + suffix, so the
            // poly-T run is g.1000..1006 (revcomp of the transcript poly-A run).
            let prefix = "N".repeat(999);
            let suffix = "N".repeat(100);
            provider.add_genomic_sequence("NC_000001.11", format!("{}TTTTTTTAC{}", prefix, suffix));
            (projector, provider)
        }

        #[test]
        fn project_to_genomic_minus_strand_deletion_renormalizes_to_genomic_3prime() {
            let (projector, provider) = make_minus_homopolymer_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // c.9del is the transcript-3'-most representation of the deletion in
            // the poly-A run. On the minus strand that maps to g.1000 — the
            // genome-5' end of the poly-T run. The spec-canonical genomic output
            // is the genome-3'-most position, g.1006del. Before #737 the projector
            // emitted the un-renormalized coordinate image (g.1000del); the
            // re-anchored output is now normalized in its own frame.
            let cds = CdsVariant {
                accession: parse_accession("NM_HOMO_MINUS.1"),
                gene_symbol: Some("HOMOGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(9)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("minus-strand c.9del should project to g.");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant, got: {}", out),
            };
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let start = g
                .loc_edit
                .location
                .start
                .inner()
                .expect("start should be concrete");
            assert_eq!(
                start.base, 1006,
                "minus-strand c.9del must re-normalize to the genome-3' end of the \
                 poly-T run (g.1006del), got g.{}del (#737)",
                start.base
            );
        }

        /// Companion to the g.-only test above, exercising the *other* #737
        /// re-anchor site: the multi-axis `.genomic` field built in
        /// `project_single_inner` (via `project_variant`). The two sites apply
        /// the identical genomic-frame renormalization and could drift silently,
        /// so assert the same minus-strand deletion lands at the genome-3' end of
        /// the poly-T run (g.1006del) through a `VariantProjection.genomic`, not
        /// just through the g.-only `project_to_genomic`.
        #[test]
        fn project_variant_minus_strand_deletion_renormalizes_genomic_axis_to_3prime() {
            let (projector, provider) = make_minus_homopolymer_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // c.9del — the transcript-3'-most representation of the deletion in
            // the poly-A run — maps to the genome-5' end of the poly-T run. The
            // spec-canonical genomic output is the genome-3'-most position,
            // g.1006del; #737 renormalizes the re-anchored `.genomic` axis here.
            let cds = CdsVariant {
                accession: parse_accession("NM_HOMO_MINUS.1"),
                gene_symbol: Some("HOMOGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(9)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let proj = vp
                .project_variant(&HgvsVariant::Cds(cds), "NM_HOMO_MINUS.1")
                .expect("minus-strand c.9del should project (multi-axis)");
            let g = match proj.genomic {
                Some(HgvsVariant::Genome(ref g)) => g,
                other => panic!("expected a Genome `.genomic` axis, got: {:?}", other),
            };
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let start = g
                .loc_edit
                .location
                .start
                .inner()
                .expect("start should be concrete");
            assert_eq!(
                start.base, 1006,
                "the multi-axis `.genomic` field must re-normalize to the genome-3' \
                 end of the poly-T run (g.1006del), got g.{}del (#737)",
                start.base
            );
        }

        // -- 4: intronic offset, plus strand -----------------------------------

        #[test]
        fn project_to_genomic_intronic_offset_plus() {
            let (projector, provider) = make_intronic_test_data();
            let vp = VariantProjector::new(projector, provider);

            // NM_INTR.1 layout (plus strand):
            //   Exon 1: genome [1000, 1010), tx [0, 10)
            //   Intron: genome [1010, 2000)
            //   Exon 2: genome [2000, 2010), tx [10, 20)
            // c.10 is the last base of exon 1 (tx_pos = 9, 0-based, genome 1009).
            // c.10+5 → 5 bases into the intron from the 5' boundary → g.1014.
            let cds = CdsVariant {
                accession: parse_accession("NM_INTR.1"),
                gene_symbol: Some("INTRGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::with_offset(10, 5)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);
            let out = vp
                .project_to_genomic(&input)
                .expect("intronic c.10+5del should project");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            // Expect a deletion at g.1014 on NC_000001.11.
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let start = g
                .loc_edit
                .location
                .start
                .inner()
                .expect("start should be concrete");
            assert_eq!(
                start.base, 1014,
                "expected g.1014 for c.10+5 on plus strand, got: {}",
                start.base
            );
        }

        // -- 5: idempotence on Genome input ------------------------------------

        #[test]
        fn project_to_genomic_idempotent_on_genome_input() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // Build a plain Genome variant.
            let g = GenomeVariant {
                accession: parse_accession("NC_000001.11"),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(1003)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let input = HgvsVariant::Genome(g.clone());
            let out = vp
                .project_to_genomic(&input)
                .expect("Genome input should pass through");
            assert_eq!(
                out.to_string(),
                input.to_string(),
                "Genome input should be idempotent under project_to_genomic"
            );
        }

        /// An `NG_`-genomic input is already in its parent frame, so
        /// `project_to_genomic` must return it unchanged. Before #702 it was fed
        /// through `reanchor_genome_output`, which mistook the parent coordinate
        /// for an `NC_` coordinate and declined (`UnsupportedProjection`: endpoint
        /// outside the placed span) instead of passing the variant through.
        #[test]
        fn project_to_genomic_idempotent_on_ng_genome_input() {
            let (projector, mut provider) = make_nc_two_transcript_setup();
            // NG_900.1 placed on NC_000001.11, plus strand: NG_ base 1 == nc 1000.
            provider.add_genomic_placement(
                "NG_900.1",
                crate::reference::GenomicPlacement {
                    nc: Accession::new("NC", "000001", Some(11)),
                    parent_start: 1,
                    nc_start: 1000,
                    nc_end: 1008,
                    strand: crate::reference::Strand::Plus,
                },
            );
            let vp = VariantProjector::new(projector, provider);

            // NG_900.1:g.4C>A is in the NG_ frame; base 4 (< nc_start 1000) would
            // be rejected by `nc_to_parent` if treated as an NC_ coordinate.
            let g = GenomeVariant {
                accession: parse_accession("NG_900.1"),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let input = HgvsVariant::Genome(g);
            let out = vp
                .project_to_genomic(&input)
                .expect("NG_ genome input should pass through unchanged");
            assert_eq!(
                out.to_string(),
                input.to_string(),
                "NG_ genome input should be idempotent under project_to_genomic"
            );
        }

        // -- 6: missing parent NG/NC -------------------------------------------

        #[test]
        fn project_to_genomic_missing_parent_returns_unsupported() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // NM_TEST.1:c.4C>A — no genomic_context attached.
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let input = HgvsVariant::Cds(cds);
            let err = vp
                .project_to_genomic(&input)
                .expect_err("missing parent should error");
            match err {
                FerroError::UnsupportedProjection { reason } => {
                    assert!(
                        reason.contains("parent reference"),
                        "expected reason to mention 'parent reference', got: {}",
                        reason
                    );
                }
                other => panic!("expected UnsupportedProjection, got: {:?}", other),
            }
        }

        // -- 7: n. (non-coding) input ------------------------------------------

        #[test]
        fn project_to_genomic_tx_input_on_noncoding() {
            let (projector, provider) = make_noncoding_test_data();
            let vp = VariantProjector::new(projector, provider);

            // NR_TEST.1:n.5A>G on plus strand → n.5 = genome 1004.
            let tx = TxVariant {
                accession: parse_accession("NR_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    TxInterval::point(TxPos::new(5)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::G,
                    },
                ),
            };
            let tx = attach_genomic_context_tx(tx, nc_parent());
            let input = HgvsVariant::Tx(tx);
            let out = vp
                .project_to_genomic(&input)
                .expect("n. → g. projection should succeed");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let s = g.to_string();
            assert!(
                s.starts_with("NC_000001.11") && s.contains(":g.1004A>G"),
                "expected NC_000001.11...:g.1004A>G, got: {}",
                s
            );
        }

        // -- 8: r. input -------------------------------------------------------

        #[test]
        fn project_to_genomic_rna_input_with_dna_bases() {
            use crate::hgvs::interval::RnaInterval;
            use crate::hgvs::location::RnaPos;
            use crate::hgvs::variant::RnaVariant;

            // Verifies r. → g. coordinate projection on the plus strand with
            // DNA-typed bases. RnaPos shape mirrors CdsPos so the coordinate
            // math is identical. (Plus-strand `Base::U` inputs are translated
            // U→T on the g. output — covered separately by
            // `project_to_genomic_rna_u_base_translated_to_t`.)
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // r.4 on NM_TEST.1 (plus strand) → g.1003 on NC_000001.11.
            let rna = RnaVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    RnaInterval::point(RnaPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let mut rna = rna;
            rna.accession = rna.accession.with_genomic_context(nc_parent());
            let input = HgvsVariant::Rna(rna);
            let out = vp
                .project_to_genomic(&input)
                .expect("r. → g. projection should succeed");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            // Pin both the interval and the bases: a coordinate regression
            // would otherwise pass `contains("C>A")` alone.
            let s = g.to_string();
            assert!(
                s.contains(":g.1003C>A"),
                "expected NC_000001.11...:g.1003C>A, got: {}",
                s
            );
            assert_eq!(g.accession.to_string(), "NC_000001.11");
        }

        /// A plus-strand `r.` input carrying `Base::U` is translated to DNA
        /// (`U`→`T`) on the g. output, so the emitted g. variant is valid DNA
        /// rather than an `…U>…` shape (#395 item 4).
        #[test]
        fn project_to_genomic_rna_u_base_translated_to_t() {
            use crate::hgvs::interval::RnaInterval;
            use crate::hgvs::location::RnaPos;
            use crate::hgvs::variant::RnaVariant;

            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // r.4u>g on NM_TEST.1 (plus strand) → g.1003; the RNA `u` reference
            // base must surface as DNA `T`, giving g.1003T>G (ref not validated
            // against the genome, so the edit's own bases are emitted).
            let mut rna = RnaVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    RnaInterval::point(RnaPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::U,
                        alternative: Base::G,
                    },
                ),
            };
            rna.accession = rna.accession.with_genomic_context(nc_parent());
            let out = vp
                .project_to_genomic(&HgvsVariant::Rna(rna))
                .expect("r. → g. projection should succeed");
            let s = out.to_string();
            assert!(
                s.contains(":g.1003T>G"),
                "RNA u> reference should emit DNA T (g.1003T>G), got: {s}"
            );
            assert!(
                !s.contains('U'),
                "g. output must not carry an RNA U base, got: {s}"
            );
        }

        // -- 8b: NG_/LRG_ re-anchor decline → error, not invalid HGVS (#655) ---

        /// An `NG_` parent with no registered placement cannot be re-anchored:
        /// cdot carries only the transcript's chromosome (`NC_`) alignment.
        /// `project_to_genomic` must decline rather than stamp the chromosome
        /// coordinate under the `NG_` accession (invalid HGVS) (#655).
        #[test]
        fn project_to_genomic_ng_parent_without_placement_declines() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("TEST", 1));
            let err = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect_err("NG_ parent without a placement must decline, not emit invalid HGVS");
            match err {
                FerroError::UnsupportedProjection { reason } => assert!(
                    reason.contains("NG_TEST.1") && reason.contains("no chromosomal placement"),
                    "expected a no-placement decline naming NG_TEST.1, got: {reason}"
                ),
                other => panic!("expected UnsupportedProjection, got: {other:?}"),
            }
        }

        /// A placement exists for the `NG_` parent, but the variant's chromosome
        /// coordinate falls outside the placed span, so it cannot be re-anchored
        /// into the parent frame. Decline rather than emit a chromosome
        /// coordinate under the `NG_` accession (#655).
        #[test]
        fn project_to_genomic_ng_parent_endpoint_outside_placement_declines() {
            let (projector, mut provider) = make_test_provider_and_projector();
            // NM_TEST.1 c.4 maps to chromosome 1003; place NG_TEST.1 on a
            // disjoint span [2000, 2010] so nc_to_parent(1003) declines.
            provider.add_genomic_placement(
                "NG_TEST.1",
                crate::reference::GenomicPlacement {
                    nc: Accession::new("NC", "000001", Some(11)),
                    parent_start: 1,
                    nc_start: 2000,
                    nc_end: 2010,
                    strand: crate::reference::Strand::Plus,
                },
            );
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("TEST", 1));
            let err = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect_err("an endpoint outside the placed span must decline");
            match err {
                FerroError::UnsupportedProjection { reason } => assert!(
                    reason.contains("outside the placed genomic span"),
                    "expected an out-of-span decline, got: {reason}"
                ),
                other => panic!("expected UnsupportedProjection, got: {other:?}"),
            }
        }

        /// De-anchor path: when the genomic axis cannot be re-anchored into the
        /// `NG_` parent frame (here the endpoint is outside the placement span),
        /// `frame_projection_owned` drops the genomic axis (reports `None`)
        /// rather than stamping the parent accession on an out-of-frame
        /// chromosome coordinate (#655). The transcript-relative coding form is
        /// still re-framed under the parent — degrading the unframable axis is
        /// preferable to dropping the whole projection.
        #[test]
        fn frame_projection_owned_drops_genomic_when_reanchor_declines() {
            let parent = ng_parent("900", 1);
            let placement = crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1010,
                strand: crate::reference::Strand::Plus,
            };
            // Genomic axis at chromosome 5000 — well outside the placed span.
            let genomic = HgvsVariant::Genome(GenomeVariant {
                accession: Accession::new("NC", "000001", Some(11)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(5000)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            });
            let coding = HgvsVariant::Cds(CdsVariant {
                accession: parse_accession("NM_900.1"),
                gene_symbol: Some("GENE900".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            });
            let proj = VariantProjection {
                genomic: Some(genomic),
                coding: Some(coding),
                noncoding: None,
                protein: None,
                rna: None,
                transcript_id: "NM_900.1".to_string(),
                gene_symbol: Some("GENE900".to_string()),
                is_frameshift: false,
                is_intronic: false,
                is_utr: false,
                affects_init: false,
            };
            let vp = VariantProjector::new(Projector::new(CdotMapper::new()), MockProvider::new());
            let framed = vp.frame_projection_owned(proj, &parent, Some(&placement));
            assert!(
                framed.genomic.is_none(),
                "genomic axis should be dropped when it cannot be re-anchored, got: {:?}",
                framed.genomic
            );
            // The coding axis is still re-framed under the NG_ parent.
            let coding_str = framed
                .coding
                .as_ref()
                .expect("coding axis retained")
                .to_string();
            assert!(
                coding_str.starts_with("NG_900.1("),
                "coding should be framed under the NG_ parent, got: {coding_str}"
            );
        }

        /// #860 Case B: under a bare `LRG_<n>` genomic parent, the fan-out
        /// coding/protein axes are relabeled to the input's LRG namespace
        /// (`LRG_<n>t<k>`/`LRG_<n>p<k>`) and rendered **bare** — not
        /// `LRG_<n>(LRG_<n>t<k>)` and not `LRG_<n>(NM_…)`.
        #[test]
        fn frame_projection_owned_relabels_lrg_namespace_bare() {
            use crate::hgvs::edit::ProteinEdit;
            use crate::hgvs::interval::ProtInterval;
            use crate::hgvs::location::{AminoAcid, ProtPos};
            use crate::hgvs::variant::ProteinVariant;

            let mut cdot = CdotMapper::new();
            // LRG file carries .1; the fan-out coding axis is .3 — base match.
            cdot.insert_lrg_mapping("LRG_24t1".to_string(), "NM_001114101.1".to_string());
            cdot.insert_lrg_mapping("LRG_24t2".to_string(), "NM_172369.2".to_string());
            let vp = VariantProjector::new(Projector::new(cdot), MockProvider::new());

            let parent = parse_accession("LRG_24"); // bare genomic LRG, is_lrg()
            let coding = HgvsVariant::Cds(CdsVariant {
                accession: parse_accession("NM_001114101.3"),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(127)),
                    NaEdit::Substitution {
                        reference: Base::G,
                        alternative: Base::A,
                    },
                ),
            });
            let protein = HgvsVariant::Protein(ProteinVariant {
                accession: parse_accession("NP_001107573.1"),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    ProtInterval::point(ProtPos::new(AminoAcid::Met, 1)),
                    ProteinEdit::whole_protein_unknown(),
                ),
            });
            let proj = VariantProjection {
                genomic: None,
                coding: Some(coding),
                noncoding: None,
                protein: Some(protein),
                rna: None,
                transcript_id: "NM_001114101.3".to_string(),
                gene_symbol: None,
                is_frameshift: false,
                is_intronic: false,
                is_utr: false,
                affects_init: false,
            };
            let framed = vp.frame_projection_owned(proj, &parent, None);

            assert_eq!(
                framed.coding.unwrap().to_string(),
                "LRG_24t1:c.127G>A",
                "coding should render bare LRG_24t1 (no NM_, no LRG_24(...) wrapper)"
            );
            assert!(
                framed.protein.unwrap().to_string().starts_with("LRG_24p1:"),
                "protein should render bare LRG_24p1"
            );
        }

        /// `reanchor_genome_to_parent` declines an endpoint with no single
        /// resolved position (a compound `(a_b)` `UncertainBoundary::Range`)
        /// with `InvalidCoordinates` rather than stamping the parent accession
        /// on an unresolved coordinate (#655). This decline is defensive — the
        /// public `project_to_genomic` pivot rejects uncertain endpoints
        /// earlier — so it is exercised here directly.
        #[test]
        fn reanchor_genome_to_parent_declines_uncertain_endpoint() {
            use crate::hgvs::interval::UncertainBoundary;
            use crate::hgvs::Mu;
            let placement = crate::reference::GenomicPlacement {
                nc: Accession::new("NC", "000001", Some(11)),
                parent_start: 1,
                nc_start: 1000,
                nc_end: 1010,
                strand: crate::reference::Strand::Plus,
            };
            // The end boundary is a compound range, so `.inner()` is `None`:
            // there is no single position to re-anchor.
            let mut interval = GenomeInterval::new(GenomePos::new(1003), GenomePos::new(1006));
            interval.end = UncertainBoundary::range(
                Mu::Certain(GenomePos::new(1004)),
                Mu::Certain(GenomePos::new(1006)),
            );
            let gv = GenomeVariant {
                accession: ng_parent("900", 1),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    interval,
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let err = VariantProjector::<MockProvider>::reanchor_genome_to_parent(gv, &placement)
                .expect_err("an uncertain endpoint must decline");
            match err {
                FerroError::InvalidCoordinates { msg } => assert!(
                    msg.contains("single resolved position"),
                    "expected an uncertain-position decline, got: {msg}"
                ),
                other => panic!("expected InvalidCoordinates, got: {other:?}"),
            }
        }

        // -- 8c: single-transcript project_variant re-anchors / degrades the
        //        genomic axis for an NG_/LRG_ parent (#702) --------------------

        /// `project_variant` re-anchors the `.genomic` axis into an `NG_`
        /// parent's own frame when a placement is known (#480/#702), rather than
        /// leaving the chromosome-frame pivot stamped with the `NG_` accession.
        #[test]
        fn project_variant_reanchors_genomic_into_ng_parent_with_placement() {
            let (projector, mut provider) = make_test_provider_and_projector();
            // NM_TEST.1 c.4 maps to chromosome 1003; place NG_TEST.1 so chromosome
            // 1000 == NG base 1, hence chromosome 1003 == NG_TEST.1:g.4.
            provider.add_genomic_placement(
                "NG_TEST.1",
                crate::reference::GenomicPlacement {
                    nc: Accession::new("NC", "000001", Some(11)),
                    parent_start: 1,
                    nc_start: 1000,
                    nc_end: 1010,
                    strand: crate::reference::Strand::Plus,
                },
            );
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("TEST", 1));
            let proj = vp
                .project_variant(&HgvsVariant::Cds(cds), "NM_TEST.1")
                .expect("project_variant should succeed");
            let g = proj
                .genomic
                .expect("genomic axis should be present and re-anchored");
            let s = g.to_string();
            assert!(
                s.starts_with("NG_TEST.1") && s.contains(":g.4C>A"),
                "expected re-anchored NG_TEST.1:g.4C>A, got: {s}"
            );
            assert!(proj.coding.is_some(), "coding axis should be present");
        }

        /// `project_variant` drops the `.genomic` axis to `None` for an `NG_`
        /// parent with no known placement rather than stamping the `NG_`
        /// accession on a chromosome coordinate (#655/#702). The coding/protein
        /// axes — still valid — are kept; unlike the g.-only `project_to_genomic`,
        /// the multi-axis path degrades rather than hard-errors.
        #[test]
        fn project_variant_drops_genomic_for_ng_parent_without_placement() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("TEST", 1));
            let proj = vp
                .project_variant(&HgvsVariant::Cds(cds), "NM_TEST.1")
                .expect("project_variant should still succeed (multi-axis)");
            assert!(
                proj.genomic.is_none(),
                "genomic axis should be dropped for an unframable NG_ parent, got: {:?}",
                proj.genomic
            );
            assert!(proj.coding.is_some(), "coding axis should be retained");
        }

        /// An `NC_` chromosome parent shares cdot's frame, so `project_variant`
        /// emits its `.genomic` unchanged (regression guard for the #702
        /// re-anchor path).
        #[test]
        fn project_variant_emits_nc_parent_genomic_unchanged() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(4)),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let proj = vp
                .project_variant(&HgvsVariant::Cds(cds), "NM_TEST.1")
                .expect("project_variant should succeed");
            let g = proj.genomic.expect("genomic axis present for NC_ parent");
            let s = g.to_string();
            assert!(
                s.starts_with("NC_000001.11") && s.contains(":g.1003C>A"),
                "expected NC_000001.11:g.1003C>A unchanged, got: {s}"
            );
        }

        // -- 8b: #887 — NG-parented pter/qter termini resolve on project_variant ---

        /// A projector whose provider knows only the NG parent's sequence, so a
        /// whole-arm terminus resolves against it. The reference's first and last
        /// bases are deliberately NOT part of a homopolymer run (`C…G`) so that
        /// `normalize(g.1del)` stays `g.1del` and `normalize(g.<len>del)` stays
        /// `g.<len>del` (an all-`A` ref would 3'-shift both to the same
        /// coordinate). Returns the projector and the transcript id to project
        /// against (`NM_TEST.1`, the input's own transcript).
        fn ng_parent_projector_with_seq(seq: &str) -> (VariantProjector<MockProvider>, String) {
            let projector = Projector::new(CdotMapper::new());
            let mut provider = MockProvider::new();
            provider.add_genomic_sequence("NG_TEST.1", seq.to_string());
            (
                VariantProjector::new(projector, provider),
                "NM_TEST.1".to_string(),
            )
        }

        #[test]
        fn project_variant_ng_parent_pter_resolves_genomic_887() {
            // len = 50, first base `C` (not shiftable), last base `G`.
            let (projector, tx_id) = ng_parent_projector_with_seq(&format!("C{}G", "A".repeat(48)));
            let v = crate::parse_hgvs("NG_TEST.1(NM_TEST.1):c.pterdel").unwrap();
            let proj = projector
                .project_variant(&v, &tx_id)
                .expect("pter must resolve, not error");
            assert_eq!(
                proj.genomic.as_ref().map(|g| g.to_string()),
                Some("NG_TEST.1:g.1del".to_string())
            );
            assert!(proj.protein.is_none());
            assert!(proj.coding.is_some(), "coding echoes the bare input form");
        }

        #[test]
        fn project_variant_ng_parent_qter_resolves_genomic_887() {
            let (projector, tx_id) = ng_parent_projector_with_seq(&format!("C{}G", "A".repeat(48)));
            let v = crate::parse_hgvs("NG_TEST.1(NM_TEST.1):c.qterdel").unwrap();
            let proj = projector
                .project_variant(&v, &tx_id)
                .expect("qter must resolve");
            assert_eq!(
                proj.genomic.as_ref().map(|g| g.to_string()),
                Some("NG_TEST.1:g.50del".to_string())
            );
        }

        #[test]
        fn project_variant_ng_parent_pter_qter_range_resolves_887() {
            let (projector, tx_id) = ng_parent_projector_with_seq(&format!("C{}G", "A".repeat(48)));
            let v = crate::parse_hgvs("NG_TEST.1(NM_TEST.1):c.pter_qterdel").unwrap();
            let proj = projector
                .project_variant(&v, &tx_id)
                .expect("pter_qter range must resolve");
            assert_eq!(
                proj.genomic.as_ref().map(|g| g.to_string()),
                Some("NG_TEST.1:g.1_50del".to_string())
            );
        }

        /// An unsupported terminus marker pair (`cen`) must decline with `None`
        /// — falling through to the normal rejection path — **even when the
        /// parent-length lookup would itself fail**. The marker-pair
        /// classification runs *before* the parent-length lookup, so an
        /// unsupported pair never leaks an incidental length-lookup `Some(Err)`
        /// (which is meaningful only for a supported `pter`/`qter` terminus) and
        /// so aborts the whole projection instead of the intended fall-through.
        #[test]
        fn unsupported_cen_terminus_declines_before_parent_length_lookup_887() {
            let (projector, _tx) = ng_parent_projector_with_seq(&format!("C{}G", "A".repeat(48)));
            // Parent `NG_MISSING.1` is unknown to the provider, so a length
            // lookup against it errors — the exact condition that previously
            // leaked a `Some(Err(..))` for the unsupported `cen` pair.
            let HgvsVariant::Cds(cds) = crate::parse_hgvs("NM_TEST.1:c.cendel").unwrap() else {
                panic!("c.cendel must parse as a Cds variant");
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("MISSING", 1));
            let variant = HgvsVariant::Cds(cds);
            assert!(
                projector.project_cds_terminus_to_parent(&variant).is_none(),
                "unsupported cen terminus must decline with None, not surface a \
                 parent-length-lookup error",
            );
        }

        #[test]
        fn project_variant_pter_transcript_id_mismatch_still_errors_887() {
            // The early terminus guard must sit AFTER a transcript_id match check,
            // so a mismatched target still errors instead of silently resolving.
            let (projector, _tx) = ng_parent_projector_with_seq(&format!("C{}G", "A".repeat(48)));
            let v = crate::parse_hgvs("NG_TEST.1(NM_TEST.1):c.pterdel").unwrap();
            assert!(projector.project_variant(&v, "NM_OTHER.9").is_err());
        }

        // -- 9: protein input → unsupported ------------------------------------

        #[test]
        fn project_to_genomic_unsupported_for_protein() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            let input = crate::parse_hgvs("NP_TEST.1:p.Arg2Cys").expect("p. parse should succeed");
            assert!(matches!(input, HgvsVariant::Protein(_)));
            let err = vp
                .project_to_genomic(&input)
                .expect_err("p. → g. should be unsupported");
            assert!(
                matches!(err, FerroError::UnsupportedProjection { .. }),
                "expected UnsupportedProjection for p., got: {:?}",
                err
            );
        }

        // -- 10: genome allele input → idempotent (#851) -----------------------

        #[test]
        fn project_to_genomic_genome_allele_is_idempotent() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            let v1 = crate::parse_hgvs("chr1:g.1003C>A").expect("parse v1");
            let v2 = crate::parse_hgvs("chr1:g.1006T>A").expect("parse v2");
            let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![v1, v2]));
            // Members are already Genome → each passes through unchanged and the
            // allele reassembles as itself (#851 routes alleles per member).
            let out = vp
                .project_to_genomic(&allele)
                .expect("genome allele projects idempotently (#851)");
            assert_eq!(out, allele);
        }

        // -- 11: ? position sentinel → unsupported -----------------------------

        #[test]
        fn project_to_genomic_unknown_position_returns_unsupported() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // Build NG_TEST.1(NM_TEST.1):c.?C>A — start/end are CdsPos::unknown(None).
            let unknown_pos = CdsPos::unknown(None);
            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1").with_genomic_context(ng_parent("TEST", 1)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(unknown_pos),
                    NaEdit::Substitution {
                        reference: Base::C,
                        alternative: Base::A,
                    },
                ),
            };
            let input = HgvsVariant::Cds(cds);
            let err = vp
                .project_to_genomic(&input)
                .expect_err("? position should be unsupported");
            assert!(
                matches!(err, FerroError::UnsupportedProjection { .. }),
                "expected UnsupportedProjection for unknown position, got: {:?}",
                err
            );
        }

        // -- 12: parsed `c.?` routes to UnsupportedProjection ------------------

        /// Pins that a *parsed* unknown-position input (`UncertainBoundary::
        /// Single(Mu::Unknown)`) routes to `UnsupportedProjection`, not
        /// `InvalidCoordinates`. Test 11 above covers the synthetic
        /// `CdsPos::unknown()` shape; this one covers the parser-produced
        /// shape, which has different AST internals and previously slipped
        /// through `.inner()` as `None` indistinguishable from a Range.
        #[test]
        fn project_to_genomic_parsed_question_position_returns_unsupported() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // `c.?_10del` parses to a Cds variant whose start endpoint is
            // `UncertainBoundary::Single(Mu::Unknown)`. We need to attach a
            // genomic_context so the variant kind is correct; reading the
            // parsed variant out and mutating its accession is sufficient.
            let parsed = crate::parse_hgvs("NM_TEST.1:c.?_10del").expect("parse should succeed");
            let cds = match parsed {
                HgvsVariant::Cds(c) => c,
                other => panic!("expected Cds variant, got: {:?}", other),
            };
            let cds = attach_genomic_context_cds(cds, ng_parent("TEST", 1));
            let input = HgvsVariant::Cds(cds);

            let err = vp
                .project_to_genomic(&input)
                .expect_err("parsed c.? should be unsupported");
            match err {
                FerroError::UnsupportedProjection { reason } => {
                    assert!(
                        reason.contains("unknown position") || reason.contains("`?`"),
                        "expected reason to mention unknown position, got: {}",
                        reason
                    );
                }
                other => panic!(
                    "expected UnsupportedProjection for parsed c.?_10del, got: {:?}",
                    other
                ),
            }
        }

        // -- 13 + 14: UTR position projection (c.*N, c.-N) ---------------------

        /// Build a fixture with 5'UTR (3 bases) + CDS (9 bases) + 3'UTR (3 bases) so
        /// c.-1, c.1, c.9, c.*1 etc. all have well-defined genomic counterparts.
        ///
        /// Layout (plus strand):
        ///   Exon 1: genome [1000, 1015), tx [0, 15)  — 15 bp
        ///   cds_start = 3 (0-based)  ⇒ c.-3 = g.1000, c.-1 = g.1002
        ///   cds_end   = 12 (exclusive) ⇒ c.1 = g.1003, c.9 = g.1011, c.*1 = g.1012
        ///   sequence  = "TTTATGCGCTAAGGG" (5'UTR + CDS + 3'UTR)
        fn make_utr_test_provider_and_projector() -> (Projector, MockProvider) {
            let mut cdot = CdotMapper::new();
            cdot.add_transcript(
                "NM_UTR.1".to_string(),
                CdotTranscript {
                    gene_name: Some("UTRGENE".to_string()),
                    contig: "chr1".to_string(),
                    strand: ProvStrand::Plus,
                    exons: vec![[1000, 1015, 0, 15]],
                    cds_start: Some(3),
                    cds_end: Some(12),
                    gene_id: None,
                    protein: Some("NP_UTR.1".to_string()),
                    exon_cigars: Vec::new(),
                },
            );
            let projector = Projector::new(cdot);

            let mut provider = MockProvider::new();
            provider.add_transcript(Transcript {
                id: "NM_UTR.1".to_string(),
                gene_symbol: Some("UTRGENE".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("TTTATGCGCTAAGGG".to_string()),
                cds_start: Some(4), // 1-based inclusive
                cds_end: Some(12),
                exons: vec![Exon::with_genomic(1, 1, 15, 1000, 1014)],
                chromosome: Some("chr1".to_string()),
                genomic_start: Some(1000),
                genomic_end: Some(1014),
                genome_build: Default::default(),
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
                protein_id: Some("NP_UTR.1".to_string()),
            });
            // 999 leading Ns so HGVS g.1000 (the exon's `genome_start`) maps to
            // 0-based index 999, matching every other fixture and the real
            // genome-coordinate convention (HGVS g.X == 0-based X-1). The exon
            // genome and transcript sequences are then identical, so the #644
            // sequence-aware correction is a no-op here.
            let prefix = "N".repeat(999);
            let suffix = "N".repeat(100);
            provider.add_genomic_sequence("chr1", format!("{}TTTATGCGCTAAGGG{}", prefix, suffix));
            (projector, provider)
        }

        #[test]
        fn project_to_genomic_3utr_star_position() {
            let (projector, provider) = make_utr_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // c.*1A>G on NM_UTR.1 → tx_pos = cds_end (12) → genome 1012.
            // (Without the utr3-aware fix in data::mapping::cds_to_genome, ferro
            // silently mapped c.*1 as c.1 and produced g.1003.)
            let cds = CdsVariant {
                accession: parse_accession("NM_UTR.1"),
                gene_symbol: Some("UTRGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos {
                        base: 1,
                        offset: None,
                        utr3: true,
                        special: None,
                    }),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::G,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);

            let out = vp
                .project_to_genomic(&input)
                .expect("c.*1A>G should project");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let s = g.to_string();
            assert!(
                s.contains(":g.1012A>G"),
                "expected NC_000001.11...:g.1012A>G for c.*1A>G, got: {}",
                s
            );
        }

        /// A transcript whose cdot exon→genome map (e.g. #790-derived) **strips**
        /// the post-transcriptional poly-A tail: the cdot exon ends at tx 12
        /// (genome 1011 inclusive / 1012 exclusive), but the transcript FASTA is
        /// 15 bases long — a 3-base poly-A tail (tx 12..15) with no exon. This
        /// mirrors NM_003002.2 (#797), where the last exon ends at tx 1364 but
        /// the FASTA is 1382 bases.
        fn make_polya_test_provider_and_projector() -> (Projector, MockProvider) {
            let mut cdot = CdotMapper::new();
            cdot.add_transcript(
                "NM_POLYA.1".to_string(),
                CdotTranscript {
                    gene_name: Some("POLYAGENE".to_string()),
                    contig: "chr1".to_string(),
                    strand: ProvStrand::Plus,
                    // exon covers only the genomic core (poly-A stripped): tx 0..12.
                    exons: vec![[1000, 1012, 0, 12]],
                    cds_start: Some(3),
                    cds_end: Some(9),
                    gene_id: None,
                    protein: Some("NP_POLYA.1".to_string()),
                    exon_cigars: Vec::new(),
                },
            );
            let projector = Projector::new(cdot);

            let mut provider = MockProvider::new();
            // FASTA carries the full 15-base transcript incl. the 3-base poly-A
            // tail ("AAA"); cds is c.1..c.6 (tx 3..9), 3'UTR core c.*1..c.*3
            // (tx 9..12), poly-A c.*4..c.*6 (tx 12..15).
            provider.add_transcript(Transcript {
                id: "NM_POLYA.1".to_string(),
                gene_symbol: Some("POLYAGENE".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("TTTATGCGCGTAAA".to_string() + "A"), // 15 bases, trailing AAA tail
                cds_start: Some(4),                                 // 1-based inclusive
                cds_end: Some(9),
                exons: vec![Exon::with_genomic(1, 1, 12, 1000, 1011)],
                chromosome: Some("chr1".to_string()),
                genomic_start: Some(1000),
                genomic_end: Some(1011),
                genome_build: Default::default(),
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
                protein_id: Some("NP_POLYA.1".to_string()),
            });
            // Genomic bases: 999 leading bases so HGVS g.1000 == 0-based index 999,
            // then the 12-base exon core (indices 999..1010), then a downstream
            // tail past the stripped exon (g.1012, g.1013, …). The poly-A walk is
            // base-agnostic, so the downstream bases don't affect the *walk*; but
            // they DO feed the genomic-frame renormalization the public
            // `project_to_genomic` runs on its output, so the flanks are
            // **non-repeating** (a repeating `ACGT`/`CAGT` cycle, no homopolymer at
            // any single position) — a single-base del at g.1012 / g.1013 then
            // stays put under 3'-shifting, exercising normalization without moving
            // the asserted coordinate. The same bases are registered under the
            // `NC_000001.11` output accession so normalization can fetch them
            // (registering only `chr1` would silently fall back to the
            // un-normalized form; #797 review).
            let genomic = polya_genomic_sequence();
            provider.add_genomic_sequence("chr1", genomic.clone());
            provider.add_genomic_sequence("NC_000001.11", genomic);
            (projector, provider)
        }

        /// Build the shared genomic backbone for the poly-A fixtures: a
        /// non-repeating `ACGT`-cycle prefix (indices 0..999), the 12-base exon
        /// core `TTTATGCGCGTA` (indices 999..1011 → HGVS g.1000..g.1011), and a
        /// non-repeating downstream `CAGT`-cycle tail (indices 1011.. → g.1012+).
        /// Non-repeating flanks keep single-base dels at the asserted poly-A
        /// coordinates stable under genomic-frame renormalization.
        fn polya_genomic_sequence() -> String {
            let cycle = |pat: &str, len: usize| -> String {
                pat.chars().cycle().take(len).collect::<String>()
            };
            let prefix = cycle("ACGT", 999);
            let exon = "TTTATGCGCGTA"; // 12 bases, indices 999..1011
            let downstream = cycle("CAGT", 100); // indices 1011..
            format!("{prefix}{exon}{downstream}")
        }

        /// A `c.*` position in the genomic CORE of the 3'UTR projects normally
        /// (unchanged by #797). c.*1 = tx 9 → genome 1009.
        #[test]
        fn project_to_genomic_polya_core_3utr_unchanged() {
            let (projector, provider) = make_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_POLYA.1"),
                gene_symbol: Some("POLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos {
                        base: 1,
                        offset: None,
                        utr3: true,
                        special: None,
                    }),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let out = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect("c.*1del should project (genomic-core 3'UTR)");
            let s = match out {
                HgvsVariant::Genome(ref g) => g.to_string(),
                _ => panic!("expected Genome variant"),
            };
            assert!(
                s.contains(":g.1009del"),
                "expected g.1009del for c.*1del, got: {s}"
            );
        }

        /// A `c.*` position in the POLY-A region (past the stripped exon end, but
        /// within the true transcript length) projects via the contiguous
        /// downstream genome walk (#797). c.*4 = tx 12 → genome 1012 (exon
        /// genome_end, exclusive → first downstream base).
        #[test]
        fn project_to_genomic_polya_region_walks_downstream() {
            let (projector, provider) = make_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_POLYA.1"),
                gene_symbol: Some("POLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos {
                        base: 4,
                        offset: None,
                        utr3: true,
                        special: None,
                    }),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let out = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect("c.*4del (poly-A region) should project via contiguous walk");
            let s = match out {
                HgvsVariant::Genome(ref g) => g.to_string(),
                _ => panic!("expected Genome variant"),
            };
            assert!(
                s.contains(":g.1012del"),
                "expected g.1012del for c.*4del (poly-A walk), got: {s}"
            );
        }

        /// A `c.*N+k` whose base is in the poly-A region folds the downstream
        /// offset linearly (#797), matching mutalyzer / the normalizer's
        /// `c.*N+k → c.*(N+k)`. c.*4+1 → effective tx 13 → genome 1013.
        #[test]
        fn project_to_genomic_polya_region_with_offset_folds_linearly() {
            let (projector, provider) = make_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_POLYA.1"),
                gene_symbol: Some("POLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos {
                        base: 4,
                        offset: Some(1),
                        utr3: true,
                        special: None,
                    }),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let out = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect("c.*4+1del (poly-A region, +offset) should project");
            let s = match out {
                HgvsVariant::Genome(ref g) => g.to_string(),
                _ => panic!("expected Genome variant"),
            };
            assert!(
                s.contains(":g.1013del"),
                "expected g.1013del for c.*4+1del (poly-A walk + offset), got: {s}"
            );
        }

        /// #797 (multi-axis): a `c.*` poly-A position routed through the full
        /// `project_variant` / `project` / `project_variant_all` public APIs must
        /// not be rejected as `TranscriptNotOverlapping`. The synthetic poly-A
        /// genome coordinate is outside the cdot exon span by construction, so the
        /// usual genome→CDS re-derivation cannot round-trip it; the projection must
        /// short-circuit, keeping the input's own coding form and the contiguous
        /// downstream genomic coordinate.
        #[test]
        fn project_variant_polya_region_multiaxis_short_circuits() {
            let (projector, provider) = make_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_POLYA.1"),
                gene_symbol: Some("POLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos {
                        base: 4,
                        offset: None,
                        utr3: true,
                        special: None,
                    }),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let proj = vp
                .project_variant(&HgvsVariant::Cds(cds), "NM_POLYA.1")
                .expect("poly-A c.*4del must project through project_variant (multi-axis)");
            // Genomic axis: the contiguous downstream walk, NC parent unchanged.
            let g = proj
                .genomic
                .as_ref()
                .expect("poly-A multi-axis projection must carry a genomic axis");
            assert!(
                g.to_string().contains(":g.1012del"),
                "expected g.1012del genomic axis, got: {g}"
            );
            // Coding axis: the input's own c.*4del form is preserved (poly-A
            // positions cannot be re-derived from the synthetic genome coordinate).
            let c = proj
                .coding
                .as_ref()
                .expect("poly-A multi-axis projection must carry a coding axis");
            assert!(
                c.to_string().contains("c.*4del"),
                "expected the input c.*4del coding axis preserved, got: {c}"
            );
            // A 3'UTR / poly-A position has no protein consequence.
            assert!(!proj.is_intronic, "poly-A 3'UTR position is not intronic");
        }

        /// Minus-strand analogue: a `c.*` poly-A position projects through the
        /// multi-axis public API without a `TranscriptNotOverlapping` rejection.
        #[test]
        fn project_variant_minus_polya_region_multiaxis_short_circuits() {
            let (projector, provider) = make_minus_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_MPOLYA.1"),
                gene_symbol: Some("MPOLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos {
                        base: 4,
                        offset: None,
                        utr3: true,
                        special: None,
                    }),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let proj = vp
                .project_variant(&HgvsVariant::Cds(cds), "NM_MPOLYA.1")
                .expect("minus-strand poly-A c.*4del must project through project_variant");
            let g = proj
                .genomic
                .as_ref()
                .expect("poly-A multi-axis projection must carry a genomic axis");
            assert!(
                g.to_string().contains(":g.999del"),
                "expected g.999del genomic axis, got: {g}"
            );
            let c = proj
                .coding
                .as_ref()
                .expect("poly-A multi-axis projection must carry a coding axis");
            assert!(
                c.to_string().contains("c.*4del"),
                "expected the input c.*4del coding axis preserved, got: {c}"
            );
        }

        /// #797 (multi-axis gate): the poly-A short-circuit must fire **only** for
        /// a `c.*` (CDS) endpoint. A bare transcript-coordinate `n.` input that
        /// lands off the cdot exon span is NOT a poly-A 3'UTR endpoint — it has no
        /// `c.*` poly-A semantics — so it must still raise the canonical
        /// `TranscriptNotOverlapping` decline rather than being smuggled through
        /// the bypass and returned in the `.coding` slot.
        ///
        /// `n.13` on `NM_POLYA.1` is tx index 12, which is off the cdot exon
        /// (`tx 0..12`) — the same genome neighborhood the poly-A walk targets for
        /// the `c.*4` case. Because the poly-A genome→genome walk fallback in
        /// `project_to_genomic_nc` self-gates on a `c.*` (`p.utr3`) endpoint, the
        /// `n.` input never synthesizes an outside-span genome coordinate: it
        /// declines cleanly during projection rather than being smuggled through
        /// the multi-axis short-circuit and returned in the `.coding` slot. With
        /// the looser `!Genome` gate, a future `project_to_genomic_nc` that did
        /// synthesize such a coordinate for an `n.`/`r.` input would have wrongly
        /// taken the bypass; the `Cds`-only gate forecloses that.
        #[test]
        fn project_variant_offexon_noncoding_input_is_not_polya_short_circuited() {
            use crate::hgvs::interval::TxInterval;
            use crate::hgvs::location::TxPos;
            use crate::hgvs::variant::TxVariant;

            let (projector, provider) = make_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let tx = TxVariant {
                accession: parse_accession("NM_POLYA.1"),
                gene_symbol: Some("POLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    TxInterval::point(TxPos::new(13)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let tx = attach_genomic_context_tx(tx, nc_parent());
            let result = vp.project_variant(&HgvsVariant::Tx(tx), "NM_POLYA.1");
            let err = result.expect_err(
                "an off-exon n. input is not a c.* poly-A endpoint and must \
                 decline, not take the poly-A multi-axis bypass",
            );
            // The decline is a clean projection failure, NOT a multi-axis bypass:
            // the input must never come back with the `n.` form parked in the
            // `.coding` slot (which is what the over-broad `!Genome` gate allowed).
            assert!(
                matches!(
                    err,
                    FerroError::TranscriptNotOverlapping { .. }
                        | FerroError::InvalidCoordinates { .. }
                ),
                "expected a clean projection decline for an off-exon n. input, \
                 got: {err:?}",
            );
        }

        /// #797 (poly-A fallback gate, c.-only): the `project_to_genomic`
        /// poly-A walk fallback must fire **only** for a `c.` (CDS) input. An
        /// off-exon `n.*` (`Tx`) endpoint is reshaped into a `CdsPos` with
        /// `utr3` set for this projection closure, so without the `c.`-only gate
        /// it would be reinterpreted as a `c.*` poly-A walk and emit a genomic
        /// coordinate — a coordinate class that has no poly-A semantics. It must
        /// instead decline. `n.13` on `NM_POLYA.1` is tx index 12, off the cdot
        /// exon (`tx 0..12`) — the genome neighborhood the poly-A walk targets.
        #[test]
        fn project_to_genomic_offexon_noncoding_n_input_declines() {
            use crate::hgvs::interval::TxInterval;
            use crate::hgvs::location::TxPos;
            use crate::hgvs::variant::TxVariant;

            let (projector, provider) = make_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let tx = TxVariant {
                accession: parse_accession("NM_POLYA.1"),
                gene_symbol: Some("POLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    TxInterval::point(TxPos::new(13)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let tx = attach_genomic_context_tx(tx, nc_parent());
            let err = vp
                .project_to_genomic(&HgvsVariant::Tx(tx))
                .expect_err("off-exon n. input must decline, not take the c.* poly-A walk");
            assert!(
                matches!(
                    err,
                    FerroError::TranscriptNotOverlapping { .. }
                        | FerroError::InvalidCoordinates { .. }
                ),
                "expected a clean decline for an off-exon n. input, got: {err:?}",
            );
        }

        /// #797 (poly-A fallback gate, c.-only): RNA (`r.`) analogue of the `n.`
        /// decline above. An off-exon `r.*` endpoint must not be reinterpreted as
        /// a `c.*` poly-A walk by `project_to_genomic`.
        #[test]
        fn project_to_genomic_offexon_rna_input_declines() {
            use crate::hgvs::interval::RnaInterval;
            use crate::hgvs::location::RnaPos;
            use crate::hgvs::variant::RnaVariant;

            let (projector, provider) = make_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            // r.*4 reshapes to a CdsPos { base: 4, utr3: true } — the same shape a
            // genuine c.*4 poly-A endpoint carries — so it directly exercises the
            // `c.`-only gate: an RNA input with this shape must still decline.
            let rna = RnaVariant {
                accession: parse_accession("NM_POLYA.1"),
                gene_symbol: Some("POLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    RnaInterval::point(RnaPos {
                        base: 4,
                        offset: None,
                        utr3: true,
                    }),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let rna = {
                let mut r = rna;
                r.accession = r.accession.with_genomic_context(nc_parent());
                r
            };
            let err = vp
                .project_to_genomic(&HgvsVariant::Rna(rna))
                .expect_err("off-exon r. input must decline, not take the c.* poly-A walk");
            assert!(
                matches!(
                    err,
                    FerroError::TranscriptNotOverlapping { .. }
                        | FerroError::InvalidCoordinates { .. }
                ),
                "expected a clean decline for an off-exon r. input, got: {err:?}",
            );
        }

        /// #797 (boundary-spanning): a `c.*` interval that spans the core→poly-A
        /// boundary — one endpoint inside the stripped terminal exon (3'UTR
        /// genomic core), the other in the poly-A tail past the exon — must
        /// project through the multi-axis short-circuit, not be rejected as
        /// `TranscriptNotOverlapping`. On `NM_POLYA.1`, `c.*3` is tx 11 (genome
        /// 1011, inside exon `tx 0..12`) and `c.*4` is tx 12 (genome 1012, poly-A
        /// tail). With the over-strict `outside(start) && outside(end)` gate the
        /// single outside (poly-A) endpoint was rejected; the `any-outside` gate
        /// short-circuits correctly.
        #[test]
        fn project_variant_polya_boundary_spanning_short_circuits() {
            let (projector, provider) = make_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_POLYA.1"),
                gene_symbol: Some("POLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::new(
                        CdsPos {
                            base: 3,
                            offset: None,
                            utr3: true,
                            special: None,
                        },
                        CdsPos {
                            base: 4,
                            offset: None,
                            utr3: true,
                            special: None,
                        },
                    ),
                    // A delins (not a plain del) so the boundary-spanning interval
                    // does not collapse into the surrounding poly-A homopolymer as a
                    // repeat-contraction during normalization — that would shift the
                    // interval and obscure the boundary-spanning shape under test.
                    NaEdit::Delins {
                        sequence: InsertedSequence::Literal("GG".parse().unwrap()),
                        deleted: None,
                        deleted_length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let proj = vp
                .project_variant(&HgvsVariant::Cds(cds), "NM_POLYA.1")
                .expect("boundary-spanning c.*3_*4delinsGG must project (multi-axis)");
            // Genomic axis: the walked interval (g.1011 inside exon, g.1012 the
            // first poly-A walk coordinate) re-anchored under the NC parent. The
            // public path renormalizes it in the genome frame, so assert the axis
            // is present on the NC contig rather than an exact (pre-normalization)
            // coordinate — the load-bearing assertion is that it short-circuits at
            // all instead of raising TranscriptNotOverlapping.
            let g = proj
                .genomic
                .as_ref()
                .expect("boundary-spanning poly-A projection must carry a genomic axis")
                .to_string();
            // Pin the exact normalized interval, not just the accession: the start
            // (g.1011) is inside the stripped terminal exon and the end (g.1012) is
            // the first poly-A walk coordinate, and the plus-strand edit keeps the
            // literal `GG` insertion. A wrong endpoint or a non-reverse-complemented
            // edit would fail here, where `starts_with("NC_000001.11")` would not.
            assert!(
                g.contains(":g.1011_1012delinsGG"),
                "expected plus-strand boundary interval g.1011_1012delinsGG, got: {g}"
            );
            // Coding axis: the input's own c.*3_*4delinsGG form is preserved.
            let c = proj
                .coding
                .as_ref()
                .expect("boundary-spanning poly-A projection must carry a coding axis");
            assert!(
                c.to_string().contains("c.*3_*4delinsGG"),
                "expected the input c.*3_*4delinsGG coding axis preserved, got: {c}"
            );
        }

        /// Minus-strand analogue of the boundary-spanning case. On the minus
        /// `NM_MPOLYA.1`, `c.*3` maps to genome 1000 (inside the exon) and `c.*4`
        /// to genome 999 (poly-A tail, one base downstream on the minus strand),
        /// so the genomic interval spans g.999_1000.
        #[test]
        fn project_variant_minus_polya_boundary_spanning_short_circuits() {
            let (projector, provider) = make_minus_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_MPOLYA.1"),
                gene_symbol: Some("MPOLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::new(
                        CdsPos {
                            base: 3,
                            offset: None,
                            utr3: true,
                            special: None,
                        },
                        CdsPos {
                            base: 4,
                            offset: None,
                            utr3: true,
                            special: None,
                        },
                    ),
                    NaEdit::Delins {
                        sequence: InsertedSequence::Literal("GG".parse().unwrap()),
                        deleted: None,
                        deleted_length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let proj = vp
                .project_variant(&HgvsVariant::Cds(cds), "NM_MPOLYA.1")
                .expect("minus-strand boundary-spanning c.*3_*4delinsGG must project");
            // As in the plus-strand case, assert the axis is present on the NC
            // contig rather than an exact coordinate (genome-frame renormalization
            // may shuffle it); the load-bearing check is that it short-circuits.
            let g = proj
                .genomic
                .as_ref()
                .expect("boundary-spanning poly-A projection must carry a genomic axis")
                .to_string();
            // Pin the exact normalized interval: on the minus strand `c.*3` maps to
            // the higher genome coordinate (g.1000, inside the exon) and `c.*4` to
            // the lower poly-A walk coordinate (g.999), so the genome-frame interval
            // is g.999_1000, and the inserted `GG` is reverse-complemented to `CC`.
            // Pinning the full interval+edit catches a wrong endpoint or a missing
            // reverse-complement that `starts_with("NC_000001.11")` would miss.
            assert!(
                g.contains(":g.999_1000delinsCC"),
                "expected minus-strand boundary interval g.999_1000delinsCC, got: {g}"
            );
            let c = proj
                .coding
                .as_ref()
                .expect("boundary-spanning poly-A projection must carry a coding axis");
            assert!(
                c.to_string().contains("c.*3_*4delinsGG"),
                "expected the input c.*3_*4delinsGG coding axis preserved, got: {c}"
            );
        }

        /// #797 (NG/LRG fan-out): an NG_-parented poly-A `c.*` input routed
        /// through `project_variant_all` must NOT be de-anchored to a synthetic
        /// (off-exon) `NC_` `Genome`. The poly-A walk produces a coordinate past
        /// the cdot exon span by construction; de-anchoring to it would seed the
        /// fan-out stab off the transcript and lose enumeration. The de-anchor
        /// step instead keeps the transcript-coordinate (`c.`) form for a confirmed
        /// poly-A endpoint, so the `c.` seed path anchors to the terminal exon,
        /// enumerates the transcript, and reaches the multi-axis short-circuit —
        /// framing the coding axis under the NG_ parent.
        #[test]
        fn project_variant_all_ng_parented_polya_enumerates_under_parent() {
            let (projector, mut provider) = make_polya_test_provider_and_projector();
            // Place NG_POLYA.1 onto the same NC_000001.11 the poly-A fixture uses:
            // NG base 1 == nc 1000, covering the 12-base exon core [1000,1012).
            provider.add_genomic_placement(
                "NG_POLYA.1",
                crate::reference::GenomicPlacement {
                    nc: Accession::new("NC", "000001", Some(11)),
                    parent_start: 1,
                    nc_start: 1000,
                    nc_end: 1012,
                    strand: crate::reference::Strand::Plus,
                },
            );
            let vp = VariantProjector::new(projector, provider);

            // NG_POLYA.1(NM_POLYA.1):c.*4del — a poly-A c.* endpoint with an NG_
            // genomic_context.
            let cds = CdsVariant {
                accession: parse_accession("NM_POLYA.1").with_genomic_context(Accession::new(
                    "NG",
                    "POLYA",
                    Some(1),
                )),
                gene_symbol: Some("POLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos {
                        base: 4,
                        offset: None,
                        utr3: true,
                        special: None,
                    }),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let results = vp
                .project_variant_all(&HgvsVariant::Cds(cds))
                .expect("NG-parented poly-A c.*4del must enumerate, not de-anchor off-exon");
            let proj = results
                .iter()
                .find(|p| p.transcript_id == "NM_POLYA.1")
                .expect("NM_POLYA.1 must be enumerated for the NG-parented poly-A input");
            // Coding axis preserved and framed under the NG_ parent.
            let c = proj
                .coding
                .as_ref()
                .expect("NG-parented poly-A projection must carry a coding axis")
                .to_string();
            assert!(
                c.starts_with("NG_POLYA.1(") && c.contains("c.*4del"),
                "expected NG_-framed c.*4del coding axis, got: {c}"
            );
        }

        /// #797 (fan-out): a `c.*` poly-A input routed through the
        /// `project_variant_all` public API must enumerate its
        /// terminal-exon-overlapping transcript rather than failing to seed
        /// fan-out. The seed-position extraction (`extract_contig_and_pos`)
        /// previously hard-failed on the stripped poly-A tail (`cds_to_genome`
        /// declines); it now seeds the stab query from the 3'-terminal exon anchor
        /// when `try_extend_polya_to_genome` confirms a genuine poly-A endpoint, so
        /// the per-transcript projection reaches the multi-axis short-circuit.
        #[test]
        fn project_variant_all_polya_region_enumerates_transcript() {
            let (projector, provider) = make_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_POLYA.1"),
                gene_symbol: Some("POLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos {
                        base: 4,
                        offset: None,
                        utr3: true,
                        special: None,
                    }),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let projections = vp
                .project_variant_all(&HgvsVariant::Cds(cds))
                .expect("poly-A c.*4del must enumerate through project_variant_all");
            assert!(
                !projections.is_empty(),
                "project_variant_all must enumerate the poly-A transcript, got none"
            );
            let proj = projections
                .iter()
                .find(|p| p.transcript_id == "NM_POLYA.1")
                .expect("NM_POLYA.1 must be among the enumerated transcripts");
            let g = proj
                .genomic
                .as_ref()
                .expect("enumerated poly-A projection must carry a genomic axis");
            assert!(
                g.to_string().contains(":g.1012del"),
                "expected g.1012del genomic axis from fan-out, got: {g}"
            );
            let c = proj
                .coding
                .as_ref()
                .expect("enumerated poly-A projection must carry a coding axis");
            assert!(
                c.to_string().contains("c.*4del"),
                "expected the input c.*4del coding axis preserved, got: {c}"
            );
        }

        /// Minus-strand analogue of `make_polya_test_provider_and_projector`:
        /// the 3'-terminal exon is the one with the smallest genome_start, and
        /// the poly-A walk descends in genome from `genome_start - 1`.
        fn make_minus_polya_test_provider_and_projector() -> (Projector, MockProvider) {
            let mut cdot = CdotMapper::new();
            cdot.add_transcript(
                "NM_MPOLYA.1".to_string(),
                CdotTranscript {
                    gene_name: Some("MPOLYAGENE".to_string()),
                    contig: "chr1".to_string(),
                    strand: ProvStrand::Minus,
                    // single 3'-terminal exon, genome [1000,1012), tx 0..12.
                    exons: vec![[1000, 1012, 0, 12]],
                    cds_start: Some(3),
                    cds_end: Some(9),
                    gene_id: None,
                    protein: Some("NP_MPOLYA.1".to_string()),
                    exon_cigars: Vec::new(),
                },
            );
            let projector = Projector::new(cdot);
            let mut provider = MockProvider::new();
            // 15-base transcript incl. 3-base poly-A tail (tx 12..15).
            provider.add_transcript(Transcript {
                id: "NM_MPOLYA.1".to_string(),
                gene_symbol: Some("MPOLYAGENE".to_string()),
                strand: TxStrand::Minus,
                sequence: Some("TTTATGCGCGTAAAA".to_string()), // 15 bases
                cds_start: Some(4),
                cds_end: Some(9),
                exons: vec![Exon::with_genomic(1, 1, 12, 1000, 1011)],
                chromosome: Some("chr1".to_string()),
                genomic_start: Some(1000),
                genomic_end: Some(1011),
                genome_build: Default::default(),
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
                protein_id: Some("NP_MPOLYA.1".to_string()),
            });
            // Same non-repeating backbone as the plus fixture, registered under
            // both `chr1` (cdot contig) and the `NC_000001.11` output accession so
            // genomic-frame renormalization can fetch bases around the minus-strand
            // walk targets (g.999 / g.998) without falling back (#797 review).
            let genomic = polya_genomic_sequence();
            provider.add_genomic_sequence("chr1", genomic.clone());
            provider.add_genomic_sequence("NC_000001.11", genomic);
            (projector, provider)
        }

        /// Minus-strand poly-A walk: c.*4 = tx 12 → genome_start(1000) - 1 = 999.
        #[test]
        fn project_to_genomic_minus_polya_region_walks_downstream() {
            let (projector, provider) = make_minus_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_MPOLYA.1"),
                gene_symbol: Some("MPOLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos {
                        base: 4,
                        offset: None,
                        utr3: true,
                        special: None,
                    }),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let out = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect("minus-strand c.*4del (poly-A region) should project");
            let s = match out {
                HgvsVariant::Genome(ref g) => g.to_string(),
                _ => panic!("expected Genome variant"),
            };
            assert!(
                s.contains(":g.999del"),
                "expected g.999del for minus-strand c.*4del, got: {s}"
            );
        }

        /// Minus-strand poly-A walk with a downstream offset: c.*4+1 = tx 13 →
        /// genome 998 (one base further downstream on the minus strand).
        #[test]
        fn project_to_genomic_minus_polya_region_with_offset_folds_linearly() {
            let (projector, provider) = make_minus_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_MPOLYA.1"),
                gene_symbol: Some("MPOLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos {
                        base: 4,
                        offset: Some(1),
                        utr3: true,
                        special: None,
                    }),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let out = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect("minus-strand c.*4+1del should project");
            let s = match out {
                HgvsVariant::Genome(ref g) => g.to_string(),
                _ => panic!("expected Genome variant"),
            };
            assert!(
                s.contains(":g.998del"),
                "expected g.998del for minus-strand c.*4+1del, got: {s}"
            );
        }

        /// A `c.*` position at/past the TRUE transcript end (beyond the poly-A
        /// tail) declines cleanly — never walks into non-transcript genome.
        /// c.*7 = tx 15 = true_tx_length → decline.
        #[test]
        fn project_to_genomic_beyond_polya_tail_declines() {
            let (projector, provider) = make_polya_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);
            let cds = CdsVariant {
                accession: parse_accession("NM_POLYA.1"),
                gene_symbol: Some("POLYAGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos {
                        base: 7,
                        offset: None,
                        utr3: true,
                        special: None,
                    }),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let err = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect_err("c.*7del is past the transcript end; must decline");
            let msg = err.to_string();
            assert!(
                msg.contains("Cannot map tx position") || msg.contains("15"),
                "expected a clean tx-mapping decline for c.*7del, got: {msg}"
            );
        }

        #[test]
        fn project_to_genomic_5utr_negative_position() {
            let (projector, provider) = make_utr_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            // c.-1T>G on NM_UTR.1 → tx_pos = cds_start - 1 = 2 → genome 1002
            // (seq[2] = 'T' in "TTTATGCGCTAAGGG").
            let cds = CdsVariant {
                accession: parse_accession("NM_UTR.1"),
                gene_symbol: Some("UTRGENE".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(-1)),
                    NaEdit::Substitution {
                        reference: Base::T,
                        alternative: Base::G,
                    },
                ),
            };
            let cds = attach_genomic_context_cds(cds, nc_parent());
            let input = HgvsVariant::Cds(cds);

            let out = vp
                .project_to_genomic(&input)
                .expect("c.-1T>G should project");
            let g = match out {
                HgvsVariant::Genome(ref g) => g,
                _ => panic!("expected Genome variant"),
            };
            assert_eq!(g.accession.to_string(), "NC_000001.11");
            let s = g.to_string();
            assert!(
                s.contains(":g.1002T>G"),
                "expected NC_000001.11...:g.1002T>G for c.-1T>G, got: {}",
                s
            );
        }

        // -- 15: intronic offset on non-coding transcript → unsupported (#332) -

        // -- special positions (pter/qter/cen) in project_to_genomic ----------

        /// #537: a `pter` CDS marker projects to the parent reference's 5'-most
        /// genomic coordinate (`g.1`), not the old special-sentinel rejection.
        /// (PR #534 keeps the *c.-axis* refusal; this covers the g. axis.) The
        /// raw projection is `g.1del`; the normalizer's 3' rule rolls it
        /// downstream.
        #[test]
        fn project_to_genomic_special_pter_projects_to_parent_terminus() {
            let (projector, mut provider) = make_test_provider_and_projector();
            // pter/qter resolve against the parent reference's length.
            provider.add_genomic_sequence("NG_TEST.1", "A".repeat(40));
            let vp = VariantProjector::new(projector, provider);

            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1").with_genomic_context(ng_parent("TEST", 1)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::pter()),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let input = HgvsVariant::Cds(cds);
            let g = vp
                .project_to_genomic(&input)
                .expect("pter must project to the parent's 5' terminus");
            assert_eq!(g.to_string(), "NG_TEST.1:g.1del");
        }

        /// A *mixed* numeric/`qter` range (`c.1_qter`) is out of #537 scope —
        /// the numeric endpoint still needs a transcript mapping — so it remains
        /// `UnsupportedProjection`. (A pure `c.qter` resolves to the parent
        /// terminus; that is covered by the `tests/it/issue_537_*` integration
        /// tests, which have a parent length to resolve against.)
        #[test]
        fn project_to_genomic_special_qter_returns_unsupported() {
            let (projector, provider) = make_test_provider_and_projector();
            let vp = VariantProjector::new(projector, provider);

            let cds = CdsVariant {
                accession: parse_accession("NM_TEST.1").with_genomic_context(ng_parent("TEST", 1)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::new(CdsPos::new(1), CdsPos::qter()),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let input = HgvsVariant::Cds(cds);
            let err = vp
                .project_to_genomic(&input)
                .expect_err("qter position must not project to genomic");
            assert!(
                matches!(err, FerroError::UnsupportedProjection { .. }),
                "expected UnsupportedProjection for qter position, got: {:?}",
                err
            );
        }

        /// Pins the deferred-to-#332 limitation: intronic offsets on `n.`
        /// (non-coding) inputs return `UnsupportedProjection` with a reason
        /// that references #332. If/when #332 lands and the runner is wired
        /// up, this test should start failing (XPASS) and be promoted.
        #[test]
        fn project_to_genomic_intronic_offset_on_noncoding_unsupported() {
            let (projector, provider) = make_noncoding_test_data();
            let vp = VariantProjector::new(projector, provider);

            // n.5+3del on NR_TEST.1 — non-coding transcript with intronic offset.
            let tx = TxVariant {
                accession: parse_accession("NR_TEST.1"),
                gene_symbol: Some("TESTGENE".to_string()),
                loc_edit: LocEdit::new(
                    TxInterval::point(TxPos::with_offset(5, 3)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            let tx = attach_genomic_context_tx(tx, ng_parent("NRTEST", 1));
            let input = HgvsVariant::Tx(tx);

            let err = vp
                .project_to_genomic(&input)
                .expect_err("intronic n. should be unsupported pending #332");
            match err {
                FerroError::UnsupportedProjection { reason } => {
                    assert!(
                        reason.contains("#332"),
                        "expected reason to reference #332, got: {}",
                        reason
                    );
                }
                other => panic!("expected UnsupportedProjection, got: {:?}", other),
            }
        }
    }

    // -------------------------------------------------------------------------
    // Multi-build cdot build-hint propagation (closes #389 item 2)
    //
    // Both `project_to_genomic` (c./n./r. → g.) and `project_single_inner`
    // (g. → c./n./r.) MUST consult the genome build implied by the
    // input's NG/NC parent (resp. the g. variant's NC_* accession) when
    // looking up the transcript and converting coordinates. Pre-#389
    // those sites used `cdot.get_transcript` / `transcript_genome_span` /
    // `cds_to_genome` / `genome_to_cds` directly, which always returned
    // the primary build's view — silently mis-projecting GRCh37 inputs
    // against a GRCh38-primary cdot (and vice versa).
    // -------------------------------------------------------------------------
    mod multi_build_projection_tests {
        use super::*;
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::{CdsInterval, GenomeInterval};
        use crate::hgvs::location::{CdsPos, GenomePos};
        use crate::hgvs::variant::{Accession, CdsVariant, GenomeVariant, HgvsVariant, LocEdit};

        /// Minimal multi-build cdot fixture: one coding transcript whose
        /// GRCh37 view sits at genome [10000, 10010) on NC_000001.10 and
        /// whose GRCh38 view sits at genome [20000, 20010) on
        /// NC_000001.11. Different genomic coords across builds is what
        /// lets us tell which view the projector consulted.
        ///
        /// `start_codon = 0`, `stop_codon = 10` ⇒ the whole exon is CDS
        /// (no UTRs), so c.1..c.10 maps to tx_pos 0..9. The fixture is
        /// coding so `project_single_inner`'s `genome_to_cds_on_build`
        /// call lands on the CDS branch — non-coding would hit a separate
        /// "Cannot convert tx position to CDS" error unrelated to #389.
        fn multi_build_cdot_json() -> &'static str {
            r#"
            {
                "transcripts": {
                    "NM_MB_TEST.1": {
                        "gene_name": "MBTEST",
                        "genome_builds": {
                            "GRCh37": {
                                "contig": "NC_000001.10",
                                "strand": "+",
                                "exons": [[10000, 10010, 1, 0, 10, "M10"]]
                            },
                            "GRCh38": {
                                "contig": "NC_000001.11",
                                "strand": "+",
                                "exons": [[20000, 20010, 1, 0, 10, "M10"]]
                            }
                        },
                        "start_codon": 0,
                        "stop_codon": 10
                    }
                }
            }
            "#
        }

        /// Build a `VariantProjector` whose cdot mapper has the two-build
        /// fixture above. `primary` selects which build the primary
        /// transcripts table is populated from (the other build lives in
        /// `alt_build_transcripts`). The `MockProvider` carries
        /// NM_MB_TEST.1 so that downstream protein prediction inside
        /// `project_single_inner` can fetch the transcript without
        /// hitting `ReferenceNotFound`.
        fn make_multi_build_projector(primary: &str) -> VariantProjector<MockProvider> {
            use crate::reference::transcript::{
                Exon, GenomeBuild as RefGenomeBuild, ManeStatus, Strand as TxStrand,
                Transcript as RefTranscript,
            };
            use std::sync::OnceLock;

            let cdot =
                CdotMapper::from_reader_with_build(multi_build_cdot_json().as_bytes(), primary)
                    .expect("multi-build cdot JSON should parse");
            let projector = Projector::new(cdot);
            let mut provider = MockProvider::new();
            provider.add_transcript(RefTranscript {
                id: "NM_MB_TEST.1".to_string(),
                gene_symbol: Some("MBTEST".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("ATGCGCTAATC".to_string()),
                cds_start: Some(1),
                cds_end: Some(10),
                exons: vec![Exon::new(1, 1, 10)],
                chromosome: Some("NC_000001.11".to_string()),
                genomic_start: Some(20000),
                genomic_end: Some(20009),
                genome_build: RefGenomeBuild::GRCh38,
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                protein_id: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
            });
            VariantProjector::new(projector, provider)
        }

        fn nc_chr1(version: u32) -> Accession {
            Accession::new("NC", "000001", Some(version))
        }

        fn ng_parent(num: &str, version: u32) -> Accession {
            Accession::new("NG", num, Some(version))
        }

        /// Smoke test the JSON parses into a usable multi-build mapper —
        /// keeps the later test failures readable if the JSON ever drifts.
        #[test]
        fn fixture_loads_both_builds() {
            let cdot =
                CdotMapper::from_reader_with_build(multi_build_cdot_json().as_bytes(), "GRCh38")
                    .expect("multi-build JSON must parse");
            assert_eq!(cdot.primary_build(), Some("GRCh38"));
            assert_eq!(cdot.transcript_count(), 1);
            let tx38 = cdot
                .get_transcript_on_build("NM_MB_TEST.1", "GRCh38")
                .expect("GRCh38 view must be populated in primary transcripts");
            assert_eq!(tx38.contig, "NC_000001.11");
            // #742: genome bounds load in the HGVS 1-based convention (+1 vs the
            // raw cdot 0-based alt_start/alt_end).
            assert_eq!(tx38.exons, vec![[20001, 20011, 0, 10]]);
            assert_eq!(tx38.cds_start, Some(0));
            assert_eq!(tx38.cds_end, Some(10));
            let tx37 = cdot
                .get_transcript_on_build("NM_MB_TEST.1", "GRCh37")
                .expect("GRCh37 view must be populated in alt_build_transcripts");
            assert_eq!(tx37.contig, "NC_000001.10");
            assert_eq!(tx37.exons, vec![[10001, 10011, 0, 10]]);
        }

        /// c.5A>T against NC_000001.10(NM_MB_TEST.1) → expect GRCh37
        /// coordinates (g.10005; cdot's internal exon[0] is HGVS-numeric,
        /// so c.1 sits at g.10001 and c.5 at g.10005). Pre-fix, this
        /// would have used the primary (GRCh38) view and produced
        /// ~g.20005, despite the parent's GRCh37 version.
        #[test]
        fn project_to_genomic_uses_grch37_view_for_nc_v10_parent() {
            let vp = make_multi_build_projector("GRCh38");
            let cds = CdsVariant {
                accession: parse_accession("NM_MB_TEST.1").with_genomic_context(nc_chr1(10)),
                gene_symbol: Some("MBTEST".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(5)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let out = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect("GRCh37 parent should project against GRCh37 exon table");
            match out {
                HgvsVariant::Genome(g) => {
                    let s = g.to_string();
                    assert!(
                        s.contains(":g.10005A>T"),
                        "expected GRCh37 coord g.10005A>T, got: {}",
                        s
                    );
                }
                other => panic!("expected Genome, got: {:?}", other),
            }
        }

        /// c.5A>T against NC_000001.11(NM_MB_TEST.1) → expect GRCh38
        /// coordinates (g.20005). Sibling to the GRCh37 case; together
        /// they prove the parent's build version is what drives the
        /// cdot view selection.
        #[test]
        fn project_to_genomic_uses_grch38_view_for_nc_v11_parent() {
            let vp = make_multi_build_projector("GRCh38");
            let cds = CdsVariant {
                accession: parse_accession("NM_MB_TEST.1").with_genomic_context(nc_chr1(11)),
                gene_symbol: Some("MBTEST".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(5)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let out = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect("GRCh38 parent should project against GRCh38 exon table");
            match out {
                HgvsVariant::Genome(g) => {
                    let s = g.to_string();
                    assert!(
                        s.contains(":g.20005A>T"),
                        "expected GRCh38 coord g.20005A>T, got: {}",
                        s
                    );
                }
                other => panic!("expected Genome, got: {:?}", other),
            }
        }

        /// g. → c.: NC_000001.10:g.10005A>T projected onto NM_MB_TEST.1
        /// → expect c.5A>T. Pre-fix, `project_single_inner` used the
        /// primary (GRCh38) `transcript_genome_span` of [20000, 20010),
        /// so position 10005 fell outside and the call returned
        /// `TranscriptNotOverlapping`. Post-fix the GRCh37 span
        /// [10000, 10010) is consulted and the projection succeeds.
        #[test]
        fn project_normalized_uses_grch37_view_for_nc_v10_input() {
            let vp = make_multi_build_projector("GRCh38");
            let g = GenomeVariant {
                accession: nc_chr1(10),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(10005)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let projection = vp
                .project_normalized(&HgvsVariant::Genome(g), "NM_MB_TEST.1")
                .expect("GRCh37 g. input must project via the alt-build cdot view");
            let coding = projection
                .coding
                .as_ref()
                .expect("coding (c.) form of the projection must be present");
            let s = coding.to_string();
            assert!(
                s.contains(":c.5A>T"),
                "expected c.5A>T for NC_000001.10 g. input, got: {}",
                s
            );
        }

        /// NG_* parents are build-agnostic per the HGVS spec — the chromosome
        /// (`NC_`) pivot falls back to the primary cdot view rather than forcing
        /// an alt-build lookup. With primary=GRCh38 the c.5 position must come
        /// back at the GRCh38 exon's position (g.20005) regardless of any GRCh37
        /// alt-build data.
        ///
        /// This is asserted on the `project_to_genomic_nc` pivot, which is where
        /// the build fallback lives. The public `project_to_genomic` re-anchors
        /// its output and now *declines* an NG_/LRG_ parent with no known
        /// placement rather than stamping a chromosome coordinate under the NG_
        /// accession (#655 — see `project_to_genomic_ng_parent_without_placement_declines`).
        #[test]
        fn project_to_genomic_with_ng_parent_uses_primary_build() {
            let vp = make_multi_build_projector("GRCh38");
            let cds = CdsVariant {
                accession: parse_accession("NM_MB_TEST.1")
                    .with_genomic_context(ng_parent("MBTEST", 1)),
                gene_symbol: Some("MBTEST".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(5)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let out = vp
                .project_to_genomic_nc(&HgvsVariant::Cds(cds))
                .expect("NG-parented variant should pivot via the primary cdot view");
            match out {
                HgvsVariant::Genome(g) => {
                    let s = g.to_string();
                    assert!(
                        s.starts_with("NG_MBTEST.1"),
                        "expected pivot to inherit NG_ parent accession, got: {}",
                        s
                    );
                    assert!(
                        s.contains(":g.20005A>T"),
                        "NG (build-agnostic) parent must use primary (GRCh38) coords, got: {}",
                        s
                    );
                }
                other => panic!("expected Genome, got: {:?}", other),
            }
        }

        /// `transcript_cache` and `ref_protein_cache` must partition by
        /// genome-build identity, not collapse to a single key per
        /// `transcript_id` (#389 follow-up).
        ///
        /// Pre-fix: the `coding` variant built inside `project_single_inner`
        /// dropped the genomic accession (`parse_accession(transcript_id)`
        /// carries no `genomic_context`), so `parent_key_for(&coding)`
        /// returned `None` for both `NC_*.10` and `NC_*.11` inputs. The
        /// two builds therefore shared cache entries — a GRCh38-built
        /// `Transcript` / `RefProteinBundle` could be served for a GRCh37
        /// projection and vice versa.
        ///
        /// Post-fix: the parallel `cache_variant` stamps the projected g.
        /// accession (`NC_*.10` vs `NC_*.11`) as `genomic_context`, so
        /// `parent_key_for(&cache_variant)` produces distinct keys for
        /// the two builds and the caches keep them apart. The test
        /// projects an indel (deletion) on each build so both branches
        /// of the protein-prediction dispatch
        /// (`cached_get_transcript_for_variant` + `cached_ref_translation`)
        /// execute.
        #[test]
        fn caches_partition_by_genome_build_for_nc_inputs() {
            let vp = make_multi_build_projector("GRCh38");

            // c.5del on GRCh37 (g.10005del on NC_000001.10).
            let g37 = GenomeVariant {
                accession: nc_chr1(10),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(10005)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            vp.project_normalized(&HgvsVariant::Genome(g37), "NM_MB_TEST.1")
                .expect("GRCh37 deletion projection must succeed");

            // c.5del on GRCh38 (g.20005del on NC_000001.11).
            let g38 = GenomeVariant {
                accession: nc_chr1(11),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(20005)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            };
            vp.project_normalized(&HgvsVariant::Genome(g38), "NM_MB_TEST.1")
                .expect("GRCh38 deletion projection must succeed");

            let tx_keys: std::collections::HashSet<_> = vp
                .transcript_cache
                .read()
                .expect("transcript cache poisoned")
                .keys()
                .cloned()
                .collect();
            assert_eq!(
                tx_keys.len(),
                2,
                "transcript_cache must keep GRCh37 and GRCh38 entries separate, \
                 got keys: {:?}",
                tx_keys
            );
            assert!(
                tx_keys.contains(&("NM_MB_TEST.1".to_string(), Some("NC_000001.10".to_string()))),
                "expected an NC_000001.10 key, got: {:?}",
                tx_keys
            );
            assert!(
                tx_keys.contains(&("NM_MB_TEST.1".to_string(), Some("NC_000001.11".to_string()))),
                "expected an NC_000001.11 key, got: {:?}",
                tx_keys
            );

            let rp_keys: std::collections::HashSet<_> = vp
                .ref_protein_cache
                .read()
                .expect("ref-protein cache poisoned")
                .keys()
                .cloned()
                .collect();
            assert_eq!(
                rp_keys.len(),
                2,
                "ref_protein_cache must keep GRCh37 and GRCh38 entries separate, \
                 got keys: {:?}",
                rp_keys
            );
        }

        /// Same as above but for NC_000001.11 (GRCh38, the primary build).
        /// Confirms the primary-build path still works alongside the new
        /// alt-build dispatch.
        #[test]
        fn project_normalized_uses_grch38_view_for_nc_v11_input() {
            let vp = make_multi_build_projector("GRCh38");
            let g = GenomeVariant {
                accession: nc_chr1(11),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(20005)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let projection = vp
                .project_normalized(&HgvsVariant::Genome(g), "NM_MB_TEST.1")
                .expect("GRCh38 g. input must project via the primary cdot view");
            let coding = projection
                .coding
                .as_ref()
                .expect("coding (c.) form of the projection must be present");
            let s = coding.to_string();
            assert!(
                s.contains(":c.5A>T"),
                "expected c.5A>T for NC_000001.11 g. input, got: {}",
                s
            );
        }

        /// End-to-end: `project_variant_all` on a GRCh37 genome input must
        /// DISCOVER the overlapping transcript (no transcript id supplied) and
        /// project it. This is the corpus's actual path: the input is a bare
        /// `NC_*.10:g.…` HGVSg string and the projector must route overlap
        /// discovery to the GRCh37 (secondary) build via the contig, then
        /// build-route the downstream genome→cds projection to the GRCh37 exon
        /// table. Pre-fix, overlap discovery only indexed the primary (GRCh38)
        /// build, so `project_variant_all` returned an empty `Vec` ("got []").
        #[test]
        fn project_variant_all_discovers_grch37_transcript_by_overlap() {
            let vp = make_multi_build_projector("GRCh38");
            // NC_000001.10 is GRCh37; the transcript id is NOT supplied — it
            // must be discovered by overlap.
            let g = GenomeVariant {
                accession: nc_chr1(10),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(10005)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let projections = vp
                .project_variant_all(&HgvsVariant::Genome(g))
                .expect("project_variant_all must succeed for a GRCh37 g. input");

            assert!(
                !projections.is_empty(),
                "project_variant_all must discover the overlapping GRCh37 transcript \
                 (pre-fix this returned [])"
            );
            let proj = projections
                .iter()
                .find(|p| p.transcript_id == "NM_MB_TEST.1")
                .expect("NM_MB_TEST.1 must be among the discovered transcripts");
            let coding = proj
                .coding
                .as_ref()
                .expect("coding (c.) form must be present on the discovered projection");
            let s = coding.to_string();
            assert!(
                s.contains(":c.5A>T"),
                "expected c.5A>T from the GRCh37-discovered projection, got: {}",
                s
            );
        }

        /// Sibling to the above on the primary (GRCh38) build: overlap
        /// discovery on NC_000001.11 must keep working through the same path.
        #[test]
        fn project_variant_all_discovers_grch38_transcript_by_overlap() {
            let vp = make_multi_build_projector("GRCh38");
            let g = GenomeVariant {
                accession: nc_chr1(11),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(20005)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let projections = vp
                .project_variant_all(&HgvsVariant::Genome(g))
                .expect("project_variant_all must succeed for a GRCh38 g. input");
            let proj = projections
                .iter()
                .find(|p| p.transcript_id == "NM_MB_TEST.1")
                .expect("NM_MB_TEST.1 must be discovered on the primary build too");
            let s = proj
                .coding
                .as_ref()
                .expect("coding form present")
                .to_string();
            assert!(
                s.contains(":c.5A>T"),
                "expected c.5A>T from the GRCh38-discovered projection, got: {}",
                s
            );
        }

        // ---------------------------------------------------------------------
        // Predicted-RNA allele build-scoping (#843)
        //
        // The predicted `r.` allele in `project_allele_inner` fetches the
        // transcript for `predict_rna` via the input. For a bare `g.` allele the
        // member `NC_*` accession is build-bearing but carries no
        // `genomic_context`, so before #843 the fetch resolved the primary-build
        // transcript regardless of build — and `predict_rna`'s deletion 3'-shift
        // (which walks the transcript *sequence*) then ran against the wrong
        // build's bases, emitting a wrong-build `r.`. These tests model a
        // transcript whose GRCh37 vs GRCh38 bases differ at a deletion-shift
        // boundary (the real divergence a genome-reconstructed transcript shows;
        // FASTA-backed transcripts serve build-identical bases) and assert the
        // allele `r.` is now build-correct.
        // ---------------------------------------------------------------------

        /// A 30 nt single-exon coding fixture whose GRCh37 and GRCh38 transcript
        /// *bases* differ only in the length of a poly-A run starting at c.5:
        /// GRCh37 runs `A` over c.5..7 (blocked by `C` at c.8); GRCh38 runs `A`
        /// over c.5..9 (blocked by `C` at c.10). A `c.5del` (deleting an `A`)
        /// therefore 3'-shifts to `r.(7del)` on GRCh37 but `r.(9del)` on GRCh38
        /// — a clean, build-dependent predicted-RNA output. The genome contigs
        /// differ per build (`NC_000001.10` [10000,10030) GRCh37 /
        /// `NC_000001.11` [20000,20030) GRCh38), so a `g.` input's build is
        /// inferred from its `NC_*` version.
        fn make_rna_divergence_projector() -> VariantProjector<MockProvider> {
            make_rna_divergence_projector_inner(true)
        }

        /// Variant of [`make_rna_divergence_projector`] that omits the GRCh37
        /// build-keyed transcript record, leaving only the GRCh38 build-keyed
        /// record plus the build-agnostic primary. A GRCh37 `g.` input therefore
        /// *misses* the build-keyed lookup (the mock surfaces `ReferenceNotFound`
        /// for the partial setup), modeling a transcript whose sequence is present
        /// only in the primary build's data.
        fn make_rna_divergence_projector_grch38_keyed_only() -> VariantProjector<MockProvider> {
            make_rna_divergence_projector_inner(false)
        }

        fn make_rna_divergence_projector_inner(
            register_grch37_keyed: bool,
        ) -> VariantProjector<MockProvider> {
            use crate::reference::transcript::{
                Exon, GenomeBuild as RefGenomeBuild, ManeStatus, Strand as TxStrand,
                Transcript as RefTranscript,
            };
            use std::sync::OnceLock;

            let cdot_json = r#"
            {
                "transcripts": {
                    "NM_RNA_DIV.1": {
                        "gene_name": "RNADIV",
                        "genome_builds": {
                            "GRCh37": {
                                "contig": "NC_000001.10",
                                "strand": "+",
                                "exons": [[10000, 10030, 1, 0, 30, "M30"]]
                            },
                            "GRCh38": {
                                "contig": "NC_000001.11",
                                "strand": "+",
                                "exons": [[20000, 20030, 1, 0, 30, "M30"]]
                            }
                        },
                        "start_codon": 0,
                        "stop_codon": 30
                    }
                }
            }
            "#;
            // c.:  1 2 3 4 5 6 7 8 9 ...                      (1-based)
            // 37:  G G C G A A A C T ...   poly-A c.5..7, C at c.8
            // 38:  G G C G A A A A A C ... poly-A c.5..9, C at c.10
            // The deletion projects to c.5 (a poly-A base on BOTH builds), so the
            // g.→c. step is build-identical; only the run length downstream — and
            // thus the predicted `r.` 3'-shift target — differs: GRCh37 shifts to
            // c.7 (`r.(7del)`), GRCh38 to c.9 (`r.(9del)`).
            let seq37 = "GGCGAAACTTTTTTTTTTTTTTTTTTTTTT"; // len 30
            let seq38full = "GGCGAAAAACTTTTTTTTTTTTTTTTTTTTTT"; // len 32 -> trim to 30
            let seq38 = &seq38full[..30];
            assert_eq!(seq37.len(), 30, "seq37 must be 30 nt");
            assert_eq!(seq38.len(), 30, "seq38 must be 30 nt");

            let cdot = CdotMapper::from_reader_with_build(cdot_json.as_bytes(), "GRCh38")
                .expect("rna-divergence cdot JSON should parse");
            let projector = Projector::new(cdot);
            let mut provider = MockProvider::new();
            let mk = |build: RefGenomeBuild, seq: &str, chrom: &str, gstart: u64| RefTranscript {
                id: "NM_RNA_DIV.1".to_string(),
                gene_symbol: Some("RNADIV".to_string()),
                strand: TxStrand::Plus,
                sequence: Some(seq.to_string()),
                cds_start: Some(1),
                cds_end: Some(30),
                exons: vec![Exon::new(1, 1, 30)],
                chromosome: Some(chrom.to_string()),
                genomic_start: Some(gstart),
                genomic_end: Some(gstart + 29),
                genome_build: build,
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                protein_id: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
            };
            // Build-keyed transcripts (the #843 divergence vector). The GRCh37
            // record is optional so a partial setup (GRCh38-keyed + primary only)
            // can model a build-keyed miss for a GRCh37 input.
            if register_grch37_keyed {
                provider.add_transcript_on_build(
                    "GRCh37",
                    mk(RefGenomeBuild::GRCh37, seq37, "NC_000001.10", 10000),
                );
            }
            provider.add_transcript_on_build(
                "GRCh38",
                mk(RefGenomeBuild::GRCh38, seq38, "NC_000001.11", 20000),
            );
            // Build-agnostic fallback = the primary (GRCh38) bases, matching how a
            // real provider serves the primary transcript when no build is keyed.
            // This is exactly what the pre-#843 build-agnostic fetch returned for
            // BOTH builds, so a GRCh37 input wrongly read the GRCh38 run length.
            provider.add_transcript(mk(RefGenomeBuild::GRCh38, seq38, "NC_000001.11", 20000));
            VariantProjector::new(projector, provider)
        }

        /// Build a single-deletion `g.` allele at `g.<pos>del` on `NC_000001.<ver>`.
        fn g_del_allele(version: u32, pos: u64) -> HgvsVariant {
            let member = HgvsVariant::Genome(GenomeVariant {
                accession: nc_chr1(version),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    GenomeInterval::point(GenomePos::new(pos)),
                    NaEdit::Deletion {
                        sequence: None,
                        length: None,
                    },
                ),
            });
            HgvsVariant::Allele(AlleleVariant::new(vec![member], AllelePhase::Cis))
        }

        /// The allele predicted `r.` must read the *input's build* transcript, so
        /// a `g.` allele on GRCh37 (`NC_000001.10`) and one on GRCh38
        /// (`NC_000001.11`) — same `c.6del` after projection — yield
        /// build-distinct 3'-shifted `r.` outputs. Pre-#843 both read the primary
        /// (GRCh38) transcript, so the GRCh37 allele emitted the GRCh38 shift
        /// (`r.(10del)`): a wrong-build `r.`. Post-fix the GRCh37 allele emits its
        /// own `r.(8del)`.
        #[test]
        fn allele_predicted_rna_is_build_scoped() {
            let vp = make_rna_divergence_projector();

            // g.10005del on NC_000001.10 (GRCh37) -> c.5del.
            let g37 = g_del_allele(10, 10005);
            let r37 = vp
                .project_normalized(&g37, "NM_RNA_DIV.1")
                .expect("GRCh37 g. allele must project")
                .rna
                .expect("GRCh37 allele must predict an r.")
                .to_string();

            // g.20005del on NC_000001.11 (GRCh38) -> c.6del.
            let g38 = g_del_allele(11, 20005);
            let r38 = vp
                .project_normalized(&g38, "NM_RNA_DIV.1")
                .expect("GRCh38 g. allele must project")
                .rna
                .expect("GRCh38 allele must predict an r.")
                .to_string();

            // The two builds 3'-shift the deletion to different positions because
            // their poly-A run lengths differ. Pre-fix r37 == r38 (both GRCh38).
            // Assert the FULL rendered string (exact, not a `contains` substring
            // match that would also accept e.g. `17del`/`19del`).
            assert_eq!(
                r37, "NM_RNA_DIV.1:r.7del",
                "GRCh37 allele r. must 3'-shift to its own run end (r.7del)"
            );
            assert_eq!(
                r38, "NM_RNA_DIV.1:r.9del",
                "GRCh38 allele r. must 3'-shift to its own run end (r.9del)"
            );
            assert_ne!(
                r37, r38,
                "GRCh37 and GRCh38 alleles must predict build-distinct r.; \
                 equal output means predict_rna read the same (wrong-build) transcript"
            );
        }

        /// The allele predicted-`r.` transcript fetch must partition the
        /// transcript cache by build (mirrors
        /// `caches_partition_by_genome_build_for_nc_inputs` for the single path).
        /// Pre-#843 the bare `g.` allele fetched under `(id, None)` for both
        /// builds (one colliding key); post-fix each build stamps its `NC_*`
        /// accession, yielding two build-distinct keys.
        #[test]
        fn allele_predicted_rna_does_not_pollute_cache_with_parentless_key() {
            // The build-distinct keys `(id, Some(NC_*.10))` / `(id, Some(NC_*.11))`
            // are inserted by each member's genome-pivot projection regardless of
            // the #843 bug, so their presence is NOT a fix signal. The bug's
            // fingerprint is the PARENTLESS `(id, None)` key the buggy allele-site
            // fetch (keyed on the unparented `g.` allele) inserted in ADDITION:
            // post-fix `predict_rna_transcript` stamps the build-bearing accession
            // first, so no `(id, None)` entry is ever created for the allele path.
            let vp = make_rna_divergence_projector();
            vp.project_normalized(&g_del_allele(10, 10005), "NM_RNA_DIV.1")
                .expect("GRCh37 allele projection");
            vp.project_normalized(&g_del_allele(11, 20005), "NM_RNA_DIV.1")
                .expect("GRCh38 allele projection");

            let tx_keys: std::collections::HashSet<_> = vp
                .transcript_cache
                .read()
                .expect("transcript cache poisoned")
                .keys()
                .cloned()
                .collect();
            // Sanity: the build-scoped entries exist (from the member pivots).
            assert!(
                tx_keys.contains(&("NM_RNA_DIV.1".to_string(), Some("NC_000001.10".to_string())))
                    && tx_keys
                        .contains(&("NM_RNA_DIV.1".to_string(), Some("NC_000001.11".to_string()))),
                "build-scoped member-pivot keys must be present, got: {tx_keys:?}"
            );
            // The actual fix signal: no parentless entry for the allele transcript.
            assert!(
                !tx_keys.contains(&("NM_RNA_DIV.1".to_string(), None)),
                "predict_rna must not fetch the allele transcript build-agnostically;                  a parentless (id, None) key means it read the primary build, got: {tx_keys:?}"
            );
        }

        /// A build-keyed miss under a stamped build context must NOT fall back to
        /// the build-agnostic primary transcript (#843 / CodeRabbit). With only
        /// the GRCh38 transcript build-keyed (plus the primary), a GRCh37 `g.`
        /// allele's build-stamped fetch misses; `predict_rna_transcript` must
        /// surface that miss instead of silently re-serving the GRCh38 primary —
        /// which would emit the wrong-build `r.9del`. The c./g. forms stay
        /// build-identical (they read cdot, not transcript bases), so the
        /// projection still succeeds; only the predicted `r.` is withheld.
        /// Pre-fix the `ReferenceNotFound` fallback masked the miss (rna = r.9del).
        #[test]
        fn allele_predicted_rna_does_not_fall_back_to_primary_on_build_miss() {
            let vp = make_rna_divergence_projector_grch38_keyed_only();
            let g37 = g_del_allele(10, 10005);
            let proj = vp
                .project_normalized(&g37, "NM_RNA_DIV.1")
                .expect("GRCh37 g. allele must still project (c./g. are build-identical)");
            // The coding form is build-identical, so it must be present — this
            // confirms the projection reached the predicted-`r.` step rather than
            // bailing earlier (which would make `rna.is_none()` vacuous).
            assert!(
                proj.coding.is_some(),
                "the build-identical c. form must project even with the GRCh37 \
                 transcript absent"
            );
            assert!(
                proj.rna.is_none(),
                "a build-keyed miss must not fall back to the primary transcript; \
                 a predicted r. here ({:?}) means the wrong-build (GRCh38) bases \
                 were served",
                proj.rna.as_ref().map(|r| r.to_string())
            );
        }

        /// The build-match guard must decline *only* on a wrong build, never
        /// over-decline the matching one. On the same partial fixture (GRCh38
        /// build-keyed + primary), a GRCh38 `g.` allele resolves its own
        /// build-keyed transcript, so the guard passes and the predicted `r.` is
        /// emitted as usual. This locks the guard's precision: a future "just
        /// suppress the fallback" regression that broke the build-identical axes
        /// (or dropped the matching-build `r.`) would fail here.
        #[test]
        fn grch38_predicted_rna_survives_build_guard_on_partial_fixture() {
            let vp = make_rna_divergence_projector_grch38_keyed_only();
            let g38 = g_del_allele(11, 20005);
            let r38 = vp
                .project_normalized(&g38, "NM_RNA_DIV.1")
                .expect("GRCh38 g. allele must project")
                .rna
                .expect("GRCh38 allele must still predict an r. (build-keyed record present)")
                .to_string();
            assert_eq!(
                r38, "NM_RNA_DIV.1:r.9del",
                "GRCh38 allele must 3'-shift to its own run end (r.9del); the build guard \
                 must not decline a matching-build transcript"
            );
        }

        // ---------------------------------------------------------------------
        // `--assembly` override for build-agnostic NG_ inputs (#715)
        //
        // Builds on the two-build cdot fixture above by giving NM_MB_TEST.1's
        // NG_ parent (NG_000999.1) per-build chromosomal placements with
        // *different* parent offsets, so the re-anchored NG_ coordinate reveals
        // which build's placement was selected: GRCh38 → g.105, GRCh37 → g.205.
        // This is the end-to-end GRCh37 NG_ path #713 deferred (only reachable
        // once the override supplies a build for a bare NG_).
        // ---------------------------------------------------------------------

        const NG_PARENT: &str = "NG_000999.1";

        /// Two-build projector whose NG_000999.1 has a GRCh38 placement
        /// (NC_000001.11, parent_start 100) and — when `with_grch37` — a GRCh37
        /// placement (NC_000001.10, parent_start 200).
        fn make_ng_assembly_projector(with_grch37: bool) -> VariantProjector<MockProvider> {
            use crate::reference::transcript::{
                Exon, GenomeBuild as RefGenomeBuild, ManeStatus, Strand as TxStrand,
                Transcript as RefTranscript,
            };
            use std::sync::OnceLock;

            let cdot =
                CdotMapper::from_reader_with_build(multi_build_cdot_json().as_bytes(), "GRCh38")
                    .expect("multi-build cdot JSON should parse");
            let projector = Projector::new(cdot);
            let mut provider = MockProvider::new();
            provider.add_transcript(RefTranscript {
                id: "NM_MB_TEST.1".to_string(),
                gene_symbol: Some("MBTEST".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("ATGCGCTAATC".to_string()),
                cds_start: Some(1),
                cds_end: Some(10),
                exons: vec![Exon::new(1, 1, 10)],
                chromosome: Some("NC_000001.11".to_string()),
                genomic_start: Some(20000),
                genomic_end: Some(20009),
                genome_build: RefGenomeBuild::GRCh38,
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                protein_id: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
            });
            provider.add_genomic_placement(
                NG_PARENT,
                crate::reference::provider::GenomicPlacement {
                    nc: parse_accession("NC_000001.11"),
                    parent_start: 100,
                    nc_start: 20000,
                    nc_end: 20010,
                    strand: crate::reference::transcript::Strand::Plus,
                },
            );
            if with_grch37 {
                provider.add_genomic_placement(
                    NG_PARENT,
                    crate::reference::provider::GenomicPlacement {
                        nc: parse_accession("NC_000001.10"),
                        parent_start: 200,
                        nc_start: 10000,
                        nc_end: 10010,
                        strand: crate::reference::transcript::Strand::Plus,
                    },
                );
            }
            VariantProjector::new(projector, provider)
        }

        /// c.5A>T on NM_MB_TEST.1 with a *bare* NG_000999.1 parent (no build).
        fn ng_cds_input() -> HgvsVariant {
            HgvsVariant::Cds(CdsVariant {
                accession: parse_accession("NM_MB_TEST.1")
                    .with_genomic_context(parse_accession(NG_PARENT)),
                gene_symbol: Some("MBTEST".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(5)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            })
        }

        /// `--assembly GRCh37` makes a bare NG_ select its GRCh37 placement and
        /// re-anchor against the GRCh37 chromosome view (g.205, not g.105).
        #[test]
        fn assembly_override_selects_grch37_placement_for_bare_ng() {
            let vp = make_ng_assembly_projector(true).with_assembly(Some("GRCh37"));
            let out = vp
                .project_to_genomic(&ng_cds_input())
                .expect("GRCh37 override should re-anchor via the GRCh37 placement");
            let s = out.to_string();
            assert!(
                s.contains(":g.205A>T"),
                "expected the GRCh37 placement (parent_start 200 -> g.205), got: {s}"
            );
        }

        /// Without the override a bare NG_ keeps the GRCh38-preferred placement
        /// (g.105) — adding a GRCh37 placement must not shift the default (#713).
        #[test]
        fn no_assembly_override_keeps_grch38_placement_for_bare_ng() {
            let vp = make_ng_assembly_projector(true);
            let out = vp
                .project_to_genomic(&ng_cds_input())
                .expect("default should re-anchor via the GRCh38 placement");
            let s = out.to_string();
            assert!(
                s.contains(":g.105A>T"),
                "expected the GRCh38 placement (parent_start 100 -> g.105), got: {s}"
            );
        }

        /// `--assembly GRCh37` with no GRCh37 placement present declines rather
        /// than mis-anchoring onto the GRCh38 placement (#655/#702 guard).
        #[test]
        fn assembly_override_declines_when_grch37_placement_absent() {
            let vp = make_ng_assembly_projector(false).with_assembly(Some("GRCh37"));
            let err = vp
                .project_to_genomic(&ng_cds_input())
                .expect_err("must decline when the requested build has no placement");
            assert!(
                matches!(err, FerroError::UnsupportedProjection { .. }),
                "expected an UnsupportedProjection decline, got: {err:?}"
            );
        }

        /// A build the input accession already encodes (NC_000001.11) wins over a
        /// contradictory `--assembly GRCh37`: the override only fills in for
        /// build-agnostic inputs, never reinterprets an explicit accession (#715).
        #[test]
        fn input_encoded_build_wins_over_assembly_override() {
            let vp = make_ng_assembly_projector(true).with_assembly(Some("GRCh37"));
            let cds = CdsVariant {
                accession: parse_accession("NM_MB_TEST.1").with_genomic_context(nc_chr1(11)),
                gene_symbol: Some("MBTEST".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(5)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            };
            let out = vp
                .project_to_genomic(&HgvsVariant::Cds(cds))
                .expect("explicit NC_*.11 parent projects on GRCh38 despite the override");
            let s = out.to_string();
            assert!(
                s.starts_with("NC_000001.11") && s.contains(":g.20005A>T"),
                "input-encoded GRCh38 (g.20005) must win over --assembly GRCh37 (g.10005), got: {s}"
            );
        }

        /// Regression (#716): the downstream cdot lookup in `project_single_inner`
        /// must reuse the build that `project_to_genomic_nc` used for the pivot,
        /// not re-derive it from the projected `NC_` accession.
        ///
        /// Fixture: a bare `NG_` parent whose GRCh37 placement re-anchors onto a
        /// contig (`NC_000001.99`) whose version the hardcoded heuristic does NOT
        /// recognize — exactly the case the data-driven inference targets, and the
        /// case where re-deriving from the accession yields `None`. Under
        /// `--assembly GRCh37` the pivot correctly selects the GRCh37 alignment,
        /// but pre-fix the protein/cdot lookup re-inferred `None` from
        /// `NC_000001.99`, fell back to cdot's primary (GRCh38) view, and the
        /// GRCh37-frame genomic position fell outside the GRCh38 exon →
        /// `TranscriptNotOverlapping`. With the build hint preserved, the full
        /// projection succeeds on GRCh37.
        #[test]
        fn assembly_override_survives_downstream_cdot_lookup_for_unrecognized_contig() {
            use crate::reference::transcript::{
                Exon, GenomeBuild as RefGenomeBuild, ManeStatus, Strand as TxStrand,
                Transcript as RefTranscript,
            };
            use std::sync::OnceLock;

            // GRCh37 alignment lives on an accession the version heuristic cannot
            // classify; GRCh38 on the canonical NC_000001.11.
            const UR_PARENT: &str = "NG_000777.1";
            let cdot_json = r#"
            {
                "transcripts": {
                    "NM_UR_TEST.1": {
                        "gene_name": "URTEST",
                        "genome_builds": {
                            "GRCh37": {
                                "contig": "NC_000001.99",
                                "strand": "+",
                                "exons": [[10000, 10010, 1, 0, 10, "M10"]]
                            },
                            "GRCh38": {
                                "contig": "NC_000001.11",
                                "strand": "+",
                                "exons": [[20000, 20010, 1, 0, 10, "M10"]]
                            }
                        },
                        "start_codon": 0,
                        "stop_codon": 10
                    }
                }
            }
            "#;

            // Sanity-check the premise: the heuristic genuinely declines on the
            // re-anchored GRCh37 contig, so a re-derive-from-accession path would
            // lose the build and this test would distinguish the two.
            assert_eq!(
                crate::liftover::aliases::infer_genome_build_from_accession(&parse_accession(
                    "NC_000001.99"
                )),
                None,
                "premise: NC_000001.99 is unrecognized by the version heuristic",
            );

            let cdot = CdotMapper::from_reader_with_build(cdot_json.as_bytes(), "GRCh38")
                .expect("multi-build cdot JSON should parse");
            let projector = Projector::new(cdot);
            let mut provider = MockProvider::new();
            provider.add_transcript(RefTranscript {
                id: "NM_UR_TEST.1".to_string(),
                gene_symbol: Some("URTEST".to_string()),
                strand: TxStrand::Plus,
                sequence: Some("ATGCGCTAATC".to_string()),
                cds_start: Some(1),
                cds_end: Some(10),
                exons: vec![Exon::new(1, 1, 10)],
                chromosome: Some("NC_000001.11".to_string()),
                genomic_start: Some(20000),
                genomic_end: Some(20009),
                genome_build: RefGenomeBuild::GRCh38,
                mane_status: ManeStatus::default(),
                refseq_match: None,
                ensembl_match: None,
                protein_id: None,
                exon_cigars: Vec::new(),
                cached_introns: OnceLock::new(),
            });
            provider.add_genomic_placement(
                UR_PARENT,
                crate::reference::provider::GenomicPlacement {
                    nc: parse_accession("NC_000001.11"),
                    parent_start: 100,
                    nc_start: 20000,
                    nc_end: 20010,
                    strand: crate::reference::transcript::Strand::Plus,
                },
            );
            provider.add_genomic_placement(
                UR_PARENT,
                crate::reference::provider::GenomicPlacement {
                    nc: parse_accession("NC_000001.99"),
                    parent_start: 200,
                    nc_start: 10000,
                    nc_end: 10010,
                    strand: crate::reference::transcript::Strand::Plus,
                },
            );
            let vp = VariantProjector::new(projector, provider).with_assembly(Some("GRCh37"));

            // c.5A>T on NM_UR_TEST.1 with a *bare* NG_ parent (build-agnostic).
            let input = HgvsVariant::Cds(CdsVariant {
                accession: parse_accession("NM_UR_TEST.1")
                    .with_genomic_context(parse_accession(UR_PARENT)),
                gene_symbol: Some("URTEST".to_string()),
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(5)),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::T,
                    },
                ),
            });

            let proj = vp.project_variant(&input, "NM_UR_TEST.1").expect(
                "GRCh37 override must drive the downstream cdot lookup; pre-fix this errored \
                 because the lookup re-inferred None from the unrecognized NC_000001.99 contig",
            );
            let coding = proj
                .coding
                .as_ref()
                .expect("coding (c.) form must be present on the GRCh37 projection")
                .to_string();
            assert!(
                coding.contains(":c.5A>T"),
                "expected c.5A>T round-tripped via the GRCh37 alignment, got: {coding}",
            );
        }
    }

    #[test]
    fn project_empty_allele_unknown_transcript_returns_reference_not_found() {
        // Regression: the empty-allele fast path used to silently return Ok(...)
        // for any transcript_id, even one not present in the cdot mapper. That
        // was inconsistent with every other projection path, which surfaces
        // ReferenceNotFound for unknown accessions.
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let empty_allele = HgvsVariant::Allele(AlleleVariant::cis(vec![]));
        let err = vp
            .project_variant(&empty_allele, "NM_DOES_NOT_EXIST.1")
            .expect_err("empty allele on unknown transcript should error");
        assert!(
            matches!(err, FerroError::ReferenceNotFound { ref id } if id == "NM_DOES_NOT_EXIST.1"),
            "expected ReferenceNotFound for unknown transcript, got: {:?}",
            err
        );
    }

    // ------------------------------------------------------------------
    // Direct c.→p. path for bare-NM_ inputs (#498)
    //
    // A bare transcript coding input (no genomic_context parent) carries
    // no genome alignment, so the c.→g.→CDS roundtrip cannot run. Protein
    // prediction only needs the CDS sequence + the 1-based CDS position,
    // which for an exonic c. variant the input already provides. These
    // tests pin that the direct path predicts protein and reports
    // `genomic = None` (there is no genomic form for a bare-NM_ input).
    // ------------------------------------------------------------------

    #[test]
    fn project_bare_nm_substitution_predicts_protein_without_genome() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_TEST.1:c.4C>A — bare, no genomic_context parent.
        let variant = make_coding_variant("NM_TEST.1", None);
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare-NM_ projection should succeed via the direct c.→p. path");
        // No alignment ingested for a bare NM_ ⇒ no genomic representation.
        assert!(
            proj.genomic.is_none(),
            "bare-NM_ projection should have genomic = None, got {:?}",
            proj.genomic
        );
        // The coding axis is the input itself.
        assert!(proj.coding.is_some(), "coding form should be present");
        // c.4C>A: codon 2 CGC(Arg) → AGC(Ser) ⇒ p.(Arg2Ser).
        let protein = proj.protein.expect("protein should be predicted");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Arg2Ser)");
    }

    #[test]
    fn project_nonaug_init_codon_reports_unknown_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_NOATG.1 starts with CTG (non-AUG). c.1C>G affects the initiation
        // codon → unknown consequence, reported as the spec-correct whole-protein
        // p.? (not p.(Met1?), since residue 1 is not Met) without translating
        // (#771). Before the fix this declined via the #625 non-ATG guard.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_NOATG.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::G,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_NOATG.1")
            .expect("projection should succeed");
        let protein = proj
            .protein
            .expect("an initiation-codon variant must report p.?");
        assert_eq!(format!("{protein}"), "NP_NOATG.1:p.?");
    }

    #[test]
    fn project_unreadable_start_init_codon_reports_whole_protein_unknown() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_NOSEQ.1 has no transcript sequence, so `read_cds_start_codon`
        // cannot observe the start codon (returns `ProteinSequenceUnavailable`).
        // An initiation-codon edit (`c.1A>G`) must NOT assume a canonical `ATG`
        // start and emit `p.(Met1?)`: with no observed `ATG`, the strict
        // form-choice falls to the conservative whole-protein unknown `p.?`.
        // (Pins the `start_is_atg` default; the downstream decline gate is
        // unaffected and keeps its permissive default.)
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_NOSEQ.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_NOSEQ.1")
            .expect("projection should succeed");
        let protein = proj
            .protein
            .expect("an initiation-codon variant must report a protein consequence");
        assert_eq!(format!("{protein}"), "NP_NOSEQ.1:p.?");
    }

    #[test]
    fn project_nonaug_downstream_variant_translates() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // A DOWNSTREAM variant on a legitimate non-AUG (CTG) transcript now
        // translates (#780): the near-cognate start keeps the frame, and the
        // start-codon identity is irrelevant to a downstream consequence.
        // NM_NOATG.1 = "CTGCGCTAA": codon 2 CGC (Arg); c.4C>A → AGC (Ser).
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_NOATG.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_NOATG.1")
            .expect("projection should succeed");
        let protein = proj
            .protein
            .expect("a downstream variant on a non-AUG transcript must now translate (#780)");
        assert_eq!(format!("{protein}"), "NP_NOATG.1:p.(Arg2Ser)");
    }

    #[test]
    fn project_drift_start_downstream_variant_declines() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // A downstream variant on a GENUINE off-the-start drift codon (AGT, not a
        // recognized initiator) must STILL decline — #780 widens the guard only to
        // recognized near-cognate initiators (CTG/GTG/TTG), preserving #625's
        // protection against an inconsistent cds_start. NM_DRIFT.1 = "AGTCGCTAA".
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_DRIFT.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_DRIFT.1")
            .expect("projection should succeed");
        assert!(
            proj.protein.is_none(),
            "a downstream variant on a genuine drift start must decline (#625 preserved)"
        );
    }

    #[test]
    fn project_nonaug_downstream_indel_translates() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // A DOWNSTREAM *indel* on a legitimate non-AUG (CTG) transcript now
        // translates (#780) via `predict_indel_protein` — the sibling of the
        // substitution path covered by `project_nonaug_downstream_variant_translates`.
        // NM_NOATG.1 = "CTG|CGC|TAA": `c.4_6del` removes codon 2 (CGC = Arg), an
        // in-frame single-residue deletion strictly downstream of the CTG start ⇒
        // p.(Arg2del). Before #780 the #625 non-ATG guard declined this.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_NOATG.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(4), CdsPos::new(6)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_NOATG.1")
            .expect("projection should succeed");
        let protein = proj
            .protein
            .expect("a downstream indel on a non-AUG transcript must now translate (#780)");
        assert_eq!(format!("{protein}"), "NP_NOATG.1:p.(Arg2del)");
    }

    #[test]
    fn project_nonaug_init_codon_indel_reports_unknown_protein() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // Companion to the indel-translates case: an indel OVERLAPPING the
        // initiation codon on the same non-AUG (CTG) transcript still yields the
        // spec-correct whole-protein `p.?` (consistent with #771/#772), without
        // translating. NM_NOATG.1 = "CTG|CGC|TAA": `c.1_3del` deletes the CTG
        // start codon, so the consequence is unknown — `p.?` (residue 1 is not
        // Met, so the `Met1?` form would be wrong). The `affects_init`
        // short-circuit returns this before the indel path is reached.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_NOATG.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(3)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_NOATG.1")
            .expect("projection should succeed");
        let protein = proj
            .protein
            .expect("an initiation-codon indel must report p.?");
        assert_eq!(format!("{protein}"), "NP_NOATG.1:p.?");
    }

    #[test]
    fn project_cis_allele_init_codon_member_nonaug_reports_unknown_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_NOATG.1 starts with CTG (non-AUG). The cis allele `c.[1C>G;4C>A]`
        // has an initiation-codon member (`c.1C>G`) and a downstream member
        // (`c.4C>A`, codon 2). Once initiation is disrupted the whole-allele
        // protein consequence is unknown, so it collapses to the spec-correct
        // `p.?` — even though the downstream member declines under the #625
        // non-ATG guard (which would otherwise make the all-or-nothing allele
        // protein `None`). Mirrors the single-variant #771 fix and the corpus
        // row `NM_024426.4:c.1_4delinsATGA` → `c.[1C>A;4C>A]` → `p.?`.
        let member_init = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_NOATG.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::G,
                },
            ),
        });
        let member_downstream = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_NOATG.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![member_init, member_downstream]));
        let proj = vp
            .project_variant(&allele, "NM_NOATG.1")
            .expect("projection should succeed");
        let protein = proj
            .protein
            .expect("a cis allele with an init-codon member must report a protein consequence");
        assert_eq!(format!("{protein}"), "NP_NOATG.1:p.?");
    }

    #[test]
    fn project_cis_allele_init_codon_member_atg_reports_met1_unknown() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_TEST.1 starts with ATG (canonical). The cis allele `c.[1A>G;4C>A]`
        // has an initiation-codon member (`c.1A>G`); on an ATG transcript the
        // unknown form is `p.(Met1?)`. The whole-allele protein collapses to it
        // (initiation unknown dominates the downstream missense member).
        let member_init = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let member_downstream = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![member_init, member_downstream]));
        let proj = vp
            .project_variant(&allele, "NM_TEST.1")
            .expect("projection should succeed");
        let protein = proj
            .protein
            .expect("a cis allele with an init-codon member must report a protein consequence");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Met1?)");
    }

    #[test]
    fn project_trans_allele_init_codon_member_does_not_collapse() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // Pins the `allele.phase == AllelePhase::Cis` guard at the init-codon
        // collapse (projector.rs ~1672). Members `c.1C>G` (initiation codon) and
        // `c.4C>A` (codon 2) on the genuine-drift NM_DRIFT.1 (a non-recognized
        // `AGT` start) are assembled here as a *trans* allele. The init-codon
        // collapse is cis-only: trans members are independent alleles, so they
        // fall through to the all-or-nothing protein rule. The downstream member
        // declines under the #625 guard (an `AGT` start is genuine drift, not a
        // recognized near-cognate initiator — #780 does not translate it), so the
        // all-or-nothing rule yields no whole-allele protein (`None`) — i.e. the
        // `p.?` collapse is NOT applied. Dropping the `== Cis` guard (or widening
        // it to `!= Trans`) would start collapsing this case and fail here.
        // (NM_NOATG.1's CTG start now translates its downstream member, so the
        // drift transcript is used to keep a declining member here.)
        //
        // The five other non-cis phases (`Unknown`, `Mosaic`, `Chimeric`,
        // `AndOr`, `Products`) take the same non-cis fall-through path; `Trans`
        // is the representative negative case.
        let member_init = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_DRIFT.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let member_downstream = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_DRIFT.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let allele =
            HgvsVariant::Allele(AlleleVariant::trans(vec![member_init, member_downstream]));
        let proj = vp
            .project_variant(&allele, "NM_DRIFT.1")
            .expect("projection should succeed");
        assert!(
            proj.protein.is_none(),
            "a trans allele must not trigger the cis init-codon collapse; with a downstream \
             member declining under the #625 non-ATG guard the all-or-nothing protein is None, \
             got {:?}",
            proj.protein
        );
    }

    /// #855 (delins.md): two in-cis substitutions affecting ADJACENT residues are
    /// described as a single combined delins, not a bracketed list. Mirrors the
    /// corpus row `c.[274G>T;278A>G]` → `p.(Asp92_Tyr93delinsTyrCys)`.
    /// CDS "ATGGATTATTAA" = Met-Asp-Tyr-Stop; `c.[4G>T;8A>G]` makes Asp2→Tyr and
    /// Tyr3→Cys (adjacent codons 2 and 3).
    #[test]
    fn project_cis_adjacent_substitutions_combine_to_delins() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_multicodon_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let m1 = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_DELINS.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::G,
                    alternative: Base::T,
                },
            ),
        });
        let m2 = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_DELINS.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(8)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![m1, m2]));
        let proj = vp
            .project_variant(&allele, "NM_DELINS.1")
            .expect("projection should succeed");
        let protein = proj
            .protein
            .expect("adjacent cis subs must report a protein");
        assert_eq!(
            format!("{protein}"),
            "NP_DELINS.1:p.(Asp2_Tyr3delinsTyrCys)"
        );
    }

    /// #855: an in-cis frameshift allele collapses to the whole-protein-unknown
    /// `p.?` (uncertain.md). Mirrors the corpus row `c.[41A>C;250del]` → `p.?`.
    /// Uses NON-adjacent members (`c.[4G>T;11del]`) so the allele is not merged
    /// into a single delins before the allele path runs; the `c.11del` member
    /// frameshifts the product. Output is bare `p.?` (single ProteinVariant).
    #[test]
    fn project_cis_frameshift_allele_collapses_to_unknown() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_multicodon_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let sub = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_SEP.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::G,
                    alternative: Base::T,
                },
            ),
        });
        let del = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_SEP.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(11)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![sub, del]));
        let proj = vp
            .project_variant(&allele, "NM_SEP.1")
            .expect("projection should succeed");
        let protein = proj
            .protein
            .expect("an in-cis frameshift allele must report a protein consequence");
        assert_eq!(format!("{protein}"), "NP_SEP.1:p.?");
    }

    /// #855 regression guard (delins.md): two in-cis substitutions SEPARATED by an
    /// unchanged residue stay individual (bracketed), NOT merged into a delins.
    /// CDS "ATGGATTATTGCTAA" = Met-Asp-Tyr-Cys-Stop; `c.[4G>T;11G>A]` makes
    /// Asp2→Tyr and Cys4→Tyr (codons 2 and 4; codon 3 Tyr is unchanged).
    #[test]
    fn project_cis_separated_substitutions_stay_bracketed() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_multicodon_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let m1 = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_SEP.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::G,
                    alternative: Base::T,
                },
            ),
        });
        let m2 = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_SEP.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(11)),
                NaEdit::Substitution {
                    reference: Base::G,
                    alternative: Base::A,
                },
            ),
        });
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![m1, m2]));
        let proj = vp
            .project_variant(&allele, "NM_SEP.1")
            .expect("projection should succeed");
        let protein = proj
            .protein
            .expect("separated cis subs must report a protein");
        // Separated subs stay an Allele; the Allele Display hoists the shared gene
        // symbol (unchanged pre-existing behavior). The point of this guard is that
        // they are NOT merged into a single delins.
        assert_eq!(
            format!("{protein}"),
            "NP_SEP.1(SEPGENE):p.[(Asp2Tyr);(Cys4Tyr)]"
        );
    }

    #[test]
    fn project_bare_nm_substitution_populates_noncoding_axis() {
        // c↔n axis: a bare-NM_ c. SNV carries an n. (transcript-relative) form
        // even with no genome alignment. The n. edit matches the c. edit; only
        // the coordinate frame reframes.
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_TEST.1:c.4C>A — bare, no genomic_context parent.
        let variant = make_coding_variant("NM_TEST.1", None);
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare-NM_ projection should succeed");
        assert!(proj.coding.is_some(), "coding (c.) form should be present");
        let n = proj
            .noncoding
            .as_ref()
            .expect("noncoding (n.) axis should be populated for a coding transcript");
        assert!(
            matches!(n, HgvsVariant::Tx(_)),
            "noncoding form should be an n. (Tx) variant, got {n:?}"
        );
        let n_str = n.to_string();
        assert!(
            n_str.starts_with("NM_TEST.1"),
            "n. form should render on the transcript accession, got {n_str}"
        );
        // NM_TEST.1's CDS begins at transcript position 1 (sequence ATGCGCTAA,
        // no 5'UTR), so c.4 → n.4 and the edit is unchanged.
        assert!(
            n_str.contains(":n.4C>A"),
            "c.4C>A should reframe to n.4C>A (cds_start at tx pos 1), got {n_str}"
        );
    }

    #[test]
    fn project_rejects_versioned_request_on_silent_substitution() {
        // Issue #785: projection normalizes first (#637), and normalization now
        // declines an EXPLICITLY-versioned c./n./r. request whose exact version
        // the reference does not serve but a sibling does (a silent
        // substitution). Every downstream axis — coding coordinates included —
        // would otherwise be computed against the sibling's frame, a mapping
        // whose stated reference and frame disagree. So the whole projection
        // declines with `TranscriptVersionNotExact` rather than emit a
        // confidently-wrong coding form. This extends the protein-only gate of
        // #505 to the coordinate axes for versioned requests.
        let (projector, mut provider) = make_test_provider_and_projector();
        // NM_TEST.2 is absent; a request for it is silently served the .1 bases.
        provider.mark_version_substitution("NM_TEST.2", "NM_TEST.1");
        let vp = VariantProjector::new(projector, provider);
        let variant = make_coding_variant("NM_TEST.2", None);
        let err = vp
            .project_variant(&variant, "NM_TEST.2")
            .expect_err("an explicitly-versioned silent substitution must decline wholesale");
        assert!(
            matches!(err, FerroError::TranscriptVersionNotExact { ref requested } if requested == "NM_TEST.2"),
            "expected TranscriptVersionNotExact for NM_TEST.2, got {err:?}"
        );
    }

    #[test]
    fn project_declines_protein_but_keeps_coding_for_bare_non_version_exact() {
        // Issue #505 (still in force for BARE requests): a bare (unversioned)
        // accession keeps lenient "latest version" resolution, so the #785
        // version-substitution gate does NOT apply — coding still projects
        // against the resolved transcript. But the protein path reads the CDS
        // bases directly, so when the resolved bases are not version-exact it
        // must still decline to predict a protein rather than attribute a
        // possibly-wrong product to the request. Protein is declined; the coding
        // axis is unaffected.
        let (projector, mut provider) = make_test_provider_and_projector();
        provider.mark_non_version_exact("NM_TEST");
        let vp = VariantProjector::new(projector, provider);
        // Bare NM_TEST (no .version) — not gated by #785.
        let variant = make_coding_variant("NM_TEST", None);
        let proj = vp
            .project_variant(&variant, "NM_TEST")
            .expect("a bare request must still project; only protein prediction declines");
        assert!(
            proj.protein.is_none(),
            "protein must be declined for a non-version-exact transcript, got {:?}",
            proj.protein
        );
        assert!(
            proj.coding.is_some(),
            "coding form should still be present for a bare request despite the protein gate"
        );
    }

    #[test]
    fn project_bare_nm_inframe_deletion_predicts_protein_without_genome() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_TEST.1:c.4_6del — delete codon 2 (CGC = Arg) from ATG|CGC|TAA,
        // an in-frame single-residue deletion ⇒ p.(Arg2del). Bare, no parent.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(4), CdsPos::new(6)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare-NM_ deletion projection should succeed via the direct path");
        assert!(
            proj.genomic.is_none(),
            "bare-NM_ projection has no genomic form"
        );
        let protein = proj.protein.expect("protein should be predicted");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Arg2del)");
    }

    #[test]
    fn project_bare_nm_intronic_has_no_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_TEST.1:c.4+5A>G — intronic offset on a bare NM_. Without a
        // genome alignment the intronic position cannot be placed, so
        // protein prediction is not attempted: protein = None, no panic.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::with_offset(4, 5)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare-NM_ intronic projection should succeed (no protein)");
        assert!(proj.genomic.is_none());
        assert!(
            proj.protein.is_none(),
            "intronic variant should not predict a protein consequence"
        );
        assert!(proj.is_intronic, "is_intronic flag should be set");
    }

    #[test]
    fn project_whole_exon_deletion_predicts_inframe() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // c.4-1_6+1del: both endpoints intronic, brackets the exonic span c.4_6
        // (codon 2, CGC = Arg). Clamp [4, 6] → in-frame single-residue deletion.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::with_offset(4, -1), CdsPos::with_offset(6, 1)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("whole-exon deletion projection should succeed");
        let protein = proj.protein.expect("protein should be predicted");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Arg2del)");
    }

    #[test]
    fn project_pure_intron_deletion_has_no_protein() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // c.4+1_4+9del: both offsets on base 4 → lo 5 > hi 4 → no exonic CDS → None.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::with_offset(4, 1), CdsPos::with_offset(4, 9)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("pure-intron deletion projection should succeed (no protein)");
        assert!(
            proj.protein.is_none(),
            "pure-intron deletion must not predict a protein"
        );
    }

    #[test]
    fn project_partial_exon_deletion_has_no_protein() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // c.4-1_5del: one exonic endpoint (c.5) → partial-exon / splice → None.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::with_offset(4, -1), CdsPos::new(5)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("partial-exon deletion projection should succeed (no protein)");
        assert!(
            proj.protein.is_none(),
            "partial-exon deletion must not predict a protein"
        );
    }

    // ------------------------------------------------------------------
    // `.genomic` is derived from inner projections, never the raw input
    // allele (#508 review). A transcript-coordinate (c./n./r.) allele has a
    // non-genomic `original`, and a bare-NM_ inner projection has no genomic
    // form at all — so the aggregated `.genomic` must come from the inner
    // projections' `.genomic` fields, dropping to `None` when any member
    // lacks one.
    // ------------------------------------------------------------------

    #[test]
    fn project_bare_nm_coding_allele_has_no_genomic() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // [c.4C>A;c.7T>A] on a bare NM_ (no genomic_context parent): each member
        // projects via the direct c.→p. path with genomic = None, so the
        // aggregated allele must report genomic = None rather than the input
        // c. allele.
        let mk = |base: i64, r: Base, a: Base| {
            HgvsVariant::Cds(CdsVariant {
                accession: parse_accession("NM_TEST.1"),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    CdsInterval::point(CdsPos::new(base)),
                    NaEdit::Substitution {
                        reference: r,
                        alternative: a,
                    },
                ),
            })
        };
        let allele = HgvsVariant::Allele(AlleleVariant::cis(vec![
            mk(4, Base::C, Base::A),
            mk(7, Base::T, Base::A),
        ]));
        let proj = vp
            .project_variant(&allele, "NM_TEST.1")
            .expect("bare-NM_ coding allele should project via the direct path");
        assert!(
            proj.genomic.is_none(),
            "bare-NM_ coding allele must report genomic = None, got {:?}",
            proj.genomic
        );
        assert!(proj.coding.is_some(), "coding allele should be present");
    }

    #[test]
    fn project_empty_allele_reports_no_genomic() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // An empty allele has no inner variant to derive a genomic form from,
        // so `.genomic` must be None (consistent with coding/protein, which are
        // already None for the empty allele) rather than the input allele.
        let empty = HgvsVariant::Allele(AlleleVariant::cis(vec![]));
        let proj = vp
            .project_variant(&empty, "NM_TEST.1")
            .expect("empty allele on a known transcript should project");
        assert!(
            proj.genomic.is_none(),
            "empty allele must report genomic = None, got {:?}",
            proj.genomic
        );
        assert!(proj.coding.is_none(), "empty allele has no coding form");
        assert!(proj.protein.is_none(), "empty allele has no protein form");
    }

    #[test]
    fn project_bare_nm_unknown_position_is_unsupported_not_invalid() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::uncertainty::Mu;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // A bare coding input whose position is the `?` sentinel
        // (`Single(Mu::Unknown)`) must classify as UnsupportedProjection — the
        // same contract `project_to_genomic` preserves — not InvalidCoordinates,
        // which is reserved for a genuinely missing/range coordinate (#508).
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::with_uncertainty(Mu::Unknown, Mu::Unknown),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        // Use `project_normalized` to exercise the boundary classification in
        // `project_coding_direct` directly, without the normalizer rewriting the
        // `?` position first.
        let err = vp
            .project_normalized(&variant, "NM_TEST.1")
            .expect_err("bare c.? must not project");
        assert!(
            matches!(err, FerroError::UnsupportedProjection { .. }),
            "bare c.? should be UnsupportedProjection, got {err:?}"
        );
    }

    #[test]
    fn project_bare_nm_concrete_unknown_sentinel_is_unsupported_not_utr() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // A programmatically-built `CdsPos::unknown(None)` carries the `base == 0`
        // sentinel inside a `Mu::Certain` (not `Mu::Unknown`), so it slips past
        // `resolve_uncertain_boundary`. Without an explicit `is_unknown()` guard
        // the direct path would treat `base <= 0` as 5'UTR and return `Ok`,
        // diverging from `project_to_genomic`, which rejects the sentinel as
        // `UnsupportedProjection` (#508 review).
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::unknown(None)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let err = vp
            .project_normalized(&variant, "NM_TEST.1")
            .expect_err("concrete `?` sentinel must not project as UTR");
        assert!(
            matches!(err, FerroError::UnsupportedProjection { .. }),
            "concrete `?` sentinel should be UnsupportedProjection, got {err:?}"
        );
    }

    // -------------------------------------------------------------------------
    // Issue #504: CDS-boundary edits that reach the initiation codon
    // -------------------------------------------------------------------------

    /// Bare-`NM_` coding transcript with a 5 nt 5'UTR ("GCACC") preceding the
    /// CDS "ATGTACTAA" (Met-Tyr-Stop). The 5'UTR lets `c.-1` exist so an edit
    /// can span the 5'UTR→CDS boundary into the translation initiation codon.
    ///
    /// Full transcript sequence: "GCACCATGTACTAA" (14 nt).
    ///   tx 1-5  = 5'UTR  G C A C C   → c.-5 c.-4 c.-3 c.-2 c.-1
    ///   tx 6-8  = ATG    (initiation codon) → c.1 c.2 c.3
    ///   tx 9-14 = TACTAA (Tyr-Stop)
    /// So `c.-1`=C(tx5), `c.1`=A(tx6), `c.2`=T(tx7), `c.3`=G(tx8).
    fn make_utr5_provider_and_projector() -> (Projector, MockProvider) {
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_UTR5.1".to_string(),
            CdotTranscript {
                gene_name: Some("TESTGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                // [genome_start(0-based), genome_end(0-based excl), tx_start(1-based), tx_end]
                exons: vec![[1000, 1014, 0, 14]],
                cds_start: Some(5), // 0-based: CDS starts at tx_pos 5 (5 nt 5'UTR)
                cds_end: Some(14),  // 0-based exclusive
                gene_id: None,
                protein: Some("NP_UTR5.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_UTR5.1".to_string(),
            gene_symbol: Some("TESTGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("GCACCATGTACTAA".to_string()),
            cds_start: Some(6), // 1-based inclusive: first CDS base (A of ATG)
            cds_end: Some(14),
            exons: vec![Exon::new(1, 1, 14)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1013),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: Some("NP_UTR5.1".to_string()),
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "GCACCATGTACTAA", suffix));
        (projector, provider)
    }

    /// A boundary-spanning duplication whose span reaches into the initiation
    /// codon (CDS 1–3) has an unpredictable protein consequence and must be
    /// reported as `p.(Met1?)` — even though its 5' end is in the 5'UTR
    /// (`c.-1`, CDS base ≤ 0). Before #504 the coarse `is_utr` gate dropped it
    /// to "no protein predicted" (`None`). Exercised via `project_normalized`
    /// so the canonical boundary-spanning form is fed straight to the gate.
    #[test]
    fn boundary_dup_reaching_initiation_codon_is_met1_unknown() {
        let (projector, provider) = make_utr5_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // c.-1_2dup duplicates c.-1,c.1,c.2 — the copy is inserted between c.2
        // and c.3, splitting the ATG initiation codon.
        let v = crate::parse_hgvs("NM_UTR5.1:c.-1_2dup").expect("parse c.-1_2dup");
        let result = vp
            .project_normalized(&v, "NM_UTR5.1")
            .expect("projection should succeed");
        let p = result
            .protein
            .as_ref()
            .map(|p| p.to_string())
            .expect("a boundary edit reaching the initiation codon must predict a protein");
        assert_eq!(p, "NP_UTR5.1:p.(Met1?)");
    }

    /// A 5'UTR edit that does **not** reach the initiation codon (both
    /// endpoints have CDS base ≤ 0) has no predictable CDS consequence, so the
    /// refined gate must still report `None` — not a spurious `p.(Met1?)`.
    /// Guards against #504 widening the gate too far.
    #[test]
    fn pure_5utr_edit_not_reaching_initiation_codon_predicts_no_protein() {
        let (projector, provider) = make_utr5_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // c.-3_-2del deletes two 5'UTR bases; the CDS is untouched.
        let v = crate::parse_hgvs("NM_UTR5.1:c.-3_-2del").expect("parse c.-3_-2del");
        let result = vp
            .project_normalized(&v, "NM_UTR5.1")
            .expect("projection should succeed");
        assert!(
            result.protein.is_none(),
            "pure-5'UTR edit must not predict a protein, got {:?}",
            result.protein.map(|p| p.to_string())
        );
    }

    /// A substitution in the initiation codon proper (CDS base 1–3, not UTR)
    /// also funnels through the gate's initiation-codon short-circuit, yielding
    /// `p.(Met1?)`. Confirms #504's refactor did not regress the in-CDS
    /// start-codon path (#512).
    #[test]
    fn in_cds_initiation_codon_substitution_is_met1_unknown() {
        let (projector, provider) = make_utr5_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // c.1A>G changes the first base of the ATG start codon.
        let v = crate::parse_hgvs("NM_UTR5.1:c.1A>G").expect("parse c.1A>G");
        let result = vp
            .project_normalized(&v, "NM_UTR5.1")
            .expect("projection should succeed");
        let p = result
            .protein
            .as_ref()
            .map(|p| p.to_string())
            .expect("a start-codon substitution must predict a protein");
        assert_eq!(p, "NP_UTR5.1:p.(Met1?)");
    }

    // -------------------------------------------------------------------------
    // Direct c.→p. path: special positions (pter/qter/cen) must be rejected
    // as UnsupportedProjection — not silently mis-classified as 5'UTR because
    // `base == 0` satisfies `cds_start.base <= 0`.
    // -------------------------------------------------------------------------

    #[test]
    fn project_bare_nm_special_pter_is_unsupported_not_5utr() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // A `pter` CDS position has `base == 0` and `special == Some(Pter)`.
        // Without the `is_special()` guard, `cds_start.base <= 0` evaluates to
        // `true` and the variant is silently classified as 5'UTR, producing a
        // wrong protein consequence. The guard must reject it as
        // `UnsupportedProjection` (#488).
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::pter()),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let err = vp
            .project_normalized(&variant, "NM_TEST.1")
            .expect_err("pter sentinel must not project as 5'UTR");
        assert!(
            matches!(err, FerroError::UnsupportedProjection { .. }),
            "pter sentinel should be UnsupportedProjection, got {err:?}"
        );
    }

    #[test]
    fn project_bare_nm_special_qter_is_unsupported() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::variant::{CdsVariant, LocEdit};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // `qter` endpoint: base == 0, special == Some(Qter). Same bug path as
        // pter — must be rejected rather than silently misclassified.
        let variant = HgvsVariant::Cds(CdsVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::qter()),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });
        let err = vp
            .project_normalized(&variant, "NM_TEST.1")
            .expect_err("qter sentinel must not project as UTR");
        assert!(
            matches!(err, FerroError::UnsupportedProjection { .. }),
            "qter sentinel should be UnsupportedProjection, got {err:?}"
        );
    }

    // ------------------------------------------------------------------
    // Direct n./r.→CDS→p. path for bare transcript inputs (#506)
    //
    // A bare non-coding (`n.`) or RNA (`r.`) input on a coding transcript
    // (no `genomic_context` parent) carries no genome alignment, so the
    // c.→g.→CDS roundtrip cannot run. But an `n.` base is a 1-based
    // transcript position: convert it to a CDS position via the
    // transcript's `cds_start` (exon/CIGAR-aware `tx_to_cds`), then run the
    // same protein-consequence prediction the coding path uses. The
    // resulting projection has `genomic = None`.
    // ------------------------------------------------------------------

    /// A non-coding transcript (no CDS) for the "no protein on n. input"
    /// guard. Sequence "ACGTACGTAC" (10 nt), exon tx 1..10, plus strand, no
    /// `cds_start`/`cds_end`.
    fn make_noncoding_provider_and_projector() -> (Projector, MockProvider) {
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NR_TEST.1".to_string(),
            CdotTranscript {
                gene_name: Some("NCGENE".to_string()),
                contig: "chr1".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1010, 0, 10]],
                cds_start: None,
                cds_end: None,
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NR_TEST.1".to_string(),
            gene_symbol: Some("NCGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("ACGTACGTAC".to_string()),
            cds_start: None,
            cds_end: None,
            exons: vec![Exon::new(1, 1, 10)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1009),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ACGTACGTAC", suffix));
        (projector, provider)
    }

    #[test]
    fn project_bare_n_substitution_predicts_protein_without_genome() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{LocEdit, TxVariant};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NM_TEST.1's CDS begins at transcript position 1 (sequence ATGCGCTAA,
        // no 5'UTR), so n.4 == c.4. n.4C>A: codon 2 CGC(Arg) → AGC(Ser).
        let variant = HgvsVariant::Tx(TxVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare n. input should project via the direct n.→CDS→p. path");
        assert!(
            proj.genomic.is_none(),
            "bare n. input should report genomic = None, got {:?}",
            proj.genomic
        );
        // The non-coding axis is the input itself; the coding axis is derived.
        assert!(
            proj.noncoding.is_some(),
            "noncoding (n.) form should be present"
        );
        let protein = proj.protein.expect("protein should be predicted");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Arg2Ser)");
    }

    #[test]
    fn project_bare_n_populates_coding_axis() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{LocEdit, TxVariant};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // n.4C>A → c.4C>A (cds_start at tx pos 1). The coding form reframes the
        // transcript-relative position to CDS-relative; the edit is unchanged.
        let variant = HgvsVariant::Tx(TxVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare n. input should project");
        let coding = proj
            .coding
            .as_ref()
            .expect("coding (c.) axis should be derived for a coding transcript");
        assert!(
            matches!(coding, HgvsVariant::Cds(_)),
            "coding form should be a c. (Cds) variant, got {coding:?}"
        );
        assert!(
            coding.to_string().contains(":c.4C>A"),
            "n.4C>A should reframe to c.4C>A, got {coding}"
        );
    }

    #[test]
    fn project_bare_r_substitution_predicts_protein_without_genome() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::location::RnaPos;
        use crate::hgvs::variant::{LocEdit, RnaVariant};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // r.4c>a mirrors n.4C>A on the sense strand. The RNA edit carries
        // `Base::U`-style RNA bases (here `c`/`a`, which are DNA-compatible);
        // the direct path translates U→T before reading the DNA CDS, so the
        // protein consequence matches the coding path: p.(Arg2Ser).
        let variant = HgvsVariant::Rna(RnaVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare r. input should project via the direct r.→CDS→p. path");
        assert!(
            proj.genomic.is_none(),
            "bare r. input should report genomic = None, got {:?}",
            proj.genomic
        );
        let protein = proj.protein.expect("protein should be predicted");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Arg2Ser)");
    }

    #[test]
    fn project_bare_r_with_u_base_predicts_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::location::RnaPos;
        use crate::hgvs::variant::{LocEdit, RnaVariant};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // r.2u>a: tx pos 2 is the `T` of ATG (RNA reads `u`). The U→T
        // translation must run so the ref base matches the DNA CDS. c.2T>A
        // changes the initiation codon ATG→AAG ⇒ p.(Met1?).
        let variant = HgvsVariant::Rna(RnaVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(2)),
                NaEdit::Substitution {
                    reference: Base::U,
                    alternative: Base::A,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare r. input with a U base should project");
        let protein = proj
            .protein
            .expect("a start-codon r. substitution should predict a protein");
        assert_eq!(format!("{protein}"), "NP_TEST.1:p.(Met1?)");
    }

    #[test]
    fn project_bare_n_intronic_has_no_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{LocEdit, TxVariant};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // n.4+5A>G — intronic offset on a bare NM_. Without a genome alignment
        // the intronic position cannot be placed, so protein prediction is not
        // attempted: protein = None, no panic.
        let variant = HgvsVariant::Tx(TxVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::with_offset(4, 5)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_TEST.1")
            .expect("bare n. intronic projection should succeed (no protein)");
        assert!(proj.genomic.is_none());
        assert!(
            proj.protein.is_none(),
            "intronic n. variant should not predict a protein consequence"
        );
        assert!(proj.is_intronic, "is_intronic flag should be set");
    }

    #[test]
    fn project_bare_n_on_noncoding_transcript_has_no_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{LocEdit, TxVariant};
        let (projector, provider) = make_noncoding_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // NR_TEST.1 has no CDS, so an n. variant on it has no protein
        // consequence — the projection succeeds with protein = None.
        let variant = HgvsVariant::Tx(TxVariant {
            accession: parse_accession("NR_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::T,
                    alternative: Base::G,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NR_TEST.1")
            .expect("bare n. on a non-coding transcript should project (no protein)");
        assert!(proj.genomic.is_none());
        assert!(
            proj.protein.is_none(),
            "n. on a non-coding transcript must not predict a protein, got {:?}",
            proj.protein
        );
        // The non-coding axis still carries the input n. form.
        assert!(proj.noncoding.is_some(), "noncoding form should be present");
    }

    #[test]
    fn project_bare_n_utr_has_no_protein() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{LocEdit, TxVariant};
        // NM_UTR5.1: 5 nt 5'UTR then CDS at tx pos 6. n.3 lies in the 5'UTR
        // (c.-3) and does not reach the initiation codon, so no protein.
        let (projector, provider) = make_utr5_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let variant = HgvsVariant::Tx(TxVariant {
            accession: parse_accession("NM_UTR5.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(3)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let proj = vp
            .project_variant(&variant, "NM_UTR5.1")
            .expect("bare n. 5'UTR projection should succeed (no protein)");
        assert!(proj.genomic.is_none());
        assert!(
            proj.is_utr,
            "is_utr flag should be set for a 5'UTR n. position"
        );
        assert!(
            proj.protein.is_none(),
            "pure-5'UTR n. position must not predict a protein, got {:?}",
            proj.protein
        );
    }

    #[test]
    fn project_bare_n_transcript_id_mismatch_is_rejected() {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::location::TxPos;
        use crate::hgvs::variant::{LocEdit, TxVariant};
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // An n. input on NM_TEST.1 projected against a different transcript
        // must be rejected, mirroring the coding path's guard.
        let variant = HgvsVariant::Tx(TxVariant {
            accession: parse_accession("NM_TEST.1"),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(4)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::A,
                },
            ),
        });
        let err = vp
            .project_normalized(&variant, "NM_OTHER.1")
            .expect_err("transcript_id mismatch must be rejected");
        assert!(
            matches!(err, FerroError::UnsupportedProjection { .. }),
            "transcript_id mismatch should be UnsupportedProjection, got {err:?}"
        );
    }

    #[test]
    fn project_variant_unrelated_transcript_id_mismatch_is_rejected() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // A bare `c.` input on NM_TEST.1 projected against an unrelated transcript
        // must still error, exercising `reconcile_self_projection_target`'s
        // fallthrough on the `project_variant` path: the variant carries no
        // genomic_context, so reconciliation returns the requested target
        // unchanged and the `project_coding_direct` guard rejects the genuine
        // mismatch (#783).
        let variant = make_coding_variant("NM_TEST.1", None);
        let err = vp
            .project_variant(&variant, "NM_NOATG.1")
            .expect_err("unrelated transcript_id mismatch must be rejected");
        assert!(
            matches!(err, FerroError::UnsupportedProjection { .. }),
            "transcript_id mismatch should be UnsupportedProjection, got {err:?}"
        );
    }

    /// Guard test: `MockProvider` (no assembly-report table) must return the same
    /// build inference as the bare hardcoded free function for an NC_ accession —
    /// i.e. routing through the provider is behavior-preserving (#716).
    #[test]
    fn infer_input_build_routes_through_provider() {
        use crate::hgvs::variant::Accession;
        let provider = crate::reference::mock::MockProvider::new();
        let acc = Accession::new("NC", "000017", Some(11));
        assert_eq!(
            provider.infer_genome_build(&acc),
            crate::liftover::aliases::infer_genome_build_from_accession(&acc)
        );
    }

    /// A provider that delegates everything to an inner [`MockProvider`] but
    /// reports a deliberately *different* build than the hardcoded heuristic for
    /// one accession. Used to prove a projector path observes the provider's
    /// `infer_genome_build` rather than the hardcoded free function (#716).
    #[derive(Clone)]
    struct BuildOverrideProvider {
        inner: MockProvider,
    }

    impl ReferenceProvider for BuildOverrideProvider {
        fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
            self.inner.get_transcript(id)
        }

        fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
            self.inner.get_sequence(id, start, end)
        }

        fn infer_genome_build(&self, accession: &Accession) -> Option<&'static str> {
            // NC_000017.11 is GRCh38 under the hardcoded heuristic; force GRCh37
            // so a test can tell the provider's answer apart from the fallback.
            if accession.full() == "NC_000017.11" {
                Some("GRCh37")
            } else {
                self.inner.infer_genome_build(accession)
            }
        }
    }

    /// Override test: the projector must consult the provider's
    /// `infer_genome_build`, not the hardcoded helper. With a provider that
    /// deliberately classifies the `NC_000017.11` parent as GRCh37 (the heuristic
    /// says GRCh38), `infer_input_build` must observe the provider's GRCh37 — if
    /// it regressed to the hardcoded helper it would return GRCh38 and this would
    /// fail (#716).
    #[test]
    fn infer_input_build_observes_provider_override() {
        let parent = Accession::new("NC", "000017", Some(11));
        // Sanity-check the premise: the hardcoded helper disagrees with the
        // override, so the assertion below genuinely distinguishes the two paths.
        assert_eq!(
            crate::liftover::aliases::infer_genome_build_from_accession(&parent),
            Some("GRCh38"),
        );
        let provider = BuildOverrideProvider {
            inner: MockProvider::new(),
        };
        let vp = VariantProjector::new(Projector::new(CdotMapper::new()), provider);
        // A bare NM (no inherent build) whose genomic_context is the NC_ parent —
        // `infer_input_build` falls through to consult the parent via the provider.
        let variant = make_coding_variant("NM_TEST.1", Some(parent));
        assert_eq!(vp.infer_input_build(&variant), Some("GRCh37"));
    }

    // ------------------------------------------------------------------
    // #693: reframe_rna_from_input
    // ------------------------------------------------------------------

    /// Build a bare predicted RNA `NM_x:r.(<...>)` carrying a synthesized gene
    /// symbol, as `predict_rna` produces on the genome-pivot / direct paths.
    fn synthesized_rna(accession: &str, gene: Option<&str>) -> HgvsVariant {
        use crate::hgvs::edit::{Base, NaEdit};
        use crate::hgvs::interval::RnaInterval;
        use crate::hgvs::location::RnaPos;
        use crate::hgvs::variant::{LocEdit, RnaVariant};
        let pos = RnaPos {
            base: 41,
            offset: None,
            utr3: false,
        };
        let edit = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::C,
        };
        HgvsVariant::Rna(RnaVariant {
            accession: parse_accession(accession),
            gene_symbol: gene.map(str::to_string),
            loc_edit: LocEdit::new_predicted(RnaInterval::new(pos, pos), edit),
        })
    }

    /// An `NG_(NM_)` input re-frames the predicted RNA under its genomic parent
    /// and drops the synthesized gene symbol — `NM_x(GENE):r.(…)` →
    /// `NG_y(NM_x):r.(…)` (#693).
    #[test]
    fn reframe_rna_carries_input_genomic_context() {
        let input = crate::parse_hgvs("NG_009299.1(NM_017668.3):c.41A>C").unwrap();
        let rna = synthesized_rna("NM_017668.3", Some("NDE1"));
        let reframed = VariantProjector::<MockProvider>::reframe_rna_from_input(Some(rna), &input)
            .expect("some");
        assert_eq!(reframed.to_string(), "NG_009299.1(NM_017668.3):r.(41a>c)");
    }

    /// A bare input (no genomic context, no input gene symbol) strips any
    /// synthesized gene symbol so the RNA renders bare — `NM_x(GENE):r.(…)` →
    /// `NM_x:r.(…)` (#693, #121: never fabricate a selector absent from input).
    #[test]
    fn reframe_rna_strips_synthesized_gene_for_bare_input() {
        let input = crate::parse_hgvs("NM_017668.3:c.41A>C").unwrap();
        let rna = synthesized_rna("NM_017668.3", Some("NDE1"));
        let reframed = VariantProjector::<MockProvider>::reframe_rna_from_input(Some(rna), &input)
            .expect("some");
        assert_eq!(reframed.to_string(), "NM_017668.3:r.(41a>c)");
    }

    /// A predicted RNA allele re-frames every member so the first-member-derived
    /// compact prefix carries the input's genomic context (#693).
    #[test]
    fn reframe_rna_reframes_allele_members() {
        let input = crate::parse_hgvs("NG_012337.1(NM_003002.2):c.[274G>T;278A>G]").unwrap();
        let members = vec![
            synthesized_rna("NM_003002.2", Some("SDHD")),
            synthesized_rna("NM_003002.2", Some("SDHD")),
        ];
        let allele = HgvsVariant::Allele(AlleleVariant::new_uncertain(members, AllelePhase::Cis));
        let reframed =
            VariantProjector::<MockProvider>::reframe_rna_from_input(Some(allele), &input)
                .expect("some");
        // First-member prefix carries NG_; predicted wrapper inside the bracket.
        assert!(
            reframed
                .to_string()
                .starts_with("NG_012337.1(NM_003002.2):r.["),
            "got {reframed}"
        );
        // The compact prefix only reflects the first member, so assert every RNA
        // member is reframed: each carries the NG_ parent and no synthesized gene
        // symbol. Guards against a regression that leaves later members bare.
        let HgvsVariant::Allele(allele) = &reframed else {
            panic!("expected RNA allele, got {reframed}");
        };
        for member in &allele.variants {
            let HgvsVariant::Rna(r) = member else {
                panic!("expected RNA member, got {member}");
            };
            assert_eq!(
                r.accession.genomic_context.as_deref().map(|a| a.full()),
                Some("NG_012337.1".to_string()),
                "member was not framed under the NG_ parent: {member}"
            );
            assert!(
                r.gene_symbol.is_none(),
                "member kept a synthesized gene symbol: {member}"
            );
        }
    }

    /// An allele input that carries a user-provided gene selector on its members
    /// (e.g. `NM_x(GENE):c.[...]`) keeps that selector on the predicted `r.`
    /// allele. `gene_symbol()` returns `None` for an `Allele`, so without
    /// deriving the symbol from the members the selector would be wrongly
    /// stripped (#693).
    #[test]
    fn reframe_rna_preserves_allele_member_gene_symbol() {
        let input = crate::parse_hgvs("NM_003002.2(SDHD):c.[274G>T;278A>G]").unwrap();
        // The synthesized RNA members deliberately carry a *different* selector
        // (`SYNTH`) than the input's `SDHD`. `reframe_rna_from_input` derives the
        // output selector from the INPUT allele members, not the predicted ones,
        // so the output must show `SDHD`. Asserting `SDHD` here (rather than the
        // synthesized `SYNTH`) pins the input-derived path: a regression that
        // read the predicted-member symbol would surface `SYNTH` and fail.
        let members = vec![
            synthesized_rna("NM_003002.2", Some("SYNTH")),
            synthesized_rna("NM_003002.2", Some("SYNTH")),
        ];
        let allele = HgvsVariant::Allele(AlleleVariant::new_uncertain(members, AllelePhase::Cis));
        let reframed =
            VariantProjector::<MockProvider>::reframe_rna_from_input(Some(allele), &input)
                .expect("some");
        let HgvsVariant::Allele(allele) = &reframed else {
            panic!("expected RNA allele, got {reframed}");
        };
        for member in &allele.variants {
            let HgvsVariant::Rna(r) = member else {
                panic!("expected RNA member, got {member}");
            };
            assert!(
                r.accession.genomic_context.is_none(),
                "member fabricated a genomic context: {member}"
            );
            assert_eq!(
                r.gene_symbol.as_deref(),
                Some("SDHD"),
                "member dropped the input gene selector: {member}"
            );
        }
    }

    /// When allele members disagree on a gene selector, none is carried onto the
    /// predicted `r.` allele rather than picking the first member's (#693).
    #[test]
    fn reframe_rna_drops_disagreeing_allele_gene_symbols() {
        // Disagreement must live on the INPUT allele members: the first parses
        // with `SDHD`, then the second's selector is forced to `OTHER`. The
        // selector is derived from these input members, so unanimity fails and
        // no symbol is carried. (The synthesized members keep a consistent
        // `SDHD`, proving the drop is driven by input disagreement, not the
        // predicted-member symbols — which `reframe_rna_from_input` ignores.)
        let mut input = crate::parse_hgvs("NM_003002.2(SDHD):c.[274G>T;278A>G]").unwrap();
        {
            let HgvsVariant::Allele(input_allele) = &mut input else {
                panic!("expected input allele, got {input}");
            };
            let [_, second] = input_allele.variants.as_mut_slice() else {
                panic!("expected two input members, got {input}");
            };
            let HgvsVariant::Cds(cds) = second else {
                panic!("expected Cds input member, got {second}");
            };
            cds.gene_symbol = Some("OTHER".to_string());
        }
        let members = vec![
            synthesized_rna("NM_003002.2", Some("SDHD")),
            synthesized_rna("NM_003002.2", Some("SDHD")),
        ];
        let allele = HgvsVariant::Allele(AlleleVariant::new_uncertain(members, AllelePhase::Cis));
        let reframed =
            VariantProjector::<MockProvider>::reframe_rna_from_input(Some(allele), &input)
                .expect("some");
        let HgvsVariant::Allele(allele) = &reframed else {
            panic!("expected RNA allele, got {reframed}");
        };
        for member in &allele.variants {
            let HgvsVariant::Rna(r) = member else {
                panic!("expected RNA member, got {member}");
            };
            assert!(
                r.gene_symbol.is_none(),
                "member kept a gene selector despite disagreement: {member}"
            );
        }
    }

    /// `None` (axis unavailable) passes through unchanged.
    #[test]
    fn reframe_rna_passes_none_through() {
        let input = crate::parse_hgvs("NM_017668.3:c.41A>C").unwrap();
        assert!(VariantProjector::<MockProvider>::reframe_rna_from_input(None, &input).is_none());
    }
}
