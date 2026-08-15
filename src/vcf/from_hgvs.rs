//! HGVS to VCF conversion
//!
//! This module provides functionality to convert HGVS variants to VCF records.
//!
//! # Coordinate System
//!
//! | Context | Basis | Notes |
//! |---------|-------|-------|
//! | HGVS positions | 1-based | Input positions from parsed variants |
//! | VCF POS field | 1-based | Output position in VCF record |
//! | Array indexing | 0-based | Internal sequence access via `- 1` |
//!
//! VCF requires an "anchor base" for insertions and deletions. This module
//! fetches the preceding reference base when converting HGVS insertions.
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::convert::mapper::CoordinateMapper;
use crate::error::FerroError;
use crate::hgvs::edit::NaEdit;
use crate::hgvs::interval::{interval_is_wraparound, GenomeInterval};
use crate::hgvs::location::{CdsPos, RnaPos, TxPos};
use crate::hgvs::variant::{CdsVariant, GenomeVariant, HgvsVariant, RnaVariant, TxVariant};
use crate::project::edit::u_to_t_edit;
use crate::reference::transcript::{GenomeBuild, Transcript};
use crate::reference::ReferenceProvider;
use crate::sequence::reverse_complement;

use super::record::VcfRecord;

/// Refuse a genomic interval whose start or end holds something VCF cannot be
/// given a coordinate for: a `+`/`-` offset, or a `pter`/`qter`/`cen` special
/// position.
///
/// The two are the same defect in two fields of one struct.
/// [`GenomePos`](crate::hgvs::location::GenomePos) can hold either because the
/// parser accepts both, and neither half of the conversion can honour them:
///
/// - an **offset** has nothing on a genomic accession to be measured against —
///   there is no exon structure, and `checklist.md:16` prohibits an offset on a
///   genomic position outright — and a VCF `POS` is a single integer with no
///   offset notation;
/// - a **special position** names a landmark of the assembled chromosome
///   (`pter`, `qter`) or of its centromere annotation (`cen`), none of which a
///   sequence accession carries. `GenomePos` stores `base: 0` beside the
///   `special` marker precisely because there is no coordinate to store.
///
/// What [`HgvsToVcfConverter::convert_genome`] did instead was read `base` and
/// drop both, which is the one answer worse than either honest option. It
/// produced the same failure the SPDI path was fixed for twice over — see
/// `spdi::convert::reject_unresolvable_genomic_position`, whose ruling this
/// mirrors, including the verdict:
///
/// - **Dropping the offset** made `g.10+2del`, `g.10-2del` and `g.10del` — three
///   descriptions `normalize` keeps distinct — emit one byte-identical VCF row
///   (`POS=9 REF=AC ALT=["A"]`), with no error at all in either build profile.
/// - **Flattening a special position** onto `base: 0` sent 0 into
///   [`HgvsToVcfConverter::build_vcf_record`]'s `start - 1`, which panicked in
///   debug and wrapped to `18446744073709551614` in release. The wrap was caught
///   only because the fixture's provider bounds-checks its fetch; a provider
///   that does not would emit a record at the wrapped coordinate.
///
/// Refusing is what the rest of this module already does with the positions it
/// cannot lower — the `c.`/`n.`/`r.` arms decline an intronic position rather
/// than pretend it sits on the exon, and every unrepresentable edit declines by
/// name. The genomic arm was the one hole in an otherwise consistent policy.
///
/// **One guard, not six.** The SPDI sibling needed a call site per genomic arm
/// and passed each its axis label; here `g.`, `m.` and `o.` all funnel through
/// `convert_genome` (the `m.`/`o.` arms of [`HgvsToVcfConverter::convert`]
/// rebuild the variant as a `GenomeVariant`), so there is one guarded funnel and
/// no per-site label to get wrong. The offending description is quoted instead,
/// which is what a caller needs to act.
fn reject_unresolvable_genomic_position(
    interval: &GenomeInterval,
    hgvs_string: &str,
) -> Result<(), FerroError> {
    // Name the side the offending endpoint came from. Both endpoints are checked, so
    // `g.4_10+2del` reports the value without saying which of the two carried it, and an
    // interval whose two ends fail differently is indistinguishable from one whose start
    // fails twice. The label is diagnostic only — the refusal itself is unchanged.
    for (side, endpoint) in [
        ("interval start", interval.start.inner()),
        ("interval end", interval.end.inner()),
    ]
    .into_iter()
    .filter_map(|(side, endpoint)| endpoint.map(|endpoint| (side, endpoint)))
    {
        if endpoint.offset.is_some() {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "genomic position {endpoint} ({side}) carries a +/- offset: a genomic \
                     position cannot have one, and a VCF POS has no offset notation to express \
                     it. Drop the offset, or describe the variant on a transcript accession \
                     where an offset is meaningful. Variant: {hgvs_string}"
                ),
            });
        }
        if let Some(special) = endpoint.special {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "genomic position `{special}` ({side}) names no numeric coordinate: `pter`, \
                     `qter` and `cen` are landmarks of the assembled chromosome, not \
                     positions on the sequence, and a VCF POS has no notation for them. \
                     Give the position a number, or keep the description in HGVS. \
                     Variant: {hgvs_string}"
                ),
            });
        }
    }
    Ok(())
}

/// The 1-based reference position an anchored edit hangs its anchor base on:
/// the base immediately before `start`.
///
/// VCF anchors a deletion — and a length-changing `delins` — on the preceding
/// reference base, so `start` must have a predecessor. This is `start - 1` on
/// `u64`, and it was written that way: at `start = 0` it panicked in debug and
/// wrapped to `u64::MAX - 1` in release, which then either escaped as a record
/// at a nonsense coordinate or came back as a bounds error blaming the provider.
///
/// The check lives here rather than only at the genomic arm because
/// [`HgvsToVcfConverter::build_vcf_record`] is the **shared trunk**:
/// `convert_cds`, `convert_tx`, `convert_rna` and `convert_genome` all reach it.
/// The `pter` guard above closes the one caller that is known to have fed it a
/// zero; this closes the arithmetic for every caller, present and future, and
/// makes the two build profiles answer identically.
fn anchor_position_before(start: u64, hgvs_string: &str) -> Result<u64, FerroError> {
    start
        .checked_sub(1)
        .ok_or_else(|| FerroError::InvalidCoordinates {
            msg: format!(
                "1-based position {start} has no anchor base before it: VCF anchors a \
             deletion or a length-changing delins on the preceding reference base, and \
             position {start} has no predecessor. Variant: {hgvs_string}"
            ),
        })
}

/// The nucleotides an inserted payload contributes to a VCF ALT allele.
///
/// A VCF allele field holds bases. `Display for InsertedSequence` renders the
/// *HGVS spelling*, which coincides with the bases for exactly one shape —
/// [`InsertedSequence::Literal`] — and for a single-part
/// [`InsertedSequence::Complex`] wrapping one. Every other shape renders as a
/// description: `212_221inv` for a range, `10` for a count, `(?)` for an
/// uncertainty. Formatting one straight into an allele produced records like
/// `ALT=C212_221inv`, which is not a VCF record at all.
///
/// So this is an **allowlist**: it spells the shapes that are already bases and
/// declines the rest. The alternative — resolving a range against the reference
/// — cannot be done safely here. [`HgvsToVcfConverter::build_vcf_record`] is the
/// shared trunk for `convert_cds`, `convert_tx`, `convert_rna` and
/// `convert_genome`, and it receives the edit in the *description's own*
/// coordinate frame while
/// [`HgvsToVcfConverter::get_reference_sequence`] reads the genomic contig by
/// chromosome name. A `c.`-frame range resolved against genomic coordinates
/// would return the wrong bases with nothing to show for it — trading a visibly
/// malformed record for a silently wrong one. Doing it correctly needs the
/// axis-aware lookup `crate::normalize::rules::fetch_position_range_bases` has
/// and this converter does not, which is a feature rather than this guard.
///
/// Declining is not a loss of reach: the range spellings are alternative
/// renderings of a payload that is equally legal as a literal, and the literal
/// converts. The error names the payload so a caller can see which one to
/// re-spell.
fn inserted_sequence_bases(
    sequence: &crate::hgvs::edit::InsertedSequence,
    hgvs_string: &str,
) -> Result<String, FerroError> {
    use crate::hgvs::edit::{InsertedPart, InsertedSequence};

    match sequence {
        InsertedSequence::Literal(seq) => Ok(seq.to_string()),
        // A single-part `Complex` renders as its part (HGVS v21 reserves
        // brackets for two or more), so an all-literal `Complex` already
        // produced correct bases in the one-part case and a bracketed
        // `[A;B]` in the rest. Concatenating is right for both.
        InsertedSequence::Complex(parts)
            if parts
                .iter()
                .all(|part| matches!(part, InsertedPart::Literal(_))) =>
        {
            Ok(parts
                .iter()
                .map(|part| match part {
                    InsertedPart::Literal(seq) => seq.to_string(),
                    _ => unreachable!("guarded by the match arm above"),
                })
                .collect())
        }
        // `ins` with no payload contributes nothing; the anchor base alone is
        // already what the surrounding arms emit for it.
        InsertedSequence::Empty => Ok(String::new()),
        other => Err(FerroError::UnsupportedVariant {
            variant_type: format!(
                "inserted payload `{other}` is an HGVS description, not nucleotides, and a VCF \
                 allele field holds only bases; ferro does not resolve a payload against the \
                 reference on this path because the edit arrives in the description's own \
                 coordinate frame. Re-spell the payload as literal bases. Variant: {hgvs_string}"
            ),
        }),
    }
}

/// Result of converting an HGVS variant to VCF
#[derive(Debug, Clone)]
pub struct VcfConversionResult {
    /// The generated VCF record
    pub record: VcfRecord,
    /// Original HGVS string (for INFO annotation)
    pub hgvs_string: String,
    /// Warning messages (e.g., if normalization was applied)
    pub warnings: Vec<String>,
}

/// Converter for HGVS variants to VCF records
pub struct HgvsToVcfConverter<'a, P: ReferenceProvider> {
    transcript: &'a Transcript,
    mapper: CoordinateMapper<'a>,
    provider: &'a P,
    genome_build: GenomeBuild,
}

impl<'a, P: ReferenceProvider> HgvsToVcfConverter<'a, P> {
    /// Create a new converter for a specific transcript
    pub fn new(transcript: &'a Transcript, provider: &'a P) -> Self {
        let mapper = CoordinateMapper::new(transcript);
        let genome_build = transcript.genome_build;
        Self {
            transcript,
            mapper,
            provider,
            genome_build,
        }
    }

    /// Set the genome build for output VCF records
    pub fn with_genome_build(mut self, build: GenomeBuild) -> Self {
        self.genome_build = build;
        self
    }

    /// Convert an HGVS variant to a VCF record
    pub fn convert(&self, variant: &HgvsVariant) -> Result<VcfConversionResult, FerroError> {
        match variant {
            HgvsVariant::Cds(cds) => self.convert_cds(cds),
            HgvsVariant::Tx(tx) => self.convert_tx(tx),
            HgvsVariant::Genome(genome) => self.convert_genome(genome),
            HgvsVariant::Rna(rna) => self.convert_rna(rna),
            HgvsVariant::Protein(_) => Err(FerroError::ConversionError {
                msg: "Protein variant conversion to VCF not supported (requires back-translation)"
                    .to_string(),
            }),
            HgvsVariant::Mt(mt) => {
                if interval_is_wraparound(&mt.loc_edit.location) {
                    return Err(FerroError::ConversionError {
                        msg: format!(
                            "Cannot convert wraparound m. variant to VCF: \
                             VCF has no representation for circular-contig records. \
                             Variant: {}",
                            mt
                        ),
                    });
                }
                // Mitochondrial variants use the same format as genomic
                let genome = GenomeVariant {
                    accession: mt.accession.clone(),
                    gene_symbol: mt.gene_symbol.clone(),
                    loc_edit: mt.loc_edit.clone(),
                };
                self.convert_genome(&genome)
            }
            HgvsVariant::Allele(allele) => {
                // For simple alleles with a single variant, convert directly
                // For complex alleles with multiple variants, convert the first and warn
                if allele.variants.is_empty() {
                    return Err(FerroError::ConversionError {
                        msg: "Empty allele variant cannot be converted to VCF".to_string(),
                    });
                }

                if allele.variants.len() == 1 {
                    // Single variant allele - convert the inner variant
                    return self.convert(&allele.variants[0]);
                }

                // Multiple variants in allele: each maps to its own VCF record,
                // so a single `VcfConversionResult` cannot represent them.
                // Direct callers to `convert_all`, which decomposes the allele
                // into one result per member, rather than silently dropping all
                // but one.
                Err(FerroError::ConversionError {
                    msg: format!(
                        "Complex allele with {} variants maps to multiple VCF records; \
                         use `convert_all` to decompose it into per-member records",
                        allele.variants.len()
                    ),
                })
            }
            HgvsVariant::NullAllele => Err(FerroError::ConversionError {
                msg: "Null allele marker [0] cannot be converted to VCF".to_string(),
            }),
            HgvsVariant::UnknownAllele => Err(FerroError::ConversionError {
                msg: "Unknown allele marker [?] cannot be converted to VCF".to_string(),
            }),
            HgvsVariant::Circular(circular) => {
                if interval_is_wraparound(&circular.loc_edit.location) {
                    return Err(FerroError::ConversionError {
                        msg: format!(
                            "Cannot convert wraparound o. variant to VCF: \
                             VCF has no representation for circular-contig records. \
                             Variant: {}",
                            circular
                        ),
                    });
                }
                // Circular variants use the same format as genomic, convert to linear coordinates
                let genome = GenomeVariant {
                    accession: circular.accession.clone(),
                    gene_symbol: circular.gene_symbol.clone(),
                    loc_edit: circular.loc_edit.clone(),
                };
                self.convert_genome(&genome)
            }
            HgvsVariant::RnaFusion(_) => Err(FerroError::ConversionError {
                msg: "RNA fusion (::) conversion to VCF not supported (complex structural variant)"
                    .to_string(),
            }),
            HgvsVariant::GenomeRing(_) => Err(FerroError::ConversionError {
                msg:
                    "genome ring (::) conversion to VCF not supported (complex structural variant)"
                        .to_string(),
            }),
            HgvsVariant::Supernumerary(_) => Err(FerroError::ConversionError {
                msg: "supernumerary (sup) conversion to VCF not supported \
                      (complex structural variant)"
                    .to_string(),
            }),
        }
    }

    /// Convert an HGVS variant to one VCF record per genomic event.
    ///
    /// This is the decomposition-aware entry point. For every variant that
    /// [`Self::convert`] supports, it returns a single-element vector wrapping
    /// the same result. For a multi-member allele (`var.[a;b;...]`) — which
    /// maps to multiple, independently-positioned VCF records and therefore
    /// cannot be expressed as one [`VcfConversionResult`] — it decomposes the
    /// allele into its member variants and converts each, returning one result
    /// per member in declaration order.
    ///
    /// Decomposition is recursive: a member that is itself an allele is
    /// flattened, so a nested `[[a;b];c]` yields `[a, b, c]`. The allele's
    /// phase relationship (cis/trans/mosaic/…) is not encoded in the per-record
    /// output — VCF has no phase column at the record level — so this lowers
    /// every member to its own record regardless of phase.
    ///
    /// A member whose shape cannot be lowered to VCF (e.g. `[0]`, an RNA
    /// whole-entity edit, or a structural arm) fails the whole call with a
    /// clear error naming the offending member's 0-based index and HGVS string,
    /// rather than silently dropping it. An empty allele is rejected.
    pub fn convert_all(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Vec<VcfConversionResult>, FerroError> {
        match variant {
            HgvsVariant::Allele(allele) => {
                if allele.variants.is_empty() {
                    return Err(FerroError::ConversionError {
                        msg: "Empty allele variant cannot be converted to VCF".to_string(),
                    });
                }

                let mut results = Vec::new();
                for (index, member) in allele.variants.iter().enumerate() {
                    let member_results =
                        self.convert_all(member)
                            .map_err(|err| FerroError::ConversionError {
                                msg: format!(
                                    "Could not convert allele member {} ({}): {}",
                                    index, member, err
                                ),
                            })?;
                    results.extend(member_results);
                }
                Ok(results)
            }
            other => Ok(vec![self.convert(other)?]),
        }
    }

    /// Convert a CDS variant to VCF
    fn convert_cds(&self, variant: &CdsVariant) -> Result<VcfConversionResult, FerroError> {
        let hgvs_string = variant.to_string();
        let mut warnings = Vec::new();

        // Get the CDS interval
        let interval = &variant.loc_edit.location;
        let edit = variant
            .loc_edit
            .edit
            .inner()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "CDS variant with unknown edit cannot be converted to VCF: {}",
                    hgvs_string
                ),
            })?;

        // Map CDS positions to genomic positions
        let start_cds = interval
            .start
            .inner()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "CDS variant with unknown start position cannot be converted to VCF: {}",
                    hgvs_string
                ),
            })?;
        let end_cds = interval
            .end
            .inner()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "CDS variant with unknown end position cannot be converted to VCF: {}",
                    hgvs_string
                ),
            })?;

        let genomic_start = self.cds_to_genomic(start_cds)?;
        let genomic_end = self.cds_to_genomic(end_cds)?;

        // Get chromosome from transcript
        let chrom = self
            .transcript
            .chromosome
            .as_ref()
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no chromosome information".to_string(),
            })?
            .clone();

        // Handle intronic positions
        if start_cds.is_intronic() || end_cds.is_intronic() {
            warnings.push("Variant includes intronic positions".to_string());
        }

        // Build VCF record from edit
        let record =
            self.build_vcf_record(&chrom, genomic_start, genomic_end, edit, &hgvs_string)?;

        Ok(VcfConversionResult {
            record,
            hgvs_string,
            warnings,
        })
    }

    /// Convert a transcript variant to VCF
    fn convert_tx(&self, variant: &TxVariant) -> Result<VcfConversionResult, FerroError> {
        let hgvs_string = variant.to_string();
        let warnings = Vec::new();

        let interval = &variant.loc_edit.location;
        let edit = variant
            .loc_edit
            .edit
            .inner()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "Transcript variant with unknown edit cannot be converted to VCF: {}",
                    hgvs_string
                ),
            })?;

        // Map transcript positions to genomic
        let start_tx = interval
            .start
            .inner()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "Transcript variant with unknown start position cannot be converted to VCF: {}",
                    hgvs_string
                ),
            })?;
        let end_tx = interval
            .end
            .inner()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "Transcript variant with unknown end position cannot be converted to VCF: {}",
                    hgvs_string
                ),
            })?;

        let genomic_start =
            self.mapper
                .tx_to_genomic(start_tx)?
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!("Could not map transcript position {} to genomic", start_tx),
                })?;
        let genomic_end =
            self.mapper
                .tx_to_genomic(end_tx)?
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!("Could not map transcript position {} to genomic", end_tx),
                })?;

        let chrom = self
            .transcript
            .chromosome
            .as_ref()
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no chromosome information".to_string(),
            })?
            .clone();

        let record =
            self.build_vcf_record(&chrom, genomic_start, genomic_end, edit, &hgvs_string)?;

        Ok(VcfConversionResult {
            record,
            hgvs_string,
            warnings,
        })
    }

    /// Convert an RNA (`r.`) variant to VCF.
    ///
    /// HGVS RNA numbering mirrors the matching DNA reference's numbering:
    /// - on a **coding** transcript ([`Transcript::is_coding`]), `r.` is
    ///   CDS-relative — identical to `c.` (negative base in the 5' UTR, `*N` in
    ///   the 3' UTR) — so each [`RnaPos`] lowers to the corresponding [`CdsPos`]
    ///   and reuses the same CDS→genomic mapping as [`Self::convert_cds`];
    /// - on a **non-coding** transcript, `r.` is transcript-relative —
    ///   identical to `n.` — so each [`RnaPos`] base lowers to a [`TxPos`] and
    ///   reuses the transcript→genomic mapping of [`Self::convert_tx`].
    ///
    /// Either way the genomic coordinates feed the shared
    /// [`Self::build_vcf_record`], so REF/ALT/anchor handling is identical to
    /// the c./n. paths.
    ///
    /// RNA-only shapes with no single-locus genomic representation are declined
    /// rather than mapped to wrong coordinates: whole-entity edits (`r.=`,
    /// `r.?`, `r.spl`, `r.0`) are caught by the [`NaEdit::is_whole_entity`]
    /// guard before any coordinate mapping (`r.=` in particular would otherwise
    /// emit a spurious no-op record at the parser's placeholder position), and
    /// an unknown edit (`Mu::Unknown`) is declined instead of panicking.
    fn convert_rna(&self, variant: &RnaVariant) -> Result<VcfConversionResult, FerroError> {
        let hgvs_string = variant.to_string();
        let mut warnings = Vec::new();

        // Decline unknown/whole-entity RNA edits up front. These carry no
        // genomic single-locus representation; mapping their (placeholder)
        // coordinates would emit a wrong record.
        let edit = variant
            .loc_edit
            .edit
            .inner()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "RNA variant with unknown edit (r.?) cannot be converted to VCF: {}",
                    hgvs_string
                ),
            })?;
        if edit.is_whole_entity() {
            return Err(FerroError::ConversionError {
                msg: format!(
                    "RNA whole-entity edit (r.=/r.?/r.spl/r.0) has no single-locus \
                     genomic VCF representation: {}",
                    hgvs_string
                ),
            });
        }

        let interval = &variant.loc_edit.location;
        let start_rna = interval
            .start
            .inner()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "RNA variant with unknown start position cannot be converted to VCF: {}",
                    hgvs_string
                ),
            })?;
        let end_rna = interval
            .end
            .inner()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "RNA variant with unknown end position cannot be converted to VCF: {}",
                    hgvs_string
                ),
            })?;

        // Map RNA positions to genomic, choosing the numbering frame by whether
        // the transcript is coding (CDS-relative, like c.) or non-coding
        // (transcript-relative, like n.).
        let (genomic_start, genomic_end) = if self.transcript.is_coding() {
            (
                self.cds_to_genomic(&rna_pos_to_cds(start_rna))?,
                self.cds_to_genomic(&rna_pos_to_cds(end_rna))?,
            )
        } else {
            let map_tx = |pos: &RnaPos| -> Result<u64, FerroError> {
                // A non-coding transcript has no CDS, so it has no `*N` 3'-UTR
                // frame and no `-N` 5'-UTR frame: `r.` numbering is plain
                // transcript-relative (== n.). A `utr3`/negative position here
                // is meaningless and `rna_pos_to_tx` would silently drop the
                // `utr3` flag, mapping `*N` as `N`. Decline instead, matching
                // the non-coding guard in `normalize::rna_pos_to_txpos`.
                if pos.utr3 || pos.base < 1 {
                    return Err(FerroError::ConversionError {
                        msg: format!(
                            "RNA UTR position {} has no meaning on a non-coding transcript \
                             (no CDS frame); cannot convert to VCF",
                            pos
                        ),
                    });
                }
                self.mapper
                    .tx_to_genomic(&rna_pos_to_tx(pos))?
                    .ok_or_else(|| FerroError::ConversionError {
                        msg: format!("Could not map RNA position {} to genomic", pos),
                    })
            };
            (map_tx(start_rna)?, map_tx(end_rna)?)
        };

        let chrom = self
            .transcript
            .chromosome
            .as_ref()
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no chromosome information".to_string(),
            })?
            .clone();

        if start_rna.is_intronic() || end_rna.is_intronic() {
            warnings.push("Variant includes intronic positions".to_string());
        }

        // The VCF record lives on the genomic (DNA) axis, so any RNA `u`/`U`
        // bases the edit carries (substitution ref/alt, inserted/deleted
        // sequences, repeat units) must be translated to `T` before they reach
        // `build_vcf_record` — otherwise they would be emitted verbatim as the
        // invalid DNA letter `U`. This mirrors the SPDI path's `apply_alphabet`
        // (`'U' => 'T'`). Bases recovered from the genomic provider are already
        // DNA and are untouched.
        let dna_edit = u_to_t_edit(edit);

        let record =
            self.build_vcf_record(&chrom, genomic_start, genomic_end, &dna_edit, &hgvs_string)?;

        Ok(VcfConversionResult {
            record,
            hgvs_string,
            warnings,
        })
    }

    /// Convert a genomic variant to VCF
    fn convert_genome(&self, variant: &GenomeVariant) -> Result<VcfConversionResult, FerroError> {
        let hgvs_string = variant.to_string();
        let warnings = Vec::new();

        let interval = &variant.loc_edit.location;
        reject_unresolvable_genomic_position(interval, &hgvs_string)?;
        let edit = variant
            .loc_edit
            .edit
            .inner()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "Genomic variant with unknown edit cannot be converted to VCF: {}",
                    hgvs_string
                ),
            })?;

        let genomic_start = interval
            .start
            .inner()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "Genomic variant with unknown start position cannot be converted to VCF: {}",
                    hgvs_string
                ),
            })?
            .base;
        let genomic_end = interval
            .end
            .inner()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "Genomic variant with unknown end position cannot be converted to VCF: {}",
                    hgvs_string
                ),
            })?
            .base;

        // For genomic variants, we need to derive chromosome from accession
        // NC_000001.11 -> chr1, etc.
        let chrom = accession_to_chromosome(&variant.accession.to_string())
            .unwrap_or_else(|| variant.accession.number.to_string());

        let record =
            self.build_vcf_record(&chrom, genomic_start, genomic_end, edit, &hgvs_string)?;

        Ok(VcfConversionResult {
            record,
            hgvs_string,
            warnings,
        })
    }

    /// Build a VCF record from genomic coordinates and edit
    fn build_vcf_record(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        edit: &NaEdit,
        hgvs_string: &str,
    ) -> Result<VcfRecord, FerroError> {
        let (pos, reference, alternate) = match edit {
            NaEdit::Substitution {
                reference: ref_base,
                alternative: alt_base,
            } => {
                // SNV: simple one-base change
                (
                    start,
                    ref_base.to_char().to_string(),
                    vec![alt_base.to_char().to_string()],
                )
            }
            NaEdit::SubstitutionNoRef {
                alternative: alt_base,
            } => {
                // Substitution without reference - need to fetch reference from provider
                let ref_base = self.get_reference_base(chrom, start)?;
                (
                    start,
                    ref_base.to_string(),
                    vec![alt_base.to_char().to_string()],
                )
            }
            NaEdit::Deletion { sequence, .. } => {
                // Deletion: need anchor base before the deletion
                // VCF format: REF includes anchor + deleted bases, ALT is just anchor
                let anchor_pos = anchor_position_before(start, hgvs_string)?;
                let anchor_base = self.get_reference_base(chrom, anchor_pos)?;

                let deleted_seq = if let Some(seq) = sequence {
                    seq.to_string()
                } else {
                    // Need to fetch from reference
                    self.get_reference_sequence(chrom, start, end)?
                };

                let ref_allele = format!("{}{}", anchor_base, deleted_seq);
                let alt_allele = anchor_base.to_string();

                (anchor_pos, ref_allele, vec![alt_allele])
            }
            NaEdit::Insertion { sequence } => {
                // Insertion: anchor base at position, ALT includes anchor + inserted sequence
                let anchor_pos = start;
                let anchor_base = self.get_reference_base(chrom, anchor_pos)?;
                let inserted_seq = inserted_sequence_bases(sequence, hgvs_string)?;

                let ref_allele = anchor_base.to_string();
                let alt_allele = format!("{}{}", anchor_base, inserted_seq);

                (anchor_pos, ref_allele, vec![alt_allele])
            }
            NaEdit::BreakpointInsertion { .. } => {
                return Err(FerroError::ConversionError {
                    msg: "Positionless breakpoint insertions cannot be converted to VCF (no anchor coordinate)"
                        .to_string(),
                });
            }
            NaEdit::Delins { sequence, .. } => {
                // Deletion-insertion: REF is deleted sequence, ALT is replacement
                let deleted_seq = self.get_reference_sequence(chrom, start, end)?;
                let alt_seq = inserted_sequence_bases(sequence, hgvs_string)?;

                // If lengths match, no anchor needed
                // If different lengths, might need anchor for proper VCF normalization
                if deleted_seq.len() == alt_seq.len() {
                    (start, deleted_seq, vec![alt_seq])
                } else {
                    // Add anchor for length-changing variants
                    let anchor_pos = anchor_position_before(start, hgvs_string)?;
                    let anchor_base = self.get_reference_base(chrom, anchor_pos)?;
                    let ref_allele = format!("{}{}", anchor_base, deleted_seq);
                    let alt_allele = format!("{}{}", anchor_base, alt_seq);
                    (anchor_pos, ref_allele, vec![alt_allele])
                }
            }
            NaEdit::Duplication { sequence, .. } => {
                // Duplication is an insertion of the duplicated sequence
                let dup_seq = if let Some(seq) = sequence {
                    seq.to_string()
                } else {
                    self.get_reference_sequence(chrom, start, end)?
                };

                // Insert at end of duplicated region
                let anchor_pos = end;
                let anchor_base = self.get_reference_base(chrom, anchor_pos)?;

                let ref_allele = anchor_base.to_string();
                let alt_allele = format!("{}{}", anchor_base, dup_seq);

                (anchor_pos, ref_allele, vec![alt_allele])
            }
            NaEdit::DupIns { sequence } => {
                // Duplication-insertion: duplication of region followed by insertion
                // VCF representation: REF is anchor base, ALT is anchor + dup_seq + inserted_seq
                let dup_seq = self.get_reference_sequence(chrom, start, end)?;
                let ins_seq = inserted_sequence_bases(sequence, hgvs_string)?;

                let anchor_pos = end;
                let anchor_base = self.get_reference_base(chrom, anchor_pos)?;

                let ref_allele = anchor_base.to_string();
                let alt_allele = format!("{}{}{}", anchor_base, dup_seq, ins_seq);

                (anchor_pos, ref_allele, vec![alt_allele])
            }
            NaEdit::Inversion { .. } => {
                // Inversion: REF is original, ALT is reverse complement
                let original_seq = self.get_reference_sequence(chrom, start, end)?;
                let inverted_seq = reverse_complement(&original_seq);
                (start, original_seq, vec![inverted_seq])
            }
            NaEdit::Repeat {
                sequence,
                count,
                additional_counts,
                trailing,
            } => {
                use crate::hgvs::edit::RepeatCount;

                // Helper to extract count from RepeatCount
                let get_count = |rc: &RepeatCount| -> Result<u64, FerroError> {
                    match rc {
                        RepeatCount::Exact(n) => Ok(*n),
                        RepeatCount::Range(min, _)
                        | RepeatCount::UncertainRange(min, _)
                        | RepeatCount::MinUncertain(min) => Ok(*min),
                        RepeatCount::MaxUncertain(_) | RepeatCount::Unknown => {
                            Err(FerroError::ConversionError {
                                msg: "Repeat variant with unknown count cannot be converted to VCF"
                                    .to_string(),
                            })
                        }
                    }
                };

                // Get repeat unit sequence (including any trailing bases)
                let repeat_unit = if let Some(seq) = sequence {
                    let mut unit = seq.to_string();
                    if let Some(trail) = trailing {
                        unit.push_str(&trail.to_string());
                    }
                    unit
                } else if let Some(trail) = trailing {
                    trail.to_string()
                } else {
                    return Err(FerroError::ConversionError {
                        msg: "Repeat variant without sequence cannot be converted to VCF"
                            .to_string(),
                    });
                };

                // Get all repeat counts (primary + any additional for genotype notation)
                let primary_count = get_count(count)?;
                let mut all_counts = vec![primary_count];
                for additional in additional_counts {
                    all_counts.push(get_count(additional)?);
                }

                // Repeat expansion is an insertion at the given position
                let anchor_pos = start;
                let anchor_base = self.get_reference_base(chrom, anchor_pos)?;
                let ref_allele = anchor_base.to_string();

                // Generate alt alleles for each count (deduplicating)
                let mut alt_alleles: Vec<String> = all_counts
                    .iter()
                    .map(|&n| {
                        let expanded_seq = repeat_unit.repeat(n as usize);
                        format!("{}{}", anchor_base, expanded_seq)
                    })
                    .collect();
                alt_alleles.sort();
                alt_alleles.dedup();

                (anchor_pos, ref_allele, alt_alleles)
            }
            NaEdit::Identity { .. } => {
                // Identity means no change - return a record with same REF and ALT
                let ref_base = self.get_reference_base(chrom, start)?;
                (start, ref_base.to_string(), vec![ref_base.to_string()])
            }
            NaEdit::NPaddedDeletion { .. } => {
                return Err(FerroError::ConversionError {
                    msg: "N-padded deletions over an uncertain range cannot be converted to VCF"
                        .to_string(),
                });
            }
            NaEdit::Conversion { .. } => {
                return Err(FerroError::ConversionError {
                    msg: "Conversion variants cannot be directly converted to VCF".to_string(),
                });
            }
            NaEdit::Unknown { .. } => {
                return Err(FerroError::ConversionError {
                    msg: "Unknown variants cannot be converted to VCF".to_string(),
                });
            }
            NaEdit::Methylation { .. } => {
                return Err(FerroError::ConversionError {
                    msg:
                        "Methylation variants cannot be converted to VCF (epigenetic, not sequence)"
                            .to_string(),
                });
            }
            NaEdit::CopyNumber { count } => {
                // Copy number variant: use VCF structural variant notation
                // REF: single base at position
                // ALT: <CNV> symbolic allele
                // INFO: SVTYPE=CNV, CN=count, END=end
                let ref_base = self.get_reference_base(chrom, start)?;

                let mut record = VcfRecord::new(
                    chrom.to_string(),
                    start,
                    ref_base.to_string(),
                    vec!["<CNV>".to_string()],
                )
                .with_genome_build(self.genome_build);

                // Add CNV-specific INFO fields
                record.info.insert(
                    "SVTYPE".to_string(),
                    super::record::InfoValue::String("CNV".to_string()),
                );
                record.info.insert(
                    "CN".to_string(),
                    super::record::InfoValue::Integer(*count as i64),
                );
                record.info.insert(
                    "END".to_string(),
                    super::record::InfoValue::Integer(end as i64),
                );
                record.info.insert(
                    "HGVS".to_string(),
                    super::record::InfoValue::String(hgvs_string.to_string()),
                );

                return Ok(record);
            }
            NaEdit::Splice { .. } => {
                return Err(FerroError::ConversionError {
                    msg: "Splice variants (r.spl/r.spl?) cannot be converted to VCF (RNA-specific effect)"
                        .to_string(),
                });
            }
            NaEdit::NoProduct => {
                return Err(FerroError::ConversionError {
                    msg:
                        "No product variants (r.0) cannot be converted to VCF (RNA-specific effect)"
                            .to_string(),
                });
            }
            NaEdit::MultiRepeat { .. } => {
                return Err(FerroError::ConversionError {
                    msg: "Multi-repeat variants cannot be directly converted to VCF yet"
                        .to_string(),
                });
            }
            NaEdit::PositionOnly => {
                return Err(FerroError::ConversionError {
                    msg: "Position-only variants cannot be converted to VCF (no edit specified)"
                        .to_string(),
                });
            }
        };

        let mut record = VcfRecord::new(chrom.to_string(), pos, reference, alternate)
            .with_genome_build(self.genome_build);

        // Add HGVS annotation to INFO
        record.info.insert(
            "HGVS".to_string(),
            super::record::InfoValue::String(hgvs_string.to_string()),
        );

        Ok(record)
    }

    /// Map CDS position to genomic position
    fn cds_to_genomic(&self, cds_pos: &CdsPos) -> Result<u64, FerroError> {
        self.mapper
            .cds_to_genomic(cds_pos)?
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!("Could not map CDS position {} to genomic", cds_pos),
            })
    }

    /// Get a single reference base at a 1-based genomic position on `chrom`.
    fn get_reference_base(&self, chrom: &str, pos: u64) -> Result<char, FerroError> {
        let bases = self.get_reference_sequence(chrom, pos, pos)?;
        bases
            .chars()
            .next()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!("Could not get reference base at {}:{}", chrom, pos),
            })
    }

    /// Fetch reference bases over a 1-based inclusive genomic range `[start,
    /// end]` on contig `chrom`.
    ///
    /// All callers in [`Self::build_vcf_record`] supply genomic coordinates
    /// (the CDS/transcript paths map to genomic before calling), so the bases
    /// are fetched from the provider's genomic contig — keyed by the same
    /// `chrom` name used for the output VCF record — rather than from the
    /// transcript-relative sequence (which is in a different coordinate frame
    /// and would yield the wrong base for genomic positions).
    ///
    /// Coordinates are validated (`start >= 1`, `end >= start`) before the
    /// 1-based → 0-based half-open conversion, and the returned slice is
    /// length-checked so a short provider read errors rather than silently
    /// truncating the reference allele. When the provider has no genomic data
    /// for `chrom`, the error propagates so del/ins anchor recovery declines
    /// cleanly instead of emitting a wrong base.
    fn get_reference_sequence(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        if start < 1 || end < start {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "invalid 1-based genomic range [{}, {}] for {}",
                    start, end, chrom
                ),
            });
        }

        // 1-based inclusive [start, end] → 0-based half-open [start - 1, end).
        let zb_start = start - 1;
        let zb_end = end;
        let expected_len = (zb_end - zb_start) as usize;

        // The genomic accession is the contig name; prefer the dedicated
        // genomic path and fall back to `get_sequence` for providers that
        // store contigs under the generic sequence map. Only the "no genomic
        // data" decline (`GenomicReferenceNotAvailable`) triggers the fallback;
        // any other genomic error (e.g. an out-of-range fetch from a provider
        // that *does* serve this contig) is propagated as-is rather than masked
        // by whatever the generic path returns.
        let bases = match self.provider.get_genomic_sequence(chrom, zb_start, zb_end) {
            Ok(seq) => seq,
            Err(FerroError::GenomicReferenceNotAvailable { .. }) => {
                self.provider.get_sequence(chrom, zb_start, zb_end)?
            }
            Err(other) => return Err(other),
        };

        if bases.len() != expected_len {
            return Err(FerroError::ConversionError {
                msg: format!(
                    "reference fetch for {}:{}-{} returned {} bases, expected {}",
                    chrom,
                    start,
                    end,
                    bases.len(),
                    expected_len
                ),
            });
        }

        Ok(bases)
    }
}

/// Convert a genomic HGVS variant directly to VCF (no transcript needed).
///
/// This standalone path has no [`ReferenceProvider`], so it cannot recover the
/// VCF anchor base that deletions and insertions require — those edits decline
/// with a clear error. To convert deletions/insertions, use
/// [`HgvsToVcfConverter`], which fetches the anchor base from the provider.
///
/// Like the provider-backed path, a `+`/`-` offset or a `pter`/`qter`/`cen`
/// special position on either endpoint is refused rather than read through to
/// `base` — see `reject_unresolvable_genomic_position`. Both entry points are
/// guarded, because both read the same field of the same struct: unguarded,
/// this one lowered `g.10+2C>T` and `g.10C>T` to the same `POS=10` row.
pub fn genomic_hgvs_to_vcf(variant: &GenomeVariant) -> Result<VcfRecord, FerroError> {
    let interval = &variant.loc_edit.location;
    reject_unresolvable_genomic_position(interval, &variant.to_string())?;
    let edit = variant
        .loc_edit
        .edit
        .inner()
        .ok_or_else(|| FerroError::ConversionError {
            msg: format!(
                "Genomic variant with unknown edit cannot be converted to VCF: {}",
                variant
            ),
        })?;

    let start = interval
        .start
        .inner()
        .ok_or_else(|| FerroError::ConversionError {
            msg: format!(
                "Genomic variant with unknown start position cannot be converted to VCF: {}",
                variant
            ),
        })?
        .base;
    let _end = interval
        .end
        .inner()
        .ok_or_else(|| FerroError::ConversionError {
            msg: format!(
                "Genomic variant with unknown end position cannot be converted to VCF: {}",
                variant
            ),
        })?
        .base;

    let chrom = accession_to_chromosome(&variant.accession.to_string())
        .unwrap_or_else(|| variant.accession.number.to_string());

    let (pos, reference, alternate) = match edit {
        NaEdit::Substitution {
            reference: ref_base,
            alternative: alt_base,
        } => (
            start,
            ref_base.to_char().to_string(),
            vec![alt_base.to_char().to_string()],
        ),
        NaEdit::Deletion { sequence, .. } => {
            // VCF deletions need the preceding (anchor) base, which requires
            // reference access this provider-less path does not have. Surface a
            // distinct message when the deleted sequence is also unknown.
            let deleted =
                sequence.as_ref().ok_or_else(|| {
                    FerroError::ConversionError {
                msg: "Cannot convert deletion to VCF without deleted sequence (no reference data)"
                    .to_string(),
            }
                })?;
            return Err(FerroError::ConversionError {
                msg: format!(
                    "Cannot convert deletion '{}' to VCF: anchor base required (no reference data)",
                    deleted
                ),
            });
        }
        NaEdit::Insertion { sequence } => {
            // VCF insertions need the anchor base; without reference access we
            // cannot recover it.
            return Err(FerroError::ConversionError {
                msg: format!(
                    "Cannot convert insertion '{}' to VCF: anchor base required (no reference data)",
                    sequence
                ),
            });
        }
        _ => {
            return Err(FerroError::ConversionError {
                msg: format!("Edit type {:?} not supported for simple conversion", edit),
            });
        }
    };

    Ok(VcfRecord::new(chrom, pos, reference, alternate))
}

/// Map a RefSeq chromosome accession to its UCSC chromosome name.
///
/// The primary-assembly `NC_` → UCSC mapping is resolved through the shared,
/// build-aware [`ContigAliases`](crate::liftover::aliases::ContigAliases)
/// reverse table ([`refseq_to_ucsc`](crate::liftover::aliases::ContigAliases::refseq_to_ucsc))
/// rather than a hand-rolled accession ladder, so the single source of truth in
/// `liftover::aliases` governs both the forward and reverse directions. The
/// table is keyed by the fully versioned accession and carries both GRCh37 and
/// GRCh38 versions of chr1–22, X, Y, and M (e.g. `NC_012920.1` → `chrM`).
///
/// Accessions the table does not describe (alternate loci `NC_06*`, bacterial
/// `NZ_`, unplaced scaffolds `NT_`/`NW_`) are returned verbatim; anything else
/// declines with `None`.
fn accession_to_chromosome(accession: &str) -> Option<String> {
    use crate::liftover::aliases::default_human_aliases;

    if let Some(ucsc) = default_human_aliases().refseq_to_ucsc(accession) {
        return Some(ucsc.to_string());
    }

    // Non-primary contigs the shared table does not cover are passed through
    // by accession (without version) so callers retain a usable contig label.
    let base = accession.split('.').next().unwrap_or(accession);
    if base.starts_with("NC_06")
        || base.starts_with("NZ_")
        || base.starts_with("NT_")
        || base.starts_with("NW_")
    {
        return Some(base.to_string());
    }

    None
}

/// Lower an [`RnaPos`] to the equivalent [`CdsPos`] for a coding transcript.
///
/// HGVS `r.` numbering on a coding transcript is identical to `c.` numbering:
/// `base` is relative to the start codon (negative in the 5' UTR), `offset` is
/// the intronic offset, and `utr3` marks `*N` 3'-UTR positions. `RnaPos` has no
/// `special` (pter/qter/cen) component, so it maps to `special: None`.
fn rna_pos_to_cds(pos: &RnaPos) -> CdsPos {
    CdsPos {
        base: pos.base,
        offset: pos.offset,
        utr3: pos.utr3,
        special: None,
    }
}

/// Lower an [`RnaPos`] to the equivalent [`TxPos`] for a non-coding transcript.
///
/// On a non-coding transcript HGVS `r.` numbering is transcript-relative
/// (identical to `n.`): `base` is the transcript coordinate and `offset` is the
/// intronic offset. There is no 3'-UTR `*N` frame on a non-coding transcript,
/// so `utr3` cannot be represented in `TxPos` (`downstream` stays `false`).
/// Callers MUST reject `utr3`/negative-`base` positions before calling this —
/// `convert_rna`'s non-coding branch does so — otherwise the dropped `utr3`
/// flag would map `*N` as the unrelated transcript position `N`.
///
/// That obligation is an ordering, not a type, so it is pinned rather than
/// assumed: `the_dropped_utr3_marker_is_shielded_by_the_non_coding_guard`
/// asserts both that the marker is really lost here and that the sole call
/// site refuses ahead of it, so moving the guard — or adding a second caller —
/// fails a test instead of silently answering a different question.
fn rna_pos_to_tx(pos: &RnaPos) -> TxPos {
    TxPos {
        base: pos.base,
        offset: pos.offset,
        downstream: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::edit::{
        Base, InsertedSequence, MethylationStatus, ProteinEdit, RepeatCount, Sequence,
    };
    use crate::hgvs::interval::{
        CdsInterval, GenomeInterval, ProtInterval, RnaInterval, TxInterval,
    };
    use crate::hgvs::location::{CdsPos, GenomePos, ProtPos, RnaPos, TxPos};
    use crate::hgvs::variant::AllelePhase;
    use crate::hgvs::variant::{
        Accession, AlleleVariant, CircularVariant, LocEdit, MtVariant, ProteinVariant,
        RnaFusionBreakpoint, RnaFusionVariant, RnaVariant,
    };
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use crate::reference::MockProvider;
    use std::str::FromStr;
    use std::sync::OnceLock;

    fn create_test_transcript() -> Transcript {
        Transcript {
            cds_start_incomplete: false,
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            exons: vec![
                Exon::with_genomic(1, 1, 100, 1000, 1099),
                Exon::with_genomic(2, 101, 200, 2000, 2099),
            ],
            cds_start: Some(50),
            cds_end: Some(150),
            sequence: Some("ATGCATGCATGC".repeat(20)),
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2099),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATGC"), "GCAT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement(""), "");
    }

    #[test]
    fn test_reverse_complement_lowercase() {
        assert_eq!(reverse_complement("atgc"), "gcat");
        assert_eq!(reverse_complement("AaTt"), "aAtT");
    }

    #[test]
    fn test_reverse_complement_unknown_bases() {
        assert_eq!(reverse_complement("ATNC"), "GNAT");
    }

    #[test]
    fn test_accession_to_chromosome() {
        assert_eq!(
            accession_to_chromosome("NC_000001.11"),
            Some("chr1".to_string())
        );
        assert_eq!(
            accession_to_chromosome("NC_000022.11"),
            Some("chr22".to_string())
        );
        assert_eq!(
            accession_to_chromosome("NC_000023.11"),
            Some("chrX".to_string())
        );
        assert_eq!(
            accession_to_chromosome("NC_000024.10"),
            Some("chrY".to_string())
        );
        assert_eq!(
            accession_to_chromosome("NC_012920.1"),
            Some("chrM".to_string())
        );
        assert_eq!(accession_to_chromosome("NM_000088.3"), None);
    }

    #[test]
    fn test_accession_to_chromosome_invalid() {
        assert_eq!(accession_to_chromosome("NC_000025.11"), None);
        assert_eq!(accession_to_chromosome(""), None);
        assert_eq!(accession_to_chromosome("invalid"), None);
    }

    /// Build a converter whose provider carries genomic sequence under the
    /// UCSC contig name `chr1` (the key `build_vcf_record` fetches against),
    /// so anchor-base recovery and coordinate mapping can be exercised
    /// end-to-end. The first base is at 1-based genomic position 1.
    fn converter_with_genomic_chr1<'a>(
        transcript: &'a Transcript,
        provider: &'a mut MockProvider,
        sequence: &str,
    ) -> HgvsToVcfConverter<'a, MockProvider> {
        provider.add_genomic_sequence("chr1", sequence);
        HgvsToVcfConverter::new(transcript, provider)
    }

    #[test]
    fn test_get_reference_base_reads_genomic_contig() {
        // Genomic sequence "ACGTACGT" under chr1; 1-based position N returns
        // the base at 0-based index N-1 (off-by-one boundary check).
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, "ACGTACGT");

        assert_eq!(converter.get_reference_base("chr1", 1).unwrap(), 'A');
        assert_eq!(converter.get_reference_base("chr1", 2).unwrap(), 'C');
        assert_eq!(converter.get_reference_base("chr1", 4).unwrap(), 'T');
        assert_eq!(converter.get_reference_base("chr1", 8).unwrap(), 'T');
    }

    #[test]
    fn test_get_reference_base_rejects_position_zero() {
        // A 1-based position of 0 is invalid and must error, not silently map
        // to index 0 or underflow.
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, "ACGT");

        assert!(converter.get_reference_base("chr1", 0).is_err());
    }

    #[test]
    fn test_get_reference_sequence_validates_range() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, "ACGTACGT");

        // Correct 1-based inclusive range maps to the right bases.
        assert_eq!(
            converter.get_reference_sequence("chr1", 2, 4).unwrap(),
            "CGT"
        );
        // start < 1 is invalid.
        assert!(converter.get_reference_sequence("chr1", 0, 4).is_err());
        // end < start is invalid.
        assert!(converter.get_reference_sequence("chr1", 4, 2).is_err());
    }

    #[test]
    fn test_genomic_deletion_recovers_anchor_from_provider() {
        // chr1 = "ACGTACGT"; delete genomic 3..4 ("GT"). Anchor is the base
        // at genomic position 2 ("C"). REF = anchor + deleted, ALT = anchor,
        // POS = anchor position (start - 1).
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, "ACGTACGT");

        let variant = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::new(GenomePos::new(3), GenomePos::new(4)),
                NaEdit::Deletion {
                    sequence: Some(Sequence::from_str("GT").unwrap()),
                    length: None,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 2);
        assert_eq!(record.reference, "CGT");
        assert_eq!(record.alternate, vec!["C"]);
    }

    #[test]
    fn test_genomic_insertion_recovers_anchor_from_provider() {
        // chr1 = "ACGTACGT"; insert "TTT" after genomic position 3. Anchor is
        // the base at genomic position 3 ("G"). REF = anchor, ALT = anchor + ins.
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, "ACGTACGT");

        let variant = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(3)),
                NaEdit::Insertion {
                    sequence: InsertedSequence::Literal(Sequence::from_str("TTT").unwrap()),
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 3);
        assert_eq!(record.reference, "G");
        assert_eq!(record.alternate, vec!["GTTT"]);
    }

    #[test]
    fn test_genomic_deletion_declines_without_reference() {
        // No genomic sequence registered: anchor recovery must fail with a
        // clear error rather than emit a silently wrong base.
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::new(GenomePos::new(3), GenomePos::new(4)),
                NaEdit::Deletion {
                    sequence: Some(Sequence::from_str("GT").unwrap()),
                    length: None,
                },
            ),
        });

        assert!(converter.convert(&variant).is_err());
    }

    #[test]
    fn test_genomic_insertion_declines_without_reference() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(3)),
                NaEdit::Insertion {
                    sequence: InsertedSequence::Literal(Sequence::from_str("TTT").unwrap()),
                },
            ),
        });

        assert!(converter.convert(&variant).is_err());
    }

    #[test]
    fn test_accession_to_chromosome_mt_and_unmapped() {
        // MT resolves via the shared ContigAliases table.
        assert_eq!(
            accession_to_chromosome("NC_012920.1"),
            Some("chrM".to_string())
        );
        // An unmapped primary contig accession declines.
        assert_eq!(accession_to_chromosome("NC_000025.11"), None);
    }

    #[test]
    fn test_accession_to_chromosome_non_primary_passthrough() {
        // Non-primary contigs the shared table does not cover are returned
        // verbatim (without version) rather than mapped to a fabricated chrN.
        assert_eq!(
            accession_to_chromosome("NW_009646201.1"),
            Some("NW_009646201".to_string())
        );
        assert_eq!(
            accession_to_chromosome("NT_167244.2"),
            Some("NT_167244".to_string())
        );
        assert_eq!(
            accession_to_chromosome("NZ_CP011113.1"),
            Some("NZ_CP011113".to_string())
        );
        assert_eq!(
            accession_to_chromosome("NC_060925.1"),
            Some("NC_060925".to_string())
        );
    }

    #[test]
    fn test_genomic_hgvs_to_vcf_snv() {
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };

        let vcf = genomic_hgvs_to_vcf(&variant).unwrap();
        assert_eq!(vcf.chrom, "chr1");
        assert_eq!(vcf.pos, 12345);
        assert_eq!(vcf.reference, "A");
        assert_eq!(vcf.alternate, vec!["G"]);
    }

    #[test]
    fn test_genomic_hgvs_to_vcf_deletion_requires_anchor() {
        // Deletion without reference data for anchor base returns an error
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::new(GenomePos::new(12345), GenomePos::new(12347)),
                NaEdit::Deletion {
                    sequence: Some(Sequence::from_str("ATG").unwrap()),
                    length: None,
                },
            ),
        };

        let result = genomic_hgvs_to_vcf(&variant);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("anchor base"));
    }

    #[test]
    fn test_genomic_hgvs_to_vcf_insertion_requires_anchor() {
        // Insertion without reference data for anchor base returns an error
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Insertion {
                    sequence: InsertedSequence::Literal(Sequence::from_str("ATG").unwrap()),
                },
            ),
        };

        let result = genomic_hgvs_to_vcf(&variant);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("anchor base"));
    }

    #[test]
    fn test_genomic_hgvs_to_vcf_unsupported_edit() {
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::new(GenomePos::new(12345), GenomePos::new(12350)),
                NaEdit::Inversion {
                    sequence: None,
                    length: None,
                },
            ),
        };

        let result = genomic_hgvs_to_vcf(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_converter_cds_snv() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };

        let result = converter.convert(&HgvsVariant::Cds(variant)).unwrap();
        assert_eq!(result.record.chrom, "chr1");
        assert!(result.record.info.contains_key("HGVS"));
    }

    #[test]
    fn test_converter_with_genome_build() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter =
            HgvsToVcfConverter::new(&transcript, &provider).with_genome_build(GenomeBuild::GRCh37);

        let variant = CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };

        let result = converter.convert(&HgvsVariant::Cds(variant)).unwrap();
        assert_eq!(result.record.genome_build, GenomeBuild::GRCh37);
    }

    #[test]
    fn test_converter_protein_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Protein(ProteinVariant {
            accession: Accession::new("NP", "000079", Some(2)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                ProtInterval::point(ProtPos::new(crate::hgvs::location::AminoAcid::Met, 1)),
                ProteinEdit::Unknown {
                    predicted: false,
                    whole_protein: false,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(matches!(err, FerroError::ConversionError { .. }));
    }

    /// On a coding transcript, `r.` numbering is CDS-relative (== c.), so
    /// `r.1` maps through CDS 1 → tx 50 → genomic 1049 (cds_start=50, exon1 tx
    /// 1→genomic 1000). The substitution carries its own REF base.
    #[test]
    fn test_converter_rna_substitution_coding() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        // chr1 sequence long enough to cover genomic position 1049.
        let seq = "ACGT".repeat(400); // 1600 bases
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1049);
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["G"]);
    }

    /// SubstitutionNoRef reads the REF base from the genomic provider at the
    /// mapped position, exercising the same provider-read path as the cds/tx
    /// arms. genomic 1049 = 0-based index 1048; seq "ACGT" repeated → index
    /// 1048 % 4 == 0 → 'A'.
    #[test]
    fn test_converter_rna_substitution_no_ref_reads_provider() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(1)),
                NaEdit::SubstitutionNoRef {
                    alternative: Base::T,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1049);
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["T"]);
    }

    /// RNA deletion reuses anchor-base recovery: `r.2_3del` maps to genomic
    /// 1050..1051; the anchor is the base at 1049. REF = anchor + deleted, ALT
    /// = anchor, POS = anchor (start - 1).
    #[test]
    fn test_converter_rna_deletion_recovers_anchor() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::new(RnaPos::new(2), RnaPos::new(3)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        // CDS 2 → genomic 1050; CDS 3 → genomic 1051. Anchor at 1049.
        assert_eq!(record.pos, 1049);
        // anchor base (index 1048 % 4 == 0 → 'A') + deleted bases (1050,1051).
        // index 1049 % 4 == 1 → 'C'; index 1050 % 4 == 2 → 'G'. So deleted = "CG".
        assert_eq!(record.reference, "ACG");
        assert_eq!(record.alternate, vec!["A"]);
    }

    /// RNA insertion reuses anchor-base recovery: `r.2_3ins` anchors at the
    /// start position's genomic base. ALT = anchor + inserted.
    #[test]
    fn test_converter_rna_insertion_recovers_anchor() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::new(RnaPos::new(2), RnaPos::new(3)),
                NaEdit::Insertion {
                    sequence: InsertedSequence::Literal(Sequence::from_str("TTT").unwrap()),
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        // Insertion anchors at start (CDS 2 → genomic 1050; index 1049 → 'C').
        assert_eq!(record.pos, 1050);
        assert_eq!(record.reference, "C");
        assert_eq!(record.alternate, vec!["CTTT"]);
    }

    /// An RNA substitution carrying a uracil ALT (`r.1A>U`) must be lowered to
    /// the DNA alphabet: the emitted VCF ALT is `T`, never the invalid `U`.
    /// (REF `A` is already DNA-valid and unchanged.)
    #[test]
    fn test_converter_rna_substitution_u_to_t() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::U,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1049);
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["T"]);
    }

    /// An RNA insertion whose inserted sequence carries uracil
    /// (`r.2_3insUUU`) must be lowered to the DNA alphabet: the inserted bases
    /// appear as `T`, so ALT is `CTTT` rather than leaking `U`.
    #[test]
    fn test_converter_rna_insertion_u_to_t() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::new(RnaPos::new(2), RnaPos::new(3)),
                NaEdit::Insertion {
                    sequence: InsertedSequence::Literal(Sequence::from_str("UUU").unwrap()),
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1050);
        assert_eq!(record.reference, "C");
        assert_eq!(record.alternate, vec!["CTTT"]);
    }

    /// An RNA delins whose inserted sequence carries uracil (`r.2_3delinsUU`)
    /// must lower the inserted bases U→T while the provider-recovered deleted
    /// (REF) bases stay DNA. Genomic 1050..1051 deletes "CG" (len 2), inserts
    /// "UU"→"TT" (len 2); equal lengths so no anchor is added.
    #[test]
    fn test_converter_rna_delins_u_to_t() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::new(RnaPos::new(2), RnaPos::new(3)),
                NaEdit::Delins {
                    sequence: InsertedSequence::Literal(Sequence::from_str("UU").unwrap()),
                    deleted: None,
                    deleted_length: None,
                    substitution_reference: None,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1050);
        // Deleted REF bases come from the provider (DNA): indexes 1049,1050 → "CG".
        assert_eq!(record.reference, "CG");
        // Inserted ALT bases are lowered U→T.
        assert_eq!(record.alternate, vec!["TT"]);
    }

    /// `r.=` (whole-entity identity) has no single-locus genomic representation
    /// and must decline rather than emit a spurious no-op record at the
    /// parser's placeholder position. This is the regression
    /// `build_vcf_record`'s Identity arm would otherwise miss.
    #[test]
    fn test_converter_rna_whole_entity_identity_declines() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(1)),
                NaEdit::Identity {
                    sequence: None,
                    whole_entity: true,
                },
            ),
        });

        let err = converter.convert(&variant).unwrap_err();
        assert!(matches!(err, FerroError::ConversionError { .. }));
    }

    /// `r.spl` and `r.0` are RNA-only effects with no genomic representation.
    #[test]
    fn test_converter_rna_splice_and_no_product_decline() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        for edit in [NaEdit::Splice { unknown: false }, NaEdit::NoProduct] {
            let variant = HgvsVariant::Rna(RnaVariant {
                accession: Accession::new("NM", "000088", Some(3)),
                gene_symbol: None,
                loc_edit: LocEdit::new(RnaInterval::point(RnaPos::new(1)), edit),
            });
            assert!(converter.convert(&variant).is_err());
        }
    }

    /// An unknown RNA edit (`Mu::Unknown`, `r.?`) must decline cleanly, not
    /// panic on the `inner()` unwrap.
    #[test]
    fn test_converter_rna_unknown_edit_declines() {
        use crate::hgvs::uncertainty::Mu;
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit {
                location: RnaInterval::point(RnaPos::new(1)),
                edit: Mu::Unknown,
            },
        });

        let err = converter.convert(&variant).unwrap_err();
        assert!(matches!(err, FerroError::ConversionError { .. }));
    }

    /// An unknown CDS *position* (`c.?_123del`) must decline cleanly, not panic
    /// on the position `inner()` unwrap (#943 — the CDS/Tx/Genome siblings of
    /// the already-hardened RNA path). This is the exact input the CLI reaches
    /// with no upstream guard.
    #[test]
    fn test_converter_cds_unknown_start_position_declines() {
        use crate::hgvs::uncertainty::Mu;
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::with_uncertainty(Mu::Unknown, Mu::Certain(CdsPos::new(123))),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });

        let err = converter.convert(&variant).unwrap_err();
        assert!(matches!(err, FerroError::ConversionError { .. }));
    }

    /// An unknown CDS *edit* (`Mu::Unknown`) must decline, not panic.
    #[test]
    fn test_converter_cds_unknown_edit_declines() {
        use crate::hgvs::uncertainty::Mu;
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit {
                location: CdsInterval::point(CdsPos::new(1)),
                edit: Mu::Unknown,
            },
        });

        let err = converter.convert(&variant).unwrap_err();
        assert!(matches!(err, FerroError::ConversionError { .. }));
    }

    /// An unknown transcript (`n.`) position must decline, not panic (#943).
    #[test]
    fn test_converter_tx_unknown_start_position_declines() {
        use crate::hgvs::uncertainty::Mu;
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Tx(TxVariant {
            accession: Accession::new("NR", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                TxInterval::with_uncertainty(Mu::Unknown, Mu::Certain(TxPos::new(123))),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        });

        let err = converter.convert(&variant).unwrap_err();
        assert!(matches!(err, FerroError::ConversionError { .. }));
    }

    /// An unknown genomic (`g.`) position must decline, not panic (#943) — both
    /// via the converter and the standalone `genomic_hgvs_to_vcf`.
    #[test]
    fn test_converter_genome_unknown_position_declines() {
        use crate::hgvs::uncertainty::Mu;
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let genome = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::with_uncertainty(Mu::Unknown, Mu::Certain(GenomePos::new(123))),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };

        let err = converter
            .convert(&HgvsVariant::Genome(genome.clone()))
            .unwrap_err();
        assert!(matches!(err, FerroError::ConversionError { .. }));

        // The standalone entry point must also decline rather than panic.
        let err = genomic_hgvs_to_vcf(&genome).unwrap_err();
        assert!(matches!(err, FerroError::ConversionError { .. }));
    }

    /// On a non-coding transcript (no CDS), `r.` numbering is
    /// transcript-relative (== n.). `r.1` maps through tx 1 → genomic 1000.
    #[test]
    fn test_converter_rna_substitution_non_coding() {
        let mut transcript = create_test_transcript();
        transcript.cds_start = None;
        transcript.cds_end = None;
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NR", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1000);
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["G"]);
    }

    /// A non-coding transcript has no CDS frame, so an `r.*N` (3'-UTR) position
    /// is meaningless there. It must decline rather than silently drop the
    /// `utr3` flag and map `*N` as transcript position `N`.
    #[test]
    fn test_converter_rna_utr3_on_non_coding_declines() {
        let mut transcript = create_test_transcript();
        transcript.cds_start = None;
        transcript.cds_end = None;
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NR", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::utr3(5)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let err = converter.convert(&variant).unwrap_err();
        assert!(matches!(err, FerroError::ConversionError { .. }));
    }

    /// `rna_pos_to_tx` hard-codes `downstream: false`, so it drops the `utr3`
    /// (`*`) marker the AST carries — the same in-band-marker drop that
    /// `CoordinateMapper::tx_to_cds` is being fixed for on the `n.` axis. Here
    /// the site is genuinely **unreachable** with `utr3 == true`, so this test
    /// pins the guard that makes it so instead of changing behaviour: making
    /// the helper refuse would be dead code today, and a refusal nothing can
    /// reach is worse than none because it reads as coverage.
    ///
    /// "Handled one layer up" is only true while the layer stays put, which is
    /// what this test is for. It has three halves, and the first is what makes
    /// the other two worth having:
    ///
    /// 1. The drop is **real** — `r.*5` and `r.5` lower to the same `TxPos`,
    ///    so if control ever reached the helper the answer would be silently
    ///    wrong rather than absent.
    /// 2. The helper's **only** call site — `convert_rna`'s non-coding `map_tx`
    ///    closure — refuses first, and refuses *for that reason*: the message
    ///    is asserted, not just `is_err()`, so a decline arriving from some
    ///    other check cannot stand in for this guard.
    /// 3. The coding branch never reaches the helper at all — it lowers
    ///    through `rna_pos_to_cds`, which carries `utr3` faithfully — so a
    ///    coding transcript must not produce the non-coding guard's message.
    ///
    /// This is the same shape as
    /// `the_delins_anchor_is_shielded_by_the_reference_fetch`: an unreachable
    /// site whose unreachability depends on an ordering nothing else pins.
    /// If the guard is moved after the call, or a second caller appears, this
    /// is what notices.
    #[test]
    fn the_dropped_utr3_marker_is_shielded_by_the_non_coding_guard() {
        // 1. The helper really does lose the marker.
        let starred = rna_pos_to_tx(&RnaPos::utr3(5));
        let plain = rna_pos_to_tx(&RnaPos::new(5));
        assert!(
            !starred.downstream,
            "the helper is expected to drop the `*` marker; if it now carries it, \
             this test has been overtaken by a real fix and should be replaced by \
             one asserting the marker survives"
        );
        assert_eq!(
            starred, plain,
            "`r.*5` and `r.5` must be indistinguishable after lowering — that \
             indistinguishability is exactly why the guard below has to exist"
        );

        // 2. The only call site refuses first, and says why.
        let mut non_coding = create_test_transcript();
        non_coding.cds_start = None;
        non_coding.cds_end = None;
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&non_coding, &mut provider, &seq);

        for pos in [RnaPos::utr3(5), RnaPos::new(0), RnaPos::new(-3)] {
            let descriptor = format!("r.{pos}");
            let variant = HgvsVariant::Rna(RnaVariant {
                accession: Accession::new("NR", "000088", Some(3)),
                gene_symbol: None,
                loc_edit: LocEdit::new(
                    RnaInterval::point(pos),
                    NaEdit::Substitution {
                        reference: Base::A,
                        alternative: Base::G,
                    },
                ),
            });

            let message = converter
                .convert(&variant)
                .expect_err("a position with no meaning on a non-coding transcript must be refused")
                .to_string();
            assert!(
                message.contains("has no meaning on a non-coding transcript"),
                "`{descriptor}` was refused by some other check, so the guard that \
                 shields `rna_pos_to_tx` is no longer what stops it. Re-read \
                 `convert_rna`'s non-coding branch before trusting this decline. \
                 Got: {message}"
            );
        }

        // 3. The coding branch does not go through the helper at all.
        let coding = create_test_transcript();
        let mut coding_provider = MockProvider::new();
        let coding_converter =
            converter_with_genomic_chr1(&coding, &mut coding_provider, &"ACGT".repeat(700));
        let coding_variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::utr3(5)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        // Asserted unconditionally rather than under `if let Err(..)`: a bare
        // conditional check says nothing at all when the conversion succeeds, so
        // it could not distinguish "the coding branch works" from "the coding
        // branch stopped being reached". Both outcomes are pinned instead — the
        // conversion must succeed, and if it ever regresses to an error the
        // message must not be the non-coding guard's.
        match coding_converter.convert(&coding_variant) {
            Ok(result) => assert_eq!(
                result.hgvs_string, "NM_000088.3:r.*5a>g",
                "the coding lowering must convert `r.*5`, not merely avoid the \
                 non-coding guard"
            ),
            Err(err) => {
                let message = err.to_string();
                assert!(
                    !message.contains("has no meaning on a non-coding transcript"),
                    "a CODING transcript reached the non-coding guard, so `is_coding()` \
                     no longer selects the `rna_pos_to_cds` lowering: {message}"
                );
                panic!("the coding `r.*5` lowering must still convert: {message}");
            }
        }
    }

    #[test]
    fn test_converter_null_allele_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let result = converter.convert(&HgvsVariant::NullAllele);
        assert!(result.is_err());
    }

    #[test]
    fn test_converter_unknown_allele_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let result = converter.convert(&HgvsVariant::UnknownAllele);
        assert!(result.is_err());
    }

    #[test]
    fn test_converter_empty_allele_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Allele(AlleleVariant::new(vec![], AllelePhase::Cis));

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_converter_single_variant_allele() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let inner = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let variant = HgvsVariant::Allele(AlleleVariant::new(vec![inner], AllelePhase::Cis));

        let result = converter.convert(&variant);
        assert!(result.is_ok());
    }

    #[test]
    fn test_converter_complex_allele_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let inner1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let inner2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(10)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let variant =
            HgvsVariant::Allele(AlleleVariant::new(vec![inner1, inner2], AllelePhase::Cis));

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    /// Build a coding-transcript CDS substitution member (carries its own REF
    /// base, so no provider read is needed). `c.1` → genomic 1049, `c.2` →
    /// 1050, etc.
    fn cds_sub_member(cds_base: i64, ref_b: Base, alt_b: Base) -> HgvsVariant {
        HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(cds_base)),
                NaEdit::Substitution {
                    reference: ref_b,
                    alternative: alt_b,
                },
            ),
        })
    }

    /// `convert_all` decomposes a 2-member allele into two records, one per
    /// member, at the members' respective positions.
    #[test]
    fn test_convert_all_two_member_allele() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Allele(AlleleVariant::new(
            vec![
                cds_sub_member(1, Base::A, Base::G),
                cds_sub_member(2, Base::C, Base::T),
            ],
            AllelePhase::Cis,
        ));

        let results = converter.convert_all(&variant).unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].record.pos, 1049); // c.1 → genomic 1049
        assert_eq!(results[0].record.alternate, vec!["G"]);
        assert_eq!(results[1].record.pos, 1050); // c.2 → genomic 1050
        assert_eq!(results[1].record.alternate, vec!["T"]);
    }

    /// `convert_all` decomposes a 3-member allele into three records.
    #[test]
    fn test_convert_all_three_member_allele() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Allele(AlleleVariant::new(
            vec![
                cds_sub_member(1, Base::A, Base::G),
                cds_sub_member(2, Base::C, Base::T),
                cds_sub_member(3, Base::G, Base::A),
            ],
            AllelePhase::Cis,
        ));

        let results = converter.convert_all(&variant).unwrap();
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].record.pos, 1049);
        assert_eq!(results[1].record.pos, 1050);
        assert_eq!(results[2].record.pos, 1051);
    }

    /// `convert_all` on a single-member allele returns one record (and the
    /// `convert` single-record path remains unregressed — covered separately by
    /// `test_converter_single_variant_allele`).
    #[test]
    fn test_convert_all_single_member_allele() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Allele(AlleleVariant::new(
            vec![cds_sub_member(1, Base::A, Base::G)],
            AllelePhase::Cis,
        ));

        let results = converter.convert_all(&variant).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].record.pos, 1049);
    }

    /// `convert_all` on a scalar (non-allele) variant returns a 1-element vec.
    #[test]
    fn test_convert_all_scalar_variant() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = cds_sub_member(1, Base::A, Base::G);
        let results = converter.convert_all(&variant).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].record.pos, 1049);
    }

    /// `convert_all` recursively flattens nested alleles into a flat record
    /// list in declaration order.
    #[test]
    fn test_convert_all_flattens_nested_allele() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let nested = HgvsVariant::Allele(AlleleVariant::new(
            vec![
                cds_sub_member(1, Base::A, Base::G),
                cds_sub_member(2, Base::C, Base::T),
            ],
            AllelePhase::Cis,
        ));
        let outer = HgvsVariant::Allele(AlleleVariant::new(
            vec![nested, cds_sub_member(3, Base::G, Base::A)],
            AllelePhase::Cis,
        ));

        let results = converter.convert_all(&outer).unwrap();
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].record.pos, 1049);
        assert_eq!(results[1].record.pos, 1050);
        assert_eq!(results[2].record.pos, 1051);
    }

    /// A member whose shape cannot be lowered to VCF (here `[0]` / NullAllele)
    /// fails the whole `convert_all` call with a clear per-member error rather
    /// than silently dropping it.
    #[test]
    fn test_convert_all_unsupported_member_errors() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Allele(AlleleVariant::new(
            vec![cds_sub_member(1, Base::A, Base::G), HgvsVariant::NullAllele],
            AllelePhase::Cis,
        ));

        let err = converter.convert_all(&variant).unwrap_err();
        match err {
            FerroError::ConversionError { msg } => {
                // Names the offending member index (1) so the failure is
                // diagnosable, not silent.
                assert!(msg.contains("allele member 1"), "got: {msg}");
            }
            other => panic!("expected ConversionError, got {other:?}"),
        }
    }

    /// `convert_all` rejects an empty allele.
    #[test]
    fn test_convert_all_empty_allele_errors() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Allele(AlleleVariant::new(vec![], AllelePhase::Cis));
        assert!(converter.convert_all(&variant).is_err());
    }

    #[test]
    fn test_converter_mt_variant() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Mt(MtVariant {
            accession: Accession::new("NC", "012920", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(3243)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().record.chrom, "chrM");
    }

    #[test]
    fn test_converter_circular_variant() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Circular(CircularVariant {
            accession: Accession::new("NC", "012920", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_ok());
    }

    #[test]
    fn test_converter_rna_fusion_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let five_prime = RnaFusionBreakpoint {
            accession: Accession::new("NM", "000001", Some(1)),
            gene_symbol: None,
            interval: RnaInterval::point(RnaPos::new(100)),
        };
        let three_prime = RnaFusionBreakpoint {
            accession: Accession::new("NM", "000002", Some(1)),
            gene_symbol: None,
            interval: RnaInterval::point(RnaPos::new(200)),
        };

        let variant = HgvsVariant::RnaFusion(RnaFusionVariant::new(five_prime, three_prime));

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_converter_tx_variant() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Tx(TxVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(50)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().record.chrom, "chr1");
    }

    #[test]
    fn test_converter_genome_variant() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_ok());
        let vcf = result.unwrap();
        assert_eq!(vcf.record.chrom, "chr1");
        assert_eq!(vcf.record.pos, 12345);
    }

    #[test]
    fn test_vcf_conversion_result_fields() {
        let result = VcfConversionResult {
            record: VcfRecord::new(
                "chr1".to_string(),
                100,
                "A".to_string(),
                vec!["G".to_string()],
            ),
            hgvs_string: "NC_000001.11:g.100A>G".to_string(),
            warnings: vec!["test warning".to_string()],
        };

        assert_eq!(result.hgvs_string, "NC_000001.11:g.100A>G");
        assert_eq!(result.warnings.len(), 1);
        assert_eq!(result.record.chrom, "chr1");
    }

    /// `c.1_3delinsTTT`: CDS 1..3 → genomic 1049..1051. The deleted span and the
    /// replacement are both 3 bases, so no anchor base is added — REF is the
    /// deleted genomic bases, ALT is the literal replacement, POS is the start.
    #[test]
    fn test_build_vcf_record_delins_same_length() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(3)),
                NaEdit::Delins {
                    sequence: InsertedSequence::Literal(Sequence::from_str("TTT").unwrap()),
                    deleted: None,
                    deleted_length: None,
                    substitution_reference: None,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        // genomic 1049..1051 = indices 1048..1050 of "ACGT".repeat → "ACG".
        assert_eq!(record.pos, 1049);
        assert_eq!(record.reference, "ACG");
        assert_eq!(record.alternate, vec!["TTT"]);
    }

    /// `c.1_3dup` of "ATG": a duplication is an insertion of the duplicated
    /// sequence anchored on the base at the end of the duplicated region (CDS 3
    /// → genomic 1051). REF is that anchor base, ALT is anchor + duplicated unit.
    #[test]
    fn test_build_vcf_record_duplication() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(3)),
                NaEdit::Duplication {
                    sequence: Some(Sequence::from_str("ATG").unwrap()),
                    length: None,
                    uncertain_extent: None,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        // anchor = genomic 1051 = index 1050 of "ACGT".repeat → 'G'.
        assert_eq!(record.pos, 1051);
        assert_eq!(record.reference, "G");
        assert_eq!(record.alternate, vec!["GATG"]);
    }

    /// `c.1_4inv`: CDS 1..4 → genomic 1049..1052. REF is the original genomic
    /// span, ALT is its reverse complement. A non-palindromic seed is used so
    /// REF != ALT (a palindrome would make the assertion vacuous).
    #[test]
    fn test_build_vcf_record_inversion() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        // "AACCGGTT" tiled: genomic 1049..1052 = indices 1048..1051 = "AACC".
        let seq = "AACCGGTT".repeat(200);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(4)),
                NaEdit::Inversion {
                    sequence: None,
                    length: None,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1049);
        assert_eq!(record.reference, "AACC");
        // reverse_complement("AACC") == "GGTT".
        assert_eq!(record.alternate, vec!["GGTT"]);
    }

    /// An inversion over a span containing an IUPAC ambiguity code must complement
    /// that code, not emit it verbatim into the ALT field (#1347).
    ///
    /// `R` (A or G) complements to `Y` (C or T). Passing it through unchanged
    /// claims the ambiguity is self-complementary, which is wrong for every code
    /// except `S`, `W` and `N` — and it writes that claim into produced VCF.
    #[test]
    fn inversion_alt_complements_iupac_ambiguity_codes() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        // Same tiling as `test_build_vcf_record_inversion`, with the third base
        // replaced by `R`: genomic 1049..1052 = "AARC".
        let seq = "AARCGGTT".repeat(200);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(4)),
                NaEdit::Inversion {
                    sequence: None,
                    length: None,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.reference, "AARC");
        // reverse_complement("AARC") == "GYTT": R -> Y, not R -> R.
        assert_eq!(record.alternate, vec!["GYTT"]);
    }

    /// `c.1 CAG[5]`: a repeat expansion is an insertion anchored on the base at
    /// the start position (CDS 1 → genomic 1049). REF is the anchor base, ALT is
    /// anchor + the unit repeated `count` times.
    #[test]
    fn test_build_vcf_record_repeat() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::Exact(5),
                    additional_counts: vec![],
                    trailing: None,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        // anchor = genomic 1049 = index 1048 of "ACGT".repeat → 'A'.
        assert_eq!(record.pos, 1049);
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec![format!("A{}", "CAG".repeat(5))]);
    }

    /// `c.1 CAG[5];[10]`: additional counts emit one ALT per count, sorted and
    /// deduped. The x5 ALT is a string prefix of the x10 ALT, so it sorts first.
    #[test]
    fn test_build_vcf_record_repeat_with_additional_counts() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::Exact(5),
                    additional_counts: vec![RepeatCount::Exact(10)],
                    trailing: None,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1049);
        assert_eq!(record.reference, "A");
        assert_eq!(
            record.alternate,
            vec![
                format!("A{}", "CAG".repeat(5)),
                format!("A{}", "CAG".repeat(10)),
            ]
        );
    }

    #[test]
    fn test_build_vcf_record_repeat_unknown_count() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::Unknown,
                    additional_counts: vec![],
                    trailing: None,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_vcf_record_repeat_no_sequence() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: None,
                    count: RepeatCount::Exact(5),
                    additional_counts: vec![],
                    trailing: None,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    /// `c.1=` (locus identity, not whole-entity): a no-op record whose REF and
    /// ALT are both the reference base at the mapped position (CDS 1 → genomic
    /// 1049). REF must equal ALT and equal the provider's base there.
    #[test]
    fn test_build_vcf_record_identity() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Identity {
                    sequence: Some(Sequence::from_str("A").unwrap()),
                    whole_entity: false,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        // genomic 1049 = index 1048 of "ACGT".repeat → 'A'; identity REF == ALT.
        assert_eq!(record.pos, 1049);
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["A"]);
    }

    #[test]
    fn test_build_vcf_record_conversion_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(10)),
                NaEdit::Conversion {
                    source: "NM_other:c.100_110".to_string(),
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_vcf_record_unknown_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Unknown {
                    whole_entity: false,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_vcf_record_methylation_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Methylation {
                    status: MethylationStatus::GainOfMethylation,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    /// `c.1_100` copy number 3: emitted as a symbolic `<CNV>` structural record.
    /// POS/REF come from the start base (CDS 1 = tx 50 in exon1 → genomic 1049);
    /// INFO carries SVTYPE=CNV, CN=count, END=genomic end. CDS 100 = tx 149,
    /// which falls in exon2 (tx 101..200 → genomic 2000..2099), so it maps to
    /// genomic 2000 + (149 - 101) = 2048 — the exon boundary is crossed.
    #[test]
    fn test_build_vcf_record_copy_number() {
        use crate::vcf::record::InfoValue;

        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(100)),
                NaEdit::CopyNumber { count: 3 },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1049);
        // genomic 1049 = index 1048 of "ACGT".repeat → 'A'.
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["<CNV>"]);
        assert_eq!(
            record.info.get("SVTYPE"),
            Some(&InfoValue::String("CNV".to_string()))
        );
        assert_eq!(record.info.get("CN"), Some(&InfoValue::Integer(3)));
        // CDS 100 = tx 149 → exon2 → genomic 2048.
        assert_eq!(record.info.get("END"), Some(&InfoValue::Integer(2048)));
    }

    #[test]
    fn test_cds_with_intronic_position_warns() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::with_offset(10, 5)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_ok());
        let vcf = result.unwrap();
        assert!(!vcf.warnings.is_empty());
        assert!(vcf.warnings[0].contains("intronic"));
    }

    /// `c.1 CAG[(5_10)]`: a bounded range collapses to its lower bound (5) for
    /// VCF, so the ALT is anchor + unit repeated 5 times.
    #[test]
    fn test_repeat_count_range() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::Range(5, 10),
                    additional_counts: vec![],
                    trailing: None,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1049);
        assert_eq!(record.reference, "A");
        // Range lower bound (5) drives the expansion.
        assert_eq!(record.alternate, vec![format!("A{}", "CAG".repeat(5))]);
    }

    /// `c.1 CAG[(5_?)]`: a min-uncertain count uses its known lower bound (5),
    /// so the ALT is anchor + unit repeated 5 times.
    #[test]
    fn test_repeat_count_min_uncertain() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::MinUncertain(5),
                    additional_counts: vec![],
                    trailing: None,
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1049);
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec![format!("A{}", "CAG".repeat(5))]);
    }

    #[test]
    fn test_repeat_count_max_uncertain() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::MaxUncertain(10),
                    additional_counts: vec![],
                    trailing: None,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    /// `c.1 CAG[5]` with a trailing "T": the trailing bases are appended to the
    /// repeat unit before expansion, so the unit is "CAGT" repeated 5 times.
    #[test]
    fn test_repeat_with_trailing() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::Exact(5),
                    additional_counts: vec![],
                    trailing: Some(Sequence::from_str("T").unwrap()),
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1049);
        assert_eq!(record.reference, "A");
        // unit = "CAG" + trailing "T" = "CAGT", repeated 5 times.
        assert_eq!(record.alternate, vec![format!("A{}", "CAGT".repeat(5))]);
    }

    /// `c.1` repeat with no primary sequence but a trailing "T": the trailing
    /// bases alone form the repeat unit, so the unit is "T" repeated 5 times.
    #[test]
    fn test_repeat_trailing_only() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let seq = "ACGT".repeat(400);
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, &seq);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: None,
                    count: RepeatCount::Exact(5),
                    additional_counts: vec![],
                    trailing: Some(Sequence::from_str("T").unwrap()),
                },
            ),
        });

        let record = converter.convert(&variant).unwrap().record;
        assert_eq!(record.pos, 1049);
        assert_eq!(record.reference, "A");
        // unit = trailing "T" only, repeated 5 times.
        assert_eq!(record.alternate, vec![format!("A{}", "T".repeat(5))]);
    }

    // ------------------------------------------------------------------
    // Unresolvable genomic positions: offsets and pter/qter/cen
    // ------------------------------------------------------------------

    /// Sixteen bases of `chr1`, enough for every descriptor below to have a
    /// real anchor and a real deleted base when it is *not* refused. The
    /// guard under test fires before any fetch, so the content matters only
    /// to the control cases that must still convert.
    const CHR1_16: &str = "ACGTACGTACGTACGT";

    /// Convert `descriptor` through the full provider-backed converter, with
    /// sixteen bases of `chr1` available.
    fn convert_descriptor(descriptor: &str) -> Result<VcfConversionResult, FerroError> {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("chr1", CHR1_16);
        provider.add_genomic_sequence("chrM", CHR1_16);
        provider.add_genomic_sequence("001416", CHR1_16);
        let converter = HgvsToVcfConverter::new(&transcript, &provider);
        let variant = crate::parse_hgvs(descriptor).expect("fixture must parse");
        converter.convert(&variant)
    }

    /// A `g.` position carrying a `+`/`-` offset must not be silently
    /// flattened onto its base. `convert_genome` read only `GenomePos::base`,
    /// so `g.10+2del` and `g.10del` emitted byte-identical VCF records — a
    /// wrong row with no error at all, in both build profiles.
    ///
    /// The verdict mirrors #1641's ruling for the SPDI path: **refuse**. A
    /// genomic accession has no exon structure to measure an offset against,
    /// and VCF is purely positional, so there is nothing to resolve and no
    /// notation to carry it. The `m.` and `o.` axes hold the same
    /// `Interval<GenomePos>` and funnel through the same `convert_genome`, so
    /// all three are exercised.
    ///
    /// The message is asserted rather than a bare `is_err()`: several of these
    /// shapes can fail for unrelated reasons (a fetch that runs off the end of
    /// a 16-base contig, for one), and "failed somehow" would be satisfied by
    /// a build with the guard removed.
    #[test]
    fn a_genomic_offset_is_refused_rather_than_dropped() {
        for descriptor in [
            "NC_000001.11:g.10+2del",
            "NC_000001.11:g.10-2del",
            "NC_000001.11:g.10+2_12del",
            "NC_000001.11:g.4_10+2del",
            "NC_000001.11:g.10+2C>T",
            "NC_012920.1:m.10+2del",
            "NC_001416.1:o.10+2del",
        ] {
            let err = match convert_descriptor(descriptor) {
                Err(err) => err,
                Ok(result) => panic!(
                    "`{descriptor}` converted to POS={} REF={} ALT={:?}; the offset was dropped",
                    result.record.pos, result.record.reference, result.record.alternate
                ),
            };
            let message = err.to_string();
            assert!(
                message.contains("carries a +/- offset"),
                "`{descriptor}` was refused for some other reason: {message}"
            );
        }
    }

    /// The invariant the drop breaks, stated directly: two descriptions that
    /// the converter maps to the same VCF record must be the same variant.
    /// `g.10+2del`, `g.10-2del` and `g.10del` are three descriptions and were
    /// one row — `POS=9 REF=AC ALT=["A"]` for all three.
    ///
    /// The denominator is asserted, because once the guard is in place the
    /// two offset-carrying rows never reach `seen` and the collision check
    /// has nothing left to compare. If that count grows, an offset is
    /// convertible again and the collision check has just gone live.
    #[test]
    fn distinct_genomic_descriptions_do_not_share_one_vcf_record() {
        /// Only the offset-free spelling converts today.
        const CONVERTIBLE_TODAY: usize = 1;

        let descriptors = [
            "NC_000001.11:g.10+2del",
            "NC_000001.11:g.10-2del",
            "NC_000001.11:g.10del",
        ];
        let mut seen: Vec<(&str, String)> = Vec::new();
        for descriptor in descriptors {
            if let Ok(result) = convert_descriptor(descriptor) {
                let row = format!(
                    "{}:{} {} -> {:?}",
                    result.record.chrom,
                    result.record.pos,
                    result.record.reference,
                    result.record.alternate
                );
                if let Some((other, _)) = seen.iter().find(|(_, r)| *r == row) {
                    panic!(
                        "`{descriptor}` and `{other}` are different variants but share \
                         the VCF row `{row}`"
                    );
                }
                seen.push((descriptor, row));
            }
        }

        assert_eq!(
            seen.len(),
            CONVERTIBLE_TODAY,
            "{} of {} descriptors converted, expected {CONVERTIBLE_TODAY}: {seen:?}. \
             If this grew, an offset-carrying description is convertible again and the \
             collision check above has just become live",
            seen.len(),
            descriptors.len()
        );
    }

    /// `pter`, `qter` and `cen` name landmarks of the assembled chromosome,
    /// not positions on the sequence, and `GenomePos` stores `base: 0` for all
    /// three. `convert_genome` never called `is_special()`, so base 0 reached
    /// `build_vcf_record`, whose Deletion arm computes `start - 1` on `u64`:
    /// a debug panic (`attempt to subtract with overflow`) and a release wrap
    /// to `18446744073709551614`.
    ///
    /// The verdict mirrors #1662's ruling for the SPDI path: **refuse**, for
    /// the same reason — a sequence accession carries no telomere or
    /// centromere annotation to resolve the landmark against.
    ///
    /// The message is asserted, and that is load-bearing here: in release the
    /// wrapped coordinate is *already* caught by the provider's bounds check,
    /// so `is_err()` passes on unguarded code. It passes for the wrong reason,
    /// and a provider that does not bounds-check its fetch would emit a record
    /// at the wrapped coordinate instead.
    #[test]
    fn a_genomic_special_position_is_refused_rather_than_flattened() {
        for descriptor in [
            "NC_000001.11:g.pter_4del",
            "NC_000001.11:g.4_qterdel",
            "NC_000001.11:g.4_cendel",
            "NC_000001.11:g.pter_4delinsACGT",
            "NC_000001.11:g.4_qterdup",
            "NC_012920.1:m.pter_4del",
            "NC_001416.1:o.pter_4del",
        ] {
            let err = match convert_descriptor(descriptor) {
                Err(err) => err,
                Ok(result) => panic!(
                    "`{descriptor}` converted to POS={} REF={} ALT={:?}; the special \
                     position was flattened onto base 0",
                    result.record.pos, result.record.reference, result.record.alternate
                ),
            };
            let message = err.to_string();
            assert!(
                message.contains("names no numeric coordinate"),
                "`{descriptor}` was refused for some other reason: {message}"
            );
        }
    }

    /// The provider-less entry point reads the same field of the same struct,
    /// so it carried the same defect and gets the same guard. It matters most
    /// for the substitution arm, which is the one edit it lowers without
    /// needing a reference: unguarded, `g.10+2C>T` and `g.10C>T` both emitted
    /// `POS=10 REF=C ALT=[T]`.
    ///
    /// `ferro`'s `hgvs2vcf` subcommand is this function's only non-test
    /// caller, so this is the CLI-visible half of the defect.
    #[test]
    fn the_provider_less_path_refuses_an_unresolvable_genomic_position() {
        for (descriptor, expected) in [
            ("NC_000001.11:g.10+2C>T", "carries a +/- offset"),
            ("NC_000001.11:g.10-2C>T", "carries a +/- offset"),
            (
                "NC_000001.11:g.pter_4delACGT",
                "names no numeric coordinate",
            ),
            (
                "NC_000001.11:g.4_qterdelACGT",
                "names no numeric coordinate",
            ),
        ] {
            let variant = match crate::parse_hgvs(descriptor).expect("fixture must parse") {
                HgvsVariant::Genome(genome) => genome,
                other => panic!("`{descriptor}` parsed as {other:?}, expected a genomic variant"),
            };
            let err = match genomic_hgvs_to_vcf(&variant) {
                Err(err) => err,
                Ok(record) => panic!(
                    "`{descriptor}` converted to POS={} REF={} ALT={:?}",
                    record.pos, record.reference, record.alternate
                ),
            };
            let message = err.to_string();
            assert!(
                message.contains(expected),
                "`{descriptor}` was refused for some other reason: {message}"
            );
        }

        // Control: the offset-free spelling still converts on this path.
        let plain = match crate::parse_hgvs("NC_000001.11:g.10C>T").unwrap() {
            HgvsVariant::Genome(genome) => genome,
            other => panic!("expected a genomic variant, got {other:?}"),
        };
        let record = genomic_hgvs_to_vcf(&plain).expect("the legal spelling must still convert");
        assert_eq!(record.pos, 10);
    }

    /// The refusal names WHICH endpoint was unresolvable, not just the value.
    ///
    /// Both endpoints are checked in one loop, so before this the two sides were
    /// indistinguishable in the message: `g.pter_4del` and `g.4_qterdel` differ only
    /// in which end carries the landmark, and both reported the landmark alone. A
    /// caller holding a two-ended interval could not tell which end to correct.
    ///
    /// The two directions are asserted against each other rather than each on its
    /// own, so a label hardcoded to one side — the obvious wrong fix, and the shape
    /// a single-sided fixture would pass — fails here.
    #[test]
    fn the_refusal_names_which_endpoint_was_unresolvable() {
        for (descriptor, expected_side, forbidden_side) in [
            (
                "NC_000001.11:g.pter_4delACGT",
                "interval start",
                "interval end",
            ),
            (
                "NC_000001.11:g.4_qterdelACGT",
                "interval end",
                "interval start",
            ),
        ] {
            let variant = match crate::parse_hgvs(descriptor).expect("fixture must parse") {
                HgvsVariant::Genome(genome) => genome,
                other => panic!("`{descriptor}` parsed as {other:?}, expected a genomic variant"),
            };
            let message = match genomic_hgvs_to_vcf(&variant) {
                Err(err) => err.to_string(),
                Ok(record) => panic!("`{descriptor}` converted to POS={}", record.pos),
            };
            assert!(
                message.contains(expected_side),
                "`{descriptor}` must name `{expected_side}`: {message}"
            );
            assert!(
                !message.contains(forbidden_side),
                "`{descriptor}` named the wrong end (`{forbidden_side}`): {message}"
            );
        }
    }

    /// The offset drop is confined to `convert_genome`, but `start - 1` lives
    /// in `build_vcf_record`, which is the shared trunk: `convert_cds`,
    /// `convert_tx`, `convert_rna` and `convert_genome` all reach it. So the
    /// anchor arithmetic is checked at the trunk as well as guarded at the
    /// genomic arm, and this test calls the trunk directly with `start = 0` —
    /// no caller is needed for the arithmetic to be wrong.
    ///
    /// Without the check the Deletion arm panics in debug (`attempt to
    /// subtract with overflow`) and wraps in release; with it, both profiles
    /// refuse identically and say why.
    #[test]
    fn the_vcf_anchor_cannot_underflow_at_position_zero() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, CHR1_16);

        let edit = NaEdit::Deletion {
            sequence: Some(Sequence::from_str("A").unwrap()),
            length: None,
        };
        let err = converter
            .build_vcf_record("chr1", 0, 1, &edit, "synthetic")
            .expect_err("an anchor before position 0 does not exist");
        let message = err.to_string();
        assert!(
            message.contains("has no anchor base before it"),
            "the deletion arm underflowed or declined for another reason: {message}"
        );
    }

    /// The **second** `start - 1` in the trunk — the length-changing `Delins`
    /// branch — is checked too, but it is not what refuses `start = 0` and
    /// this test pins that rather than letting a bare `is_err()` imply
    /// otherwise. That arm fetches the deleted bases *before* it computes an
    /// anchor, and the fetch's own `start >= 1` validation fires first, so the
    /// subtraction is unreachable at zero.
    ///
    /// This is exactly the shape the SPDI siblings warn about: a decline that
    /// looks like coverage but comes from a check asking a different question.
    /// The arithmetic is still made checked, because "unreachable today"
    /// depends on an evaluation order nothing else pins — and this assertion
    /// is what will notice if that order changes.
    #[test]
    fn the_delins_anchor_is_shielded_by_the_reference_fetch() {
        let transcript = create_test_transcript();
        let mut provider = MockProvider::new();
        let converter = converter_with_genomic_chr1(&transcript, &mut provider, CHR1_16);

        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::from_str("ACGT").unwrap()),
            deleted: None,
            deleted_length: None,
            substitution_reference: None,
        };
        let err = converter
            .build_vcf_record("chr1", 0, 1, &edit, "synthetic")
            .expect_err("position 0 is not a valid 1-based coordinate");
        let message = err.to_string();
        assert!(
            message.contains("invalid 1-based genomic range [0, 1]"),
            "the delins arm no longer declines at the fetch; if it now declines at the \
             anchor check instead, that is an improvement — re-read the arm's order and \
             re-pin this. Got: {message}"
        );
    }
}
