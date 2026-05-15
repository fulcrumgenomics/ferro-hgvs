//! `VariantProjector` orchestrator.

use crate::data::mapping::MappingInfo;
use crate::data::projection::Projector;
use crate::error::FerroError;
use crate::hgvs::edit::NaEdit;
use crate::hgvs::interval::{CdsInterval, TxInterval};
use crate::hgvs::location::{CdsPos, GenomePos, TxPos};
use crate::hgvs::variant::{is_frameshift, CdsVariant, HgvsVariant, LocEdit, TxVariant};
use crate::normalize::{NormalizeConfig, Normalizer};
use crate::project::accession::parse_accession;
use crate::project::edit::transform_edit_for_strand;
use crate::project::protein::{predict_indel_protein, predict_substitution_protein};
use crate::project::result::VariantProjection;
use crate::reference::ReferenceProvider;

pub struct VariantProjector<P: ReferenceProvider + Clone> {
    projector: Projector,
    provider: P,
    normalizer: Normalizer<P>,
}

impl<P: ReferenceProvider + Clone> VariantProjector<P> {
    pub fn new(projector: Projector, provider: P) -> Self {
        let normalizer = Normalizer::with_config(provider.clone(), NormalizeConfig::default());
        Self {
            projector,
            provider,
            normalizer,
        }
    }

    pub fn with_normalize_config(mut self, config: NormalizeConfig) -> Self {
        self.normalizer = Normalizer::with_config(self.provider.clone(), config);
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
    pub fn project_variant(
        &self,
        variant: &HgvsVariant,
        transcript_id: &str,
    ) -> Result<VariantProjection, FerroError> {
        // 1. Normalize the genomic variant. The normalizer is built once at
        // construction time so we don't clone the (potentially heavy) provider
        // on every call.
        let normalized = self.normalizer.normalize(variant)?;

        // 2. Require a g. variant for now.
        let genome_variant = match &normalized {
            HgvsVariant::Genome(g) => g.clone(),
            _ => {
                return Err(FerroError::UnsupportedProjection {
                    reason: "VariantProjector currently only accepts g. variants".to_string(),
                });
            }
        };

        let edit = genome_variant
            .loc_edit
            .edit
            .inner()
            .cloned()
            .ok_or_else(|| FerroError::UnsupportedProjection {
                reason: "g. variant has no concrete edit".to_string(),
            })?;

        // 3. Look up the transcript in the cdot mapper.
        let cdot_tx = self
            .projector
            .mapper()
            .cdot()
            .get_transcript(transcript_id)
            .ok_or_else(|| FerroError::ReferenceNotFound {
                id: transcript_id.to_string(),
            })?;
        let strand = cdot_tx.strand;
        let gene_symbol = cdot_tx.gene_name.clone();
        let cdot_protein = cdot_tx.protein.clone();
        let is_coding = cdot_tx.cds_start.is_some();

        // 4. Extract start and end genomic positions from the variant interval.
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
        let normalized_str = normalized.to_string();

        // Compute the overall genomic extent of the transcript from its exons.
        // Exon format: [genome_start(0-based), genome_end(0-based excl), tx_start, tx_end].
        // The transcript's genomic extent is [min(genome_start), max(genome_end)) — 0-based half-open.
        let (tx_genome_start, tx_genome_end) = {
            let exons = &cdot_tx.exons;
            if exons.is_empty() {
                return Err(FerroError::ReferenceNotFound {
                    id: transcript_id.to_string(),
                });
            }
            let starts = exons.iter().map(|e| e[0]).min().unwrap();
            let ends = exons.iter().map(|e| e[1]).max().unwrap();
            (starts, ends)
        };

        // We map start and end positions independently rather than calling
        // CoordinateMapper::genome_interval_to_cds because the latter discards
        // the end position's MappingInfo (only start_info is propagated for the
        // is_intronic / in_5utr / in_3utr flags). For a range variant that crosses
        // an exon-intron boundary, we want is_intronic=true when EITHER end is
        // intronic, which requires both MappingInfos. The strand swap is otherwise
        // equivalent to what genome_interval_to_cds does internally.
        // Helper: map one GenomePos → CdsPos, converting out-of-range errors to
        // TranscriptNotOverlapping.
        let map_position = |gp: &GenomePos| -> Result<(CdsPos, MappingInfo), FerroError> {
            // Check whether the position falls within the transcript's overall genomic
            // extent before calling genome_to_cds.  genome_to_cds uses intronic
            // arithmetic to handle positions between exons, so it will happily return
            // a large offset for a position that is completely outside the transcript.
            // Exon bounds are 0-based half-open [genome_start, genome_end), so a
            // position is outside when it is before the first exon's start or at/after
            // the last exon's end.
            if gp.base < tx_genome_start || gp.base >= tx_genome_end {
                return Err(FerroError::TranscriptNotOverlapping {
                    variant: normalized_str.clone(),
                    transcript_id: transcript_id.to_string(),
                });
            }
            match mapper.genome_to_cds(transcript_id, gp) {
                Ok(res) => Ok((res.variant, res.info)),
                Err(FerroError::InvalidCoordinates { .. })
                | Err(FerroError::ConversionError { .. }) => {
                    Err(FerroError::TranscriptNotOverlapping {
                        variant: normalized_str.clone(),
                        transcript_id: transcript_id.to_string(),
                    })
                }
                Err(other) => Err(other),
            }
        };

        let (cds_start_raw, info_start) = map_position(&g_start)?;
        let (cds_end_raw, info_end) = map_position(&g_end)?;

        // On minus strand the start and end of the c. interval are swapped
        // relative to the g. interval (transcript reads in the opposite direction).
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

        // 5. Transform the edit for the transcript strand.
        let c_edit = transform_edit_for_strand(&edit, strand);

        // 6. Build the c./n. HGVS variant.
        let coding = if is_coding {
            let interval = CdsInterval::new(cds_start, cds_end);
            HgvsVariant::Cds(CdsVariant {
                accession: parse_accession(transcript_id),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, c_edit.clone()),
            })
        } else {
            // Non-coding transcript: report as n. using transcript positions.
            // CdsPos.base is signed; for non-coding it acts as the tx position.
            let tx_start = TxPos::new(cds_start.base);
            let tx_end = TxPos::new(cds_end.base);
            HgvsVariant::Tx(TxVariant {
                accession: parse_accession(transcript_id),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(TxInterval::new(tx_start, tx_end), c_edit.clone()),
            })
        };

        // 7. Derive per-position flags from MappingInfo.
        let is_intronic = info_start.is_intronic || info_end.is_intronic;
        let is_utr = !is_intronic
            && (info_start.in_5utr || info_start.in_3utr || info_end.in_5utr || info_end.in_3utr);

        // 8. Predict protein consequence for CDS variants.
        let mut protein = None;
        if !is_intronic && !is_utr && is_coding {
            // Prefer the explicit cdot.protein accession; otherwise infer NP_/XP_
            // from NM_/XM_ by stripping the prefix (not by substring replace, which
            // would mangle accessions whose suffix happens to contain "NM_"). If
            // we cannot infer a protein accession, skip protein prediction rather
            // than fabricate a bogus one.
            let prot_acc = match cdot_protein {
                Some(p) => Some(p),
                None => transcript_id
                    .strip_prefix("NM_")
                    .map(|rest| format!("NP_{rest}"))
                    .or_else(|| {
                        transcript_id
                            .strip_prefix("XM_")
                            .map(|rest| format!("XP_{rest}"))
                    }),
            };
            if let Some(prot_acc) = prot_acc {
                match &c_edit {
                    NaEdit::Substitution { .. } => {
                        let tx_for_codon = self.provider.get_transcript(transcript_id)?;
                        protein = Some(predict_substitution_protein(
                            &tx_for_codon,
                            cds_start.base,
                            &c_edit,
                            &prot_acc,
                        )?);
                    }
                    NaEdit::Deletion { .. }
                    | NaEdit::Insertion { .. }
                    | NaEdit::Duplication { .. }
                    | NaEdit::Delins { .. }
                    | NaEdit::Inversion { .. }
                        // Only predict when both CDS positions are concrete exonic positions
                        // (no intronic offsets — already guarded above).
                        if cds_start.offset.is_none()
                            && cds_end.offset.is_none()
                            && cds_start.base > 0
                            && cds_end.base > 0 =>
                    {
                        let tx_for_codon = self.provider.get_transcript(transcript_id)?;
                        match predict_indel_protein(
                            &tx_for_codon,
                            cds_start.base,
                            cds_end.base,
                            &c_edit,
                            &prot_acc,
                        ) {
                            Ok(pv) => protein = Some(pv),
                            // Non-fatal: unsupported edits or missing sequence → leave protein=None.
                            Err(FerroError::UnsupportedProjection { .. })
                            | Err(FerroError::ProteinSequenceUnavailable { .. }) => {}
                            Err(other) => return Err(other),
                        }
                    }
                    _ => {}
                }
            }
        }

        let frameshift = is_frameshift(&coding);
        Ok(VariantProjection {
            genomic: normalized.clone(),
            coding: Some(coding),
            protein,
            transcript_id: transcript_id.to_string(),
            gene_symbol,
            is_frameshift: frameshift,
            is_intronic,
            is_utr,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::cdot::{CdotMapper, CdotTranscript};
    use crate::data::projection::Projector;
    use crate::reference::mock::MockProvider;
    use crate::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
    use crate::reference::Strand as ProvStrand;
    use std::sync::OnceLock;

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
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
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
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        // Genomic sequence: 999 N's + "ATGCGCTAA" + 100 N's.
        // At 0-based index 999 = 'A' (g.1000), ..., 0-based index 1002 = 'C' (g.1003).
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ATGCGCTAA", suffix));
        (projector, provider)
    }

    fn make_minus_strand_provider_and_projector() -> (Projector, MockProvider) {
        // Same 9bp CDS on chr1, but transcript is on the minus strand.
        //
        // Plus-strand chr1:1000..1008 reads "TTAGCGCAT".  The minus-strand
        // transcript reads the complement strand from g.1008 backwards to
        // g.1000, giving CDS "ATGCGCTAA" (Met-Arg-Stop).
        //
        // cdot coordinate convention (matches make_test_provider_and_projector):
        //   exons:     [genome_start(0-based), genome_end(0-based excl), tx_start, tx_end]
        //   cds_start: 0-based tx position of first CDS base
        //   cds_end:   0-based tx position one-past the last CDS base
        //
        // Minus-strand mapping (genome_to_tx for Minus):
        //   offset = genome_end - 1 - genome_pos
        //   tx_pos = tx_start + offset
        //
        //   g.1008 → offset = 1009-1-1008 = 0 → tx_pos 0 (c.1, 'A')
        //   g.1007 → offset = 1 → tx_pos 1 (c.2, 'T')
        //   g.1005 → offset = 3 → tx_pos 3 → tx_to_cds(3, cds_start=0) = 4 → c.4
        //
        // Variant: g.1005G>A on the plus strand.
        //   The plus-strand base at g.1005 is 'G' (stated; actual byte is 'C'
        //   — but for substitutions the normalizer trusts the stated ref).
        //   After transform_edit_for_strand(Minus): G>A → revcomp → C>T.
        //   c.4 of "ATGCGCTAA" is 'C' (first base of codon 2 "CGC" = Arg).
        //   C>T → "TGC" = Cys.  Protein: p.(Arg2Cys).
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
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        // Plus-strand sequence on chr1 (not needed for substitution normalization,
        // provided for completeness).
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}TTAGCGCAT{}", prefix, suffix));
        (projector, provider)
    }

    #[test]
    fn project_substitution_minus_strand_revcomps_ref_alt() {
        let (projector, provider) = make_minus_strand_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // g.1005G>A: plus-strand G at position 1005.
        // Minus-strand transcript maps this to c.4 and revcomps G>A → C>T.
        // Codon 2 "CGC" (Arg) → "TGC" (Cys) = missense.
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
        assert_eq!(p, "NP_TEST_MINUS.1(TESTGENE):p.(Arg2Cys)");
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
        // g.1015 is 6 positions into the intron after exon 1 (which ends at genome
        // position 1009 inclusive / 1010 exclusive).
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
        // chr1:5000 is far outside the transcript at 1000..1009.
        let err = vp
            .project("NC_000001.11:g.5000A>G", "NM_TEST.1")
            .expect_err("should fail to project outside the transcript");
        assert!(
            matches!(err, FerroError::TranscriptNotOverlapping { .. }),
            "expected TranscriptNotOverlapping, got: {:?}",
            err
        );
    }

    #[test]
    fn project_substitution_plus_strand_missense() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);

        // g.1003 is the 4th base of the exon (0-based index 3 in "ATGCGCTAA" = 'C').
        // HGVS: g.1003C>A (cdot 0-based genome coord 1003 → tx_pos 3 → c.4).
        // CDS codon 2: CGC (Arg). Frame 0 substitution C>A → AGC (Ser) = missense.
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
        assert_eq!(p, "NP_TEST.1(TESTGENE):p.(Arg2Ser)");

        assert!(!result.is_frameshift);
        assert!(!result.is_intronic);
        assert!(!result.is_utr);
    }

    #[test]
    fn project_single_base_deletion_is_frameshift() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // g.1004del deletes a single base (c.5, the 'G' in the Arg codon CGC) — net=-1 → frameshift.
        // Mutated CDS: "ATGCCTAA" → Met-Pro-Stop → alt=[Met,Pro].
        // first_diff=1 (Arg≠Pro) → p.(Arg2Profs).
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
        // g.1003_1005del deletes 3 bases (c.4_6 = entire Arg codon CGC) — in-frame.
        // Mutated CDS: "ATGTAA" → [Met] → p.(Arg2del).
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
        // g.1003_1004insA inserts one base — frameshift (net=+1).
        let result = vp
            .project("NC_000001.11:g.1003_1004insA", "NM_TEST.1")
            .expect("insertion should project");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains("ins"), "expected ins notation, got: {}", c);
        // p. should be produced (frameshift) — exact form depends on insertion position.
        assert!(result.protein.is_some(), "p. expected for CDS insertion");
        assert!(result.is_frameshift, "1-base insertion is frameshift");
    }

    #[test]
    fn project_duplication() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // g.1004dup duplicates a single base (c.5 'G') — net=+1 → frameshift.
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
        // g.1003_1004delinsAT — replace 2 bases (c.4_5, first 2 bases of Arg codon CGC)
        // with AT → new codon "ATC" = Ile. Mutated CDS: "ATGATCTAA" → [Met, Ile].
        // In-frame (net=0), straddles codon → build_inframe_delins.
        // first_diff=1 (Arg≠Ile) → p.(Arg2delinsIle).
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
        assert!(p.contains("delins"), "expected delins in p.: {}", p);
        assert!(
            !result.is_frameshift,
            "delins of equal length is not a frameshift"
        );
    }

    #[test]
    fn project_inversion() {
        let (projector, provider) = make_test_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        // g.1003_1005inv — inversion of c.4_6 (CGC → GCG = Ala). Preserves length.
        // Mutated CDS: "ATGGCGTAA" → [Met, Ala] → p.(Arg2delinsAla).
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
        assert!(p.contains("Arg2"), "expected Arg2 in p.: {}", p);
        assert!(p.contains("delins"), "expected delins in p.: {}", p);
        assert!(p.contains("Ala"), "expected Ala in p.: {}", p);
        assert!(!result.is_frameshift, "inversion is not a frameshift");
    }

    fn make_ensembl_provider_and_projector() -> (Projector, MockProvider) {
        // Same CDS as make_test_provider_and_projector but the transcript ID is
        // an Ensembl accession (no NM_/XM_ prefix) and cdot.protein is absent.
        // The projector must NOT fabricate a bogus protein accession via
        // substring substitution; it should skip protein prediction instead.
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
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence("chr1", format!("{}{}{}", prefix, "ATGCGCTAA", suffix));
        (projector, provider)
    }

    #[test]
    fn project_substitution_no_protein_accession_skips_protein() {
        let (projector, provider) = make_ensembl_provider_and_projector();
        let vp = VariantProjector::new(projector, provider);
        let result = vp
            .project("chr1:g.1003C>A", "ENST00000000001.1")
            .expect("projection should succeed");
        let c = result.coding.as_ref().expect("c. expected").to_string();
        assert!(c.contains(":c.4C>A"), "expected c.4C>A, got: {}", c);
        // Neither NM_/XM_ prefix nor cdot.protein: must skip protein prediction.
        assert!(
            result.protein.is_none(),
            "expected no protein for accession with no NM_/XM_ prefix and no cdot.protein, got: {:?}",
            result.protein
        );
    }
}
