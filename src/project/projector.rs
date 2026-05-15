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
use crate::project::protein::predict_substitution_protein;
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

        // 8. Predict protein consequence for CDS substitutions.
        let mut protein = None;
        if !is_intronic && !is_utr && is_coding && matches!(c_edit, NaEdit::Substitution { .. }) {
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
                let tx_for_codon = self.provider.get_transcript(transcript_id)?;
                protein = Some(predict_substitution_protein(
                    &tx_for_codon,
                    cds_start.base,
                    &c_edit,
                    &prot_acc,
                )?);
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
}
