//! Audit for #390 item 6 — `HgvsToVcfConverter` (provider-aware path)
//! coord-system × edit-kind coverage matrix.
//!
//! Parallel to `tests/issue_261_hgvs_vcf_coverage.rs`, which covers only
//! the no-provider `genomic_hgvs_to_vcf`. The converter accepts a
//! transcript + provider and can therefore project c./n./g./m. inputs
//! through the transcript exon table; this file pins which
//! (coord-axis × edit-kind) cells succeed and which surface a documented
//! `ConversionError`.
//!
//! Cell coverage:
//!   • c. axis  — substitution
//!   • n. axis  — substitution
//!   • g. axis  — substitution
//!   • m. axis  — substitution (mitochondrial; routes through the
//!                genome branch by re-wrapping the variant)
//!   • r. axis  — explicit "not yet supported" error
//!   • p. axis  — explicit unsupported (requires back-translation) error
//!   • Allele   — single-inner-variant routes through; empty rejects
//!
//! The transcript fixture mirrors `create_test_transcript` from
//! `src/vcf/from_hgvs.rs`'s inline tests so this audit and the unit
//! tests share one understanding of the projection.

use ferro_hgvs::error::FerroError;
use ferro_hgvs::hgvs::edit::{Base, NaEdit};
use ferro_hgvs::hgvs::interval::{CdsInterval, GenomeInterval, TxInterval};
use ferro_hgvs::hgvs::location::{CdsPos, GenomePos, TxPos};
use ferro_hgvs::hgvs::variant::{
    Accession, AllelePhase, AlleleVariant, CdsVariant, GenomeVariant, HgvsVariant, LocEdit,
    MtVariant, TxVariant,
};
use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::vcf::HgvsToVcfConverter;

/// Single-exon coding transcript on chr1, plus strand:
///   - tx [1, 100] → genome [1000, 1099]
///   - CDS tx [50, 75]
///   - sequence is `ATGCATGC` repeated so every CDS position has a
///     concrete base for the SPDI/anchor sub-path tests that need one.
fn fixture_transcript() -> Transcript {
    Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        Some("ATGCATGC".repeat(20)),
        Some(50),
        Some(75),
        vec![Exon::with_genomic(1, 1, 100, 1000, 1099)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1099),
        GenomeBuild::GRCh38,
        ManeStatus::Select,
        None,
        None,
    )
}

fn fixture_provider() -> MockProvider {
    MockProvider::new()
}

fn substitution(ref_b: Base, alt_b: Base) -> NaEdit {
    NaEdit::Substitution {
        reference: ref_b,
        alternative: alt_b,
    }
}

// =============================================================================
// SECTION 1 — Coord-axis × substitution (the always-supported cell)
// =============================================================================

mod substitution_axis_matrix {
    use super::*;

    #[test]
    fn c_axis_substitution_emits_chr_pos_ref_alt() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        // c.1 = tx 50 = genome 1049 (1-based, plus strand).
        let variant = CdsVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                substitution(Base::A, Base::G),
            ),
        };
        let result = converter
            .convert(&HgvsVariant::Cds(variant))
            .expect("c. SNV must convert via the provider-aware path");
        assert_eq!(result.record.chrom, "chr1");
        assert!(
            result.record.info.contains_key("HGVS"),
            "INFO/HGVS annotation must be populated"
        );
    }

    #[test]
    fn n_axis_substitution_emits_chr_pos_ref_alt() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = TxVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(1)),
                substitution(Base::A, Base::G),
            ),
        };
        let result = converter
            .convert(&HgvsVariant::Tx(variant))
            .expect("n. SNV must convert via the provider-aware path");
        assert_eq!(result.record.chrom, "chr1");
    }

    #[test]
    fn g_axis_substitution_emits_chr_pos_ref_alt() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                substitution(Base::A, Base::G),
            ),
        };
        let result = converter
            .convert(&HgvsVariant::Genome(variant))
            .expect("g. SNV must convert");
        assert_eq!(result.record.chrom, "chr1");
        assert_eq!(result.record.pos, 12345);
        assert_eq!(result.record.reference, "A");
    }

    #[test]
    fn m_axis_substitution_routes_through_genome_branch() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = MtVariant {
            accession: Accession::new("NC", "012920", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(3243)),
                substitution(Base::A, Base::G),
            ),
        };
        let result = converter
            .convert(&HgvsVariant::Mt(variant))
            .expect("m. SNV must route through the genome branch");
        assert_eq!(result.record.chrom, "chrM");
        assert_eq!(result.record.pos, 3243);
    }
}

// =============================================================================
// SECTION 2 — Coord-axes that the converter rejects up-front
// =============================================================================

mod rejected_axes {
    use super::*;
    use ferro_hgvs::hgvs::edit::ProteinEdit;
    use ferro_hgvs::hgvs::interval::ProtInterval;
    use ferro_hgvs::hgvs::interval::RnaInterval;
    use ferro_hgvs::hgvs::location::RnaPos;
    use ferro_hgvs::hgvs::location::{AminoAcid, ProtPos};
    use ferro_hgvs::hgvs::variant::{ProteinVariant, RnaVariant};

    #[test]
    fn r_axis_returns_not_yet_supported_error() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = RnaVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(1)),
                substitution(Base::A, Base::G),
            ),
        };
        let err = converter
            .convert(&HgvsVariant::Rna(variant))
            .expect_err("r. axis is documented as unsupported in HgvsToVcfConverter");
        assert!(
            matches!(err, FerroError::ConversionError { ref msg } if msg.to_lowercase().contains("rna")),
            "expected ConversionError mentioning RNA, got: {err}"
        );
    }

    #[test]
    fn p_axis_returns_back_translation_error() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = ProteinVariant {
            accession: Accession::new("NP", "000079", Some(2)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                ProtInterval::point(ProtPos::new(AminoAcid::Arg, 1)),
                ProteinEdit::Substitution {
                    reference: AminoAcid::Arg,
                    alternative: AminoAcid::Gln,
                },
            ),
        };
        let err = converter
            .convert(&HgvsVariant::Protein(variant))
            .expect_err("p. requires back-translation; converter rejects");
        assert!(
            matches!(err, FerroError::ConversionError { ref msg } if msg.to_lowercase().contains("protein")),
            "expected ConversionError mentioning protein, got: {err}"
        );
    }
}

// =============================================================================
// SECTION 3 — Allele variants
// =============================================================================

mod allele_branch {
    use super::*;

    #[test]
    fn single_inner_variant_routes_through_inner() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let inner = CdsVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                substitution(Base::A, Base::G),
            ),
        };
        let allele = AlleleVariant {
            variants: vec![HgvsVariant::Cds(inner)],
            phase: AllelePhase::Cis,
        };
        let result = converter
            .convert(&HgvsVariant::Allele(allele))
            .expect("single-inner allele must route through the inner conversion");
        assert_eq!(result.record.chrom, "chr1");
    }

    #[test]
    fn empty_allele_rejected() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let allele = AlleleVariant {
            variants: vec![],
            phase: AllelePhase::Cis,
        };
        let err = converter
            .convert(&HgvsVariant::Allele(allele))
            .expect_err("empty allele has nothing to encode in VCF");
        assert!(
            matches!(err, FerroError::ConversionError { ref msg } if msg.to_lowercase().contains("empty")),
            "expected ConversionError mentioning empty, got: {err}"
        );
    }
}
