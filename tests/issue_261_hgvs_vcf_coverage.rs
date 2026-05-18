//! Audit for #81 § K2 — HGVS → VCF conversion surface.
//!
//! `genomic_hgvs_to_vcf` is the no-provider path. It supports
//! substitutions outright and rejects edits that need an anchor base
//! (deletion, insertion, delins) or reference data (inversion,
//! duplication). The full path uses `HgvsToVcfConverter` with a
//! provider; pinning that here would require a transcript fixture.
//!
//! This audit pins the no-provider surface plus the accession →
//! chromosome mapping that both paths share.

use ferro_hgvs::parse_hgvs;
use ferro_hgvs::vcf::genomic_hgvs_to_vcf;
use ferro_hgvs::HgvsVariant;

fn genome_variant(s: &str) -> ferro_hgvs::hgvs::variant::GenomeVariant {
    match parse_hgvs(s).expect("parse") {
        HgvsVariant::Genome(g) => g,
        other => panic!("expected GenomeVariant, got {other:?}"),
    }
}

// =============================================================================
// SECTION 1 — Substitution converts without reference data
// =============================================================================

mod substitution_no_anchor {
    use super::*;

    #[test]
    fn chr1_sub_emits_chr_pos_ref_alt() {
        let v = genome_variant("NC_000001.11:g.12345A>G");
        let vcf = genomic_hgvs_to_vcf(&v).expect("convert");
        assert_eq!(vcf.chrom, "chr1");
        assert_eq!(vcf.pos, 12345);
        assert_eq!(vcf.reference, "A");
        assert_eq!(vcf.alternate, vec!["G".to_string()]);
    }

    #[test]
    fn chr17_sub_emits_chr17() {
        let v = genome_variant("NC_000017.11:g.7676160C>T");
        let vcf = genomic_hgvs_to_vcf(&v).expect("convert");
        assert_eq!(vcf.chrom, "chr17");
        assert_eq!(vcf.reference, "C");
        assert_eq!(vcf.alternate, vec!["T".to_string()]);
    }

    #[test]
    fn chrx_sub_emits_chrx() {
        let v = genome_variant("NC_000023.11:g.100A>G");
        let vcf = genomic_hgvs_to_vcf(&v).expect("convert");
        assert_eq!(vcf.chrom, "chrX");
    }

    #[test]
    fn chry_sub_emits_chry() {
        let v = genome_variant("NC_000024.10:g.100A>G");
        let vcf = genomic_hgvs_to_vcf(&v).expect("convert");
        assert_eq!(vcf.chrom, "chrY");
    }

    #[test]
    fn chrm_sub_emits_chrm() {
        let v = genome_variant("NC_012920.1:g.100A>G");
        let vcf = genomic_hgvs_to_vcf(&v).expect("convert");
        assert_eq!(vcf.chrom, "chrM");
    }
}

// =============================================================================
// SECTION 2 — Edits requiring anchor base error without reference data
// =============================================================================

mod anchor_required {
    use super::*;

    #[test]
    fn del_with_seq_requires_anchor() {
        let v = genome_variant("NC_000001.11:g.100_102delATG");
        let err = genomic_hgvs_to_vcf(&v).expect_err("must error");
        assert!(
            err.to_string().to_lowercase().contains("anchor"),
            "expected anchor-base error, got: {err}"
        );
    }

    #[test]
    fn ins_requires_anchor() {
        let v = genome_variant("NC_000001.11:g.100_101insATG");
        let err = genomic_hgvs_to_vcf(&v).expect_err("must error");
        assert!(err.to_string().to_lowercase().contains("anchor"));
    }
}

// =============================================================================
// SECTION 3 — Other edit types unsupported by simple path
// =============================================================================

mod unsupported_simple_path {
    use super::*;

    #[test]
    fn inv_unsupported() {
        let v = genome_variant("NC_000001.11:g.100_105invTAGCA");
        assert!(genomic_hgvs_to_vcf(&v).is_err());
    }

    #[test]
    fn dup_unsupported() {
        let v = genome_variant("NC_000001.11:g.100_102dupATG");
        assert!(genomic_hgvs_to_vcf(&v).is_err());
    }

    #[test]
    fn delins_unsupported() {
        let v = genome_variant("NC_000001.11:g.100_102delinsTTCC");
        assert!(genomic_hgvs_to_vcf(&v).is_err());
    }
}

// =============================================================================
// SECTION 4 — Non-genomic coords cannot use the simple path
// =============================================================================
//
// The simple path takes a `GenomeVariant` only — `c.`/`n.`/`r.`/`p.`
// variants can't be passed in. We just pin the type-shape contract by
// matching the parsed variant.

mod non_genomic_not_genome_variant {
    use super::*;

    #[test]
    fn cds_variant_is_not_a_genome_variant() {
        assert!(matches!(
            parse_hgvs("NM_000088.3:c.100A>G").unwrap(),
            HgvsVariant::Cds(_)
        ));
    }

    #[test]
    fn protein_variant_is_not_a_genome_variant() {
        assert!(matches!(
            parse_hgvs("NP_003997.1:p.Arg8Gln").unwrap(),
            HgvsVariant::Protein(_)
        ));
    }

    #[test]
    fn rna_variant_is_not_a_genome_variant() {
        assert!(matches!(
            parse_hgvs("NM_000088.3:r.100a>g").unwrap(),
            HgvsVariant::Rna(_)
        ));
    }

    #[test]
    fn mt_variant_is_not_a_genome_variant() {
        assert!(matches!(
            parse_hgvs("NC_012920.1:m.100A>G").unwrap(),
            HgvsVariant::Mt(_)
        ));
    }
}

// =============================================================================
// SECTION 5 — Accession → chromosome mapping edge cases
// =============================================================================

mod accession_chromosome_mapping {
    use super::*;

    /// chr22 — last numbered autosome.
    #[test]
    fn chr22_emits_chr22() {
        let v = genome_variant("NC_000022.11:g.100A>G");
        assert_eq!(genomic_hgvs_to_vcf(&v).expect("convert").chrom, "chr22");
    }

    /// Older GRCh37 chr1 accession version still maps to chr1.
    #[test]
    fn grch37_chr1_emits_chr1() {
        let v = genome_variant("NC_000001.10:g.100A>G");
        assert_eq!(genomic_hgvs_to_vcf(&v).expect("convert").chrom, "chr1");
    }
}
