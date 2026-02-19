//! Parsing utilities for CLI operations

use crate::error::FerroError;
use crate::normalize::ShuffleDirection;
use crate::reference::transcript::GenomeBuild;
use crate::vcf::VcfRecord;

/// Parse a genome build string into a GenomeBuild enum
///
/// Accepts various common aliases:
/// - GRCh37, grch37, hg19, HG19 -> GRCh37
/// - GRCh38, grch38, hg38, HG38 -> GRCh38
/// - Other values default to GRCh38
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::parse_genome_build;
/// use ferro_hgvs::reference::transcript::GenomeBuild;
///
/// assert!(matches!(parse_genome_build("GRCh38"), GenomeBuild::GRCh38));
/// assert!(matches!(parse_genome_build("hg19"), GenomeBuild::GRCh37));
/// assert!(matches!(parse_genome_build("unknown"), GenomeBuild::GRCh38));
/// ```
pub fn parse_genome_build(build: &str) -> GenomeBuild {
    match build.to_uppercase().as_str() {
        "GRCH37" | "HG19" => GenomeBuild::GRCh37,
        "GRCH38" | "HG38" => GenomeBuild::GRCh38,
        _ => GenomeBuild::GRCh38,
    }
}

/// Parse a shuffle direction string into a ShuffleDirection enum
///
/// Accepts various common formats:
/// - 3prime, 3', 3 -> ThreePrime
/// - 5prime, 5', 5 -> FivePrime
/// - Other values default to ThreePrime
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::parse_shuffle_direction;
/// use ferro_hgvs::ShuffleDirection;
///
/// assert!(matches!(parse_shuffle_direction("3prime"), ShuffleDirection::ThreePrime));
/// assert!(matches!(parse_shuffle_direction("5'"), ShuffleDirection::FivePrime));
/// assert!(matches!(parse_shuffle_direction("unknown"), ShuffleDirection::ThreePrime));
/// ```
pub fn parse_shuffle_direction(direction: &str) -> ShuffleDirection {
    match direction.to_lowercase().as_str() {
        "5prime" | "5'" | "5" => ShuffleDirection::FivePrime,
        _ => ShuffleDirection::ThreePrime,
    }
}

/// Parse a VCF line into a VcfRecord
///
/// Parses a tab-separated VCF data line (not header lines starting with #).
/// Requires at least 5 fields: CHROM, POS, ID, REF, ALT.
///
/// **Note**: This function defaults to GRCh38 for the genome build. Use
/// [`parse_vcf_line_with_build`] if you need to specify a different build.
///
/// # Arguments
///
/// * `line` - A tab-separated VCF line
///
/// # Returns
///
/// A VcfRecord on success, or a FerroError if parsing fails
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::parse_vcf_line;
///
/// let record = parse_vcf_line("chr1\t12345\t.\tA\tG\t.\t.\t.").unwrap();
/// assert_eq!(record.chrom, "chr1");
/// assert_eq!(record.pos, 12345);
/// assert_eq!(record.reference, "A");
/// assert_eq!(record.alternate, vec!["G"]);
/// ```
pub fn parse_vcf_line(line: &str) -> Result<VcfRecord, FerroError> {
    parse_vcf_line_with_build(line, GenomeBuild::GRCh38)
}

/// Parse a VCF line into a VcfRecord with a specified genome build
///
/// Like [`parse_vcf_line`] but allows specifying the genome build.
///
/// # Arguments
///
/// * `line` - A tab-separated VCF line
/// * `build` - The genome build to use (GRCh37 or GRCh38)
///
/// # Returns
///
/// A VcfRecord on success, or a FerroError if parsing fails
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::parse_vcf_line_with_build;
/// use ferro_hgvs::reference::transcript::GenomeBuild;
///
/// let record = parse_vcf_line_with_build("chr1\t12345\t.\tA\tG", GenomeBuild::GRCh37).unwrap();
/// assert!(matches!(record.genome_build, GenomeBuild::GRCh37));
/// ```
pub fn parse_vcf_line_with_build(line: &str, build: GenomeBuild) -> Result<VcfRecord, FerroError> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 5 {
        return Err(FerroError::Parse {
            msg: format!(
                "Invalid VCF line: expected at least 5 fields, got {}",
                fields.len()
            ),
            pos: 0,
            diagnostic: None,
        });
    }

    let chrom = fields[0].to_string();
    let pos: u64 = fields[1].parse().map_err(|_| FerroError::Parse {
        msg: format!("Invalid position '{}': not a valid integer", fields[1]),
        pos: 0,
        diagnostic: None,
    })?;
    let id = if fields[2] == "." {
        None
    } else {
        Some(fields[2].to_string())
    };
    let reference = fields[3].to_string();
    let alternate: Vec<String> = fields[4].split(',').map(|s| s.to_string()).collect();

    Ok(VcfRecord {
        chrom,
        pos,
        id,
        reference,
        alternate,
        quality: None,
        filter: None,
        info: Default::default(),
        format: None,
        samples: Vec::new(),
        genome_build: build,
    })
}

/// Parse a VCF line with extended fields (QUAL, FILTER, INFO)
///
/// Like `parse_vcf_line` but also parses the optional QUAL, FILTER, and INFO fields.
///
/// # Arguments
///
/// * `line` - A tab-separated VCF line
///
/// # Returns
///
/// A VcfRecord on success, or a FerroError if parsing fails
pub fn parse_vcf_line_extended(line: &str) -> Result<VcfRecord, FerroError> {
    let mut record = parse_vcf_line(line)?;
    let fields: Vec<&str> = line.split('\t').collect();

    // Parse QUAL if present
    if fields.len() > 5 && fields[5] != "." {
        if let Ok(qual) = fields[5].parse::<f32>() {
            record.quality = Some(qual);
        }
    }

    // Parse FILTER if present
    if fields.len() > 6 && fields[6] != "." {
        record.filter = Some(fields[6].to_string());
    }

    // INFO field parsing would go here if needed

    Ok(record)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ===== Genome Build Parsing Tests =====

    #[test]
    fn test_parse_genome_build_grch38() {
        assert!(matches!(parse_genome_build("GRCh38"), GenomeBuild::GRCh38));
        assert!(matches!(parse_genome_build("grch38"), GenomeBuild::GRCh38));
        assert!(matches!(parse_genome_build("GRCH38"), GenomeBuild::GRCh38));
    }

    #[test]
    fn test_parse_genome_build_hg38() {
        assert!(matches!(parse_genome_build("hg38"), GenomeBuild::GRCh38));
        assert!(matches!(parse_genome_build("HG38"), GenomeBuild::GRCh38));
    }

    #[test]
    fn test_parse_genome_build_grch37() {
        assert!(matches!(parse_genome_build("GRCh37"), GenomeBuild::GRCh37));
        assert!(matches!(parse_genome_build("grch37"), GenomeBuild::GRCh37));
    }

    #[test]
    fn test_parse_genome_build_hg19() {
        assert!(matches!(parse_genome_build("hg19"), GenomeBuild::GRCh37));
        assert!(matches!(parse_genome_build("HG19"), GenomeBuild::GRCh37));
    }

    #[test]
    fn test_parse_genome_build_default() {
        // Unknown values should default to GRCh38
        assert!(matches!(parse_genome_build("unknown"), GenomeBuild::GRCh38));
        assert!(matches!(parse_genome_build(""), GenomeBuild::GRCh38));
        assert!(matches!(parse_genome_build("hg20"), GenomeBuild::GRCh38));
    }

    // ===== Shuffle Direction Parsing Tests =====

    #[test]
    fn test_parse_shuffle_direction_three_prime() {
        assert!(matches!(
            parse_shuffle_direction("3prime"),
            ShuffleDirection::ThreePrime
        ));
        assert!(matches!(
            parse_shuffle_direction("3'"),
            ShuffleDirection::ThreePrime
        ));
        assert!(matches!(
            parse_shuffle_direction("3"),
            ShuffleDirection::ThreePrime
        ));
    }

    #[test]
    fn test_parse_shuffle_direction_five_prime() {
        assert!(matches!(
            parse_shuffle_direction("5prime"),
            ShuffleDirection::FivePrime
        ));
        assert!(matches!(
            parse_shuffle_direction("5'"),
            ShuffleDirection::FivePrime
        ));
        assert!(matches!(
            parse_shuffle_direction("5"),
            ShuffleDirection::FivePrime
        ));
    }

    #[test]
    fn test_parse_shuffle_direction_case_insensitive() {
        assert!(matches!(
            parse_shuffle_direction("5PRIME"),
            ShuffleDirection::FivePrime
        ));
        assert!(matches!(
            parse_shuffle_direction("5Prime"),
            ShuffleDirection::FivePrime
        ));
    }

    #[test]
    fn test_parse_shuffle_direction_default() {
        // Unknown values should default to ThreePrime
        assert!(matches!(
            parse_shuffle_direction("unknown"),
            ShuffleDirection::ThreePrime
        ));
        assert!(matches!(
            parse_shuffle_direction(""),
            ShuffleDirection::ThreePrime
        ));
    }

    // ===== VCF Line Parsing Tests =====

    #[test]
    fn test_parse_vcf_line_basic() {
        let record = parse_vcf_line("chr1\t12345\t.\tA\tG").unwrap();
        assert_eq!(record.chrom, "chr1");
        assert_eq!(record.pos, 12345);
        assert!(record.id.is_none());
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["G"]);
    }

    #[test]
    fn test_parse_vcf_line_with_id() {
        let record = parse_vcf_line("chr1\t12345\trs12345\tA\tG").unwrap();
        assert_eq!(record.id, Some("rs12345".to_string()));
    }

    #[test]
    fn test_parse_vcf_line_multiple_alts() {
        let record = parse_vcf_line("chr1\t12345\t.\tA\tG,C,T").unwrap();
        assert_eq!(record.alternate, vec!["G", "C", "T"]);
    }

    #[test]
    fn test_parse_vcf_line_deletion() {
        let record = parse_vcf_line("chr1\t12345\t.\tATG\tA").unwrap();
        assert_eq!(record.reference, "ATG");
        assert_eq!(record.alternate, vec!["A"]);
    }

    #[test]
    fn test_parse_vcf_line_insertion() {
        let record = parse_vcf_line("chr1\t12345\t.\tA\tATGC").unwrap();
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["ATGC"]);
    }

    #[test]
    fn test_parse_vcf_line_complex_variant() {
        let record = parse_vcf_line("chr1\t12345\t.\tATG\tGCA").unwrap();
        assert_eq!(record.reference, "ATG");
        assert_eq!(record.alternate, vec!["GCA"]);
    }

    #[test]
    fn test_parse_vcf_line_extra_fields() {
        // VCF lines can have more than 5 fields
        let record =
            parse_vcf_line("chr1\t12345\t.\tA\tG\t30\tPASS\tDP=100;AF=0.5\tGT:AD\t0/1:50,50")
                .unwrap();
        assert_eq!(record.chrom, "chr1");
        assert_eq!(record.pos, 12345);
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["G"]);
    }

    #[test]
    fn test_parse_vcf_line_chr_prefix() {
        let record = parse_vcf_line("chr22\t12345\t.\tA\tG").unwrap();
        assert_eq!(record.chrom, "chr22");
    }

    #[test]
    fn test_parse_vcf_line_no_chr_prefix() {
        let record = parse_vcf_line("22\t12345\t.\tA\tG").unwrap();
        assert_eq!(record.chrom, "22");
    }

    #[test]
    fn test_parse_vcf_line_mitochondrial() {
        let record = parse_vcf_line("chrM\t1234\t.\tA\tG").unwrap();
        assert_eq!(record.chrom, "chrM");
    }

    #[test]
    fn test_parse_vcf_line_x_chromosome() {
        let record = parse_vcf_line("chrX\t12345\t.\tA\tG").unwrap();
        assert_eq!(record.chrom, "chrX");
    }

    #[test]
    fn test_parse_vcf_line_y_chromosome() {
        let record = parse_vcf_line("chrY\t12345\t.\tA\tG").unwrap();
        assert_eq!(record.chrom, "chrY");
    }

    #[test]
    fn test_parse_vcf_line_large_position() {
        let record = parse_vcf_line("chr1\t248956422\t.\tA\tG").unwrap();
        assert_eq!(record.pos, 248956422);
    }

    #[test]
    fn test_parse_vcf_line_position_one() {
        let record = parse_vcf_line("chr1\t1\t.\tA\tG").unwrap();
        assert_eq!(record.pos, 1);
    }

    #[test]
    fn test_parse_vcf_line_too_few_fields() {
        let result = parse_vcf_line("chr1\t12345\t.\tA");
        assert!(result.is_err());
        if let Err(FerroError::Parse { msg, .. }) = result {
            assert!(msg.contains("expected at least 5 fields"));
        }
    }

    #[test]
    fn test_parse_vcf_line_invalid_position() {
        let result = parse_vcf_line("chr1\tnotanumber\t.\tA\tG");
        assert!(result.is_err());
        if let Err(FerroError::Parse { msg, .. }) = result {
            assert!(msg.contains("Invalid position"));
        }
    }

    #[test]
    fn test_parse_vcf_line_empty() {
        let result = parse_vcf_line("");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_vcf_line_symbolic_alt() {
        // Symbolic alleles like <DEL>, <DUP>, <INS>
        let record = parse_vcf_line("chr1\t12345\t.\tA\t<DEL>").unwrap();
        assert_eq!(record.alternate, vec!["<DEL>"]);
    }

    #[test]
    fn test_parse_vcf_line_star_allele() {
        // Star allele indicates upstream deletion
        let record = parse_vcf_line("chr1\t12345\t.\tA\t*").unwrap();
        assert_eq!(record.alternate, vec!["*"]);
    }

    // ===== Extended VCF Line Parsing Tests =====

    #[test]
    fn test_parse_vcf_line_extended_with_qual() {
        let record = parse_vcf_line_extended("chr1\t12345\t.\tA\tG\t30.5").unwrap();
        assert_eq!(record.quality, Some(30.5));
    }

    #[test]
    fn test_parse_vcf_line_extended_with_filter() {
        let record = parse_vcf_line_extended("chr1\t12345\t.\tA\tG\t30\tPASS").unwrap();
        assert_eq!(record.filter, Some("PASS".to_string()));
    }

    #[test]
    fn test_parse_vcf_line_extended_missing_qual() {
        let record = parse_vcf_line_extended("chr1\t12345\t.\tA\tG\t.").unwrap();
        assert!(record.quality.is_none());
    }

    #[test]
    fn test_parse_vcf_line_extended_missing_filter() {
        let record = parse_vcf_line_extended("chr1\t12345\t.\tA\tG\t30\t.").unwrap();
        assert!(record.filter.is_none());
    }

    #[test]
    fn test_parse_vcf_line_default_genome_build() {
        let record = parse_vcf_line("chr1\t12345\t.\tA\tG").unwrap();
        assert!(matches!(record.genome_build, GenomeBuild::GRCh38));
    }
}
