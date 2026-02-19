//! Liftover endpoint - convert genomic coordinates between GRCh37 ↔ GRCh38

use axum::{extract::State, http::StatusCode, response::Json};
use std::time::Instant;

use crate::reference::transcript::GenomeBuild as ReferenceGenomeBuild;
use crate::service::{
    server::AppState,
    types::{ErrorResponse, GenomeBuild, LiftoverRequest, LiftoverResponse},
};

/// Liftover genomic coordinates between genome builds
///
/// This endpoint converts genomic positions between GRCh37 (hg19) and GRCh38 (hg38).
/// Input can be:
/// - Simple format: "chr7:117120148"
/// - HGVS format: "NC_000007.13:g.117120148"
pub async fn liftover(
    State(state): State<AppState>,
    Json(request): Json<LiftoverRequest>,
) -> Result<Json<LiftoverResponse>, (StatusCode, Json<ErrorResponse>)> {
    let start = Instant::now();

    // Validate that source and target builds are different
    if request.from_build == request.to_build {
        return Ok(Json(LiftoverResponse {
            input: request.position.clone(),
            from_build: request.from_build.as_str().to_string(),
            to_build: request.to_build.as_str().to_string(),
            converted: Some(request.position.clone()),
            hgvs_g: None,
            chain_region: None,
            error: None,
            processing_time_ms: start.elapsed().as_millis() as u64,
        }));
    }

    // Parse the input position
    let parsed = parse_genomic_position(&request.position);

    match parsed {
        Ok((chrom, pos)) => {
            // Check if liftover engine is available
            if let Some(liftover) = &state.liftover {
                // Perform actual liftover using chain files
                let from = convert_genome_build(&request.from_build);
                let to = convert_genome_build(&request.to_build);

                match liftover.lift(from, to, &chrom, pos) {
                    Ok(result) => {
                        let elapsed_ms = start.elapsed().as_millis() as u64;
                        let converted = format!("{}:{}", result.target_contig, result.target_pos);
                        let hgvs_g = get_refseq_accession(&result.target_contig, &request.to_build)
                            .map(|acc| format!("{}:g.{}", acc, result.target_pos));

                        let chain_region = Some(format!(
                            "chain_id={}, score={}{}",
                            result.chain_id,
                            result.chain_score,
                            if result.in_gap { " (interpolated)" } else { "" }
                        ));

                        let error = if result.multiple_mappings {
                            Some(
                                "Warning: Multiple chain mappings exist for this position"
                                    .to_string(),
                            )
                        } else {
                            None
                        };

                        Ok(Json(LiftoverResponse {
                            input: request.position,
                            from_build: request.from_build.as_str().to_string(),
                            to_build: request.to_build.as_str().to_string(),
                            converted: Some(converted),
                            hgvs_g,
                            chain_region,
                            error,
                            processing_time_ms: elapsed_ms,
                        }))
                    }
                    Err(e) => {
                        let elapsed_ms = start.elapsed().as_millis() as u64;
                        Ok(Json(LiftoverResponse {
                            input: request.position,
                            from_build: request.from_build.as_str().to_string(),
                            to_build: request.to_build.as_str().to_string(),
                            converted: None,
                            hgvs_g: None,
                            chain_region: None,
                            error: Some(format!("Liftover failed: {}", e)),
                            processing_time_ms: elapsed_ms,
                        }))
                    }
                }
            } else {
                // No liftover engine available - return helpful message
                let elapsed_ms = start.elapsed().as_millis() as u64;
                Ok(Json(LiftoverResponse {
                    input: request.position,
                    from_build: request.from_build.as_str().to_string(),
                    to_build: request.to_build.as_str().to_string(),
                    converted: None,
                    hgvs_g: None,
                    chain_region: None,
                    error: Some(
                        "Liftover requires chain files. Configure 'data.liftover.grch37_to_38' \
                         and 'data.liftover.grch38_to_37' in the service configuration. \
                         Chain files can be downloaded from UCSC: \
                         https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/"
                            .to_string(),
                    ),
                    processing_time_ms: elapsed_ms,
                }))
            }
        }
        Err(e) => {
            let elapsed_ms = start.elapsed().as_millis() as u64;
            Ok(Json(LiftoverResponse {
                input: request.position,
                from_build: request.from_build.as_str().to_string(),
                to_build: request.to_build.as_str().to_string(),
                converted: None,
                hgvs_g: None,
                chain_region: None,
                error: Some(e),
                processing_time_ms: elapsed_ms,
            }))
        }
    }
}

/// Convert service GenomeBuild to reference GenomeBuild
fn convert_genome_build(build: &GenomeBuild) -> ReferenceGenomeBuild {
    match build {
        GenomeBuild::GRCh37 => ReferenceGenomeBuild::GRCh37,
        GenomeBuild::GRCh38 => ReferenceGenomeBuild::GRCh38,
    }
}

/// Parse genomic position from various formats
fn parse_genomic_position(position: &str) -> Result<(String, u64), String> {
    // Try simple format: chr7:117120148 or 7:117120148
    if let Some((chrom_part, pos_part)) = position.split_once(':') {
        // Handle HGVS format: NC_000007.13:g.117120148
        if pos_part.starts_with("g.") {
            let pos_str = pos_part.trim_start_matches("g.");
            // Extract just the position number (before any variant notation)
            let pos_num = pos_str
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .collect::<String>();
            let pos = pos_num
                .parse::<u64>()
                .map_err(|_| format!("Invalid position: {}", pos_str))?;
            let chrom = extract_chromosome_from_accession(chrom_part);
            return Ok((chrom, pos));
        }

        // Simple format: chr7:117120148
        let pos = pos_part
            .parse::<u64>()
            .map_err(|_| format!("Invalid position: {}", pos_part))?;
        let chrom = normalize_chromosome(chrom_part);
        return Ok((chrom, pos));
    }

    Err(format!(
        "Invalid position format: {}. Expected 'chr7:117120148' or 'NC_000007.13:g.117120148'",
        position
    ))
}

/// Extract chromosome from RefSeq accession
fn extract_chromosome_from_accession(accession: &str) -> String {
    // Map NC_ accessions to chromosome names
    // NC_000001.10 (GRCh37) or NC_000001.11 (GRCh38) -> chr1
    if accession.starts_with("NC_") {
        if let Some(num_part) = accession.strip_prefix("NC_") {
            if let Some(chrom_num) = num_part.split('.').next() {
                if let Ok(num) = chrom_num.parse::<u32>() {
                    return match num {
                        1..=22 => format!("chr{}", num),
                        23 => "chrX".to_string(),
                        24 => "chrY".to_string(),
                        12920 => "chrM".to_string(),
                        _ => format!("chr{}", num),
                    };
                }
            }
        }
    }
    normalize_chromosome(accession)
}

/// Normalize chromosome name to chr-prefixed format
fn normalize_chromosome(chrom: &str) -> String {
    let chrom = chrom.trim();
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

/// Get RefSeq accession for a chromosome and build
fn get_refseq_accession(chrom: &str, build: &GenomeBuild) -> Option<String> {
    let chrom_normalized = chrom.trim_start_matches("chr");
    let chrom_num = chrom_normalized
        .parse::<u32>()
        .ok()
        .or(match chrom_normalized {
            "X" => Some(23),
            "Y" => Some(24),
            "M" | "MT" => Some(12920),
            _ => None,
        })?;

    let version = match build {
        GenomeBuild::GRCh37 => match chrom_num {
            1 => "10",
            2 => "11",
            3 => "11",
            4 => "11",
            5 => "9",
            6 => "11",
            7 => "13",
            8 => "10",
            9 => "11",
            10 => "10",
            11 => "9",
            12 => "11",
            13 => "10",
            14 => "8",
            15 => "9",
            16 => "9",
            17 => "10",
            18 => "9",
            19 => "9",
            20 => "10",
            21 => "8",
            22 => "10",
            23 => "10",
            24 => "9",
            12920 => "1",
            _ => return None,
        },
        GenomeBuild::GRCh38 => match chrom_num {
            1 => "11",
            2 => "12",
            3 => "12",
            4 => "12",
            5 => "10",
            6 => "12",
            7 => "14",
            8 => "11",
            9 => "12",
            10 => "11",
            11 => "10",
            12 => "12",
            13 => "11",
            14 => "9",
            15 => "10",
            16 => "10",
            17 => "11",
            18 => "10",
            19 => "10",
            20 => "11",
            21 => "9",
            22 => "11",
            23 => "11",
            24 => "10",
            12920 => "1",
            _ => return None,
        },
    };

    if chrom_num == 12920 {
        Some(format!("NC_012920.{}", version))
    } else {
        Some(format!("NC_{:06}.{}", chrom_num, version))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_position() {
        let (chrom, pos) = parse_genomic_position("chr7:117120148").unwrap();
        assert_eq!(chrom, "chr7");
        assert_eq!(pos, 117120148);
    }

    #[test]
    fn test_parse_position_without_chr_prefix() {
        let (chrom, pos) = parse_genomic_position("7:117120148").unwrap();
        assert_eq!(chrom, "chr7");
        assert_eq!(pos, 117120148);
    }

    #[test]
    fn test_parse_hgvs_genomic() {
        let (chrom, pos) = parse_genomic_position("NC_000007.13:g.117120148").unwrap();
        assert_eq!(chrom, "chr7");
        assert_eq!(pos, 117120148);
    }

    #[test]
    fn test_parse_hgvs_with_variant() {
        // Should extract just the position even if there's variant notation after
        let (chrom, pos) = parse_genomic_position("NC_000007.14:g.117559593").unwrap();
        assert_eq!(chrom, "chr7");
        assert_eq!(pos, 117559593);
    }

    #[test]
    fn test_extract_chromosome_from_accession() {
        assert_eq!(extract_chromosome_from_accession("NC_000001.10"), "chr1");
        assert_eq!(extract_chromosome_from_accession("NC_000007.13"), "chr7");
        assert_eq!(extract_chromosome_from_accession("NC_000023.10"), "chrX");
        assert_eq!(extract_chromosome_from_accession("NC_000024.09"), "chrY");
    }

    #[test]
    fn test_get_refseq_accession() {
        assert_eq!(
            get_refseq_accession("chr1", &GenomeBuild::GRCh37),
            Some("NC_000001.10".to_string())
        );
        assert_eq!(
            get_refseq_accession("chr1", &GenomeBuild::GRCh38),
            Some("NC_000001.11".to_string())
        );
        assert_eq!(
            get_refseq_accession("chrX", &GenomeBuild::GRCh38),
            Some("NC_000023.11".to_string())
        );
    }

    #[test]
    fn test_invalid_position() {
        assert!(parse_genomic_position("invalid").is_err());
        assert!(parse_genomic_position("chr7:abc").is_err());
    }
}
