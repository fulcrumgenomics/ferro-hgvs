//! VCF conversion endpoints - bidirectional VCF ↔ HGVS conversion

use axum::{extract::State, http::StatusCode, response::Json};
use std::sync::Arc;
use std::time::Instant;

use crate::data::cdot::{CdotMapper, CdsPosition};
use crate::reference::Strand;
use crate::service::{
    server::AppState,
    types::{
        ErrorResponse, GenomeBuild, HgvsToVcfRequest, HgvsToVcfResponse, ServiceError, VcfRecord,
        VcfToHgvsRequest, VcfToHgvsResponse,
    },
    validation::validate_hgvs,
};

/// Convert VCF record to HGVS notation
///
/// This endpoint converts a VCF-style variant representation (CHROM, POS, REF, ALT)
/// to HGVS notation. Optionally provide a transcript for c. notation.
pub async fn vcf_to_hgvs(
    State(state): State<AppState>,
    Json(request): Json<VcfToHgvsRequest>,
) -> Result<Json<VcfToHgvsResponse>, (StatusCode, Json<ErrorResponse>)> {
    let start = Instant::now();

    // Validate inputs
    if request.ref_allele.is_empty() {
        let error = ServiceError::BadRequest("ref allele cannot be empty".to_string());
        return Err((
            StatusCode::from_u16(error.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR),
            Json(error.to_response()),
        ));
    }
    if request.alt.is_empty() {
        let error = ServiceError::BadRequest("alt allele cannot be empty".to_string());
        return Err((
            StatusCode::from_u16(error.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR),
            Json(error.to_response()),
        ));
    }

    let vcf = VcfRecord {
        chrom: request.chrom.clone(),
        pos: request.pos,
        ref_allele: request.ref_allele.clone(),
        alt: request.alt.clone(),
        build: request.build.as_str().to_string(),
    };

    // Convert VCF to HGVS
    let (hgvs_g, hgvs_c, hgvs_p, error) = convert_vcf_to_hgvs(
        &request.chrom,
        request.pos,
        &request.ref_allele,
        &request.alt,
        &request.build,
        request.transcript.as_deref(),
        state.cdot.as_ref(),
    );

    let elapsed_ms = start.elapsed().as_millis() as u64;

    Ok(Json(VcfToHgvsResponse {
        vcf,
        hgvs_g,
        hgvs_c,
        hgvs_p,
        error,
        processing_time_ms: elapsed_ms,
    }))
}

/// Convert HGVS notation to VCF record
///
/// This endpoint converts an HGVS variant to VCF-style representation.
pub async fn hgvs_to_vcf(
    State(state): State<AppState>,
    Json(request): Json<HgvsToVcfRequest>,
) -> Result<Json<HgvsToVcfResponse>, (StatusCode, Json<ErrorResponse>)> {
    let start = Instant::now();

    // Validate input HGVS
    if let Err(validation_error) = validate_hgvs(&request.hgvs) {
        let error = ServiceError::InvalidHgvs(validation_error.to_string());
        return Err((
            StatusCode::from_u16(error.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR),
            Json(error.to_response()),
        ));
    }

    // Parse the input HGVS
    let hgvs_str = request.hgvs.clone();
    let parse_result =
        tokio::task::spawn_blocking(move || crate::hgvs::parser::parse_hgvs_lenient(&hgvs_str))
            .await
            .map_err(|e| {
                let error = ServiceError::InternalError(format!("Task error: {}", e));
                (StatusCode::INTERNAL_SERVER_ERROR, Json(error.to_response()))
            })?;

    let elapsed_ms = start.elapsed().as_millis() as u64;

    match parse_result {
        Ok(result) => {
            let (vcf, error) =
                convert_hgvs_to_vcf(&result.result, &request.build, state.cdot.as_ref());

            Ok(Json(HgvsToVcfResponse {
                input: request.hgvs,
                vcf,
                error,
                processing_time_ms: elapsed_ms,
            }))
        }
        Err(e) => Ok(Json(HgvsToVcfResponse {
            input: request.hgvs,
            vcf: None,
            error: Some(format!("Failed to parse input: {}", e)),
            processing_time_ms: elapsed_ms,
        })),
    }
}

/// Convert VCF fields to HGVS notation
fn convert_vcf_to_hgvs(
    chrom: &str,
    pos: u64,
    ref_allele: &str,
    alt: &str,
    build: &GenomeBuild,
    transcript: Option<&str>,
    cdot: Option<&Arc<CdotMapper>>,
) -> (
    Option<String>,
    Option<String>,
    Option<String>,
    Option<String>,
) {
    // Get RefSeq accession for the chromosome
    let accession = match get_refseq_accession(chrom, build) {
        Some(acc) => acc,
        None => {
            return (
                None,
                None,
                None,
                Some(format!("Unknown chromosome: {}", chrom)),
            );
        }
    };

    // Determine variant type and generate genomic HGVS
    let (hgvs_g, variant_type, _variant_pos) =
        generate_genomic_hgvs(&accession, pos, ref_allele, alt);

    // c. and p. notation require transcript mapping
    let (hgvs_c, hgvs_p, tx_error) = if let Some(tx_id) = transcript {
        convert_to_transcript_notation(pos, ref_allele, alt, tx_id, cdot, &variant_type)
    } else {
        (None, None, None)
    };

    // Combine errors
    let error = if transcript.is_some() && hgvs_c.is_none() && tx_error.is_some() {
        tx_error
    } else {
        None
    };

    (hgvs_g, hgvs_c, hgvs_p, error)
}

/// Variant type for internal use
#[derive(Debug, Clone)]
enum VariantType {
    Substitution,
    Deletion,
    Insertion,
    Delins,
    Unknown,
}

/// Generate genomic HGVS notation from VCF fields
fn generate_genomic_hgvs(
    accession: &str,
    pos: u64,
    ref_allele: &str,
    alt: &str,
) -> (Option<String>, VariantType, u64) {
    if ref_allele.len() == 1 && alt.len() == 1 && ref_allele != alt {
        // Simple substitution
        let hgvs = format!("{}:g.{}{}>{}", accession, pos, ref_allele, alt);
        (Some(hgvs), VariantType::Substitution, pos)
    } else if ref_allele.len() > alt.len() && alt.len() == 1 && ref_allele.starts_with(alt) {
        // Deletion (VCF style with padding base)
        let del_start = pos + 1;
        let del_end = pos + (ref_allele.len() - 1) as u64;
        let hgvs = if del_start == del_end {
            format!("{}:g.{}del", accession, del_start)
        } else {
            format!("{}:g.{}_{}del", accession, del_start, del_end)
        };
        (Some(hgvs), VariantType::Deletion, del_start)
    } else if alt.len() > ref_allele.len() && ref_allele.len() == 1 && alt.starts_with(ref_allele) {
        // Insertion (VCF style with padding base)
        let inserted = &alt[1..];
        let hgvs = format!("{}:g.{}_{}ins{}", accession, pos, pos + 1, inserted);
        (Some(hgvs), VariantType::Insertion, pos)
    } else if ref_allele != alt {
        // Delins
        let hgvs = format!("{}:g.{}delins{}", accession, pos, alt);
        (Some(hgvs), VariantType::Delins, pos)
    } else {
        (None, VariantType::Unknown, pos)
    }
}

/// Convert genomic position to transcript (c.) and protein (p.) notation
fn convert_to_transcript_notation(
    genomic_pos: u64,
    ref_allele: &str,
    alt: &str,
    transcript_id: &str,
    cdot: Option<&Arc<CdotMapper>>,
    variant_type: &VariantType,
) -> (Option<String>, Option<String>, Option<String>) {
    // Need cdot data for transcript lookups
    let cdot = match cdot {
        Some(c) => c,
        None => {
            return (
                None,
                None,
                Some(
                    "Transcript mapping requires cdot data. Configure 'data.cdot_path' in settings."
                        .to_string(),
                ),
            );
        }
    };

    // Get the transcript
    let cdot_tx = match cdot.get_transcript(transcript_id) {
        Some(tx) => tx,
        None => {
            return (
                None,
                None,
                Some(format!(
                    "Transcript {} not found in cdot data",
                    transcript_id
                )),
            );
        }
    };

    // Try to convert genomic position to transcript position
    let tx_pos = match cdot_tx.genome_to_tx(genomic_pos) {
        Some(pos) => pos,
        None => {
            return (
                None,
                None,
                Some(format!(
                    "Genomic position {} not found in transcript {} exons",
                    genomic_pos, transcript_id
                )),
            );
        }
    };

    // Convert to CDS position if available
    let cds_pos = cdot_tx.tx_to_cds(tx_pos);

    // Generate c. notation
    let hgvs_c = generate_cds_hgvs(
        transcript_id,
        cds_pos.as_ref(),
        variant_type,
        ref_allele,
        alt,
    );

    // Generate p. notation for coding variants (if applicable)
    let hgvs_p = match &cds_pos {
        Some(CdsPosition::Cds(base)) if *base > 0 => {
            // Calculate protein position: (cds_pos - 1) / 3 + 1
            let prot_position = ((*base - 1) / 3 + 1) as u64;
            let prot_acc = transcript_id.replace("NM_", "NP_").replace("XM_", "XP_");

            if prot_acc.starts_with("NP_") || prot_acc.starts_with("XP_") {
                Some(format!("{}:p.{}", prot_acc, prot_position))
            } else {
                None // Non-coding transcript
            }
        }
        _ => None, // UTR or intronic variant
    };

    (hgvs_c, hgvs_p, None)
}

/// Generate c. HGVS notation from CDS position
fn generate_cds_hgvs(
    transcript_id: &str,
    cds_pos: Option<&CdsPosition>,
    variant_type: &VariantType,
    ref_allele: &str,
    alt: &str,
) -> Option<String> {
    let pos_str = match cds_pos {
        Some(CdsPosition::Cds(base)) => format!("{}", base),
        Some(CdsPosition::FivePrimeUtr(offset)) => format!("-{}", offset),
        Some(CdsPosition::ThreePrimeUtr(offset)) => format!("*{}", offset),
        None => return None, // Can't generate c. without CDS info
    };

    let edit_str = match variant_type {
        VariantType::Substitution => {
            format!("{}>{}", ref_allele, alt)
        }
        VariantType::Deletion => "del".to_string(),
        VariantType::Insertion => {
            let inserted = if alt.len() > 1 { &alt[1..] } else { alt };
            format!("ins{}", inserted)
        }
        VariantType::Delins => format!("delins{}", alt),
        VariantType::Unknown => return None,
    };

    Some(format!("{}:c.{}{}", transcript_id, pos_str, edit_str))
}

/// Convert parsed HGVS to VCF record
fn convert_hgvs_to_vcf(
    variant: &crate::hgvs::variant::HgvsVariant,
    build: &GenomeBuild,
    cdot: Option<&Arc<CdotMapper>>,
) -> (Option<VcfRecord>, Option<String>) {
    use crate::hgvs::variant::HgvsVariant;

    match variant {
        HgvsVariant::Genome(v) => {
            // Extract chromosome from accession
            let chrom = extract_chromosome_from_accession(&v.accession.to_string());

            // Extract position and alleles from the edit
            if let Some(edit) = v.loc_edit.edit.inner() {
                let (pos, ref_allele, alt) = extract_vcf_fields(edit, &v.loc_edit.location);

                if let (Some(pos), Some(ref_a), Some(alt_a)) = (pos, ref_allele, alt) {
                    return (
                        Some(VcfRecord {
                            chrom,
                            pos,
                            ref_allele: ref_a,
                            alt: alt_a,
                            build: build.as_str().to_string(),
                        }),
                        None,
                    );
                }
            }

            (
                None,
                Some("Could not extract VCF fields from variant".to_string()),
            )
        }
        HgvsVariant::Cds(v) => {
            // c. variant - need coordinate conversion to g.
            let accession = v.accession.to_string();
            convert_cds_to_vcf(&accession, &v.loc_edit, build, cdot)
        }
        HgvsVariant::Tx(v) => {
            // n. variant - need coordinate conversion to g.
            let accession = v.accession.to_string();
            convert_tx_to_vcf(&accession, &v.loc_edit, build, cdot)
        }
        _ => (
            None,
            Some("VCF conversion is only supported for g., c., and n. variants".to_string()),
        ),
    }
}

/// Convert c. variant to VCF via coordinate mapping
///
/// Handles both exonic and intronic variants:
/// - Exonic: c.350C>T - position within CDS
/// - Intronic: c.117-2del - position 2 bases into intron before exon position 117
fn convert_cds_to_vcf(
    transcript_id: &str,
    loc_edit: &crate::hgvs::variant::LocEdit<
        crate::hgvs::interval::CdsInterval,
        crate::hgvs::edit::NaEdit,
    >,
    build: &GenomeBuild,
    cdot: Option<&Arc<CdotMapper>>,
) -> (Option<VcfRecord>, Option<String>) {
    // Need cdot for coordinate conversion
    let cdot = match cdot {
        Some(c) => c,
        None => {
            return (
                None,
                Some(
                    "Converting c. to VCF requires cdot data. Configure 'data.cdot_path' in settings."
                        .to_string(),
                ),
            );
        }
    };

    // Get transcript
    let cdot_tx = match cdot.get_transcript(transcript_id) {
        Some(tx) => tx,
        None => {
            return (
                None,
                Some(format!("Transcript {} not found", transcript_id)),
            );
        }
    };

    // Extract CDS position and intronic offset
    let (cds_base, offset) = match loc_edit.location.start.inner() {
        Some(p) => (p.base, p.offset),
        None => (1, None),
    };

    // Convert CDS to transcript position
    let tx_pos = match cdot_tx.cds_to_tx(cds_base) {
        Some(pos) => pos,
        None => {
            return (
                None,
                Some("Could not convert CDS to transcript position".to_string()),
            );
        }
    };

    // Handle intronic offset if present (e.g., c.117-2del has offset -2)
    let genomic_pos = if let Some(intron_offset) = offset {
        // This is an intronic variant
        // Get the exon boundary genomic position first
        let exon_boundary_pos = match cdot_tx.tx_to_genome(tx_pos) {
            Some(pos) => pos,
            None => {
                return (None, Some("Exon boundary position not found".to_string()));
            }
        };

        // Apply intron offset based on strand orientation
        // Negative offset (e.g., -2) means position is before the exon in the intron
        // Positive offset (e.g., +5) means position is after the exon in the intron
        match cdot_tx.strand {
            Strand::Plus => {
                if intron_offset < 0 {
                    // Position is before exon (5' splice site)
                    // e.g., c.117-2 is 2 bases before position 117 in the transcript
                    exon_boundary_pos.saturating_sub((-intron_offset) as u64)
                } else {
                    // Position is after exon (3' splice site)
                    exon_boundary_pos.saturating_add(intron_offset as u64)
                }
            }
            Strand::Minus => {
                if intron_offset < 0 {
                    // On minus strand, negative offset goes in positive genomic direction
                    exon_boundary_pos.saturating_add((-intron_offset) as u64)
                } else {
                    // Positive offset goes in negative genomic direction
                    exon_boundary_pos.saturating_sub(intron_offset as u64)
                }
            }
        }
    } else {
        // Exonic variant - direct conversion
        match cdot_tx.tx_to_genome(tx_pos) {
            Some(pos) => pos,
            None => {
                return (None, Some("Position not found in exons".to_string()));
            }
        }
    };

    // Get chromosome
    let chrom = cdot_tx.contig.clone();

    // Extract alleles from edit
    if let Some(edit) = loc_edit.edit.inner() {
        let (ref_allele, alt) = extract_alleles_from_edit(edit);
        if let (Some(ref_a), Some(alt_a)) = (ref_allele, alt) {
            // Add warning for intronic variants using placeholder bases
            let warning = if offset.is_some() {
                Some("Intronic variant: VCF position calculated from intron offset. Reference allele may need verification.".to_string())
            } else {
                None
            };

            return (
                Some(VcfRecord {
                    chrom,
                    pos: genomic_pos,
                    ref_allele: ref_a,
                    alt: alt_a,
                    build: build.as_str().to_string(),
                }),
                warning,
            );
        }
    }

    (
        None,
        Some("Could not extract alleles from variant".to_string()),
    )
}

/// Convert n. variant to VCF via coordinate mapping
fn convert_tx_to_vcf(
    transcript_id: &str,
    loc_edit: &crate::hgvs::variant::LocEdit<
        crate::hgvs::interval::TxInterval,
        crate::hgvs::edit::NaEdit,
    >,
    build: &GenomeBuild,
    cdot: Option<&Arc<CdotMapper>>,
) -> (Option<VcfRecord>, Option<String>) {
    // Need cdot for coordinate conversion
    let cdot = match cdot {
        Some(c) => c,
        None => {
            return (
                None,
                Some(
                    "Converting n. to VCF requires cdot data. Configure 'data.cdot_path' in settings."
                        .to_string(),
                ),
            );
        }
    };

    // Get transcript
    let cdot_tx = match cdot.get_transcript(transcript_id) {
        Some(tx) => tx,
        None => {
            return (
                None,
                Some(format!("Transcript {} not found", transcript_id)),
            );
        }
    };

    // Extract transcript position
    let tx_pos = loc_edit
        .location
        .start
        .inner()
        .map(|p| p.base as u64)
        .unwrap_or(1);

    // Convert to genomic
    let genomic_pos = match cdot_tx.tx_to_genome(tx_pos) {
        Some(pos) => pos,
        None => {
            return (None, Some("Position not found in exons".to_string()));
        }
    };

    // Get chromosome
    let chrom = cdot_tx.contig.clone();

    // Extract alleles from edit
    if let Some(edit) = loc_edit.edit.inner() {
        let (ref_allele, alt) = extract_alleles_from_edit(edit);
        if let (Some(ref_a), Some(alt_a)) = (ref_allele, alt) {
            return (
                Some(VcfRecord {
                    chrom,
                    pos: genomic_pos,
                    ref_allele: ref_a,
                    alt: alt_a,
                    build: build.as_str().to_string(),
                }),
                None,
            );
        }
    }

    (
        None,
        Some("Could not extract alleles from variant".to_string()),
    )
}

/// Extract alleles from HGVS edit
///
/// VCF format requires a "padding" base for indels. When sequence data is not available,
/// "N" is used as a placeholder. The returned alleles are suitable for VCF but may need
/// reference lookup for accurate representation.
///
/// Returns (ref_allele, alt_allele) where either may be None if extraction fails.
fn extract_alleles_from_edit(edit: &crate::hgvs::edit::NaEdit) -> (Option<String>, Option<String>) {
    use crate::hgvs::edit::NaEdit;

    match edit {
        NaEdit::Substitution {
            reference,
            alternative,
        } => (Some(reference.to_string()), Some(alternative.to_string())),
        NaEdit::SubstitutionNoRef { alternative } => {
            // Reference base not specified in HGVS - use "N" placeholder
            // Accurate VCF requires fetching reference sequence
            (Some("N".to_string()), Some(alternative.to_string()))
        }
        NaEdit::Deletion { sequence, .. } => {
            if let Some(seq) = sequence {
                // VCF requires padding base before deletion
                // "N" is placeholder - accurate VCF needs reference lookup
                (Some(format!("N{}", seq)), Some("N".to_string()))
            } else {
                // Deletion without explicit sequence cannot be converted without reference
                (None, None)
            }
        }
        NaEdit::Insertion { sequence } => {
            // VCF requires padding base before insertion
            // "N" is placeholder - accurate VCF needs reference lookup
            (Some("N".to_string()), Some(format!("N{}", sequence)))
        }
        NaEdit::Delins { sequence } => {
            // Delins without original sequence - use "N" placeholder for ref
            (Some("N".to_string()), Some(sequence.to_string()))
        }
        NaEdit::Duplication { sequence, .. } => {
            if let Some(seq) = sequence {
                (Some(seq.to_string()), Some(format!("{}{}", seq, seq)))
            } else {
                (None, None)
            }
        }
        _ => (None, None),
    }
}

/// Extract VCF fields from HGVS edit
fn extract_vcf_fields(
    edit: &crate::hgvs::edit::NaEdit,
    location: &crate::hgvs::interval::GenomeInterval,
) -> (Option<u64>, Option<String>, Option<String>) {
    let pos = location.start.inner().map(|p| p.base);
    let (ref_allele, alt) = extract_alleles_from_edit(edit);
    (pos, ref_allele, alt)
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

/// Extract chromosome from RefSeq accession
fn extract_chromosome_from_accession(accession: &str) -> String {
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
    accession.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_refseq_accession_grch38() {
        assert_eq!(
            get_refseq_accession("chr7", &GenomeBuild::GRCh38),
            Some("NC_000007.14".to_string())
        );
        assert_eq!(
            get_refseq_accession("7", &GenomeBuild::GRCh38),
            Some("NC_000007.14".to_string())
        );
        assert_eq!(
            get_refseq_accession("chrX", &GenomeBuild::GRCh38),
            Some("NC_000023.11".to_string())
        );
    }

    #[test]
    fn test_get_refseq_accession_grch37() {
        assert_eq!(
            get_refseq_accession("chr7", &GenomeBuild::GRCh37),
            Some("NC_000007.13".to_string())
        );
        assert_eq!(
            get_refseq_accession("chr1", &GenomeBuild::GRCh37),
            Some("NC_000001.10".to_string())
        );
    }

    #[test]
    fn test_extract_chromosome_from_accession() {
        assert_eq!(extract_chromosome_from_accession("NC_000007.14"), "chr7");
        assert_eq!(extract_chromosome_from_accession("NC_000001.11"), "chr1");
        assert_eq!(extract_chromosome_from_accession("NC_000023.11"), "chrX");
        assert_eq!(extract_chromosome_from_accession("NC_000024.10"), "chrY");
        assert_eq!(extract_chromosome_from_accession("NC_012920.1"), "chrM");
    }

    #[test]
    fn test_generate_genomic_hgvs_substitution() {
        let (hgvs, var_type, pos) = generate_genomic_hgvs("NC_000007.14", 117559593, "G", "A");
        assert_eq!(hgvs, Some("NC_000007.14:g.117559593G>A".to_string()));
        assert!(matches!(var_type, VariantType::Substitution));
        assert_eq!(pos, 117559593);
    }

    #[test]
    fn test_generate_genomic_hgvs_deletion() {
        let (hgvs, var_type, _) = generate_genomic_hgvs("NC_000007.14", 117559592, "AG", "A");
        assert_eq!(hgvs, Some("NC_000007.14:g.117559593del".to_string()));
        assert!(matches!(var_type, VariantType::Deletion));
    }

    #[test]
    fn test_generate_genomic_hgvs_multi_deletion() {
        let (hgvs, var_type, _) = generate_genomic_hgvs("NC_000007.14", 117559592, "AGT", "A");
        assert_eq!(
            hgvs,
            Some("NC_000007.14:g.117559593_117559594del".to_string())
        );
        assert!(matches!(var_type, VariantType::Deletion));
    }

    #[test]
    fn test_generate_genomic_hgvs_insertion() {
        let (hgvs, var_type, _) = generate_genomic_hgvs("NC_000007.14", 117559593, "G", "GA");
        assert_eq!(
            hgvs,
            Some("NC_000007.14:g.117559593_117559594insA".to_string())
        );
        assert!(matches!(var_type, VariantType::Insertion));
    }

    #[test]
    fn test_unknown_chromosome() {
        assert_eq!(get_refseq_accession("chrW", &GenomeBuild::GRCh38), None);
        assert_eq!(get_refseq_accession("unknown", &GenomeBuild::GRCh37), None);
    }

    #[test]
    fn test_vcf_to_hgvs_without_cdot() {
        // Without cdot, should still produce g. notation but not c./p.
        let (hgvs_g, hgvs_c, hgvs_p, error) = convert_vcf_to_hgvs(
            "chr7",
            117559593,
            "G",
            "A",
            &GenomeBuild::GRCh38,
            Some("NM_000249.4"),
            None,
        );

        assert!(hgvs_g.is_some());
        assert_eq!(hgvs_g.unwrap(), "NC_000007.14:g.117559593G>A");
        assert!(hgvs_c.is_none());
        assert!(hgvs_p.is_none());
        assert!(error.is_some());
        assert!(error.unwrap().contains("cdot"));
    }
}
