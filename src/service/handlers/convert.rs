//! Coordinate conversion endpoint - convert HGVS between coordinate systems (c. ↔ g. ↔ p. ↔ n.)

use axum::{extract::State, http::StatusCode, response::Json};
use std::time::Instant;

use crate::data::cdot::CdsPosition;
use crate::service::handlers::tx_position::resolve_tx_start;
use crate::service::{
    server::AppState,
    types::{
        ConversionResult, ConvertRequest, ConvertResponse, CoordinateSystem, ErrorResponse,
        ServiceError,
    },
    validation::validate_hgvs,
};

/// Convert HGVS variant between coordinate systems
///
/// This endpoint converts an HGVS variant from one coordinate system to another.
/// Supported conversions:
/// - c. (coding) ↔ g. (genomic)
/// - c. (coding) → p. (protein)
/// - n. (non-coding) ↔ g. (genomic)
pub async fn convert(
    State(state): State<AppState>,
    Json(request): Json<ConvertRequest>,
) -> Result<Json<ConvertResponse>, (StatusCode, Json<ErrorResponse>)> {
    let start = Instant::now();

    // Validate input HGVS
    if let Err(validation_error) = validate_hgvs(&request.hgvs) {
        let error = ServiceError::InvalidHgvs(validation_error.to_string());
        return Err((
            StatusCode::from_u16(error.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR),
            Json(error.to_response()),
        ));
    }

    // Parse the input HGVS to detect source coordinate system
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
            // Detect source coordinate system and get transcript accession
            let (source_system, accession) = detect_coordinate_system_and_accession(&result.result);
            let target = request.target_system.as_str();

            // If same system, return as-is
            if source_system == target {
                return Ok(Json(ConvertResponse {
                    input: request.hgvs.clone(),
                    source_system,
                    target_system: target.to_string(),
                    converted: Some(request.hgvs),
                    all_conversions: None,
                    error: None,
                    processing_time_ms: elapsed_ms,
                }));
            }

            // Try to perform the conversion
            let (converted, all_conversions, error) = perform_conversion(
                &state,
                &result.result,
                &accession,
                &source_system,
                &request.target_system,
                request.include_all,
            );

            Ok(Json(ConvertResponse {
                input: request.hgvs,
                source_system,
                target_system: target.to_string(),
                converted,
                all_conversions,
                error,
                processing_time_ms: elapsed_ms,
            }))
        }
        Err(e) => Ok(Json(ConvertResponse {
            input: request.hgvs,
            source_system: "unknown".to_string(),
            target_system: request.target_system.as_str().to_string(),
            converted: None,
            all_conversions: None,
            error: Some(format!("Failed to parse input: {}", e)),
            processing_time_ms: elapsed_ms,
        })),
    }
}

/// Detect the coordinate system and accession from a parsed HGVS variant
fn detect_coordinate_system_and_accession(
    variant: &crate::hgvs::variant::HgvsVariant,
) -> (String, String) {
    use crate::hgvs::variant::HgvsVariant;

    match variant {
        HgvsVariant::Cds(v) => ("c".to_string(), v.accession.to_string()),
        HgvsVariant::Genome(v) => ("g".to_string(), v.accession.to_string()),
        HgvsVariant::Tx(v) => ("n".to_string(), v.accession.to_string()),
        HgvsVariant::Protein(v) => ("p".to_string(), v.accession.to_string()),
        HgvsVariant::Mt(v) => ("m".to_string(), v.accession.to_string()),
        _ => ("unknown".to_string(), String::new()),
    }
}

/// Perform the coordinate conversion
fn perform_conversion(
    state: &AppState,
    variant: &crate::hgvs::variant::HgvsVariant,
    accession: &str,
    source: &str,
    target: &CoordinateSystem,
    include_all: bool,
) -> (
    Option<String>,
    Option<Vec<ConversionResult>>,
    Option<String>,
) {
    use crate::hgvs::variant::HgvsVariant;

    // Check if we have cdot data for transcript lookups
    let cdot = match &state.cdot {
        Some(c) => c,
        None => {
            return (
                None,
                None,
                Some(
                    "Coordinate conversion requires transcript data. Configure 'data.cdot_path' \
                     in the service configuration with a path to a cdot JSON file."
                        .to_string(),
                ),
            );
        }
    };

    match variant {
        HgvsVariant::Cds(v) => {
            // Get transcript from cdot
            let cdot_tx = match cdot.get_transcript(accession) {
                Some(tx) => tx,
                None => {
                    return (
                        None,
                        None,
                        Some(format!("Transcript {} not found in cdot data", accession)),
                    );
                }
            };

            match target {
                CoordinateSystem::G => {
                    // c. -> g. conversion
                    let cds_base = v
                        .loc_edit
                        .location
                        .start
                        .inner()
                        .map(|p| p.base)
                        .unwrap_or(1);

                    // Convert CDS to transcript position, then to genomic
                    if let Some(tx_pos) = cdot_tx.cds_to_tx(cds_base) {
                        if let Some(genomic_pos) = cdot_tx.tx_to_genome(tx_pos) {
                            let chrom = &cdot_tx.contig;
                            let converted = format!("{}:g.{}", chrom, genomic_pos);

                            let all = if include_all {
                                Some(vec![ConversionResult {
                                    system: "g".to_string(),
                                    hgvs: converted.clone(),
                                    reference: chrom.to_string(),
                                }])
                            } else {
                                None
                            };

                            return (Some(converted), all, None);
                        }
                    }
                    (
                        None,
                        None,
                        Some("Position not found in transcript".to_string()),
                    )
                }
                CoordinateSystem::P => {
                    // c. -> p. conversion
                    let cds_base = v
                        .loc_edit
                        .location
                        .start
                        .inner()
                        .map(|p| p.base)
                        .unwrap_or(1);

                    // Calculate protein position: (cds_pos - 1) / 3 + 1
                    if cds_base > 0 {
                        let codon_num = ((cds_base - 1) / 3 + 1) as u64;
                        // Protein accession: the authoritative cdot value if
                        // present, else the transcript accession itself. We do
                        // NOT infer `NP_*`/`XP_*` from `NM_*`/`XM_*` by
                        // preserving the number — RefSeq does not guarantee the
                        // NM and NP numbers match, so that is frequently wrong
                        // (#808).
                        let prot_acc = cdot_tx
                            .protein
                            .clone()
                            .unwrap_or_else(|| accession.to_string());
                        let converted = format!("{}:p.{}", prot_acc, codon_num);

                        let all = if include_all {
                            Some(vec![ConversionResult {
                                system: "p".to_string(),
                                hgvs: converted.clone(),
                                reference: prot_acc,
                            }])
                        } else {
                            None
                        };

                        (Some(converted), all, None)
                    } else {
                        (
                            None,
                            None,
                            Some("UTR positions cannot be converted to protein".to_string()),
                        )
                    }
                }
                CoordinateSystem::N => {
                    // c. -> n. conversion
                    let cds_base = v
                        .loc_edit
                        .location
                        .start
                        .inner()
                        .map(|p| p.base)
                        .unwrap_or(1);

                    if let Some(tx_pos) = cdot_tx.cds_to_tx(cds_base) {
                        let converted = format!("{}:n.{}", accession, tx_pos);

                        let all = if include_all {
                            Some(vec![ConversionResult {
                                system: "n".to_string(),
                                hgvs: converted.clone(),
                                reference: accession.to_string(),
                            }])
                        } else {
                            None
                        };

                        (Some(converted), all, None)
                    } else {
                        (
                            None,
                            None,
                            Some("Could not convert CDS to transcript position".to_string()),
                        )
                    }
                }
                CoordinateSystem::C => {
                    // Already c., shouldn't happen due to same-system check
                    (Some(format!("{}", v)), None, None)
                }
            }
        }
        HgvsVariant::Genome(v) => {
            // g. -> c./n. requires finding overlapping transcripts
            let genome_pos = v
                .loc_edit
                .location
                .start
                .inner()
                .map(|p| p.base)
                .unwrap_or(0);
            let chrom = v.accession.to_string();

            // Extract chromosome name from accession (e.g., NC_000007.14 -> chr7)
            let chrom_name = extract_chromosome_from_accession(&chrom);

            // Find transcripts at this position
            let transcripts = cdot.transcripts_at_position(&chrom_name, genome_pos);

            if transcripts.is_empty() {
                return (
                    None,
                    None,
                    Some(format!(
                        "No transcripts found at {}:{}",
                        chrom_name, genome_pos
                    )),
                );
            }

            match target {
                CoordinateSystem::C | CoordinateSystem::N => {
                    // Try to convert using the first transcript found
                    let (tx_id, cdot_tx) = &transcripts[0];

                    if let Some(tx_pos) = cdot_tx.genome_to_tx(genome_pos) {
                        let converted = match target {
                            CoordinateSystem::C => {
                                if let Some(cds_pos) = cdot_tx.tx_to_cds(tx_pos) {
                                    match cds_pos {
                                        CdsPosition::Cds(pos) => {
                                            format!("{}:c.{}", tx_id, pos)
                                        }
                                        CdsPosition::FivePrimeUtr(offset) => {
                                            format!("{}:c.-{}", tx_id, offset)
                                        }
                                        CdsPosition::ThreePrimeUtr(offset) => {
                                            format!("{}:c.*{}", tx_id, offset)
                                        }
                                    }
                                } else {
                                    return (None, None, Some("Transcript has no CDS".to_string()));
                                }
                            }
                            CoordinateSystem::N => format!("{}:n.{}", tx_id, tx_pos),
                            _ => unreachable!(),
                        };

                        let all = if include_all && transcripts.len() > 1 {
                            Some(
                                transcripts
                                    .iter()
                                    .filter_map(|(id, tx)| {
                                        tx.genome_to_tx(genome_pos).map(|pos| {
                                            let hgvs = match target {
                                                CoordinateSystem::C => tx
                                                    .tx_to_cds(pos)
                                                    .map(|cds| match cds {
                                                        CdsPosition::Cds(p) => {
                                                            format!("{}:c.{}", id, p)
                                                        }
                                                        CdsPosition::FivePrimeUtr(o) => {
                                                            format!("{}:c.-{}", id, o)
                                                        }
                                                        CdsPosition::ThreePrimeUtr(o) => {
                                                            format!("{}:c.*{}", id, o)
                                                        }
                                                    })
                                                    .unwrap_or_else(|| format!("{}:n.{}", id, pos)),
                                                CoordinateSystem::N => format!("{}:n.{}", id, pos),
                                                _ => format!("{}:{}", id, pos),
                                            };
                                            ConversionResult {
                                                system: target.as_str().to_string(),
                                                hgvs,
                                                reference: id.to_string(),
                                            }
                                        })
                                    })
                                    .collect(),
                            )
                        } else {
                            None
                        };

                        (Some(converted), all, None)
                    } else {
                        (
                            None,
                            None,
                            Some("Position not in transcript exons".to_string()),
                        )
                    }
                }
                _ => (
                    None,
                    None,
                    Some(format!("Conversion from g. to {} is not supported", target)),
                ),
            }
        }
        HgvsVariant::Tx(v) => {
            // n. -> c./g. conversion
            let cdot_tx = match cdot.get_transcript(accession) {
                Some(tx) => tx,
                None => {
                    return (
                        None,
                        None,
                        Some(format!("Transcript {} not found in cdot data", accession)),
                    );
                }
            };

            convert_tx_position(cdot_tx, &v.loc_edit.location, accession, target)
        }
        _ => (
            None,
            None,
            Some(format!(
                "Conversion from {} to {} is not supported",
                source, target
            )),
        ),
    }
}

/// Convert an `n.` position to the requested target system.
///
/// Split out of [`perform_conversion`]'s `Tx` arm so the conversion can be
/// exercised against a hand-built [`CdotTranscript`] without standing up an
/// `AppState`: the arm's only use of the service state is the transcript
/// lookup, which stays with the caller.
///
/// The start position is resolved through [`resolve_tx_start`], which refuses
/// `n.*N`. Reading `base` alone — as this arm did — made `n.5` and `n.*5`
/// resolve to the same coordinate, so the handler reported a `c.` or `g.`
/// position for a different nucleotide with no diagnostic.
///
/// [`CdotTranscript`]: crate::data::cdot::CdotTranscript
fn convert_tx_position(
    cdot_tx: &crate::data::cdot::CdotTranscript,
    location: &crate::hgvs::interval::TxInterval,
    accession: &str,
    target: &CoordinateSystem,
) -> (
    Option<String>,
    Option<Vec<ConversionResult>>,
    Option<String>,
) {
    let tx_pos = match resolve_tx_start(location) {
        Ok(pos) => pos,
        Err(msg) => return (None, None, Some(msg)),
    };

    match target {
        CoordinateSystem::C => {
            // n. -> c. conversion
            if let Some(cds_pos) = cdot_tx.tx_to_cds(tx_pos) {
                let converted = match cds_pos {
                    CdsPosition::Cds(pos) => format!("{}:c.{}", accession, pos),
                    CdsPosition::FivePrimeUtr(offset) => {
                        format!("{}:c.-{}", accession, offset)
                    }
                    CdsPosition::ThreePrimeUtr(offset) => {
                        format!("{}:c.*{}", accession, offset)
                    }
                };
                (Some(converted), None, None)
            } else {
                (None, None, Some("Transcript has no CDS".to_string()))
            }
        }
        CoordinateSystem::G => {
            // n. -> g. conversion
            if let Some(genomic_pos) = cdot_tx.tx_to_genome(tx_pos) {
                let chrom = &cdot_tx.contig;
                let converted = format!("{}:g.{}", chrom, genomic_pos);
                (Some(converted), None, None)
            } else {
                (None, None, Some("Position not in exons".to_string()))
            }
        }
        _ => (
            None,
            None,
            Some(format!("Conversion from n. to {} is not supported", target)),
        ),
    }
}

/// Extract chromosome name from RefSeq accession
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
    use crate::data::cdot::CdotTranscript;
    use crate::hgvs::interval::TxInterval;
    use crate::hgvs::location::TxPos;
    use crate::reference::Strand;

    /// One 9-base exon at genome 1000..1009, CDS spanning the whole of it, so
    /// every `n.` base in range converts on both target axes and a divergence
    /// cannot be an artifact of one of them declining.
    fn single_exon_transcript() -> CdotTranscript {
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            exons: vec![[1000, 1009, 0, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_TEST.1".to_string()),
            exon_cigars: Vec::new(),
        }
    }

    /// The drop route, at the handler seam: `n.5` and `n.*5` are different
    /// nucleotides, so the conversion must not answer the same string for both.
    ///
    /// Before the fix this arm read `base` alone, so both sides converted as
    /// in-transcript position 5 — a *collapse* — on both target axes. Asserted
    /// as a divergence rather than against a pinned message so it cannot be
    /// satisfied by a wording change.
    #[test]
    fn n_to_c_and_n_to_g_do_not_answer_a_downstream_position_as_its_twin() {
        let tx = single_exon_transcript();

        for target in [CoordinateSystem::C, CoordinateSystem::G] {
            let convert = |pos: TxPos| {
                let (converted, _all, error) =
                    convert_tx_position(&tx, &TxInterval::point(pos), "NM_TEST.1", &target);
                (converted, error)
            };
            let (plain, plain_err) = convert(TxPos::new(5));
            let (downstream, downstream_err) = convert(TxPos::downstream(5));

            assert!(
                plain.is_some() && plain_err.is_none(),
                "control: a plain in-transcript n.5 must still convert to {target}, \
                 got {plain:?} / {plain_err:?}"
            );
            assert_ne!(
                plain, downstream,
                "n.5 and n.*5 are different nucleotides and must not convert to \
                 the same {target} answer"
            );
            assert!(
                downstream.is_none(),
                "n.*5 names a nucleotide past the transcript's last base and has \
                 no {target} position; the handler answered {downstream:?} instead \
                 of refusing"
            );
            assert!(
                downstream_err.is_some_and(|m| m.contains("n.*")),
                "the decline must name the notation it refused"
            );
        }
    }

    #[test]
    fn test_detect_coordinate_system() {
        // Test with a simple variant string
        let result = crate::hgvs::parser::parse_hgvs_lenient("NM_000249.4:c.350C>T").unwrap();
        let (system, acc) = detect_coordinate_system_and_accession(&result.result);
        assert_eq!(system, "c");
        assert_eq!(acc, "NM_000249.4");
    }

    #[test]
    fn test_detect_genomic_system() {
        let result =
            crate::hgvs::parser::parse_hgvs_lenient("NC_000007.14:g.117559593G>A").unwrap();
        let (system, acc) = detect_coordinate_system_and_accession(&result.result);
        assert_eq!(system, "g");
        assert_eq!(acc, "NC_000007.14");
    }

    #[test]
    fn test_extract_chromosome_from_accession() {
        assert_eq!(extract_chromosome_from_accession("NC_000001.10"), "chr1");
        assert_eq!(extract_chromosome_from_accession("NC_000007.13"), "chr7");
        assert_eq!(extract_chromosome_from_accession("NC_000023.10"), "chrX");
        assert_eq!(extract_chromosome_from_accession("NC_000024.09"), "chrY");
        assert_eq!(extract_chromosome_from_accession("NC_012920.1"), "chrM");
    }

    #[test]
    fn test_detect_tx_system() {
        let result = crate::hgvs::parser::parse_hgvs_lenient("NM_000249.4:n.500A>G").unwrap();
        let (system, acc) = detect_coordinate_system_and_accession(&result.result);
        assert_eq!(system, "n");
        assert_eq!(acc, "NM_000249.4");
    }
}
