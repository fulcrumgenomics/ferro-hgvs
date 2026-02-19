//! UTA database loading from cdot transcript alignments.
//!
//! This module loads missing transcript alignments from cdot JSON into
//! the UTA PostgreSQL database using the PostgreSQL COPY protocol for
//! high-performance bulk loading.
//!
//! ## Performance
//! - Old approach (docker exec per INSERT): ~20 rows/min = 12 days for 344K transcripts
//! - New approach (tokio-postgres COPY): ~50K rows/sec = 2-5 minutes for 344K transcripts

use crate::data::cdot::CdotMapper;
use crate::reference::Strand;
use crate::FerroError;
use bytes::{BufMut, Bytes, BytesMut};
use futures_util::SinkExt;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::process::Command;
use std::time::Instant;
use tokio_postgres::{Client, CopyInSink, NoTls};

/// Configuration for loading cdot alignments into UTA
#[derive(Debug, Clone)]
pub struct UtaLoadConfig {
    /// UTA schema name (e.g., "uta_20210129b")
    pub uta_schema: String,
    /// Docker container name
    pub container_name: String,
    /// Origin name for new transcripts
    pub origin_name: String,
    /// Batch size for progress reporting
    pub batch_size: usize,
}

impl Default for UtaLoadConfig {
    fn default() -> Self {
        Self {
            uta_schema: "uta_20210129b".to_string(),
            container_name: "ferro-uta".to_string(),
            origin_name: "cdot-ferro".to_string(),
            batch_size: 10000, // Larger batches for COPY
        }
    }
}

/// Result of loading cdot alignments into UTA
#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct UtaLoadResult {
    /// Number of transcripts loaded
    pub transcripts_loaded: usize,
    /// Number of transcripts skipped (already exist)
    pub transcripts_skipped: usize,
    /// Number of exon sets created
    pub exon_sets_created: usize,
    /// Number of exons created
    pub exons_created: usize,
    /// Number of exon alignments created
    pub exon_alns_created: usize,
    /// Errors encountered (non-fatal)
    pub errors: Vec<String>,
}

/// Load cdot transcript alignments into UTA database
///
/// This is a blocking wrapper around the async implementation.
/// Creates a new tokio runtime for this operation.
pub fn load_cdot_alignments_to_uta<P: AsRef<Path>>(
    cdot_path: P,
    config: &UtaLoadConfig,
) -> Result<UtaLoadResult, FerroError> {
    // Create a new tokio runtime for this operation
    let rt = tokio::runtime::Runtime::new().map_err(|e| FerroError::Io {
        msg: format!("Failed to create tokio runtime: {}", e),
    })?;

    rt.block_on(load_cdot_alignments_async(cdot_path.as_ref(), config))
}

/// Async implementation using tokio-postgres COPY protocol
async fn load_cdot_alignments_async(
    cdot_path: &Path,
    config: &UtaLoadConfig,
) -> Result<UtaLoadResult, FerroError> {
    let mut result = UtaLoadResult::default();

    // Phase 0: Load cdot and determine what needs loading
    eprintln!("Loading cdot transcript data...");
    let cdot = CdotMapper::from_json_file(cdot_path)?;

    let cdot_accessions: Vec<String> = cdot
        .transcript_ids()
        .filter(|ac| {
            ac.starts_with("NM_")
                || ac.starts_with("NR_")
                || ac.starts_with("XM_")
                || ac.starts_with("XR_")
        })
        .map(String::from)
        .collect();
    eprintln!("Found {} transcripts in cdot", cdot_accessions.len());

    // Use docker exec for initial query (simpler, one-time)
    let existing = get_existing_transcripts(&config.container_name, &config.uta_schema)?;
    eprintln!("Found {} transcripts already in UTA", existing.len());

    let missing: Vec<&String> = cdot_accessions
        .iter()
        .filter(|ac| !existing.contains(*ac))
        .collect();

    result.transcripts_skipped = cdot_accessions.len() - missing.len();
    eprintln!("Need to load {} missing transcripts", missing.len());

    if missing.is_empty() {
        eprintln!("Nothing to load - all transcripts already in UTA");
        return Ok(result);
    }

    // Connect to PostgreSQL via container IP
    let client = connect_to_uta(&config.container_name).await?;

    // Get or create origin_id
    let origin_id =
        get_or_create_origin_async(&client, &config.uta_schema, &config.origin_name).await?;

    let schema = &config.uta_schema;
    let total_start = Instant::now();

    // Phase 1: COPY transcripts
    eprintln!("\nPhase 1/6: Loading transcripts via COPY...");
    let phase_start = Instant::now();
    let tx_count = copy_transcripts(&client, schema, origin_id, &cdot, &missing).await?;
    result.transcripts_loaded = tx_count;
    eprintln!(
        "  Loaded {} transcripts in {:.1}s ({:.0}/s)",
        tx_count,
        phase_start.elapsed().as_secs_f64(),
        tx_count as f64 / phase_start.elapsed().as_secs_f64()
    );

    // Phase 2: COPY exon_sets (2 per transcript)
    eprintln!("Phase 2/6: Loading exon_sets via COPY...");
    let phase_start = Instant::now();
    let es_count = copy_exon_sets(&client, schema, &cdot, &missing).await?;
    result.exon_sets_created = es_count;
    eprintln!(
        "  Created {} exon_sets in {:.1}s ({:.0}/s)",
        es_count,
        phase_start.elapsed().as_secs_f64(),
        es_count as f64 / phase_start.elapsed().as_secs_f64()
    );

    // Phase 3: Query back exon_set_ids
    eprintln!("Phase 3/6: Building exon_set ID lookup...");
    let phase_start = Instant::now();
    let exon_set_map = build_exon_set_map(&client, schema, &missing).await?;
    eprintln!(
        "  Mapped {} exon_set_ids in {:.1}s",
        exon_set_map.len(),
        phase_start.elapsed().as_secs_f64()
    );

    // Phase 4: COPY exons (2 per exon per transcript)
    eprintln!("Phase 4/6: Loading exons via COPY...");
    let phase_start = Instant::now();
    let exon_count = copy_exons(&client, schema, &cdot, &missing, &exon_set_map).await?;
    result.exons_created = exon_count;
    eprintln!(
        "  Created {} exons in {:.1}s ({:.0}/s)",
        exon_count,
        phase_start.elapsed().as_secs_f64(),
        exon_count as f64 / phase_start.elapsed().as_secs_f64()
    );

    // Phase 5: Query back exon_ids
    eprintln!("Phase 5/6: Building exon ID lookup...");
    let phase_start = Instant::now();
    let exon_map = build_exon_map(&client, schema, &exon_set_map).await?;
    eprintln!(
        "  Mapped {} exon_ids in {:.1}s",
        exon_map.len(),
        phase_start.elapsed().as_secs_f64()
    );

    // Phase 6: COPY exon_alns
    eprintln!("Phase 6/6: Loading exon alignments via COPY...");
    let phase_start = Instant::now();
    let aln_count =
        copy_exon_alns(&client, schema, &cdot, &missing, &exon_set_map, &exon_map).await?;
    result.exon_alns_created = aln_count;
    eprintln!(
        "  Created {} exon alignments in {:.1}s ({:.0}/s)",
        aln_count,
        phase_start.elapsed().as_secs_f64(),
        aln_count as f64 / phase_start.elapsed().as_secs_f64()
    );

    // Refresh materialized views
    eprintln!("\nRefreshing UTA materialized views...");
    let phase_start = Instant::now();
    refresh_views_async(&client, schema).await?;
    eprintln!(
        "  Materialized views refreshed in {:.1}s",
        phase_start.elapsed().as_secs_f64()
    );

    eprintln!(
        "\nUTA loading complete in {:.1}s total",
        total_start.elapsed().as_secs_f64()
    );

    Ok(result)
}

/// Connect to PostgreSQL via localhost (requires port mapping)
async fn connect_to_uta(container: &str) -> Result<Client, FerroError> {
    // Check that container is running and has port 5432 mapped
    let output = Command::new("docker")
        .args(["port", container, "5432"])
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Docker port failed: {}", e),
        })?;

    if !output.status.success() {
        return Err(FerroError::Io {
            msg: format!(
                "Container {} doesn't have port 5432 mapped. Start with: docker run -p 5432:5432 ...",
                container
            ),
        });
    }

    // Parse the port mapping (e.g., "0.0.0.0:5432" -> 5432)
    let port_output = String::from_utf8_lossy(&output.stdout);
    let port: u16 = port_output
        .lines()
        .next()
        .and_then(|line| line.rsplit(':').next())
        .and_then(|p| p.trim().parse().ok())
        .unwrap_or(5432);

    eprintln!("Connecting to PostgreSQL at localhost:{}...", port);

    let conn_str = format!("host=localhost port={} user=postgres dbname=uta", port);

    let (client, connection) = tokio_postgres::connect(&conn_str, NoTls)
        .await
        .map_err(|e| FerroError::Io {
            msg: format!(
                "PostgreSQL connect failed: {}. Is the container running?",
                e
            ),
        })?;

    // Spawn the connection handler
    tokio::spawn(async move {
        if let Err(e) = connection.await {
            eprintln!("PostgreSQL connection error: {}", e);
        }
    });

    eprintln!("  Connected to UTA database");
    Ok(client)
}

/// COPY transcripts using text protocol
///
/// Text COPY format: tab-separated, \N for NULL, newline-terminated
async fn copy_transcripts(
    client: &Client,
    schema: &str,
    origin_id: i32,
    cdot: &CdotMapper,
    missing: &[&String],
) -> Result<usize, FerroError> {
    let copy_stmt = format!(
        "COPY {}.transcript (ac, origin_id, hgnc, cds_start_i, cds_end_i, added) FROM STDIN WITH (FORMAT text)",
        schema
    );

    let sink: CopyInSink<Bytes> = client
        .copy_in(&copy_stmt)
        .await
        .map_err(|e| FerroError::Io {
            msg: format!("COPY IN failed: {}", e),
        })?;

    futures_util::pin_mut!(sink);
    let mut count = 0;
    let mut buffer = BytesMut::with_capacity(64 * 1024); // 64KB buffer

    for accession in missing {
        if let Some(tx) = cdot.get_transcript(accession) {
            // Escape tabs and newlines in gene name, use \N for NULL
            let gene = tx
                .gene_name
                .as_ref()
                .map(|g| escape_copy_text(g))
                .unwrap_or_else(|| "\\N".to_string());

            // cds_start_i and cds_end_i can be NULL for non-coding transcripts
            let cds_start = tx
                .cds_start
                .map(|v| v.to_string())
                .unwrap_or_else(|| "\\N".to_string());
            let cds_end = tx
                .cds_end
                .map(|v| v.to_string())
                .unwrap_or_else(|| "\\N".to_string());

            // Format: ac \t origin_id \t hgnc \t cds_start_i \t cds_end_i \t added
            let line = format!(
                "{}\t{}\t{}\t{}\t{}\tnow\n",
                accession, origin_id, gene, cds_start, cds_end
            );
            buffer.put_slice(line.as_bytes());
            count += 1;

            // Flush buffer periodically (every 50K rows)
            if count % 50000 == 0 {
                sink.send(buffer.split().freeze())
                    .await
                    .map_err(|e| FerroError::Io {
                        msg: format!("COPY send failed: {}", e),
                    })?;
                eprintln!("  Progress: {}/{} transcripts...", count, missing.len());
            }
        }
    }

    // Send remaining data
    if !buffer.is_empty() {
        sink.send(buffer.freeze())
            .await
            .map_err(|e| FerroError::Io {
                msg: format!("COPY send failed: {}", e),
            })?;
    }

    // Finish the COPY operation
    let rows = sink.finish().await.map_err(|e| FerroError::Io {
        msg: format!("COPY finish failed: {}", e),
    })?;

    Ok(rows as usize)
}

/// Escape special characters for PostgreSQL text COPY format
fn escape_copy_text(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
}

/// COPY exon_sets (2 per transcript: self-alignment + genomic alignment)
async fn copy_exon_sets(
    client: &Client,
    schema: &str,
    cdot: &CdotMapper,
    missing: &[&String],
) -> Result<usize, FerroError> {
    let copy_stmt = format!(
        "COPY {}.exon_set (tx_ac, alt_ac, alt_strand, alt_aln_method, added) FROM STDIN WITH (FORMAT text)",
        schema
    );

    let sink: CopyInSink<Bytes> = client
        .copy_in(&copy_stmt)
        .await
        .map_err(|e| FerroError::Io {
            msg: format!("COPY IN failed: {}", e),
        })?;

    futures_util::pin_mut!(sink);
    let mut count = 0;
    let mut buffer = BytesMut::with_capacity(64 * 1024);

    for accession in missing {
        if let Some(tx) = cdot.get_transcript(accession) {
            // 1. Transcript exon_set (transcript aligned to itself)
            let line1 = format!("{}\t{}\t1\ttranscript\tnow\n", accession, accession);
            buffer.put_slice(line1.as_bytes());

            // 2. Genomic exon_set (transcript aligned to chromosome)
            // Use 'splign' as alt_aln_method for hgvs-rs compatibility
            let strand: i16 = match tx.strand {
                Strand::Plus => 1,
                Strand::Minus => -1,
            };
            let line2 = format!("{}\t{}\t{}\tsplign\tnow\n", accession, tx.contig, strand);
            buffer.put_slice(line2.as_bytes());

            count += 2;

            if count % 100000 == 0 {
                sink.send(buffer.split().freeze())
                    .await
                    .map_err(|e| FerroError::Io {
                        msg: format!("COPY send failed: {}", e),
                    })?;
                eprintln!("  Progress: {}/{} exon_sets...", count, missing.len() * 2);
            }
        }
    }

    if !buffer.is_empty() {
        sink.send(buffer.freeze())
            .await
            .map_err(|e| FerroError::Io {
                msg: format!("COPY send failed: {}", e),
            })?;
    }

    sink.finish().await.map_err(|e| FerroError::Io {
        msg: format!("COPY finish failed: {}", e),
    })?;

    Ok(count)
}

/// Build map of (tx_ac, alt_ac, method) -> exon_set_id
async fn build_exon_set_map(
    client: &Client,
    schema: &str,
    missing: &[&String],
) -> Result<HashMap<(String, String, String), i32>, FerroError> {
    let mut map = HashMap::new();
    let total_chunks = missing.len().div_ceil(1000);
    let mut chunk_num = 0;

    // Query in batches to avoid huge IN clauses
    for chunk in missing.chunks(1000) {
        chunk_num += 1;
        let placeholders: Vec<String> = chunk.iter().map(|ac| format!("'{}'", ac)).collect();

        let sql = format!(
            "SELECT exon_set_id, tx_ac, alt_ac, alt_aln_method
             FROM {}.exon_set
             WHERE tx_ac IN ({})",
            schema,
            placeholders.join(",")
        );

        let rows = client.query(&sql, &[]).await.map_err(|e| FerroError::Io {
            msg: format!("Query exon_set failed: {}", e),
        })?;

        for row in rows {
            let id: i32 = row.get(0);
            let tx_ac: String = row.get(1);
            let alt_ac: String = row.get(2);
            let method: String = row.get(3);
            map.insert((tx_ac, alt_ac, method), id);
        }

        if chunk_num % 50 == 0 {
            eprintln!(
                "  Query batch {}/{} ({} IDs so far)...",
                chunk_num,
                total_chunks,
                map.len()
            );
        }
    }

    Ok(map)
}

/// COPY exons (2 per exon: transcript coords + genomic coords)
async fn copy_exons(
    client: &Client,
    schema: &str,
    cdot: &CdotMapper,
    missing: &[&String],
    exon_set_map: &HashMap<(String, String, String), i32>,
) -> Result<usize, FerroError> {
    let copy_stmt = format!(
        "COPY {}.exon (exon_set_id, start_i, end_i, ord) FROM STDIN WITH (FORMAT text)",
        schema
    );

    let sink: CopyInSink<Bytes> = client
        .copy_in(&copy_stmt)
        .await
        .map_err(|e| FerroError::Io {
            msg: format!("COPY IN failed: {}", e),
        })?;

    futures_util::pin_mut!(sink);
    let mut count = 0;
    let mut buffer = BytesMut::with_capacity(64 * 1024);

    for accession in missing {
        if let Some(tx) = cdot.get_transcript(accession) {
            let tx_es_key = (
                accession.to_string(),
                accession.to_string(),
                "transcript".to_string(),
            );
            let alt_es_key = (
                accession.to_string(),
                tx.contig.clone(),
                "splign".to_string(),
            );

            let tx_es_id = exon_set_map.get(&tx_es_key).ok_or_else(|| FerroError::Io {
                msg: format!("Missing transcript exon_set for {}", accession),
            })?;
            let alt_es_id = exon_set_map
                .get(&alt_es_key)
                .ok_or_else(|| FerroError::Io {
                    msg: format!(
                        "Missing genomic exon_set for {} on {}",
                        accession, tx.contig
                    ),
                })?;

            for (ord, exon) in tx.get_exons().iter().enumerate() {
                // Transcript exon (convert 1-based tx_start to 0-based)
                let tx_start = (exon.tx_start - 1) as i32;
                let tx_end = exon.tx_end as i32;
                let line1 = format!("{}\t{}\t{}\t{}\n", tx_es_id, tx_start, tx_end, ord);
                buffer.put_slice(line1.as_bytes());

                // Genomic exon (already 0-based)
                let line2 = format!(
                    "{}\t{}\t{}\t{}\n",
                    alt_es_id, exon.genome_start, exon.genome_end, ord
                );
                buffer.put_slice(line2.as_bytes());

                count += 2;
            }

            if count % 500000 == 0 {
                sink.send(buffer.split().freeze())
                    .await
                    .map_err(|e| FerroError::Io {
                        msg: format!("COPY send failed: {}", e),
                    })?;
                eprintln!("  Progress: {} exons...", count);
            }
        }
    }

    if !buffer.is_empty() {
        sink.send(buffer.freeze())
            .await
            .map_err(|e| FerroError::Io {
                msg: format!("COPY send failed: {}", e),
            })?;
    }

    sink.finish().await.map_err(|e| FerroError::Io {
        msg: format!("COPY finish failed: {}", e),
    })?;

    Ok(count)
}

/// Build map of (exon_set_id, start_i) -> exon_id
async fn build_exon_map(
    client: &Client,
    schema: &str,
    exon_set_map: &HashMap<(String, String, String), i32>,
) -> Result<HashMap<(i32, i64), i32>, FerroError> {
    let mut map = HashMap::new();
    let es_ids: Vec<i32> = exon_set_map.values().copied().collect();
    let total_chunks = es_ids.len().div_ceil(2000);
    let mut chunk_num = 0;

    // Query in batches
    for chunk in es_ids.chunks(2000) {
        chunk_num += 1;
        let ids: Vec<String> = chunk.iter().map(|id| id.to_string()).collect();

        let sql = format!(
            "SELECT exon_id, exon_set_id, start_i
             FROM {}.exon
             WHERE exon_set_id IN ({})",
            schema,
            ids.join(",")
        );

        let rows = client.query(&sql, &[]).await.map_err(|e| FerroError::Io {
            msg: format!("Query exon failed: {}", e),
        })?;

        for row in rows {
            let exon_id: i32 = row.get(0);
            let es_id: i32 = row.get(1);
            let start_i: i32 = row.get(2);
            map.insert((es_id, start_i as i64), exon_id);
        }

        if chunk_num % 50 == 0 {
            eprintln!(
                "  Query batch {}/{} ({} IDs so far)...",
                chunk_num,
                total_chunks,
                map.len()
            );
        }
    }

    Ok(map)
}

/// COPY exon_alns (1 per exon: links tx_exon to alt_exon)
async fn copy_exon_alns(
    client: &Client,
    schema: &str,
    cdot: &CdotMapper,
    missing: &[&String],
    exon_set_map: &HashMap<(String, String, String), i32>,
    exon_map: &HashMap<(i32, i64), i32>,
) -> Result<usize, FerroError> {
    let copy_stmt = format!(
        "COPY {}.exon_aln (tx_exon_id, alt_exon_id, cigar, added) FROM STDIN WITH (FORMAT text)",
        schema
    );

    let sink: CopyInSink<Bytes> = client
        .copy_in(&copy_stmt)
        .await
        .map_err(|e| FerroError::Io {
            msg: format!("COPY IN failed: {}", e),
        })?;

    futures_util::pin_mut!(sink);
    let mut count = 0;
    let mut buffer = BytesMut::with_capacity(64 * 1024);

    for accession in missing {
        if let Some(tx) = cdot.get_transcript(accession) {
            let tx_es_key = (
                accession.to_string(),
                accession.to_string(),
                "transcript".to_string(),
            );
            let alt_es_key = (
                accession.to_string(),
                tx.contig.clone(),
                "splign".to_string(),
            );

            let tx_es_id = *exon_set_map.get(&tx_es_key).ok_or_else(|| FerroError::Io {
                msg: format!(
                    "Missing exon set for transcript {} (key: {:?})",
                    accession, tx_es_key
                ),
            })?;
            let alt_es_id = *exon_set_map
                .get(&alt_es_key)
                .ok_or_else(|| FerroError::Io {
                    msg: format!(
                        "Missing exon set for alignment {} on {} (key: {:?})",
                        accession, tx.contig, alt_es_key
                    ),
                })?;

            for exon in tx.get_exons().iter() {
                let tx_start = (exon.tx_start - 1) as i64;

                let tx_exon_id =
                    exon_map
                        .get(&(tx_es_id, tx_start))
                        .ok_or_else(|| FerroError::Io {
                            msg: format!(
                                "Missing tx exon for {} at start_i={}",
                                accession, tx_start
                            ),
                        })?;

                let alt_exon_id = exon_map
                    .get(&(alt_es_id, exon.genome_start as i64))
                    .ok_or_else(|| FerroError::Io {
                        msg: format!(
                            "Missing alt exon for {} at start_i={}",
                            accession, exon.genome_start
                        ),
                    })?;

                // CIGAR: simple match for the exon length
                let cigar = format!("{}=", exon.genome_end - exon.genome_start);
                let line = format!("{}\t{}\t{}\tnow\n", tx_exon_id, alt_exon_id, cigar);
                buffer.put_slice(line.as_bytes());

                count += 1;
            }

            if count % 250000 == 0 {
                sink.send(buffer.split().freeze())
                    .await
                    .map_err(|e| FerroError::Io {
                        msg: format!("COPY send failed: {}", e),
                    })?;
                eprintln!("  Progress: {} alignments...", count);
            }
        }
    }

    if !buffer.is_empty() {
        sink.send(buffer.freeze())
            .await
            .map_err(|e| FerroError::Io {
                msg: format!("COPY send failed: {}", e),
            })?;
    }

    sink.finish().await.map_err(|e| FerroError::Io {
        msg: format!("COPY finish failed: {}", e),
    })?;

    Ok(count)
}

/// Refresh UTA materialized views (required for hgvs-rs to find new transcripts)
async fn refresh_views_async(client: &Client, schema: &str) -> Result<(), FerroError> {
    // Must refresh in dependency order
    let views = [
        "exon_set_exons_fp_mv",
        "tx_exon_set_summary_mv",
        "tx_def_summary_mv",
    ];

    for view in views {
        let sql = format!("REFRESH MATERIALIZED VIEW {}.{}", schema, view);
        client
            .execute(&sql, &[])
            .await
            .map_err(|e| FerroError::Io {
                msg: format!("Refresh {} failed: {}", view, e),
            })?;
        eprintln!("  Refreshed {}", view);
    }

    Ok(())
}

/// Get or create origin_id for cdot-loaded transcripts
async fn get_or_create_origin_async(
    client: &Client,
    schema: &str,
    name: &str,
) -> Result<i32, FerroError> {
    // Try to get existing
    let sql = format!("SELECT origin_id FROM {}.origin WHERE name = $1", schema);
    let rows = client
        .query(&sql, &[&name])
        .await
        .map_err(|e| FerroError::Io {
            msg: format!("Query origin failed: {}", e),
        })?;

    if let Some(row) = rows.first() {
        return Ok(row.get(0));
    }

    // Create new origin
    let sql = format!(
        "INSERT INTO {}.origin (name, descr, url, updated) \
         VALUES ($1, 'Loaded from cdot by ferro-benchmark', 'https://github.com/SACGF/cdot', NOW()) \
         RETURNING origin_id",
        schema
    );
    let row = client
        .query_one(&sql, &[&name])
        .await
        .map_err(|e| FerroError::Io {
            msg: format!("Insert origin failed: {}", e),
        })?;

    Ok(row.get(0))
}

/// Get set of transcripts that already exist in UTA (either in transcript table or with alignments)
fn get_existing_transcripts(container: &str, schema: &str) -> Result<HashSet<String>, FerroError> {
    // Query both transcript table and exon_set to find all existing accessions
    // We need to skip transcripts that exist in the transcript table (PRIMARY KEY constraint)
    // AND transcripts that already have splign alignments in exon_set
    let sql = format!(
        "SELECT ac FROM {}.transcript UNION SELECT DISTINCT tx_ac FROM {}.exon_set WHERE alt_aln_method LIKE 'splign%'",
        schema, schema
    );

    let output = Command::new("docker")
        .args([
            "exec", container, "psql", "-U", "postgres", "-d", "uta", "-t", "-c", &sql,
        ])
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Docker exec failed: {}", e),
        })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(FerroError::Io {
            msg: format!("SQL query failed: {}", stderr),
        });
    }

    Ok(String::from_utf8_lossy(&output.stdout)
        .lines()
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_escape_copy_text() {
        assert_eq!(escape_copy_text("simple"), "simple");
        assert_eq!(escape_copy_text("with\ttab"), "with\\ttab");
        assert_eq!(escape_copy_text("with\nnewline"), "with\\nnewline");
        assert_eq!(escape_copy_text("back\\slash"), "back\\\\slash");
        assert_eq!(escape_copy_text("all\t\n\\three"), "all\\t\\n\\\\three");
    }

    #[test]
    fn test_config_default() {
        let config = UtaLoadConfig::default();
        assert_eq!(config.uta_schema, "uta_20210129b");
        assert_eq!(config.container_name, "ferro-uta");
        assert_eq!(config.origin_name, "cdot-ferro");
    }
}
