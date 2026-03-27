//! In-memory hgvs-rs `Provider` backed by data pre-loaded from UTA PostgreSQL
//! and SeqRepo.
//!
//! The goal is to eliminate all PostgreSQL I/O from benchmark measurements by
//! loading every table/view that the `Provider` trait requires into HashMaps
//! at construction time.

use std::collections::HashMap;
use std::time::Instant;

use indexmap::IndexMap;
use postgres::{Client, NoTls};

use biocommons_bioutils::assemblies::{Assembly, ASSEMBLY_INFOS};
use hgvs::data::error::Error as HgvsDataError;
use hgvs::data::interface::{
    GeneInfoRecord, Provider, TxExonsRecord, TxForRegionRecord, TxIdentityInfo, TxInfoRecord,
    TxMappingOptionsRecord, TxSimilarityRecord,
};
use seqrepo::{AliasOrSeqId, Interface as SeqRepoInterface, Query as SeqRepoQuery, SeqRepo};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration for building an [`InMemoryProvider`].
#[derive(Debug, Clone)]
pub struct InMemoryProviderConfig {
    /// PostgreSQL connection URL for UTA (e.g.
    /// `"postgresql://anonymous:anonymous@localhost:5432/uta"`).
    pub uta_db_url: String,
    /// UTA schema name, doubling as the data-version string
    /// (e.g. `"uta_20210129b"`).
    pub uta_db_schema: String,
    /// Path to a SeqRepo instance.  Not used in this task but reserved for
    /// Task 3 when sequence loading is added.
    pub seqrepo_path: Option<String>,
}

impl Default for InMemoryProviderConfig {
    fn default() -> Self {
        Self {
            uta_db_url: "postgresql://anonymous:anonymous@localhost:5432/uta".to_string(),
            uta_db_schema: "uta_20210129b".to_string(),
            seqrepo_path: None,
        }
    }
}

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

/// An in-memory record for the `tx_for_region` interval index.
///
/// Each entry represents one (tx_ac, alt_ac, alt_strand, alt_aln_method)
/// combination together with its aggregate exon span (min start, max end).
#[derive(Debug, Clone)]
pub struct TxRegionRecord {
    pub tx_ac: String,
    pub alt_ac: String,
    pub alt_strand: i16,
    pub alt_aln_method: String,
    /// Minimum `start_i` across all exons in this alignment set.
    pub start_i: i32,
    /// Maximum `end_i` across all exons in this alignment set.
    pub end_i: i32,
}

/// All UTA metadata held entirely in memory.
///
/// The field types mirror the return types of the corresponding `Provider`
/// trait methods so that the future `impl Provider` can simply clone from
/// the maps.
pub struct InMemoryProvider {
    // -- bookkeeping -------------------------------------------------------
    /// Schema/data version (e.g. `"uta_20210129b"`).
    pub data_version: String,
    /// Schema version read from `{schema}.meta`.
    pub schema_version: String,

    // -- assembly maps (pre-built from biocommons-bioutils) -----------------
    /// `Assembly -> IndexMap<refseq_ac, chromosome_name>`.
    pub assembly_maps: HashMap<Assembly, IndexMap<String, String>>,

    // -- gene info ---------------------------------------------------------
    /// `hgnc -> GeneInfoRecord`.
    pub gene_info: HashMap<String, GeneInfoRecord>,

    // -- transcript-protein associations -----------------------------------
    /// `tx_ac -> pro_ac` (first protein accession, descending order).
    pub tx_to_pro: HashMap<String, String>,

    // -- sequences ---------------------------------------------------------
    /// `accession -> full sequence`.
    pub sequences: HashMap<String, String>,

    // -- seq_anno: accession <-> seq_id ------------------------------------
    /// `seq_id -> Vec<accession>` (for `get_acs_for_protein_seq`).
    pub seq_id_to_acs: HashMap<String, Vec<String>>,

    // -- similar transcripts -----------------------------------------------
    /// `tx_ac1 -> Vec<TxSimilarityRecord>`.
    pub similar_transcripts: HashMap<String, Vec<TxSimilarityRecord>>,

    // -- transcript exon alignments ----------------------------------------
    /// `(tx_ac, alt_ac, alt_aln_method) -> Vec<TxExonsRecord>` ordered by
    /// `alt_start_i`.
    pub tx_exons: HashMap<(String, String, String), Vec<TxExonsRecord>>,

    // -- transcript info by gene -------------------------------------------
    /// `gene (hgnc) -> Vec<TxInfoRecord>`.
    pub tx_for_gene: HashMap<String, Vec<TxInfoRecord>>,

    // -- transcript info by (tx_ac, alt_ac, alt_aln_method) ----------------
    /// `(tx_ac, alt_ac, alt_aln_method) -> TxInfoRecord`.
    pub tx_info: HashMap<(String, String, String), TxInfoRecord>,

    // -- transcript identity -----------------------------------------------
    /// `tx_ac -> TxIdentityInfo`.
    pub tx_identity_info: HashMap<String, TxIdentityInfo>,

    // -- transcript mapping options ----------------------------------------
    /// `tx_ac -> Vec<TxMappingOptionsRecord>`.
    pub tx_mapping_options: HashMap<String, Vec<TxMappingOptionsRecord>>,

    // -- tx_for_region interval data ---------------------------------------
    /// `(alt_ac, alt_aln_method) -> Vec<TxRegionRecord>`.
    ///
    /// Each entry stores the aggregate span for one transcript alignment.
    /// The `get_tx_for_region` query filters these by overlap at lookup time.
    pub tx_region_index: HashMap<(String, String), Vec<TxRegionRecord>>,
}

// ---------------------------------------------------------------------------
// HGNC symbols of selenoproteins (mirrors hgvs-rs constant).
// ---------------------------------------------------------------------------

const SELENOPROTEIN_SYMBOLS: [&str; 25] = [
    "DIO1", "DIO3", "GPX1", "GPX2", "GPX3", "GPX4", "GPX6", "SELENOF", "SELENOH", "SELENOI",
    "SELENOK", "SELENOM", "SELENON", "SELENOO", "SELENOP", "MSRB1", "SELENOS", "SELENOT",
    "SELENOV", "SELENOW", "DIO2", "SEPHS2", "TXNRD1", "TXNRD2", "TXNRD3",
];

// ---------------------------------------------------------------------------
// Type aliases for complex return types
// ---------------------------------------------------------------------------

/// Sequences keyed by accession, plus a reverse index from seq_id to accessions.
type SequenceData = (HashMap<String, String>, HashMap<String, Vec<String>>);

/// Tx exon alignments keyed by (tx_ac, alt_ac, alt_aln_method), plus mapping options by tx_ac.
type TxExonData = (
    HashMap<(String, String, String), Vec<TxExonsRecord>>,
    HashMap<String, Vec<TxMappingOptionsRecord>>,
);

/// Tx info keyed by (tx_ac, alt_ac, alt_aln_method), plus tx info grouped by gene.
type TxInfoData = (
    HashMap<(String, String, String), TxInfoRecord>,
    HashMap<String, Vec<TxInfoRecord>>,
);

// ---------------------------------------------------------------------------
// Construction — load everything from UTA
// ---------------------------------------------------------------------------

impl InMemoryProvider {
    /// Connect to UTA PostgreSQL and load **all** metadata into memory.
    ///
    /// This is intentionally slow (minutes) — it only runs once before
    /// benchmarking begins.
    pub fn new(config: &InMemoryProviderConfig) -> Result<Self, String> {
        let total_start = Instant::now();

        // -- connect -------------------------------------------------------
        eprintln!("Connecting to UTA at {} ...", config.uta_db_url);
        let mut conn = Client::connect(&config.uta_db_url, NoTls)
            .map_err(|e| format!("Failed to connect to UTA: {e}"))?;

        let schema = &config.uta_db_schema;

        // -- schema version ------------------------------------------------
        eprintln!("Loading schema version ...");
        let schema_version = {
            let sql = format!("SELECT value FROM {schema}.meta WHERE key = 'schema_version'");
            let row = conn
                .query_one(&sql, &[])
                .map_err(|e| format!("Failed to query schema_version: {e}"))?;
            row.get::<_, String>("value")
        };
        eprintln!("  schema_version = {schema_version}");

        // -- assembly maps (from biocommons-bioutils, no DB needed) ---------
        eprintln!("Building assembly maps ...");
        let assembly_maps = Self::build_assembly_maps();

        // -- gene info -----------------------------------------------------
        let gene_info = Self::load_gene_info(&mut conn, schema)?;

        // -- associated accessions (tx -> protein) -------------------------
        let tx_to_pro = Self::load_tx_to_pro(&mut conn, schema)?;

        // -- sequences (seq_anno + seq) ------------------------------------
        let (sequences, seq_id_to_acs) = Self::load_sequences(&mut conn, schema)?;

        // -- similar transcripts (tx_similarity_v) -------------------------
        let similar_transcripts = Self::load_similar_transcripts(&mut conn, schema)?;

        // -- tx identity info (tx_def_summary_v) ---------------------------
        let tx_identity_info = Self::load_tx_identity_info(&mut conn, schema)?;

        // -- tx exon alignments (tx_exon_aln_v) ----------------------------
        // Also derives: tx_mapping_options, tx_exons
        let (tx_exons, tx_mapping_options) = Self::load_tx_exon_alignments(&mut conn, schema)?;

        // -- transcript info (transcript + exon_set) -----------------------
        // Derives: tx_info, tx_for_gene
        let (tx_info, tx_for_gene) = Self::load_transcript_info(&mut conn, schema)?;

        // -- tx_for_region index (exon_set + exon) -------------------------
        let tx_region_index = Self::load_tx_region_index(&mut conn, schema)?;

        // -- SeqRepo sequences (optional) ------------------------------------
        let mut sequences = sequences;
        if let Some(ref seqrepo_path) = config.seqrepo_path {
            Self::load_seqrepo_sequences(seqrepo_path, &mut sequences)?;
        }

        let elapsed = total_start.elapsed();
        eprintln!(
            "InMemoryProvider loaded in {:.1}s  (gene_info={}, sequences={}, \
             tx_identity={}, tx_exons_keys={}, tx_info_keys={}, region_index_keys={})",
            elapsed.as_secs_f64(),
            gene_info.len(),
            sequences.len(),
            tx_identity_info.len(),
            tx_exons.len(),
            tx_info.len(),
            tx_region_index.len(),
        );

        Ok(Self {
            data_version: schema.clone(),
            schema_version,
            assembly_maps,
            gene_info,
            tx_to_pro,
            sequences,
            seq_id_to_acs,
            similar_transcripts,
            tx_exons,
            tx_for_gene,
            tx_info,
            tx_identity_info,
            tx_mapping_options,
            tx_region_index,
        })
    }

    // -----------------------------------------------------------------------
    // Private helpers — each loads one logical dataset
    // -----------------------------------------------------------------------

    fn build_assembly_maps() -> HashMap<Assembly, IndexMap<String, String>> {
        let mut maps = HashMap::new();
        for assembly in [Assembly::Grch37, Assembly::Grch37p10, Assembly::Grch38] {
            let map = IndexMap::from_iter(
                ASSEMBLY_INFOS[assembly]
                    .sequences
                    .iter()
                    .map(|record| (record.refseq_ac.clone(), record.name.clone())),
            );
            maps.insert(assembly, map);
        }
        maps
    }

    /// Load `{schema}.gene` into `hgnc -> GeneInfoRecord`.
    fn load_gene_info(
        conn: &mut Client,
        schema: &str,
    ) -> Result<HashMap<String, GeneInfoRecord>, String> {
        let t = Instant::now();
        eprintln!("Loading gene info ...");

        let sql = format!("SELECT hgnc, maploc, descr, summary, aliases, added FROM {schema}.gene");
        let rows = conn
            .query(&sql, &[])
            .map_err(|e| format!("Failed to query gene table: {e}"))?;

        let mut map = HashMap::with_capacity(rows.len());
        for row in &rows {
            let hgnc: String = row.get("hgnc");
            let aliases_str: String = row.get("aliases");
            let aliases: Vec<String> = aliases_str.split(',').map(|s| s.to_owned()).collect();
            let record = GeneInfoRecord {
                hgnc: hgnc.clone(),
                maploc: row.get("maploc"),
                descr: row.get("descr"),
                summary: row.get("summary"),
                aliases,
                added: row.get("added"),
            };
            map.insert(hgnc, record);
        }

        eprintln!(
            "  gene info: {} records in {:.1}s",
            map.len(),
            t.elapsed().as_secs_f64()
        );
        Ok(map)
    }

    /// Load `{schema}.associated_accessions` into `tx_ac -> pro_ac`.
    ///
    /// Matches the UTA query: `SELECT pro_ac ... WHERE tx_ac = $1 ORDER BY pro_ac DESC`
    /// — we keep only the first (highest) pro_ac per tx_ac.
    fn load_tx_to_pro(conn: &mut Client, schema: &str) -> Result<HashMap<String, String>, String> {
        let t = Instant::now();
        eprintln!("Loading associated accessions ...");

        let sql = format!(
            "SELECT tx_ac, pro_ac FROM {schema}.associated_accessions \
             ORDER BY tx_ac, pro_ac DESC"
        );
        let rows = conn
            .query(&sql, &[])
            .map_err(|e| format!("Failed to query associated_accessions: {e}"))?;

        let mut map: HashMap<String, String> = HashMap::new();
        for row in &rows {
            let tx_ac: String = row.get("tx_ac");
            let pro_ac: String = row.get("pro_ac");
            // Keep only the first (highest by DESC order) pro_ac per tx_ac.
            map.entry(tx_ac).or_insert(pro_ac);
        }

        eprintln!(
            "  associated accessions: {} mappings in {:.1}s",
            map.len(),
            t.elapsed().as_secs_f64()
        );
        Ok(map)
    }

    /// Load sequences via `{schema}.seq_anno` + `{schema}.seq`.
    ///
    /// Returns:
    /// - `accession -> full_sequence` (for `get_seq` / `get_seq_part`)
    /// - `seq_id -> Vec<accession>` (for `get_acs_for_protein_seq`)
    fn load_sequences(conn: &mut Client, schema: &str) -> Result<SequenceData, String> {
        let t = Instant::now();
        eprintln!("Loading seq_anno -> seq_id mapping ...");

        // Step 1: load seq_anno (ac -> seq_id) and build reverse index.
        let sql = format!("SELECT ac, seq_id FROM {schema}.seq_anno");
        let rows = conn
            .query(&sql, &[])
            .map_err(|e| format!("Failed to query seq_anno: {e}"))?;

        let mut ac_to_seq_id: HashMap<String, String> = HashMap::with_capacity(rows.len());
        let mut seq_id_to_acs: HashMap<String, Vec<String>> = HashMap::new();
        for row in &rows {
            let ac: String = row.get("ac");
            let seq_id: String = row.get("seq_id");
            seq_id_to_acs
                .entry(seq_id.clone())
                .or_default()
                .push(ac.clone());
            ac_to_seq_id.insert(ac, seq_id);
        }
        eprintln!(
            "  seq_anno: {} accessions, {} unique seq_ids in {:.1}s",
            ac_to_seq_id.len(),
            seq_id_to_acs.len(),
            t.elapsed().as_secs_f64(),
        );

        // Step 2: load seq (seq_id -> sequence).
        let t2 = Instant::now();
        eprintln!("Loading sequences ...");
        let sql = format!("SELECT seq_id, seq FROM {schema}.seq");
        let rows = conn
            .query(&sql, &[])
            .map_err(|e| format!("Failed to query seq: {e}"))?;

        let mut seq_id_to_seq: HashMap<String, String> = HashMap::with_capacity(rows.len());
        for row in &rows {
            let seq_id: String = row.get("seq_id");
            let seq: Option<String> = row.get("seq");
            if let Some(seq) = seq {
                seq_id_to_seq.insert(seq_id, seq);
            }
        }
        eprintln!(
            "  seq: {} sequences in {:.1}s",
            seq_id_to_seq.len(),
            t2.elapsed().as_secs_f64(),
        );

        // Step 3: build ac -> sequence map.
        let t3 = Instant::now();
        eprintln!("Building accession -> sequence map ...");
        let mut sequences: HashMap<String, String> = HashMap::with_capacity(ac_to_seq_id.len());
        let mut missing = 0usize;
        for (ac, seq_id) in &ac_to_seq_id {
            if let Some(seq) = seq_id_to_seq.get(seq_id) {
                sequences.insert(ac.clone(), seq.clone());
            } else {
                missing += 1;
            }
        }
        if missing > 0 {
            eprintln!("  WARNING: {missing} accessions had no matching sequence");
        }
        eprintln!(
            "  accession -> sequence: {} entries in {:.1}s",
            sequences.len(),
            t3.elapsed().as_secs_f64(),
        );

        Ok((sequences, seq_id_to_acs))
    }

    /// Load `{schema}.tx_similarity_v` into `tx_ac1 -> Vec<TxSimilarityRecord>`.
    fn load_similar_transcripts(
        conn: &mut Client,
        schema: &str,
    ) -> Result<HashMap<String, Vec<TxSimilarityRecord>>, String> {
        let t = Instant::now();
        eprintln!("Loading similar transcripts (tx_similarity_v) ...");

        let sql = format!(
            "SELECT tx_ac1, tx_ac2, hgnc_eq, cds_eq, es_fp_eq, \
                    cds_es_fp_eq, cds_exon_lengths_fp_eq \
             FROM {schema}.tx_similarity_v \
             ORDER BY tx_ac1, tx_ac2"
        );
        let rows = conn
            .query(&sql, &[])
            .map_err(|e| format!("Failed to query tx_similarity_v: {e}"))?;

        let mut map: HashMap<String, Vec<TxSimilarityRecord>> = HashMap::new();
        for row in &rows {
            let record = TxSimilarityRecord {
                tx_ac1: row.get("tx_ac1"),
                tx_ac2: row.get("tx_ac2"),
                hgnc_eq: row.try_get("hgnc_eq").unwrap_or(false),
                cds_eq: row.try_get("cds_eq").unwrap_or(false),
                es_fp_eq: row.try_get("es_fp_eq").unwrap_or(false),
                cds_es_fp_eq: row.try_get("cds_es_fp_eq").unwrap_or(false),
                cds_exon_lengths_fp_eq: row.try_get("cds_exon_lengths_fp_eq").unwrap_or(false),
            };
            let key: String = row.get("tx_ac1");
            map.entry(key).or_default().push(record);
        }

        eprintln!(
            "  similar transcripts: {} keys, {} total records in {:.1}s",
            map.len(),
            map.values().map(|v| v.len()).sum::<usize>(),
            t.elapsed().as_secs_f64(),
        );
        Ok(map)
    }

    /// Load `{schema}.tx_def_summary_v` into `tx_ac -> TxIdentityInfo`.
    ///
    /// The UTA query is:
    /// ```sql
    /// SELECT DISTINCT(tx_ac), alt_ac, alt_aln_method, cds_start_i,
    ///        cds_end_i, lengths, hgnc
    /// FROM {schema}.tx_def_summary_v
    /// ORDER BY tx_ac, alt_ac, ...
    /// ```
    /// We keep only the first row per tx_ac (matching `query_one` semantics).
    fn load_tx_identity_info(
        conn: &mut Client,
        schema: &str,
    ) -> Result<HashMap<String, TxIdentityInfo>, String> {
        let t = Instant::now();
        eprintln!("Loading tx identity info (tx_def_summary_v) — this may be slow ...");

        let sql = format!(
            "SELECT DISTINCT tx_ac, alt_ac, alt_aln_method, cds_start_i, \
                    cds_end_i, lengths, hgnc \
             FROM {schema}.tx_def_summary_v \
             ORDER BY tx_ac, alt_ac, alt_aln_method, cds_start_i, cds_end_i, lengths, hgnc"
        );
        let rows = conn
            .query(&sql, &[])
            .map_err(|e| format!("Failed to query tx_def_summary_v: {e}"))?;

        let mut map: HashMap<String, TxIdentityInfo> = HashMap::new();
        for row in &rows {
            let tx_ac: String = row.get("tx_ac");
            // Keep only the first row per tx_ac (matches query_one behavior).
            if map.contains_key(&tx_ac) {
                continue;
            }
            let hgnc: &str = row.get("hgnc");
            let is_selenoprotein = SELENOPROTEIN_SYMBOLS.contains(&hgnc);
            let record = TxIdentityInfo {
                tx_ac: tx_ac.clone(),
                alt_ac: row.get("alt_ac"),
                alt_aln_method: row.get("alt_aln_method"),
                cds_start_i: row.get("cds_start_i"),
                cds_end_i: row.get("cds_end_i"),
                lengths: row.get("lengths"),
                hgnc: hgnc.to_string(),
                translation_table: if is_selenoprotein {
                    hgvs::sequences::TranslationTable::Selenocysteine
                } else {
                    Default::default()
                },
            };
            map.insert(tx_ac, record);
        }

        eprintln!(
            "  tx identity info: {} records in {:.1}s",
            map.len(),
            t.elapsed().as_secs_f64(),
        );
        Ok(map)
    }

    /// Load `{schema}.tx_exon_aln_v` to build:
    /// - `(tx_ac, alt_ac, alt_aln_method) -> Vec<TxExonsRecord>` (ordered by `alt_start_i`)
    /// - `tx_ac -> Vec<TxMappingOptionsRecord>` (distinct triples where `exon_aln_id IS NOT NULL`)
    fn load_tx_exon_alignments(conn: &mut Client, schema: &str) -> Result<TxExonData, String> {
        let t = Instant::now();
        eprintln!("Loading tx exon alignments (tx_exon_aln_v) — this may be slow ...");

        let sql = format!(
            "SELECT hgnc, tx_ac, alt_ac, alt_aln_method, alt_strand, ord, \
                    tx_start_i, tx_end_i, alt_start_i, alt_end_i, cigar, \
                    tx_aseq, alt_aseq, tx_exon_set_id, alt_exon_set_id, \
                    tx_exon_id, alt_exon_id, exon_aln_id \
             FROM {schema}.tx_exon_aln_v \
             ORDER BY tx_ac, alt_ac, alt_aln_method, alt_start_i"
        );
        let rows = conn
            .query(&sql, &[])
            .map_err(|e| format!("Failed to query tx_exon_aln_v: {e}"))?;

        let mut exons_map: HashMap<(String, String, String), Vec<TxExonsRecord>> = HashMap::new();
        // Track distinct (tx_ac, alt_ac, alt_aln_method) triples that have
        // a non-null exon_aln_id (for mapping options).
        let mut mapping_triples: HashMap<String, Vec<(String, String, String)>> = HashMap::new();
        // Use a set to deduplicate mapping triples.
        let mut mapping_seen: std::collections::HashSet<(String, String, String)> =
            std::collections::HashSet::new();

        for row in &rows {
            let tx_ac: String = row.get("tx_ac");
            let alt_ac: String = row.get("alt_ac");
            let alt_aln_method: String = row.get("alt_aln_method");
            let exon_aln_id: Option<i32> = row.get("exon_aln_id");

            let record = TxExonsRecord {
                hgnc: row.get("hgnc"),
                tx_ac: tx_ac.clone(),
                alt_ac: alt_ac.clone(),
                alt_aln_method: alt_aln_method.clone(),
                alt_strand: row.get("alt_strand"),
                ord: row.get("ord"),
                tx_start_i: row.get("tx_start_i"),
                tx_end_i: row.get("tx_end_i"),
                alt_start_i: row.get("alt_start_i"),
                alt_end_i: row.get("alt_end_i"),
                cigar: row.get("cigar"),
                tx_aseq: row.get("tx_aseq"),
                alt_aseq: row.get("alt_aseq"),
                tx_exon_set_id: row.get("tx_exon_set_id"),
                alt_exon_set_id: row.get("alt_exon_set_id"),
                tx_exon_id: row.get("tx_exon_id"),
                alt_exon_id: row.get("alt_exon_id"),
                exon_aln_id: exon_aln_id.unwrap_or(0),
            };

            let key = (tx_ac.clone(), alt_ac.clone(), alt_aln_method.clone());
            exons_map.entry(key.clone()).or_default().push(record);

            // Track mapping options (where exon_aln_id IS NOT NULL).
            if exon_aln_id.is_some() && mapping_seen.insert(key) {
                mapping_triples.entry(tx_ac).or_default().push((
                    row.get("tx_ac"),
                    alt_ac,
                    alt_aln_method,
                ));
            }
        }

        // Convert mapping triples to TxMappingOptionsRecord.
        let mut tx_mapping_options: HashMap<String, Vec<TxMappingOptionsRecord>> = HashMap::new();
        for (tx_ac, triples) in mapping_triples {
            let mut records: Vec<TxMappingOptionsRecord> = triples
                .into_iter()
                .map(|(tx, alt, method)| TxMappingOptionsRecord {
                    tx_ac: tx,
                    alt_ac: alt,
                    alt_aln_method: method,
                })
                .collect();
            records.sort_by(|a, b| {
                a.tx_ac
                    .cmp(&b.tx_ac)
                    .then(a.alt_ac.cmp(&b.alt_ac))
                    .then(a.alt_aln_method.cmp(&b.alt_aln_method))
            });
            tx_mapping_options.insert(tx_ac, records);
        }

        eprintln!(
            "  tx exon alignments: {} alignment keys, {} mapping-option keys in {:.1}s",
            exons_map.len(),
            tx_mapping_options.len(),
            t.elapsed().as_secs_f64(),
        );
        Ok((exons_map, tx_mapping_options))
    }

    /// Load transcript info from `{schema}.transcript` joined with
    /// `{schema}.exon_set` to build:
    /// - `(tx_ac, alt_ac, alt_aln_method) -> TxInfoRecord`
    /// - `gene (hgnc) -> Vec<TxInfoRecord>` (for `get_tx_for_gene`, excludes
    ///   `alt_aln_method = 'transcript'`)
    fn load_transcript_info(conn: &mut Client, schema: &str) -> Result<TxInfoData, String> {
        let t = Instant::now();
        eprintln!("Loading transcript info (transcript + exon_set) ...");

        let sql = format!(
            "SELECT hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method \
             FROM {schema}.transcript T \
             JOIN {schema}.exon_set ES ON T.ac = ES.tx_ac \
             ORDER BY hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method"
        );
        let rows = conn
            .query(&sql, &[])
            .map_err(|e| format!("Failed to query transcript+exon_set: {e}"))?;

        let mut info_map: HashMap<(String, String, String), TxInfoRecord> = HashMap::new();
        let mut gene_map: HashMap<String, Vec<TxInfoRecord>> = HashMap::new();

        for row in &rows {
            let hgnc: String = row.get("hgnc");
            let tx_ac: String = row.get("tx_ac");
            let alt_ac: String = row.get("alt_ac");
            let alt_aln_method: String = row.get("alt_aln_method");

            let record = TxInfoRecord {
                hgnc: hgnc.clone(),
                cds_start_i: row.get("cds_start_i"),
                cds_end_i: row.get("cds_end_i"),
                tx_ac: tx_ac.clone(),
                alt_ac: alt_ac.clone(),
                alt_aln_method: alt_aln_method.clone(),
            };

            let key = (tx_ac, alt_ac, alt_aln_method.clone());
            info_map.insert(key, record.clone());

            // get_tx_for_gene excludes alt_aln_method = 'transcript'
            if alt_aln_method != "transcript" {
                gene_map.entry(hgnc).or_default().push(record);
            }
        }

        eprintln!(
            "  transcript info: {} info keys, {} gene keys in {:.1}s",
            info_map.len(),
            gene_map.len(),
            t.elapsed().as_secs_f64(),
        );
        Ok((info_map, gene_map))
    }

    /// Build the `tx_for_region` interval index from `{schema}.exon_set`
    /// joined with `{schema}.exon`.
    ///
    /// The UTA query groups by (tx_ac, alt_ac, alt_strand, alt_aln_method)
    /// and computes MIN(start_i) and MAX(end_i).  We store these aggregates
    /// keyed by (alt_ac, alt_aln_method) for fast region-overlap filtering.
    fn load_tx_region_index(
        conn: &mut Client,
        schema: &str,
    ) -> Result<HashMap<(String, String), Vec<TxRegionRecord>>, String> {
        let t = Instant::now();
        eprintln!("Loading tx region index (exon_set + exon) ...");

        let sql = format!(
            "SELECT tx_ac, alt_ac, alt_strand, alt_aln_method, \
                    MIN(start_i) AS start_i, MAX(end_i) AS end_i \
             FROM {schema}.exon_set es \
             JOIN {schema}.exon e ON es.exon_set_id = e.exon_set_id \
             GROUP BY tx_ac, alt_ac, alt_strand, alt_aln_method \
             ORDER BY tx_ac, alt_ac, alt_strand, alt_aln_method"
        );
        let rows = conn
            .query(&sql, &[])
            .map_err(|e| format!("Failed to query exon_set+exon for region index: {e}"))?;

        let mut map: HashMap<(String, String), Vec<TxRegionRecord>> = HashMap::new();
        for row in &rows {
            let alt_ac: String = row.get("alt_ac");
            let alt_aln_method: String = row.get("alt_aln_method");
            let record = TxRegionRecord {
                tx_ac: row.get("tx_ac"),
                alt_ac: alt_ac.clone(),
                alt_strand: row.get("alt_strand"),
                alt_aln_method: alt_aln_method.clone(),
                start_i: row.get("start_i"),
                end_i: row.get("end_i"),
            };
            map.entry((alt_ac, alt_aln_method))
                .or_default()
                .push(record);
        }

        eprintln!(
            "  tx region index: {} keys, {} total records in {:.1}s",
            map.len(),
            map.values().map(|v| v.len()).sum::<usize>(),
            t.elapsed().as_secs_f64(),
        );
        Ok(map)
    }

    /// Load sequences from a SeqRepo instance into the existing `sequences` map.
    ///
    /// The `seqrepo_path` should be the full path to the SeqRepo instance
    /// directory (e.g. `/path/to/seqrepo/2021-01-29`).  The last component is
    /// used as the instance name and the parent as the root directory.
    fn load_seqrepo_sequences(
        seqrepo_path: &str,
        sequences: &mut HashMap<String, String>,
    ) -> Result<(), String> {
        let t = Instant::now();
        eprintln!("Loading sequences from SeqRepo at {seqrepo_path} ...");

        let path = std::path::Path::new(seqrepo_path);
        let instance = path
            .file_name()
            .ok_or_else(|| format!("Could not extract instance name from path: {seqrepo_path}"))?
            .to_str()
            .ok_or_else(|| format!("Non-UTF8 path component in: {seqrepo_path}"))?;
        let parent = path
            .parent()
            .ok_or_else(|| format!("Could not extract parent from path: {seqrepo_path}"))?;

        let repo = SeqRepo::new(parent, instance)
            .map_err(|e| format!("Failed to open SeqRepo at {seqrepo_path}: {e}"))?;

        // Step 1: Collect all aliases, grouped by seq_id.
        let query = SeqRepoQuery::default();
        let mut seqid_to_aliases: HashMap<String, Vec<String>> = HashMap::new();
        repo.alias_db()
            .find(&query, |record_result| {
                if let Ok(record) = record_result {
                    seqid_to_aliases
                        .entry(record.seqid.clone())
                        .or_default()
                        .push(record.alias);
                }
            })
            .map_err(|e| format!("Failed to query SeqRepo aliases: {e}"))?;

        eprintln!(
            "  SeqRepo aliases: {} unique seq_ids, loading sequences ...",
            seqid_to_aliases.len(),
        );

        // Step 2: For each unique seq_id, check if any alias is already loaded.
        //         If not, fetch the sequence and store under all aliases.
        let mut loaded = 0usize;
        let mut skipped = 0usize;
        for (seqid, aliases) in &seqid_to_aliases {
            // Skip if all aliases for this seq_id are already present.
            if aliases.iter().all(|a| sequences.contains_key(a)) {
                skipped += aliases.len();
                continue;
            }

            // Fetch the full sequence by seq_id.
            let seq =
                match repo.fetch_sequence_part(&AliasOrSeqId::SeqId(seqid.clone()), None, None) {
                    Ok(s) => s,
                    Err(e) => {
                        eprintln!("  WARNING: failed to fetch sequence for seq_id={seqid}: {e}");
                        continue;
                    }
                };

            // Store under each alias that is not already present.
            for alias in aliases {
                if !sequences.contains_key(alias) {
                    sequences.insert(alias.clone(), seq.clone());
                    loaded += 1;
                }
            }
        }

        eprintln!(
            "  SeqRepo: loaded {} new alias->sequence mappings, skipped {} already present in {:.1}s",
            loaded,
            skipped,
            t.elapsed().as_secs_f64(),
        );
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Provider trait implementation
// ---------------------------------------------------------------------------

impl Provider for InMemoryProvider {
    fn data_version(&self) -> &str {
        &self.data_version
    }

    fn schema_version(&self) -> &str {
        &self.schema_version
    }

    fn get_assembly_map(&self, assembly: Assembly) -> IndexMap<String, String> {
        self.assembly_maps
            .get(&assembly)
            .cloned()
            .unwrap_or_default()
    }

    fn get_gene_info(&self, hgnc: &str) -> Result<GeneInfoRecord, HgvsDataError> {
        self.gene_info
            .get(hgnc)
            .cloned()
            .ok_or_else(|| HgvsDataError::NoGeneFound(hgnc.to_string()))
    }

    fn get_pro_ac_for_tx_ac(&self, tx_ac: &str) -> Result<Option<String>, HgvsDataError> {
        Ok(self.tx_to_pro.get(tx_ac).cloned())
    }

    fn get_seq_part(
        &self,
        ac: &str,
        begin: Option<usize>,
        end: Option<usize>,
    ) -> Result<String, HgvsDataError> {
        let seq = self
            .sequences
            .get(ac)
            .ok_or_else(|| HgvsDataError::NoSequenceRecord(ac.to_string()))?;
        let start = begin.unwrap_or(0);
        let stop = end.unwrap_or(seq.len());
        Ok(seq[start..stop].to_string())
    }

    fn get_acs_for_protein_seq(&self, seq: &str) -> Result<Vec<String>, HgvsDataError> {
        let md5 = hgvs::sequences::seq_md5(seq, true)?;
        let mut acs = self.seq_id_to_acs.get(&md5).cloned().unwrap_or_default();
        acs.push(format!("MD5_{md5}"));
        Ok(acs)
    }

    fn get_similar_transcripts(
        &self,
        tx_ac: &str,
    ) -> Result<Vec<TxSimilarityRecord>, HgvsDataError> {
        Ok(self
            .similar_transcripts
            .get(tx_ac)
            .cloned()
            .unwrap_or_default())
    }

    fn get_tx_exons(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<Vec<TxExonsRecord>, HgvsDataError> {
        let key = (
            tx_ac.to_string(),
            alt_ac.to_string(),
            alt_aln_method.to_string(),
        );
        self.tx_exons.get(&key).cloned().ok_or_else(|| {
            HgvsDataError::NoTxExons(
                tx_ac.to_string(),
                alt_ac.to_string(),
                alt_aln_method.to_string(),
            )
        })
    }

    fn get_tx_for_gene(&self, gene: &str) -> Result<Vec<TxInfoRecord>, HgvsDataError> {
        Ok(self.tx_for_gene.get(gene).cloned().unwrap_or_default())
    }

    fn get_tx_for_region(
        &self,
        alt_ac: &str,
        alt_aln_method: &str,
        start_i: i32,
        end_i: i32,
    ) -> Result<Vec<TxForRegionRecord>, HgvsDataError> {
        let key = (alt_ac.to_string(), alt_aln_method.to_string());
        let records = self.tx_region_index.get(&key);
        let result = match records {
            Some(entries) => entries
                .iter()
                .filter(|r| r.start_i < end_i && start_i <= r.end_i)
                .map(|r| TxForRegionRecord {
                    tx_ac: r.tx_ac.clone(),
                    alt_ac: r.alt_ac.clone(),
                    alt_strand: r.alt_strand,
                    alt_aln_method: r.alt_aln_method.clone(),
                    start_i: r.start_i,
                    end_i: r.end_i,
                })
                .collect(),
            None => Vec::new(),
        };
        Ok(result)
    }

    fn get_tx_identity_info(&self, tx_ac: &str) -> Result<TxIdentityInfo, HgvsDataError> {
        self.tx_identity_info
            .get(tx_ac)
            .cloned()
            .ok_or_else(|| HgvsDataError::NoTranscriptFound(tx_ac.to_string()))
    }

    fn get_tx_info(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<TxInfoRecord, HgvsDataError> {
        let key = (
            tx_ac.to_string(),
            alt_ac.to_string(),
            alt_aln_method.to_string(),
        );
        self.tx_info
            .get(&key)
            .cloned()
            .ok_or_else(|| HgvsDataError::NoAlignmentFound(tx_ac.to_string(), alt_ac.to_string()))
    }

    fn get_tx_mapping_options(
        &self,
        tx_ac: &str,
    ) -> Result<Vec<TxMappingOptionsRecord>, HgvsDataError> {
        Ok(self
            .tx_mapping_options
            .get(tx_ac)
            .cloned()
            .unwrap_or_default())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Integration test that requires a running UTA PostgreSQL instance.
    ///
    /// Run with:
    /// ```bash
    /// cargo test --features hgvs-rs -- --ignored inmemory_provider
    /// ```
    #[test]
    #[ignore]
    fn test_load_from_uta() {
        let config = InMemoryProviderConfig::default();
        let provider = InMemoryProvider::new(&config)
            .expect("Failed to create InMemoryProvider — is UTA running?");

        // Basic sanity checks.
        assert!(
            !provider.schema_version.is_empty(),
            "schema_version should not be empty"
        );
        assert!(
            !provider.gene_info.is_empty(),
            "gene_info should not be empty"
        );
        assert!(
            !provider.sequences.is_empty(),
            "sequences should not be empty"
        );
        assert!(
            !provider.tx_identity_info.is_empty(),
            "tx_identity_info should not be empty"
        );
        assert!(
            !provider.tx_exons.is_empty(),
            "tx_exons should not be empty"
        );
        assert!(
            !provider.tx_mapping_options.is_empty(),
            "tx_mapping_options should not be empty"
        );
        assert!(!provider.tx_info.is_empty(), "tx_info should not be empty");
        assert!(
            !provider.tx_for_gene.is_empty(),
            "tx_for_gene should not be empty"
        );
        assert!(
            !provider.tx_region_index.is_empty(),
            "tx_region_index should not be empty"
        );
        assert_eq!(
            provider.assembly_maps.len(),
            3,
            "should have 3 assembly maps"
        );

        eprintln!("All InMemoryProvider sanity checks passed.");
    }

    /// Helper to build a provider with both UTA and SeqRepo loaded.
    fn make_provider_with_seqrepo() -> InMemoryProvider {
        let config = InMemoryProviderConfig {
            seqrepo_path: Some("manuscript/benchmark/output/data/seqrepo/2021-01-29".to_string()),
            ..Default::default()
        };
        InMemoryProvider::new(&config)
            .expect("Failed to create InMemoryProvider — is UTA + SeqRepo available?")
    }

    /// Verify that SeqRepo loading adds genomic sequences not present in UTA.
    #[test]
    #[ignore]
    fn test_seqrepo_loads_genomic_sequences() {
        let provider = make_provider_with_seqrepo();

        // A transcript sequence should be available from UTA.
        assert!(
            provider.sequences.contains_key("NM_000051.3"),
            "NM_000051.3 (transcript) should be loaded from UTA"
        );

        // A genomic contig should be available from SeqRepo (not stored in UTA seq).
        assert!(
            provider.sequences.contains_key("NC_000011.10"),
            "NC_000011.10 (genomic) should be loaded from SeqRepo"
        );
    }

    /// Verify `get_seq_part` returns correct-length results.
    #[test]
    #[ignore]
    fn test_get_seq_part() {
        let provider = make_provider_with_seqrepo();

        // Full sequence for a transcript.
        let full = provider.get_seq("NM_000051.3").expect("get_seq failed");
        assert!(
            full.len() > 1000,
            "NM_000051.3 sequence should be > 1000 bp"
        );

        // Partial sequence.
        let part = provider
            .get_seq_part("NM_000051.3", Some(100), Some(200))
            .expect("get_seq_part failed");
        assert_eq!(part.len(), 100, "slice [100..200] should be 100 bp");
        assert_eq!(&full[100..200], part, "slice content should match");
    }

    /// Verify `get_tx_mapping_options` returns non-empty for NM_000051.3.
    #[test]
    #[ignore]
    fn test_get_tx_mapping_options() {
        let provider = make_provider_with_seqrepo();

        let options = provider
            .get_tx_mapping_options("NM_000051.3")
            .expect("get_tx_mapping_options failed");
        assert!(
            !options.is_empty(),
            "NM_000051.3 should have mapping options"
        );
    }

    /// Verify `get_tx_identity_info` returns correct HGNC for NM_000051.3.
    #[test]
    #[ignore]
    fn test_get_tx_identity_info() {
        let provider = make_provider_with_seqrepo();

        let info = provider
            .get_tx_identity_info("NM_000051.3")
            .expect("get_tx_identity_info failed");
        assert_eq!(info.hgnc, "ATM", "NM_000051.3 should map to ATM gene");
    }

    /// Verify `get_tx_exons` returns non-empty for a known transcript+contig pair.
    #[test]
    #[ignore]
    fn test_get_tx_exons() {
        let provider = make_provider_with_seqrepo();

        // First find a valid alt_ac for NM_000051.3 from mapping options.
        let options = provider
            .get_tx_mapping_options("NM_000051.3")
            .expect("get_tx_mapping_options failed");
        let opt = options
            .iter()
            .find(|o| o.alt_aln_method == "splign")
            .expect("should have a splign alignment for NM_000051.3");

        let exons = provider
            .get_tx_exons("NM_000051.3", &opt.alt_ac, "splign")
            .expect("get_tx_exons failed");
        assert!(
            !exons.is_empty(),
            "NM_000051.3 should have exon records for {}",
            opt.alt_ac,
        );
    }
}
