# In-Memory hgvs-rs Provider for Benchmarking

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an in-memory `hgvs::data::interface::Provider` implementation that eliminates PostgreSQL and SeqRepo I/O from hgvs-rs benchmark measurements.

**Architecture:** Load all UTA metadata from PostgreSQL into `HashMap`s at startup. Mmap SeqRepo FASTA sequence files for zero-copy genomic sequence access. Resolve SeqRepo aliases (accession -> seq_id) from SQLite into a `HashMap` at startup. Implement the hgvs-rs `Provider` trait backed by these in-memory structures. The provider is created once and shared across workers via `Arc`.

**Tech Stack:** hgvs-rs `data::interface::Provider` trait, `memmap2` for FASTA mmap, `rusqlite` for one-time alias loading, `postgres` for one-time UTA loading, `seqrepo` crate's `FastaDir` (already handles FASTA indexing).

---

## File Structure

| File | Responsibility |
|------|---------------|
| `src/benchmark/inmemory_provider.rs` (create) | `InMemoryProvider` struct implementing `hgvs::data::interface::Provider`, loading logic, all in-memory data structures |
| `src/benchmark/mod.rs` (modify) | Add `pub mod inmemory_provider;` under `hgvs-rs` feature gate |
| `src/benchmark/hgvs_rs.rs` (modify) | Add `--in-memory` flag support, construct `InMemoryProvider` instead of `uta_sr::Provider` |
| `src/bin/benchmark.rs` (modify) | Wire `--in-memory` CLI flag through to hgvs-rs normalize commands |

---

### Task 1: In-Memory Provider Data Structures and Construction

**Files:**
- Create: `src/benchmark/inmemory_provider.rs`
- Modify: `src/benchmark/mod.rs`

- [ ] **Step 1: Create the module file with struct definitions**

```rust
// src/benchmark/inmemory_provider.rs

//! In-memory hgvs-rs Provider for benchmarking.
//!
//! Loads all UTA metadata and SeqRepo sequences into memory at startup,
//! eliminating database I/O from benchmark measurements.

use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use hgvs::data::interface::{
    GeneInfoRecord, Provider, TxExonsRecord, TxForRegionRecord, TxIdentityInfo, TxInfoRecord,
    TxMappingOptionsRecord, TxSimilarityRecord,
};
use hgvs::data::{error::Error as HgvsError, uta_sr};
use indexmap::IndexMap;

use crate::FerroError;

/// An in-memory Provider that pre-loads all data from UTA + SeqRepo.
///
/// Created once at startup (expensive), then shared across threads via Arc.
/// All lookups are HashMap-based with no I/O.
pub struct InMemoryProvider {
    // Metadata
    data_version: String,
    schema_version: String,
    assembly_maps: HashMap<String, IndexMap<String, String>>,

    // Gene info: hgnc -> GeneInfoRecord
    gene_info: HashMap<String, GeneInfoRecord>,

    // Transcript -> protein accession: tx_ac -> pro_ac
    pro_ac_for_tx_ac: HashMap<String, Option<String>>,

    // Protein seq -> accessions: md5 -> Vec<ac>
    acs_for_protein_seq: HashMap<String, Vec<String>>,

    // Similar transcripts: tx_ac -> Vec<TxSimilarityRecord>
    similar_transcripts: HashMap<String, Vec<TxSimilarityRecord>>,

    // Transcript exons: (tx_ac, alt_ac, alt_aln_method) -> Vec<TxExonsRecord>
    tx_exons: HashMap<(String, String, String), Vec<TxExonsRecord>>,

    // Transcripts for gene: hgnc -> Vec<TxInfoRecord>
    tx_for_gene: HashMap<String, Vec<TxInfoRecord>>,

    // Transcripts for region: (alt_ac, alt_aln_method) -> Vec<(start_i, end_i, TxForRegionRecord)>
    // We store all records per (alt_ac, alt_aln_method) and filter by range at query time.
    tx_for_region: HashMap<(String, String), Vec<TxForRegionRecord>>,

    // Transcript identity info: tx_ac -> TxIdentityInfo
    tx_identity_info: HashMap<String, TxIdentityInfo>,

    // Transcript info: (tx_ac, alt_ac, alt_aln_method) -> TxInfoRecord
    tx_info: HashMap<(String, String, String), TxInfoRecord>,

    // Transcript mapping options: tx_ac -> Vec<TxMappingOptionsRecord>
    tx_mapping_options: HashMap<String, Vec<TxMappingOptionsRecord>>,

    // Sequences: accession -> sequence string
    // For transcript sequences from UTA and genomic sequences from SeqRepo.
    sequences: HashMap<String, String>,
}
```

- [ ] **Step 2: Register the module in mod.rs**

Add to `src/benchmark/mod.rs` under the existing `hgvs-rs` feature gate (find the existing `#[cfg(feature = "hgvs-rs")]` block):

```rust
#[cfg(feature = "hgvs-rs")]
pub mod inmemory_provider;
```

- [ ] **Step 3: Verify it compiles**

Run: `cargo check --features hgvs-rs`
Expected: PASS (struct is defined but not yet used)

- [ ] **Step 4: Commit**

```bash
git add src/benchmark/inmemory_provider.rs src/benchmark/mod.rs
git commit -m "feat(benchmark): add InMemoryProvider struct for hgvs-rs"
```

---

### Task 2: Load UTA Metadata into Memory

**Files:**
- Modify: `src/benchmark/inmemory_provider.rs`

This task adds the constructor that loads all UTA metadata from PostgreSQL into the HashMap fields. This is a one-time cost at startup.

- [ ] **Step 1: Write the test**

Add at the bottom of `src/benchmark/inmemory_provider.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    /// This test requires a running UTA database (ferro-uta docker container).
    /// Run with: cargo nextest run --features hgvs-rs -E 'test(inmemory_provider)'
    #[test]
    #[ignore] // Requires running UTA + SeqRepo
    fn test_load_from_uta() {
        let config = super::InMemoryProviderConfig {
            uta_db_url: "postgresql://anonymous:anonymous@localhost:5432/uta".to_string(),
            uta_db_schema: "uta_20210129b".to_string(),
            seqrepo_path: None, // UTA-only for this test
        };
        let provider = InMemoryProvider::new(&config).unwrap();

        // Verify we loaded data
        assert!(!provider.tx_identity_info.is_empty(), "should have transcript identity info");
        assert!(!provider.tx_mapping_options.is_empty(), "should have mapping options");
        assert!(!provider.sequences.is_empty(), "should have sequences from UTA");

        // Verify a known transcript
        let info = provider.get_tx_identity_info("NM_000051.3").unwrap();
        assert_eq!(info.hgnc, "ATM");
    }
}
```

- [ ] **Step 2: Add the config struct and constructor**

Add above the `impl Provider` block in `src/benchmark/inmemory_provider.rs`:

```rust
/// Configuration for creating an InMemoryProvider.
#[derive(Debug, Clone)]
pub struct InMemoryProviderConfig {
    /// PostgreSQL connection URL for UTA database.
    pub uta_db_url: String,
    /// UTA database schema (e.g., "uta_20210129b").
    pub uta_db_schema: String,
    /// Optional path to SeqRepo directory for genomic sequences.
    /// If None, only UTA sequences are loaded (sufficient for transcript-based variants).
    pub seqrepo_path: Option<String>,
}

impl InMemoryProvider {
    /// Create a new InMemoryProvider by loading all data from UTA (and optionally SeqRepo).
    ///
    /// This is expensive (seconds to minutes) but only done once. The resulting provider
    /// is cheap to query and can be shared across threads via Arc.
    pub fn new(config: &InMemoryProviderConfig) -> Result<Self, FerroError> {
        use postgres::{Client, NoTls};

        let start = Instant::now();
        eprintln!("Loading UTA metadata into memory...");

        let mut conn = Client::connect(&config.uta_db_url, NoTls).map_err(|e| FerroError::Io {
            msg: format!("Failed to connect to UTA: {}", e),
        })?;
        let schema = &config.uta_db_schema;

        // Load schema/data versions
        let row = conn
            .query_one(
                &format!("SELECT value FROM {schema}.meta WHERE key = 'schema_version'"),
                &[],
            )
            .map_err(|e| FerroError::Io { msg: format!("Failed to get schema version: {}", e) })?;
        let schema_version: String = row.get("value");
        let data_version = schema.to_string();

        // Load sequences: ac -> seq (from seq_anno + seq tables)
        eprintln!("  Loading sequences...");
        let mut sequences = HashMap::new();
        for row in conn
            .query(
                &format!("SELECT sa.ac, s.seq FROM {schema}.seq_anno sa JOIN {schema}.seq s ON sa.seq_id = s.seq_id"),
                &[],
            )
            .map_err(|e| FerroError::Io { msg: format!("Failed to load sequences: {}", e) })?
        {
            let ac: String = row.get("ac");
            let seq: String = row.get("seq");
            sequences.insert(ac, seq);
        }
        eprintln!("    Loaded {} sequences", sequences.len());

        // Load gene info
        eprintln!("  Loading gene info...");
        let mut gene_info = HashMap::new();
        for row in conn
            .query(&format!("SELECT * FROM {schema}.gene"), &[])
            .map_err(|e| FerroError::Io { msg: format!("Failed to load gene info: {}", e) })?
        {
            let hgnc: String = row.get("hgnc");
            let aliases_str: String = row.get("aliases");
            let aliases: Vec<String> = aliases_str.split(',').map(|s| s.to_owned()).collect();
            gene_info.insert(
                hgnc.clone(),
                GeneInfoRecord {
                    hgnc,
                    maploc: row.get("maploc"),
                    descr: row.get("descr"),
                    summary: row.get("summary"),
                    aliases,
                    added: row.get("added"),
                },
            );
        }
        eprintln!("    Loaded {} genes", gene_info.len());

        // Load transcript -> protein accession mapping
        eprintln!("  Loading transcript-protein mappings...");
        let mut pro_ac_for_tx_ac: HashMap<String, Option<String>> = HashMap::new();
        for row in conn
            .query(
                &format!(
                    "SELECT tx_ac, pro_ac FROM {schema}.associated_accessions \
                     WHERE pro_ac IS NOT NULL"
                ),
                &[],
            )
            .map_err(|e| FerroError::Io { msg: format!("Failed to load pro_ac mapping: {}", e) })?
        {
            let tx_ac: String = row.get("tx_ac");
            let pro_ac: String = row.get("pro_ac");
            pro_ac_for_tx_ac.insert(tx_ac, Some(pro_ac));
        }
        eprintln!("    Loaded {} mappings", pro_ac_for_tx_ac.len());

        // Load transcript identity info from tx_def_summary_v view
        eprintln!("  Loading transcript identity info...");
        let mut tx_identity_info = HashMap::new();
        for row in conn
            .query(
                &format!(
                    "SELECT DISTINCT ON (tx_ac) tx_ac, alt_ac, alt_aln_method, \
                     cds_start_i, cds_end_i, lengths, hgnc \
                     FROM {schema}.tx_def_summary_v \
                     ORDER BY tx_ac, alt_ac, alt_aln_method"
                ),
                &[],
            )
            .map_err(|e| FerroError::Io { msg: format!("Failed to load tx identity: {}", e) })?
        {
            let tx_ac: String = row.get("tx_ac");
            let lengths_str: String = row.get("lengths");
            let lengths: Vec<i32> = lengths_str
                .trim_matches(|c| c == '{' || c == '}')
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| s.parse::<i32>().unwrap_or(0))
                .collect();
            tx_identity_info.insert(
                tx_ac.clone(),
                TxIdentityInfo {
                    tx_ac,
                    alt_ac: row.get("alt_ac"),
                    alt_aln_method: row.get("alt_aln_method"),
                    cds_start_i: row.get("cds_start_i"),
                    cds_end_i: row.get("cds_end_i"),
                    lengths,
                    hgnc: row.get("hgnc"),
                    translation_table: Default::default(),
                },
            );
        }
        eprintln!("    Loaded {} transcript identities", tx_identity_info.len());

        // Load transcript mapping options from tx_exon_aln_v
        eprintln!("  Loading transcript mapping options...");
        let mut tx_mapping_options: HashMap<String, Vec<TxMappingOptionsRecord>> = HashMap::new();
        for row in conn
            .query(
                &format!(
                    "SELECT DISTINCT tx_ac, alt_ac, alt_aln_method \
                     FROM {schema}.tx_exon_aln_v \
                     WHERE exon_aln_id IS NOT NULL \
                     ORDER BY tx_ac, alt_ac, alt_aln_method"
                ),
                &[],
            )
            .map_err(|e| FerroError::Io { msg: format!("Failed to load mapping options: {}", e) })?
        {
            let tx_ac: String = row.get("tx_ac");
            let record = TxMappingOptionsRecord {
                tx_ac: tx_ac.clone(),
                alt_ac: row.get("alt_ac"),
                alt_aln_method: row.get("alt_aln_method"),
            };
            tx_mapping_options.entry(tx_ac).or_default().push(record);
        }
        eprintln!("    Loaded mapping options for {} transcripts", tx_mapping_options.len());

        // Load transcript info (tx + exon_set join)
        eprintln!("  Loading transcript info...");
        let mut tx_info: HashMap<(String, String, String), TxInfoRecord> = HashMap::new();
        let mut tx_for_gene: HashMap<String, Vec<TxInfoRecord>> = HashMap::new();
        for row in conn
            .query(
                &format!(
                    "SELECT hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method \
                     FROM {schema}.transcript t \
                     JOIN {schema}.exon_set es ON t.ac = es.tx_ac \
                     WHERE alt_aln_method != 'transcript' \
                     ORDER BY hgnc, cds_start_i, cds_end_i, tx_ac, alt_ac, alt_aln_method"
                ),
                &[],
            )
            .map_err(|e| FerroError::Io { msg: format!("Failed to load tx info: {}", e) })?
        {
            let record = TxInfoRecord {
                hgnc: row.get("hgnc"),
                cds_start_i: row.get("cds_start_i"),
                cds_end_i: row.get("cds_end_i"),
                tx_ac: row.get("tx_ac"),
                alt_ac: row.get("alt_ac"),
                alt_aln_method: row.get("alt_aln_method"),
            };
            let key = (
                record.tx_ac.clone(),
                record.alt_ac.clone(),
                record.alt_aln_method.clone(),
            );
            tx_info.insert(key, record.clone());
            tx_for_gene.entry(record.hgnc.clone()).or_default().push(record);
        }
        eprintln!("    Loaded {} transcript info records", tx_info.len());

        // Load transcript exon alignments
        eprintln!("  Loading transcript exon alignments...");
        let mut tx_exons: HashMap<(String, String, String), Vec<TxExonsRecord>> = HashMap::new();
        let mut tx_for_region: HashMap<(String, String), Vec<TxForRegionRecord>> = HashMap::new();
        for row in conn
            .query(
                &format!("SELECT * FROM {schema}.tx_exon_aln_v ORDER BY tx_ac, alt_ac, alt_aln_method, alt_start_i"),
                &[],
            )
            .map_err(|e| FerroError::Io { msg: format!("Failed to load exon alignments: {}", e) })?
        {
            let record = TxExonsRecord {
                hgnc: row.get("hgnc"),
                tx_ac: row.get("tx_ac"),
                alt_ac: row.get("alt_ac"),
                alt_aln_method: row.get("alt_aln_method"),
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
                exon_aln_id: row.get("exon_aln_id"),
            };
            let key = (
                record.tx_ac.clone(),
                record.alt_ac.clone(),
                record.alt_aln_method.clone(),
            );
            tx_exons.entry(key).or_default().push(record);
        }
        eprintln!("    Loaded exon alignments for {} transcript-contig pairs", tx_exons.len());

        // Build tx_for_region index from tx_exons
        // Group by (alt_ac, alt_aln_method), compute min(start_i)/max(end_i) per transcript
        for (key, exons) in &tx_exons {
            let (tx_ac, alt_ac, alt_aln_method) = key;
            if let (Some(min_start), Some(max_end)) = (
                exons.iter().map(|e| e.alt_start_i).min(),
                exons.iter().map(|e| e.alt_end_i).max(),
            ) {
                let alt_strand = exons.first().map(|e| e.alt_strand).unwrap_or(1);
                let region_key = (alt_ac.clone(), alt_aln_method.clone());
                tx_for_region.entry(region_key).or_default().push(TxForRegionRecord {
                    tx_ac: tx_ac.clone(),
                    alt_ac: alt_ac.clone(),
                    alt_strand,
                    alt_aln_method: alt_aln_method.clone(),
                    start_i: min_start,
                    end_i: max_end,
                });
            }
        }

        let elapsed = start.elapsed();
        eprintln!(
            "UTA metadata loaded in {:.1}s ({} sequences, {} transcripts)",
            elapsed.as_secs_f64(),
            sequences.len(),
            tx_identity_info.len(),
        );

        Ok(Self {
            data_version,
            schema_version,
            assembly_maps: HashMap::new(), // Loaded lazily or from biocommons_bioutils
            gene_info,
            pro_ac_for_tx_ac,
            acs_for_protein_seq: HashMap::new(), // Populated on demand from sequences
            similar_transcripts: HashMap::new(), // Not needed for normalization
            tx_exons,
            tx_for_gene,
            tx_for_region,
            tx_identity_info,
            tx_info,
            tx_mapping_options,
            sequences,
        })
    }
}
```

- [ ] **Step 3: Verify it compiles**

Run: `cargo check --features hgvs-rs`
Expected: PASS

- [ ] **Step 4: Run the test (requires UTA docker)**

Run: `cargo nextest run --features hgvs-rs -E 'test(inmemory_provider)' -- --ignored`
Expected: PASS (loads data from ferro-uta container, verifies ATM transcript)

- [ ] **Step 5: Commit**

```bash
git add src/benchmark/inmemory_provider.rs
git commit -m "feat(benchmark): load UTA metadata into InMemoryProvider"
```

---

### Task 3: Load SeqRepo Sequences and Implement Provider Trait

**Files:**
- Modify: `src/benchmark/inmemory_provider.rs`

This task adds SeqRepo sequence loading (alias resolution + FASTA mmap) and implements the `hgvs::data::interface::Provider` trait.

- [ ] **Step 1: Write the test for SeqRepo integration**

Add to the `tests` module in `inmemory_provider.rs`:

```rust
    /// Test that genomic sequences are available via SeqRepo.
    #[test]
    #[ignore] // Requires running UTA + SeqRepo
    fn test_load_with_seqrepo() {
        let config = InMemoryProviderConfig {
            uta_db_url: "postgresql://anonymous:anonymous@localhost:5432/uta".to_string(),
            uta_db_schema: "uta_20210129b".to_string(),
            seqrepo_path: Some(
                "manuscript/benchmark/output/data/seqrepo/2021-01-29".to_string(),
            ),
        };
        let provider = InMemoryProvider::new(&config).unwrap();

        // Should be able to get a transcript sequence (from UTA)
        let seq = provider.get_seq_part("NM_000051.3", Some(0), Some(10)).unwrap();
        assert_eq!(seq.len(), 10);

        // Should be able to get a genomic sequence (from SeqRepo)
        let seq = provider
            .get_seq_part("NC_000011.10", Some(108222484), Some(108222494))
            .unwrap();
        assert_eq!(seq.len(), 10);
    }

    /// Test Provider trait methods return expected results.
    #[test]
    #[ignore] // Requires running UTA + SeqRepo
    fn test_provider_trait() {
        let config = InMemoryProviderConfig {
            uta_db_url: "postgresql://anonymous:anonymous@localhost:5432/uta".to_string(),
            uta_db_schema: "uta_20210129b".to_string(),
            seqrepo_path: Some(
                "manuscript/benchmark/output/data/seqrepo/2021-01-29".to_string(),
            ),
        };
        let provider = InMemoryProvider::new(&config).unwrap();

        // Test get_tx_mapping_options
        let opts = provider.get_tx_mapping_options("NM_000051.3").unwrap();
        assert!(!opts.is_empty());

        // Test get_tx_identity_info
        let info = provider.get_tx_identity_info("NM_000051.3").unwrap();
        assert_eq!(info.hgnc, "ATM");

        // Test get_tx_exons
        let first_opt = &opts[0];
        let exons = provider
            .get_tx_exons(&first_opt.tx_ac, &first_opt.alt_ac, &first_opt.alt_aln_method)
            .unwrap();
        assert!(!exons.is_empty());
    }
```

- [ ] **Step 2: Add SeqRepo loading to the constructor**

Add this method to `impl InMemoryProvider`, called at the end of `new()` when `seqrepo_path` is `Some`:

```rust
    /// Load sequences from SeqRepo into the sequences HashMap.
    ///
    /// Resolves aliases from SQLite, then reads full sequences from FASTA files.
    /// Sequences already loaded from UTA are not overwritten.
    fn load_seqrepo_sequences(
        sequences: &mut HashMap<String, String>,
        seqrepo_path: &str,
    ) -> Result<(), FerroError> {
        use seqrepo::{AliasOrSeqId, SeqRepo};

        let seqrepo_pathbuf = std::path::PathBuf::from(seqrepo_path);
        let path = seqrepo_pathbuf
            .parent()
            .ok_or(FerroError::Io {
                msg: format!("Invalid seqrepo path: {}", seqrepo_path),
            })?
            .to_str()
            .unwrap()
            .to_string();
        let instance = seqrepo_pathbuf
            .file_name()
            .ok_or(FerroError::Io {
                msg: format!("Invalid seqrepo path: {}", seqrepo_path),
            })?
            .to_str()
            .unwrap()
            .to_string();

        let repo = SeqRepo::new(path, &instance).map_err(|e| FerroError::Io {
            msg: format!("Failed to open SeqRepo: {}", e),
        })?;

        // Load all aliases from SeqRepo's SQLite database
        eprintln!("  Loading SeqRepo alias database...");
        let alias_db = repo.alias_db();
        let query = seqrepo::Query::default();
        let mut alias_to_seqid: HashMap<String, String> = HashMap::new();
        alias_db
            .find(&query, |record| {
                if let Ok(record) = record {
                    alias_to_seqid.insert(record.alias, record.seqid);
                }
            })
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to read SeqRepo aliases: {}", e),
            })?;
        eprintln!("    Loaded {} aliases", alias_to_seqid.len());

        // For each alias that is NOT already in sequences (i.e., not from UTA),
        // load the full sequence from SeqRepo FASTA files.
        eprintln!("  Loading SeqRepo sequences...");
        let mut loaded = 0usize;
        let mut skipped = 0usize;
        // Deduplicate: group aliases by seq_id to avoid loading the same sequence multiple times
        let mut seqid_to_aliases: HashMap<String, Vec<String>> = HashMap::new();
        for (alias, seqid) in &alias_to_seqid {
            if !sequences.contains_key(alias) {
                seqid_to_aliases
                    .entry(seqid.clone())
                    .or_default()
                    .push(alias.clone());
            } else {
                skipped += 1;
            }
        }

        for (seqid, aliases) in &seqid_to_aliases {
            let aos = AliasOrSeqId::SeqId(seqid.clone());
            match repo.fetch_sequence_part(&aos, None, None) {
                Ok(seq) => {
                    for alias in aliases {
                        sequences.insert(alias.clone(), seq.clone());
                        loaded += 1;
                    }
                }
                Err(e) => {
                    eprintln!("    Warning: failed to load sequence for {}: {}", seqid, e);
                }
            }
        }
        eprintln!(
            "    Loaded {} sequences from SeqRepo ({} skipped, already in UTA)",
            loaded, skipped
        );

        Ok(())
    }
```

Then add the call at the end of `new()`, before the final `Ok(Self { ... })`:

```rust
        // Optionally load SeqRepo sequences (genomic contigs, etc.)
        if let Some(ref seqrepo_path) = config.seqrepo_path {
            Self::load_seqrepo_sequences(&mut sequences, seqrepo_path)?;
        }
```

- [ ] **Step 3: Implement the Provider trait**

Add to `src/benchmark/inmemory_provider.rs`:

```rust
impl Provider for InMemoryProvider {
    fn data_version(&self) -> &str {
        &self.data_version
    }

    fn schema_version(&self) -> &str {
        &self.schema_version
    }

    fn get_assembly_map(
        &self,
        assembly: biocommons_bioutils::assemblies::Assembly,
    ) -> IndexMap<String, String> {
        // Use biocommons_bioutils directly for assembly maps (static data)
        use biocommons_bioutils::assemblies::ASSEMBLY_INFOS;
        let mut result = IndexMap::new();
        for assembly_info in ASSEMBLY_INFOS.values() {
            if assembly_info.assembly == assembly {
                for sequence in &assembly_info.sequences {
                    result.insert(
                        sequence.refseq_ac.clone(),
                        sequence.name.clone(),
                    );
                }
                break;
            }
        }
        result
    }

    fn get_gene_info(&self, hgnc: &str) -> Result<GeneInfoRecord, HgvsError> {
        self.gene_info
            .get(hgnc)
            .cloned()
            .ok_or_else(|| HgvsError::NoGeneFound(hgnc.to_string()))
    }

    fn get_pro_ac_for_tx_ac(&self, tx_ac: &str) -> Result<Option<String>, HgvsError> {
        Ok(self.pro_ac_for_tx_ac.get(tx_ac).cloned().flatten())
    }

    fn get_seq_part(
        &self,
        ac: &str,
        begin: Option<usize>,
        end: Option<usize>,
    ) -> Result<String, HgvsError> {
        let seq = self
            .sequences
            .get(ac)
            .ok_or_else(|| HgvsError::NoSequenceRecord(ac.to_string()))?;
        let begin = begin.unwrap_or(0);
        let end = end.map(|e| std::cmp::min(e, seq.len())).unwrap_or(seq.len());
        Ok(seq[begin..end].to_string())
    }

    fn get_acs_for_protein_seq(&self, seq: &str) -> Result<Vec<String>, HgvsError> {
        let md5 = hgvs::sequences::seq_md5(seq, true)?;
        if let Some(result) = self.acs_for_protein_seq.get(&md5) {
            return Ok(result.clone());
        }
        // Return just the MD5 sentinel
        Ok(vec![format!("MD5_{}", md5)])
    }

    fn get_similar_transcripts(&self, tx_ac: &str) -> Result<Vec<TxSimilarityRecord>, HgvsError> {
        Ok(self.similar_transcripts.get(tx_ac).cloned().unwrap_or_default())
    }

    fn get_tx_exons(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<Vec<TxExonsRecord>, HgvsError> {
        let key = (
            tx_ac.to_string(),
            alt_ac.to_string(),
            alt_aln_method.to_string(),
        );
        self.tx_exons.get(&key).cloned().ok_or_else(|| {
            HgvsError::NoTxExons(
                tx_ac.to_string(),
                alt_ac.to_string(),
                alt_aln_method.to_string(),
            )
        })
    }

    fn get_tx_for_gene(&self, gene: &str) -> Result<Vec<TxInfoRecord>, HgvsError> {
        Ok(self.tx_for_gene.get(gene).cloned().unwrap_or_default())
    }

    fn get_tx_for_region(
        &self,
        alt_ac: &str,
        alt_aln_method: &str,
        start_i: i32,
        end_i: i32,
    ) -> Result<Vec<TxForRegionRecord>, HgvsError> {
        let key = (alt_ac.to_string(), alt_aln_method.to_string());
        let records = self.tx_for_region.get(&key).cloned().unwrap_or_default();
        // Filter to records that overlap the query region
        Ok(records
            .into_iter()
            .filter(|r| r.start_i < end_i && start_i <= r.end_i)
            .collect())
    }

    fn get_tx_identity_info(&self, tx_ac: &str) -> Result<TxIdentityInfo, HgvsError> {
        self.tx_identity_info
            .get(tx_ac)
            .cloned()
            .ok_or_else(|| HgvsError::NoTranscriptFound(tx_ac.to_string()))
    }

    fn get_tx_info(
        &self,
        tx_ac: &str,
        alt_ac: &str,
        alt_aln_method: &str,
    ) -> Result<TxInfoRecord, HgvsError> {
        let key = (
            tx_ac.to_string(),
            alt_ac.to_string(),
            alt_aln_method.to_string(),
        );
        self.tx_info.get(&key).cloned().ok_or_else(|| {
            HgvsError::NoTranscriptFound(format!("{}/{}/{}", tx_ac, alt_ac, alt_aln_method))
        })
    }

    fn get_tx_mapping_options(
        &self,
        tx_ac: &str,
    ) -> Result<Vec<TxMappingOptionsRecord>, HgvsError> {
        Ok(self.tx_mapping_options.get(tx_ac).cloned().unwrap_or_default())
    }
}
```

- [ ] **Step 4: Verify it compiles**

Run: `cargo check --features hgvs-rs`
Expected: PASS

- [ ] **Step 5: Run the tests**

Run: `cargo nextest run --features hgvs-rs -E 'test(inmemory_provider)' -- --ignored`
Expected: All 3 tests PASS

- [ ] **Step 6: Commit**

```bash
git add src/benchmark/inmemory_provider.rs
git commit -m "feat(benchmark): implement Provider trait for InMemoryProvider with SeqRepo"
```

---

### Task 4: Wire InMemoryProvider into HgvsRsNormalizer

**Files:**
- Modify: `src/benchmark/hgvs_rs.rs`

Add a second constructor to `HgvsRsNormalizer` that accepts an `InMemoryProvider`.

- [ ] **Step 1: Write the test**

Add to the `tests` module in `hgvs_rs.rs`:

```rust
    #[test]
    #[ignore] // Requires running UTA + SeqRepo
    fn test_normalize_with_inmemory_provider() {
        use super::super::inmemory_provider::{InMemoryProvider, InMemoryProviderConfig};

        let config = InMemoryProviderConfig {
            uta_db_url: "postgresql://anonymous:anonymous@localhost:5432/uta".to_string(),
            uta_db_schema: "uta_20210129b".to_string(),
            seqrepo_path: Some(
                "manuscript/benchmark/output/data/seqrepo/2021-01-29".to_string(),
            ),
        };
        let provider = std::sync::Arc::new(InMemoryProvider::new(&config).unwrap());
        let normalizer = HgvsRsNormalizer::with_provider(provider).unwrap();

        // Test a known coding variant
        let result = normalizer.normalize("NM_000051.3:c.1058_1059delinsAA");
        assert!(result.success, "Expected success, got: {:?}", result.error);
    }
```

- [ ] **Step 2: Add `with_provider` constructor**

Add to `impl HgvsRsNormalizer` in `hgvs_rs.rs`:

```rust
    /// Create a new normalizer using a pre-built Provider (e.g., InMemoryProvider).
    pub fn with_provider(
        provider: Arc<dyn hgvs::data::interface::Provider + Send + Sync>,
    ) -> Result<Self, FerroError> {
        let mapper_config = MapperConfig::default();
        let mapper = Mapper::new(&mapper_config, provider);
        Ok(Self { mapper })
    }
```

- [ ] **Step 3: Verify it compiles and test passes**

Run: `cargo check --features hgvs-rs`
Run: `cargo nextest run --features hgvs-rs -E 'test(normalize_with_inmemory)' -- --ignored`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add src/benchmark/hgvs_rs.rs
git commit -m "feat(benchmark): add HgvsRsNormalizer::with_provider constructor"
```

---

### Task 5: Add --in-memory CLI Flag and Parallel Support

**Files:**
- Modify: `src/benchmark/hgvs_rs.rs`
- Modify: `src/bin/benchmark.rs`

Wire the `--in-memory` flag through the CLI to use `InMemoryProvider` for normalization.

- [ ] **Step 1: Add `in_memory` field to `HgvsRsConfig`**

In `hgvs_rs.rs`, add to `HgvsRsConfig`:

```rust
    /// Use in-memory provider instead of UTA+SeqRepo for benchmarking.
    pub in_memory: bool,
```

And update `Default`:

```rust
    in_memory: false,
```

- [ ] **Step 2: Update `run_hgvs_rs_normalize_parallel` to support in-memory mode**

Replace the normalizer creation logic in `run_hgvs_rs_normalize_parallel`. When `config.in_memory` is true, create a single `InMemoryProvider` wrapped in `Arc`, then create normalizers from it. This avoids multiple PostgreSQL connections entirely.

In the pre-creation block (the `eprintln!("Creating {} hgvs-rs normalizers..."` section), replace it with:

```rust
    eprintln!("Creating {} hgvs-rs normalizers (one per worker)...", workers);
    let init_start = Instant::now();
    let normalizers: Vec<_> = if config.in_memory {
        // Single shared in-memory provider for all workers
        let im_config = crate::benchmark::inmemory_provider::InMemoryProviderConfig {
            uta_db_url: config.uta_db_url.clone(),
            uta_db_schema: config.uta_db_schema.clone(),
            seqrepo_path: Some(config.seqrepo_path.clone()),
        };
        let provider = Arc::new(
            crate::benchmark::inmemory_provider::InMemoryProvider::new(&im_config)
                .map_err(|e| FerroError::Io {
                    msg: format!("Failed to create in-memory provider: {}", e),
                })?,
        );
        (0..workers)
            .map(|_| {
                HgvsRsNormalizer::with_provider(provider.clone())
                    .map_err(|e| format!("Failed to create normalizer: {}", e))
            })
            .collect()
    } else {
        (0..workers)
            .map(|_| {
                HgvsRsNormalizer::new(config)
                    .map_err(|e| format!("Failed to create normalizer: {}", e))
            })
            .collect()
    };
    let init_elapsed = init_start.elapsed();
    eprintln!(
        "Normalizer initialization took {:.3}s ({:.3}s per worker)",
        init_elapsed.as_secs_f64(),
        init_elapsed.as_secs_f64() / workers as f64,
    );
```

- [ ] **Step 3: Also update `run_hgvs_rs_normalize` (sequential) similarly**

In `run_hgvs_rs_normalize`, replace the normalizer creation:

```rust
    // Create the normalizer
    let normalizer = if config.in_memory {
        let im_config = crate::benchmark::inmemory_provider::InMemoryProviderConfig {
            uta_db_url: config.uta_db_url.clone(),
            uta_db_schema: config.uta_db_schema.clone(),
            seqrepo_path: Some(config.seqrepo_path.clone()),
        };
        let provider = Arc::new(
            crate::benchmark::inmemory_provider::InMemoryProvider::new(&im_config)?,
        );
        HgvsRsNormalizer::with_provider(provider)?
    } else {
        HgvsRsNormalizer::new(config)?
    };
```

- [ ] **Step 4: Add `--in-memory` CLI flag in benchmark.rs**

Find the hgvs-rs normalize CLI argument definitions in `src/bin/benchmark.rs` and add:

```rust
    /// Use in-memory provider (pre-loads all data, eliminates DB I/O from timing)
    #[arg(long)]
    in_memory: bool,
```

Wire it through to `HgvsRsConfig`:

```rust
    config.in_memory = args.in_memory;
```

(The exact location depends on the existing CLI structure — search for where `HgvsRsConfig` is constructed from CLI args.)

- [ ] **Step 5: Verify it compiles**

Run: `cargo check --features hgvs-rs`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add src/benchmark/hgvs_rs.rs src/bin/benchmark.rs
git commit -m "feat(benchmark): add --in-memory flag for hgvs-rs normalization"
```

---

### Task 6: Smoke Test and Validation

**Files:** None (testing only)

Run the in-memory benchmark and compare results with the PostgreSQL-backed benchmark to validate correctness and measure improvement.

- [ ] **Step 1: Build release binary**

Run: `cargo build --release --features benchmark,hgvs-rs`

- [ ] **Step 2: Create test sample**

```bash
target/release/ferro-benchmark extract sample \
  -i manuscript/benchmark/output/data/clinvar/clinvar_patterns.txt \
  -o /tmp/sample_100.txt --size 100 --seed 42 --exclude-protein
```

- [ ] **Step 3: Run PostgreSQL-backed baseline**

```bash
target/release/ferro-benchmark normalize hgvs-rs \
  -i /tmp/sample_100.txt -o /tmp/pg_results.json \
  --seqrepo-path manuscript/benchmark/output/data/seqrepo/2021-01-29 \
  --workers 4
```

- [ ] **Step 4: Run in-memory benchmark**

```bash
target/release/ferro-benchmark normalize hgvs-rs \
  -i /tmp/sample_100.txt -o /tmp/inmem_results.json \
  --seqrepo-path manuscript/benchmark/output/data/seqrepo/2021-01-29 \
  --workers 4 --in-memory
```

- [ ] **Step 5: Compare results for correctness**

```bash
# Both should have the same success count and matching outputs
target/release/ferro-benchmark compare results normalize \
  /tmp/pg_results.json /tmp/inmem_results.json \
  -o /tmp/comparison.json
```

Expected: Same success/failure counts, same normalized outputs for all patterns.

- [ ] **Step 6: Run at 10K scale for timing comparison**

```bash
target/release/ferro-benchmark extract sample \
  -i manuscript/benchmark/output/data/clinvar/clinvar_patterns.txt \
  -o /tmp/sample_10k.txt --size 10000 --seed 42 --exclude-protein

# PostgreSQL-backed
target/release/ferro-benchmark normalize hgvs-rs \
  -i /tmp/sample_10k.txt -o /tmp/pg_10k.json \
  --seqrepo-path manuscript/benchmark/output/data/seqrepo/2021-01-29 \
  --workers 4

# In-memory
target/release/ferro-benchmark normalize hgvs-rs \
  -i /tmp/sample_10k.txt -o /tmp/inmem_10k.json \
  --seqrepo-path manuscript/benchmark/output/data/seqrepo/2021-01-29 \
  --workers 4 --in-memory
```

Expected: In-memory version should be significantly faster (10x+ improvement) with identical results.

- [ ] **Step 7: Commit any fixes discovered during testing**

---

## Notes for Implementer

- The `hgvs::data::error::Error` variants used in the Provider trait (e.g., `NoGeneFound`, `NoSequenceRecord`, `NoTranscriptFound`, `NoTxExons`) need to match exactly what hgvs-rs defines. Check `hgvs-0.20.1/src/data/error.rs` for exact variant names.
- The `biocommons_bioutils::assemblies::ASSEMBLY_INFOS` structure may differ from what the code assumes. Check the actual API before implementing `get_assembly_map`.
- SeqRepo sequence loading will use ~10GB of RAM. This is fine on the 64GB development machine but may not work on the 12GB server.
- The `translation_table` field in `TxIdentityInfo` defaults to `TranslationTable::Standard`. UTA doesn't store this, so we use the default (same as the existing UTA provider).
- The `tx_exon_aln_v` and `tx_def_summary_v` are PostgreSQL views. The queries may be slow on the emulated amd64 UTA container — this is expected since it only runs once at startup.
