//! biocommons/hgvs integration for normalization benchmarks.
//!
//! biocommons/hgvs is the canonical Python HGVS implementation.
//! See: https://github.com/biocommons/hgvs
//!
//! # Local Setup
//!
//! For high-performance normalization, biocommons/hgvs requires two local components:
//!
//! 1. **UTA (Universal Transcript Archive)** - PostgreSQL database with transcript data
//! 2. **SeqRepo** - Local sequence repository for nucleotide/protein sequences
//!
//! Use the CLI commands `setup-uta` and `setup-seqrepo` to configure these.

#![allow(clippy::type_complexity)]

use crate::benchmark::cache::SupplementalMetadata;
use crate::benchmark::translate::translate_cds_to_protein;
use crate::benchmark::types::ParseResult;
use crate::data::cdot::CdotMapper;
use crate::FerroError;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// ============================================================================
// Configuration Types
// ============================================================================

/// Configuration for local biocommons/hgvs setup.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiocommonsLocalConfig {
    /// Docker container name for UTA.
    pub uta_container_name: String,
    /// UTA Docker image tag (e.g., "uta_20210129b").
    pub uta_image_tag: String,
    /// PostgreSQL port for UTA.
    pub uta_port: u16,
    /// Directory for SeqRepo data.
    pub seqrepo_dir: PathBuf,
    /// SeqRepo instance version (e.g., "2021-01-29").
    pub seqrepo_instance: String,
}

impl Default for BiocommonsLocalConfig {
    fn default() -> Self {
        Self {
            uta_container_name: "ferro-uta".to_string(),
            uta_image_tag: "uta_20210129b".to_string(),
            uta_port: 5432,
            seqrepo_dir: PathBuf::new(), // Must be specified by user
            seqrepo_instance: "2021-01-29".to_string(),
        }
    }
}

/// Settings loaded from a biocommons configuration file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiocommonsSettings {
    /// UTA database URL (e.g., postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b).
    pub uta_db_url: String,
    /// Path to SeqRepo data directory.
    pub seqrepo_dir: PathBuf,
}

/// Result of UTA setup operation.
#[derive(Debug)]
pub struct UtaSetupResult {
    /// Whether the Docker image was pulled (false if already existed).
    pub image_pulled: bool,
    /// Whether a new container was created (false if already existed).
    pub container_created: bool,
    /// The container name.
    pub container_name: String,
    /// The PostgreSQL port.
    pub port: u16,
    /// The UTA database URL.
    pub uta_db_url: String,
}

/// Result of SeqRepo setup operation.
#[derive(Debug)]
pub struct SeqRepoSetupResult {
    /// Path to the SeqRepo directory.
    pub seqrepo_dir: PathBuf,
    /// SeqRepo instance that was downloaded.
    pub instance: String,
    /// Whether data was downloaded (false if already existed).
    pub downloaded: bool,
}

// ============================================================================
// Docker and System Checks
// ============================================================================

/// Check if Docker is available and running.
pub fn check_docker_available() -> bool {
    Command::new("docker")
        .args(["version"])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Check if a Docker container exists (running or stopped).
pub fn check_container_exists(container_name: &str) -> bool {
    Command::new("docker")
        .args(["container", "inspect", container_name])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Check if a Docker container is running.
pub fn check_container_running(container_name: &str) -> bool {
    let output = Command::new("docker")
        .args([
            "container",
            "inspect",
            "-f",
            "{{.State.Running}}",
            container_name,
        ])
        .output();

    match output {
        Ok(o) if o.status.success() => String::from_utf8_lossy(&o.stdout).trim() == "true",
        _ => false,
    }
}

/// Check if we can connect to a local UTA PostgreSQL database.
///
/// This attempts a lightweight connection test to verify the database is accessible.
/// First tries a direct psql check (faster, no Python needed), then falls back to
/// the Python hgvs library check if available.
pub fn check_uta_local(port: u16) -> bool {
    // First try direct PostgreSQL check (works without Python)
    if check_uta_local_direct(port) {
        return true;
    }

    // Fall back to Python-based check
    let uta_url = format!(
        "postgresql://anonymous:anonymous@localhost:{}/uta/uta_20210129b",
        port
    );
    check_uta_connection(Some(&uta_url))
}

/// Direct PostgreSQL check for local UTA database (no Python required).
fn check_uta_local_direct(_port: u16) -> bool {
    // Try to query the transcript table as the anonymous user
    let output = Command::new("docker")
        .args([
            "exec",
            "ferro-uta",
            "psql",
            "-U",
            "anonymous",
            "-d",
            "uta",
            "-c",
            "SELECT 1 FROM uta_20210129b.transcript LIMIT 1;",
        ])
        .output();

    match output {
        Ok(o) => o.status.success(),
        Err(_) => false,
    }
}

/// Check if SeqRepo is available at the specified directory.
///
/// This verifies that the SeqRepo directory exists and contains valid data.
pub fn check_seqrepo(seqrepo_dir: &Path) -> bool {
    // Check if the directory exists
    if !seqrepo_dir.exists() {
        return false;
    }

    // Check for the sequences subdirectory (indicates valid SeqRepo structure)
    let sequences_dir = seqrepo_dir.join("sequences");
    if !sequences_dir.exists() {
        return false;
    }

    // Try to use Python to verify SeqRepo is accessible
    let check_script = format!(
        r#"
import os
os.environ['HGVS_SEQREPO_DIR'] = '{}'
from biocommons.seqrepo import SeqRepo
sr = SeqRepo('{}')
# Try to access a common sequence to verify it works
print("OK")
"#,
        seqrepo_dir.display(),
        seqrepo_dir.display()
    );

    Command::new("python3")
        .args(["-c", &check_script])
        .output()
        .map(|o| o.status.success() && String::from_utf8_lossy(&o.stdout).contains("OK"))
        .unwrap_or(false)
}

// ============================================================================
// Setup Functions
// ============================================================================

/// Set up a local UTA database using Docker.
///
/// This pulls the biocommons/uta Docker image and starts a PostgreSQL container.
///
/// # Arguments
/// * `config` - Configuration for the UTA setup
/// * `force` - If true, recreate container even if it exists
///
/// # Returns
/// A result with details about the setup operation.
pub fn setup_uta(
    config: &BiocommonsLocalConfig,
    force: bool,
    uta_dump: Option<&Path>,
) -> Result<UtaSetupResult, FerroError> {
    // Check Docker is available
    if !check_docker_available() {
        return Err(FerroError::Io {
            msg: "Docker is not available. Please install Docker and ensure it is running."
                .to_string(),
        });
    }

    let mut image_pulled = false;
    let mut container_created = false;

    // Check if container already exists
    let container_exists = check_container_exists(&config.uta_container_name);

    if container_exists && force {
        // Remove existing container
        eprintln!("Removing existing container: {}", config.uta_container_name);
        let _ = Command::new("docker")
            .args(["rm", "-f", &config.uta_container_name])
            .output();
    }

    // If a pre-downloaded dump file is provided, use custom setup with postgres base image
    if let Some(dump_path) = uta_dump {
        return setup_uta_from_dump(config, dump_path, force, container_exists);
    }

    // Standard flow: try to pull the biocommons/uta image
    let image_name = format!("biocommons/uta:{}", config.uta_image_tag);

    // Check if image exists
    let image_exists = Command::new("docker")
        .args(["image", "inspect", &image_name])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);

    // Pull image if needed or forced
    if !image_exists || force {
        eprintln!(
            "Pulling Docker image: {} (this may take a while)...",
            image_name
        );
        let pull_output = Command::new("docker")
            .args(["pull", &image_name])
            .output()
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to pull Docker image: {}", e),
            })?;

        if !pull_output.status.success() {
            let stderr = String::from_utf8_lossy(&pull_output.stderr);
            // Provide helpful message about manual download option
            eprintln!("\nFailed to pull Docker image: {}", stderr);
            eprintln!(
                "\nThe biocommons UTA image may require manual download due to human verification."
            );
            eprintln!("You can download the database dump manually and use --uta-dump:\n");
            eprintln!(
                "  1. Download: https://dl.biocommons.org/uta/{}.pgd.gz",
                config.uta_image_tag
            );
            eprintln!("     (Complete any human verification in your browser)");
            eprintln!(
                "  2. Run: ferro-benchmark setup-uta --uta-dump /path/to/{}.pgd.gz\n",
                config.uta_image_tag
            );
            return Err(FerroError::Io {
                msg: "Failed to pull Docker image. See above for manual download instructions."
                    .to_string(),
            });
        }
        image_pulled = true;
        eprintln!("Docker image pulled successfully.");
    }

    if !container_exists || force {
        // Create and start container
        eprintln!("Creating container: {}", config.uta_container_name);
        let create_output = Command::new("docker")
            .args([
                "run",
                "-d",
                "--name",
                &config.uta_container_name,
                "-p",
                &format!("{}:5432", config.uta_port),
                "-e",
                "POSTGRES_HOST_AUTH_METHOD=trust",
                &image_name,
            ])
            .output()
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to create Docker container: {}", e),
            })?;

        if !create_output.status.success() {
            let stderr = String::from_utf8_lossy(&create_output.stderr);
            return Err(FerroError::Io {
                msg: format!("Failed to create Docker container: {}", stderr),
            });
        }
        container_created = true;
        eprintln!("Container created successfully.");
    } else if !check_container_running(&config.uta_container_name) {
        // Container exists but is stopped - start it
        eprintln!("Starting existing container: {}", config.uta_container_name);
        start_uta(&config.uta_container_name)?;
    }

    // Wait for PostgreSQL to be ready
    wait_for_uta_ready(config.uta_port, 30)?;

    let uta_db_url = format!(
        "postgresql://anonymous:anonymous@localhost:{}/uta/{}",
        config.uta_port, config.uta_image_tag
    );

    Ok(UtaSetupResult {
        image_pulled,
        container_created,
        container_name: config.uta_container_name.clone(),
        port: config.uta_port,
        uta_db_url,
    })
}

/// Set up UTA from a pre-downloaded database dump file.
///
/// This creates a PostgreSQL container and loads the UTA database from the dump file.
fn setup_uta_from_dump(
    config: &BiocommonsLocalConfig,
    dump_path: &Path,
    force: bool,
    container_existed: bool,
) -> Result<UtaSetupResult, FerroError> {
    // Validate the dump file exists and has reasonable size
    if !dump_path.exists() {
        return Err(FerroError::Io {
            msg: format!(
                "UTA dump file not found: {}\n\
                 Download from: https://dl.biocommons.org/uta/{}.pgd.gz",
                dump_path.display(),
                config.uta_image_tag
            ),
        });
    }

    let metadata = std::fs::metadata(dump_path).map_err(|e| FerroError::Io {
        msg: format!("Cannot read dump file metadata: {}", e),
    })?;

    // UTA dumps are typically 1-3 GB compressed
    const MIN_DUMP_SIZE: u64 = 100 * 1024 * 1024; // 100 MB minimum
    if metadata.len() < MIN_DUMP_SIZE {
        return Err(FerroError::Io {
            msg: format!(
                "UTA dump file seems too small ({} bytes). Expected at least {} MB.\n\
                 This may indicate an incomplete download or wrong file.\n\
                 Download from: https://dl.biocommons.org/uta/{}.pgd.gz",
                metadata.len(),
                MIN_DUMP_SIZE / (1024 * 1024),
                config.uta_image_tag
            ),
        });
    }

    eprintln!(
        "Using pre-downloaded UTA dump: {} ({:.1} GB)",
        dump_path.display(),
        metadata.len() as f64 / (1024.0 * 1024.0 * 1024.0)
    );

    // Use postgres:14 as the base image (matching biocommons/uta)
    let image_name = "postgres:14";
    let mut image_pulled = false;

    // Check if postgres image exists
    let image_exists = Command::new("docker")
        .args(["image", "inspect", image_name])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);

    if !image_exists {
        eprintln!("Pulling PostgreSQL base image...");
        let pull_output = Command::new("docker")
            .args(["pull", image_name])
            .output()
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to pull postgres image: {}", e),
            })?;

        if !pull_output.status.success() {
            let stderr = String::from_utf8_lossy(&pull_output.stderr);
            return Err(FerroError::Io {
                msg: format!("Failed to pull postgres image: {}", stderr),
            });
        }
        image_pulled = true;
    }

    // Get absolute path for the dump file
    let dump_abs_path = dump_path.canonicalize().map_err(|e| FerroError::Io {
        msg: format!("Cannot resolve dump file path: {}", e),
    })?;

    let container_created = !container_existed || force;

    if container_created {
        // Create container with the dump file mounted
        eprintln!("Creating container: {}", config.uta_container_name);

        let create_output = Command::new("docker")
            .args([
                "run",
                "-d",
                "--name",
                &config.uta_container_name,
                "-p",
                &format!("{}:5432", config.uta_port),
                "-e",
                "POSTGRES_HOST_AUTH_METHOD=trust",
                "-e",
                "POSTGRES_DB=uta",
                "-v",
                &format!(
                    "{}:/docker-entrypoint-initdb.d/uta.pgd.gz:ro",
                    dump_abs_path.display()
                ),
                image_name,
            ])
            .output()
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to create Docker container: {}", e),
            })?;

        if !create_output.status.success() {
            let stderr = String::from_utf8_lossy(&create_output.stderr);
            return Err(FerroError::Io {
                msg: format!("Failed to create Docker container: {}", stderr),
            });
        }
        eprintln!("Container created successfully.");

        // Wait for PostgreSQL to start
        eprintln!("Waiting for PostgreSQL to start...");
        wait_for_postgres_ready(&config.uta_container_name, 30)?;

        // Load the UTA database dump
        eprintln!("Loading UTA database from dump (this may take several minutes)...");
        load_uta_dump(&config.uta_container_name, &config.uta_image_tag)?;
    } else if !check_container_running(&config.uta_container_name) {
        eprintln!("Starting existing container: {}", config.uta_container_name);
        start_uta(&config.uta_container_name)?;
    }

    // Wait for UTA to be fully ready
    wait_for_uta_ready(config.uta_port, 60)?;

    let uta_db_url = format!(
        "postgresql://anonymous:anonymous@localhost:{}/uta/{}",
        config.uta_port, config.uta_image_tag
    );

    Ok(UtaSetupResult {
        image_pulled,
        container_created,
        container_name: config.uta_container_name.clone(),
        port: config.uta_port,
        uta_db_url,
    })
}

/// Wait for PostgreSQL to be ready (basic connectivity).
fn wait_for_postgres_ready(container_name: &str, max_attempts: u32) -> Result<(), FerroError> {
    for attempt in 1..=max_attempts {
        let output = Command::new("docker")
            .args(["exec", container_name, "pg_isready", "-U", "postgres"])
            .output();

        if let Ok(o) = output {
            if o.status.success() {
                eprintln!("PostgreSQL is accepting connections.");
                return Ok(());
            }
        }

        if attempt == max_attempts {
            return Err(FerroError::Io {
                msg: "Timed out waiting for PostgreSQL to start".to_string(),
            });
        }
        std::thread::sleep(std::time::Duration::from_secs(2));
        eprint!(".");
    }
    Ok(())
}

/// Load the UTA database from the mounted dump file.
fn load_uta_dump(container_name: &str, schema_name: &str) -> Result<(), FerroError> {
    // Create the UTA database (may already exist from POSTGRES_DB env var)
    let _ = Command::new("docker")
        .args([
            "exec",
            container_name,
            "psql",
            "-U",
            "postgres",
            "-c",
            "CREATE DATABASE uta;",
        ])
        .output();
    // Ignore errors - database may already exist

    // Load the dump file using gzip and psql
    // The dump file is mounted at /docker-entrypoint-initdb.d/uta.pgd.gz
    eprintln!("Decompressing and loading dump file (this takes 5-15 minutes)...");

    let load_output = Command::new("docker")
        .args([
            "exec",
            container_name,
            "bash",
            "-c",
            "gzip -dc /docker-entrypoint-initdb.d/uta.pgd.gz | psql -U postgres -d uta",
        ])
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to load UTA dump: {}", e),
        })?;

    if !load_output.status.success() {
        let stderr = String::from_utf8_lossy(&load_output.stderr);
        // Some errors are expected (like "role already exists"), check for critical failures
        if stderr.contains("FATAL") || stderr.contains("could not") {
            return Err(FerroError::Io {
                msg: format!("Failed to load UTA dump: {}", stderr),
            });
        }
        eprintln!("Note: Some warnings during load (usually harmless)");
    }

    eprintln!("Database dump loaded, setting up permissions...");

    // Set up anonymous role and grant permissions on the loaded schema
    // This must happen AFTER the dump is loaded since the dump creates the schema and tables
    let setup_sql = format!(
        "DO $$ BEGIN CREATE ROLE anonymous WITH LOGIN; EXCEPTION WHEN duplicate_object THEN NULL; END $$; \
         GRANT USAGE ON SCHEMA {} TO anonymous; \
         GRANT SELECT ON ALL TABLES IN SCHEMA {} TO anonymous;",
        schema_name, schema_name
    );

    let output = Command::new("docker")
        .args([
            "exec",
            container_name,
            "psql",
            "-U",
            "postgres",
            "-d",
            "uta",
            "-c",
            &setup_sql,
        ])
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to set up UTA permissions: {}", e),
        })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        eprintln!("Warning: Permission setup returned: {}", stderr);
    }

    eprintln!("UTA database loaded successfully.");
    Ok(())
}

/// Wait for UTA database to be ready with the schema.
fn wait_for_uta_ready(port: u16, max_attempts: u32) -> Result<(), FerroError> {
    eprintln!("Waiting for UTA database to be ready...");
    for attempt in 1..=max_attempts {
        if check_uta_local(port) {
            eprintln!("UTA database is ready.");
            return Ok(());
        }
        if attempt == max_attempts {
            return Err(FerroError::Io {
                msg: "Timed out waiting for UTA database to be ready".to_string(),
            });
        }
        std::thread::sleep(std::time::Duration::from_secs(2));
        eprint!(".");
    }
    Ok(())
}

/// Find GNU rsync executable (not openrsync which macOS uses by default).
///
/// SeqRepo requires GNU rsync, but macOS ships with openrsync.
/// This function checks common locations for GNU rsync.
fn find_gnu_rsync() -> Option<String> {
    // Common locations for GNU rsync
    let candidates = [
        // Homebrew locations
        "/opt/homebrew/bin/rsync",
        "/usr/local/bin/rsync",
        // Conda/pixi environments (relative to current dir)
        ".pixi/envs/default/bin/rsync",
    ];

    for candidate in candidates {
        let path = if candidate.starts_with('.') {
            // Relative path - resolve from current directory
            std::env::current_dir().ok().map(|cwd| cwd.join(candidate))
        } else {
            Some(PathBuf::from(candidate))
        };

        if let Some(path) = path {
            if path.exists() {
                // Verify it's GNU rsync by checking version output
                if let Ok(output) = Command::new(&path).arg("--version").output() {
                    let version = String::from_utf8_lossy(&output.stdout);
                    // GNU rsync has "rsync  version X.X.X" format
                    // openrsync has "openrsync: protocol version X"
                    if version.contains("rsync  version") && !version.contains("openrsync") {
                        return Some(path.to_string_lossy().to_string());
                    }
                }
            }
        }
    }

    None
}

/// Add LRG genomic aliases to SeqRepo.
///
/// LRG FASTA files use headers like `>LRG_1g (genomic sequence)` but HGVS patterns
/// use `LRG_1:g.5000A>T` where the accession is `LRG_1`, not `LRG_1g`.
///
/// This function adds aliases from `LRG_X` -> same seq_id as `LRG_Xg` so that
/// biocommons/hgvs can find the sequences when normalizing.
fn add_lrg_genomic_aliases(seqrepo_dir: &Path) -> Result<(), FerroError> {
    // Python script to add aliases
    let python_script = r#"
import sys
from biocommons.seqrepo import SeqRepo

seqrepo_path = sys.argv[1]
sr = SeqRepo(seqrepo_path, writeable=True)

# Find all LRG_Xg aliases and add LRG_X
added = 0
for row in sr.aliases._db.execute("SELECT seq_id, alias FROM seqalias WHERE namespace='LRG' AND alias LIKE '%g'"):
    seq_id, alias = row[0], row[1]
    # Only process genomic aliases (end with 'g' but not 't1g', 't2g', etc.)
    if alias.endswith('g') and not any(alias.endswith(f't{i}g') for i in range(1, 10)):
        base_alias = alias[:-1]  # LRG_1g -> LRG_1
        # Check if alias already exists
        existing = list(sr.aliases._db.execute(
            "SELECT 1 FROM seqalias WHERE namespace='LRG' AND alias=?", (base_alias,)))
        if not existing:
            sr.aliases.store_alias(seq_id, "LRG", base_alias)
            added += 1

sr.commit()
print(f"Added {added} LRG genomic aliases")
"#;

    // Try pixi run first, then fall back to plain python3
    let output = Command::new("pixi")
        .arg("run")
        .arg("python3")
        .arg("-c")
        .arg(python_script)
        .arg(seqrepo_dir)
        .output();

    let output = match output {
        Ok(o) if o.status.success() => o,
        _ => {
            // Fall back to plain python3
            Command::new("python3")
                .arg("-c")
                .arg(python_script)
                .arg(seqrepo_dir)
                .output()
                .map_err(|e| FerroError::Io {
                    msg: format!("Failed to run Python script for LRG aliases: {}", e),
                })?
        }
    };

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(FerroError::Io {
            msg: format!("Failed to add LRG aliases: {}", stderr),
        });
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    eprintln!("{}", stdout.trim());

    Ok(())
}

/// Read transcript sequences from a FASTA file into a HashMap.
fn read_fasta_to_map(fasta_path: &Path) -> Result<HashMap<String, String>, FerroError> {
    let file = File::open(fasta_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", fasta_path.display(), e),
    })?;
    let reader = BufReader::new(file);

    let mut sequences = HashMap::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;

        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                if !current_seq.is_empty() {
                    sequences.insert(id, std::mem::take(&mut current_seq));
                }
            }
            current_id = header.split_whitespace().next().map(|s| s.to_string());
        } else if current_id.is_some() {
            current_seq.push_str(line.trim());
        }
    }

    if let Some(id) = current_id {
        if !current_seq.is_empty() {
            sequences.insert(id, current_seq);
        }
    }

    Ok(sequences)
}

/// Load additional sequences into SeqRepo from ferro reference data.
///
/// This loads genomic (NC_), RefSeqGene (NG_), and LRG sequences into an existing
/// SeqRepo instance, enabling biocommons/hgvs and hgvs-rs to handle these pattern
/// types without network requests.
///
/// # Arguments
/// * `seqrepo_dir` - Path to SeqRepo directory (e.g., data/seqrepo/2021-01-29)
/// * `ferro_ref` - Path to ferro reference directory containing genome/, refseqgene/, lrg/
///
/// # Sequences loaded
/// * `genome/GRCh38.fna` - GRCh38 chromosome sequences (NC_ accessions)
/// * `genome/GRCh37.fna` - GRCh37 chromosome sequences (NC_ accessions)
/// * `refseqgene/*.fna` - RefSeqGene sequences (NG_ accessions)
/// * `lrg/*.fasta` - LRG sequences (LRG_ accessions)
/// * `transcripts/*.fna` - All transcripts (NM_, NR_, XM_, XR_) if load_transcripts is true
///
/// # Arguments
/// * `seqrepo_dir` - Path to SeqRepo instance directory
/// * `ferro_ref` - Path to ferro reference directory
/// * `load_transcripts` - If true, also load all transcript FASTAs (~270K sequences)
pub fn load_sequences_to_seqrepo(
    seqrepo_dir: &Path,
    ferro_ref: &Path,
    load_transcripts: bool,
) -> Result<(), FerroError> {
    // Verify SeqRepo exists
    if !seqrepo_dir.exists() {
        return Err(FerroError::Io {
            msg: format!("SeqRepo directory not found: {}", seqrepo_dir.display()),
        });
    }

    // Make SeqRepo writable (it's typically a read-only snapshot)
    eprintln!("Making SeqRepo writable...");
    let chmod_output = Command::new("chmod")
        .args(["-R", "u+w"])
        .arg(seqrepo_dir)
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to chmod SeqRepo: {}", e),
        })?;

    if !chmod_output.status.success() {
        let stderr = String::from_utf8_lossy(&chmod_output.stderr);
        return Err(FerroError::Io {
            msg: format!("Failed to make SeqRepo writable: {}", stderr),
        });
    }

    // Find the instance directory (e.g., 2021-01-29)
    let instance_name = seqrepo_dir
        .file_name()
        .and_then(|n| n.to_str())
        .ok_or_else(|| FerroError::Io {
            msg: "Cannot determine SeqRepo instance name".to_string(),
        })?;

    // Get parent directory for --root-directory
    let root_dir = seqrepo_dir.parent().ok_or_else(|| FerroError::Io {
        msg: "Cannot determine SeqRepo root directory".to_string(),
    })?;

    let mut total_loaded = 0;
    let mut total_skipped = 0;

    // Helper to check if a marker sequence exists in SeqRepo
    let seqrepo_has = |accession: &str, namespace: &str| -> bool {
        let output = Command::new("seqrepo")
            .args([
                "--root-directory",
                &root_dir.to_string_lossy(),
                "export",
                "--instance-name",
                instance_name,
                "--namespace",
                namespace,
                accession,
            ])
            .output();
        matches!(output, Ok(o) if o.status.success() && !o.stdout.is_empty())
    };

    // Helper to run seqrepo load
    // We try pixi run seqrepo first (to get bgzip in PATH), then fall back to plain seqrepo
    let load_fasta = |files: &[PathBuf], namespace: &str| -> Result<usize, FerroError> {
        if files.is_empty() {
            return Ok(0);
        }

        let seqrepo_args = vec![
            "seqrepo".to_string(),
            "--root-directory".to_string(),
            root_dir.to_string_lossy().to_string(),
            "load".to_string(),
            "--instance-name".to_string(),
            instance_name.to_string(),
            "--namespace".to_string(),
            namespace.to_string(),
        ];

        let mut all_args = seqrepo_args.clone();
        for file in files {
            all_args.push(file.to_string_lossy().to_string());
        }

        // Try pixi run first (provides bgzip in PATH)
        let output = Command::new("pixi").arg("run").args(&all_args).output();

        let output = match output {
            Ok(o) if o.status.success() => o,
            _ => {
                // Fall back to plain seqrepo
                let plain_args: Vec<String> = all_args[1..].to_vec(); // Skip "seqrepo"
                Command::new("seqrepo")
                    .args(&plain_args)
                    .output()
                    .map_err(|e| FerroError::Io {
                        msg: format!("Failed to run seqrepo load: {}", e),
                    })?
            }
        };

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            // Some warnings are OK, check for fatal errors
            if stderr.contains("FATAL") || stderr.contains("RuntimeError") {
                return Err(FerroError::Io {
                    msg: format!("seqrepo load failed: {}", stderr),
                });
            }
        }

        Ok(files.len())
    };

    // 0. Load all transcripts if requested (makes biocommons fully self-contained)
    if load_transcripts {
        let transcripts_dir = ferro_ref.join("transcripts");
        if transcripts_dir.exists() {
            // Check if transcripts are already loaded using a marker sequence
            // NR_072979.2 is an NR_ transcript that was missing in previous runs
            if seqrepo_has("NR_072979.2", "NCBI") {
                eprintln!("Transcripts already loaded in SeqRepo (skipping)");
                total_skipped += 1;
            } else {
                let transcript_files: Vec<PathBuf> = std::fs::read_dir(&transcripts_dir)
                    .map_err(|e| FerroError::Io {
                        msg: format!("Failed to read transcripts directory: {}", e),
                    })?
                    .filter_map(|e| e.ok())
                    .map(|e| e.path())
                    .filter(|p| {
                        p.extension()
                            .map(|ext| ext == "fna" || ext == "fasta")
                            .unwrap_or(false)
                    })
                    .collect();

                if !transcript_files.is_empty() {
                    eprintln!(
                        "Loading {} transcript files (~270K sequences, this may take a few minutes)...",
                        transcript_files.len()
                    );
                    // Load in batches to show progress
                    for (i, chunk) in transcript_files.chunks(3).enumerate() {
                        eprint!(
                            "  Batch {}/{}...",
                            i + 1,
                            transcript_files.len().div_ceil(3)
                        );
                        load_fasta(chunk, "NCBI")?;
                        eprintln!(" done");
                    }
                    total_loaded += transcript_files.len();
                }
            }
        }
    }

    // 1. Load GRCh38 genome
    let grch38_path = ferro_ref.join("genome/GRCh38.fna");
    if grch38_path.exists() {
        // Check if GRCh38 already loaded using chr1 as marker
        if seqrepo_has("NC_000001.11", "NCBI") {
            eprintln!("GRCh38 genome already loaded in SeqRepo (skipping)");
            total_skipped += 1;
        } else {
            eprintln!("Loading GRCh38 genome sequences...");
            total_loaded += load_fasta(&[grch38_path], "NCBI")?;
        }
    }

    // 2. Load GRCh37 genome
    let grch37_path = ferro_ref.join("genome/GRCh37.fna");
    if grch37_path.exists() {
        // Check if GRCh37 already loaded using chr1 as marker
        if seqrepo_has("NC_000001.10", "NCBI") {
            eprintln!("GRCh37 genome already loaded in SeqRepo (skipping)");
            total_skipped += 1;
        } else {
            eprintln!("Loading GRCh37 genome sequences...");
            total_loaded += load_fasta(&[grch37_path], "NCBI")?;
        }
    }

    // 3. Load RefSeqGene files
    let refseqgene_dir = ferro_ref.join("refseqgene");
    if refseqgene_dir.exists() {
        // Check if RefSeqGene already loaded using NG_007400.1 (BRCA1) as marker
        if seqrepo_has("NG_007400.1", "NCBI") {
            eprintln!("RefSeqGene files already loaded in SeqRepo (skipping)");
            total_skipped += 1;
        } else {
            let refseqgene_files: Vec<PathBuf> = std::fs::read_dir(&refseqgene_dir)
                .map_err(|e| FerroError::Io {
                    msg: format!("Failed to read refseqgene directory: {}", e),
                })?
                .filter_map(|e| e.ok())
                .map(|e| e.path())
                .filter(|p| {
                    p.extension()
                        .map(|ext| ext == "fna" || ext == "fasta")
                        .unwrap_or(false)
                })
                .collect();

            if !refseqgene_files.is_empty() {
                eprintln!("Loading {} RefSeqGene files...", refseqgene_files.len());
                total_loaded += load_fasta(&refseqgene_files, "NCBI")?;
            }
        }
    }

    // 4. Load LRG files
    let lrg_dir = ferro_ref.join("lrg");
    let mut lrg_loaded = false;
    if lrg_dir.exists() {
        // Check if LRG already loaded using LRG_1 as marker
        if seqrepo_has("LRG_1", "LRG") {
            eprintln!("LRG files already loaded in SeqRepo (skipping)");
            total_skipped += 1;
        } else {
            let lrg_files: Vec<PathBuf> = std::fs::read_dir(&lrg_dir)
                .map_err(|e| FerroError::Io {
                    msg: format!("Failed to read lrg directory: {}", e),
                })?
                .filter_map(|e| e.ok())
                .map(|e| e.path())
                .filter(|p| {
                    p.extension()
                        .map(|ext| ext == "fna" || ext == "fasta")
                        .unwrap_or(false)
                })
                .collect();

            if !lrg_files.is_empty() {
                eprintln!("Loading {} LRG files...", lrg_files.len());
                // Load in batches to avoid command line length limits
                for (i, chunk) in lrg_files.chunks(100).enumerate() {
                    eprint!("  Batch {}/{}...", i + 1, lrg_files.len().div_ceil(100));
                    load_fasta(chunk, "LRG")?;
                    eprintln!(" done");
                }
                total_loaded += lrg_files.len();
                lrg_loaded = true;
            }
        }
    }

    if total_skipped > 0 {
        eprintln!(
            "Loaded {} sequence files into SeqRepo ({} categories skipped - already present)",
            total_loaded, total_skipped
        );
    } else {
        eprintln!("Loaded {} sequence files into SeqRepo", total_loaded);
    }

    // 5. Add LRG genomic aliases (LRG_Xg -> LRG_X)
    // LRG FASTA files use LRG_1g for genomic sequences, but HGVS uses LRG_1:g.
    // We need to add aliases so SeqRepo can find LRG_1 when biocommons queries for it.
    if lrg_loaded {
        eprintln!("Adding LRG genomic aliases...");
        add_lrg_genomic_aliases(seqrepo_dir)?;
    }

    // 6. Load supplemental transcripts (missing ClinVar transcripts fetched from NCBI)
    let supplemental_fasta = ferro_ref.join("supplemental/patterns_transcripts.fna");
    if supplemental_fasta.exists() {
        eprintln!(
            "Loading supplemental transcripts from {}...",
            supplemental_fasta.display()
        );
        let supplemental_files = vec![supplemental_fasta.clone()];
        match load_fasta(&supplemental_files, "NCBI") {
            Ok(n) => {
                eprintln!("  Loaded {} supplemental transcript file(s)", n);
            }
            Err(e) => {
                eprintln!("  Warning: Failed to load supplemental transcripts: {}", e);
            }
        }
    }

    // 7. Load legacy transcript versions (older RefSeq versions)
    let legacy_transcripts_fasta = ferro_ref.join("supplemental/legacy_transcripts.fna");
    if legacy_transcripts_fasta.exists() {
        eprintln!("Loading legacy transcript versions...");
        let legacy_files = vec![legacy_transcripts_fasta];
        match load_fasta(&legacy_files, "NCBI") {
            Ok(n) => {
                eprintln!("  Loaded {} legacy transcript file(s)", n);
            }
            Err(e) => {
                eprintln!("  Warning: Failed to load legacy transcripts: {}", e);
            }
        }
    }

    // 8. Load legacy GenBank sequences (non-RefSeq historical accessions)
    let legacy_genbank_fasta = ferro_ref.join("supplemental/legacy_genbank.fna");
    if legacy_genbank_fasta.exists() {
        eprintln!("Loading legacy GenBank sequences...");
        let genbank_files = vec![legacy_genbank_fasta];
        match load_fasta(&genbank_files, "NCBI") {
            Ok(n) => {
                eprintln!("  Loaded {} legacy GenBank file(s)", n);
            }
            Err(e) => {
                eprintln!("  Warning: Failed to load legacy GenBank sequences: {}", e);
            }
        }
    }

    // 9. Derive protein sequences from coding transcripts
    eprintln!("Deriving protein sequences from transcripts...");

    // Load cdot for CDS coordinates
    let cdot_path = ferro_ref.join("cdot");
    let cdot: Option<CdotMapper> = if cdot_path.exists() {
        // Find the cdot JSON file
        let cdot_files: Vec<_> = std::fs::read_dir(&cdot_path)
            .ok()
            .map(|entries| {
                entries
                    .filter_map(|e| e.ok())
                    .map(|e| e.path())
                    .filter(|p| p.extension().map(|e| e == "json").unwrap_or(false))
                    .collect()
            })
            .unwrap_or_default();

        if let Some(cdot_file) = cdot_files.first() {
            eprintln!("  Loading cdot from {}...", cdot_file.display());
            CdotMapper::from_json_file(cdot_file).ok()
        } else {
            None
        }
    } else {
        None
    };

    // Load supplemental metadata for CDS coordinates
    let supplemental_metadata_path =
        ferro_ref.join("supplemental/patterns_transcripts.metadata.json");
    let supplemental_cds: HashMap<String, (Option<u64>, Option<u64>, Option<String>)> =
        if supplemental_metadata_path.exists() {
            let file = File::open(&supplemental_metadata_path).map_err(|e| FerroError::Io {
                msg: format!(
                    "Failed to open {}: {}",
                    supplemental_metadata_path.display(),
                    e
                ),
            })?;
            let metadata: SupplementalMetadata =
                serde_json::from_reader(file).map_err(|e| FerroError::Json {
                    msg: format!(
                        "Failed to parse {}: {}",
                        supplemental_metadata_path.display(),
                        e
                    ),
                })?;
            metadata
                .transcripts
                .into_iter()
                .map(|(k, v)| (k, (v.cds_start, v.cds_end, v.gene_symbol)))
                .collect()
        } else {
            HashMap::new()
        };

    // Read transcripts from base FASTA files and supplemental
    let mut transcript_sequences: HashMap<String, String> = HashMap::new();

    // Load from transcripts directory
    let transcripts_dir = ferro_ref.join("transcripts");
    if transcripts_dir.exists() {
        for entry in std::fs::read_dir(&transcripts_dir)
            .ok()
            .into_iter()
            .flatten()
            .filter_map(|e| e.ok())
        {
            let path = entry.path();
            if path.extension().map(|e| e == "fna").unwrap_or(false) {
                if let Ok(seqs) = read_fasta_to_map(&path) {
                    transcript_sequences.extend(seqs);
                }
            }
        }
    }

    // Load from supplemental
    if supplemental_fasta.exists() {
        if let Ok(seqs) = read_fasta_to_map(&supplemental_fasta) {
            transcript_sequences.extend(seqs);
        }
    }

    // Derive proteins and write to temp FASTA
    let mut proteins_derived = 0usize;
    let temp_dir = std::env::temp_dir();
    let protein_fasta_path = temp_dir.join("derived_proteins.fna");

    {
        let mut protein_file = File::create(&protein_fasta_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create temp protein FASTA: {}", e),
        })?;

        for (tx_id, seq) in &transcript_sequences {
            // Skip non-coding transcripts
            if !tx_id.starts_with("NM_") && !tx_id.starts_with("XM_") {
                continue;
            }

            // Get CDS coordinates from cdot or supplemental metadata
            let (cds_start, cds_end, protein_id) = if let Some(ref cdot) = cdot {
                if let Some(tx) = cdot.get_transcript(tx_id) {
                    let start = tx.cds_start.unwrap_or(0) as usize;
                    let end = tx.cds_end.unwrap_or(seq.len() as u64) as usize;
                    let prot_id = tx
                        .protein
                        .clone()
                        .unwrap_or_else(|| tx_id.replace("NM_", "NP_").replace("XM_", "XP_"));
                    (start, end, prot_id)
                } else if let Some((start, end, _gene)) = supplemental_cds.get(tx_id) {
                    let start = start.map(|s| (s - 1) as usize).unwrap_or(0);
                    let end = end.map(|e| e as usize).unwrap_or(seq.len());
                    let prot_id = tx_id.replace("NM_", "NP_").replace("XM_", "XP_");
                    (start, end, prot_id)
                } else {
                    continue;
                }
            } else if let Some((start, end, _gene)) = supplemental_cds.get(tx_id) {
                let start = start.map(|s| (s - 1) as usize).unwrap_or(0);
                let end = end.map(|e| e as usize).unwrap_or(seq.len());
                let prot_id = tx_id.replace("NM_", "NP_").replace("XM_", "XP_");
                (start, end, prot_id)
            } else {
                continue;
            };

            // Validate coordinates
            if cds_end > seq.len() || cds_start >= cds_end {
                continue;
            }

            // Derive protein sequence
            let cds = &seq[cds_start..cds_end];
            if let Some(protein_seq) = translate_cds_to_protein(cds) {
                writeln!(protein_file, ">{}", protein_id).map_err(|e| FerroError::Io {
                    msg: format!("Failed to write protein FASTA: {}", e),
                })?;
                writeln!(protein_file, "{}", protein_seq).map_err(|e| FerroError::Io {
                    msg: format!("Failed to write protein FASTA: {}", e),
                })?;
                proteins_derived += 1;
            }
        }
    }

    // Load derived proteins into SeqRepo
    if proteins_derived > 0 {
        eprintln!(
            "  Loading {} derived protein sequences...",
            proteins_derived
        );
        match load_fasta(std::slice::from_ref(&protein_fasta_path), "NCBI") {
            Ok(_) => {
                eprintln!("  Loaded {} derived protein sequences", proteins_derived);
            }
            Err(e) => {
                eprintln!("  Warning: Failed to load derived proteins: {}", e);
            }
        }
        // Clean up temp file
        let _ = std::fs::remove_file(&protein_fasta_path);
    } else {
        eprintln!("  No protein sequences derived (no CDS coordinates available)");
    }

    Ok(())
}

/// Fetch missing accessions from a patterns file and load into SeqRepo.
///
/// This extracts accessions from the patterns file, fetches any that are missing
/// from SeqRepo via NCBI E-utilities, and loads them into SeqRepo.
///
/// # Arguments
/// * `patterns_path` - Path to file containing HGVS patterns (one per line)
/// * `seqrepo_dir` - Path to SeqRepo instance directory
/// * `ferro_reference` - Optional path to ferro reference directory (sequences here won't need fetching)
///
/// # Returns
/// The number of accessions fetched and loaded.
pub fn fetch_and_load_missing_accessions(
    patterns_path: &Path,
    seqrepo_dir: &Path,
    ferro_reference: Option<&Path>,
) -> Result<usize, FerroError> {
    use crate::benchmark::accessions::AccessionSources;
    use crate::benchmark::cache::extract_all_accessions_from_file;
    use tempfile::NamedTempFile;

    // 1. Extract accessions from patterns
    let accessions = extract_all_accessions_from_file(patterns_path)?;

    // Filter to nucleotide accessions (NG_, NM_, NR_, NC_)
    let nucleotide: Vec<String> = accessions
        .into_iter()
        .filter(|a| {
            a.starts_with("NG_")
                || a.starts_with("NM_")
                || a.starts_with("NR_")
                || a.starts_with("NC_")
        })
        .collect();

    if nucleotide.is_empty() {
        eprintln!("No nucleotide accessions found in patterns");
        return Ok(0);
    }
    eprintln!(
        "Found {} unique nucleotide accessions in patterns",
        nucleotide.len()
    );

    // Build set of accessions from ferro reference that are LOADED into SeqRepo
    // Note: We only check supplemental, genomic, refseqgene, and lrg directories
    // because those are what load_sequences_to_seqrepo() actually loads.
    // The transcripts/ directory is NOT loaded (SeqRepo already has most from the snapshot).
    let ferro_accessions = if let Some(ferro_ref) = ferro_reference {
        match AccessionSources::from_ferro_reference(ferro_ref) {
            Ok(sources) => {
                // Count only what's actually loaded to SeqRepo
                let loaded_count = sources.supplemental.len()
                    + sources.genomic.len()
                    + sources.refseqgene.len()
                    + sources.lrg.len();
                if loaded_count > 0 {
                    eprintln!(
                        "  Found {} accessions in ferro reference (loaded to SeqRepo)",
                        loaded_count
                    );
                }
                sources
            }
            Err(_) => AccessionSources::default(),
        }
    } else {
        AccessionSources::default()
    };

    // 2. Check which are missing from both SeqRepo AND ferro reference
    // First filter out accessions from ferro directories loaded to SeqRepo
    let need_seqrepo_check: Vec<_> = nucleotide
        .iter()
        .filter(|acc| {
            // Skip if in ferro directories that are loaded to SeqRepo
            // (supplemental, genomic, refseqgene, lrg - NOT transcripts)
            !ferro_accessions.supplemental.contains(*acc)
                && !ferro_accessions.genomic.contains(*acc)
                && !ferro_accessions.refseqgene.contains(*acc)
                && !ferro_accessions.lrg.contains(*acc)
        })
        .cloned()
        .collect();

    // Use batch SQLite query to check SeqRepo (much faster than spawning CLI per accession)
    let seqrepo_existing = check_seqrepo_has_accessions_batch(seqrepo_dir, &need_seqrepo_check)?;
    let missing: Vec<_> = need_seqrepo_check
        .into_iter()
        .filter(|acc| !seqrepo_existing.contains(acc))
        .collect();

    if missing.is_empty() {
        eprintln!("All accessions available (ferro + SeqRepo)");
        return Ok(0);
    }
    eprintln!("Found {} missing accessions to fetch", missing.len());

    // 3. Fetch from NCBI to temp file
    let temp_fasta = NamedTempFile::new().map_err(|e| FerroError::Io {
        msg: format!("Failed to create temp file: {}", e),
    })?;
    let temp_path = temp_fasta.path();

    let fetched = fetch_nucleotide_to_fasta(&missing, temp_path)?;

    if fetched == 0 {
        return Ok(0);
    }

    // 4. Load into SeqRepo
    eprintln!("Loading {} fetched accessions into SeqRepo...", fetched);
    load_fasta_to_seqrepo(seqrepo_dir, temp_path, "NCBI")?;

    Ok(fetched)
}

/// Check which accessions exist in SeqRepo using direct SQLite query (batch)
/// Returns the set of accessions that ARE present in SeqRepo
fn check_seqrepo_has_accessions_batch(
    seqrepo_dir: &Path,
    accessions: &[String],
) -> Result<std::collections::HashSet<String>, FerroError> {
    use rusqlite::Connection;
    use std::collections::HashSet;

    if accessions.is_empty() {
        return Ok(HashSet::new());
    }

    let db_path = seqrepo_dir.join("aliases.sqlite3");
    if !db_path.exists() {
        return Err(FerroError::Io {
            msg: format!("SeqRepo aliases database not found: {}", db_path.display()),
        });
    }

    let conn = Connection::open(&db_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open SeqRepo database: {}", e),
    })?;

    // Query all existing aliases in batches to avoid SQLite parameter limits
    // We query once for all aliases and filter in Rust
    eprintln!(
        "  Checking {} accessions against SeqRepo...",
        accessions.len()
    );

    let mut existing = HashSet::new();

    // SQLite has a limit of ~1000 parameters in IN clauses, so we batch
    const BATCH_SIZE: usize = 500;
    let total_batches = accessions.len().div_ceil(BATCH_SIZE);

    for (i, batch) in accessions.chunks(BATCH_SIZE).enumerate() {
        // Build parameterized query
        let placeholders: Vec<_> = (1..=batch.len()).map(|i| format!("?{}", i)).collect();
        let query = format!(
            "SELECT DISTINCT alias FROM seqalias WHERE alias IN ({}) AND is_current = 1",
            placeholders.join(", ")
        );

        let mut stmt = conn.prepare(&query).map_err(|e| FerroError::Io {
            msg: format!("Failed to prepare SeqRepo query: {}", e),
        })?;

        // Bind parameters
        let params: Vec<&dyn rusqlite::ToSql> =
            batch.iter().map(|s| s as &dyn rusqlite::ToSql).collect();
        let rows = stmt
            .query_map(params.as_slice(), |row| row.get::<_, String>(0))
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to query SeqRepo: {}", e),
            })?;

        for alias in rows.flatten() {
            existing.insert(alias);
        }

        if (i + 1) % 100 == 0 || i + 1 == total_batches {
            eprintln!(
                "  Checked {}/{} batches ({} accessions found in SeqRepo)...",
                i + 1,
                total_batches,
                existing.len()
            );
        }
    }

    eprintln!(
        "  Found {} of {} accessions in SeqRepo",
        existing.len(),
        accessions.len()
    );
    Ok(existing)
}

/// Fetch nucleotide accessions from NCBI and write to FASTA file
fn fetch_nucleotide_to_fasta(accessions: &[String], output: &Path) -> Result<usize, FerroError> {
    use std::io::Write;
    use std::process::Command;

    const BATCH_SIZE: usize = 200;

    if accessions.is_empty() {
        return Ok(0);
    }

    let mut fetched = 0;
    let mut file = std::fs::File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create output file: {}", e),
    })?;

    let total_batches = accessions.len().div_ceil(BATCH_SIZE);
    eprintln!(
        "  Fetching {} accessions in {} batches...",
        accessions.len(),
        total_batches
    );

    // Batch fetch from NCBI
    for (batch_idx, batch) in accessions.chunks(BATCH_SIZE).enumerate() {
        let ids = batch.join(",");
        let url = format!(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text",
            ids
        );

        // Add API key if available
        let url = if let Ok(api_key) = std::env::var("NCBI_API_KEY") {
            format!("{}&api_key={}", url, api_key)
        } else {
            url
        };

        let output_result = Command::new("curl").args(["-s", "-L", &url]).output();

        match output_result {
            Ok(o) if o.status.success() => {
                let text = String::from_utf8_lossy(&o.stdout);
                if !text.is_empty() && text.starts_with('>') {
                    file.write_all(text.as_bytes())
                        .map_err(|e| FerroError::Io {
                            msg: format!("Failed to write to file: {}", e),
                        })?;
                    // Count sequences by counting '>' at start of lines
                    let count = text.lines().filter(|l| l.starts_with('>')).count();
                    fetched += count;
                }
            }
            _ => {
                eprintln!(
                    "  Warning: Failed to fetch batch {}/{}",
                    batch_idx + 1,
                    total_batches
                );
            }
        }

        // Rate limit
        std::thread::sleep(std::time::Duration::from_millis(350));

        eprintln!(
            "  Batch {}/{}: fetched {} sequences so far",
            batch_idx + 1,
            total_batches,
            fetched
        );
    }

    Ok(fetched)
}

/// Load a FASTA file into SeqRepo
fn load_fasta_to_seqrepo(
    seqrepo_dir: &Path,
    fasta_path: &Path,
    namespace: &str,
) -> Result<(), FerroError> {
    use std::process::Command;

    // Get root and instance from seqrepo path
    let root_dir = seqrepo_dir.parent().ok_or_else(|| FerroError::Io {
        msg: "SeqRepo directory has no parent".to_string(),
    })?;
    let instance_name = seqrepo_dir
        .file_name()
        .and_then(|n| n.to_str())
        .ok_or_else(|| FerroError::Io {
            msg: "SeqRepo directory has no name".to_string(),
        })?;

    // Try pixi run first (provides bgzip in PATH), then fall back to plain seqrepo
    let output = Command::new("pixi")
        .args([
            "run",
            "seqrepo",
            "--root-directory",
            &root_dir.to_string_lossy(),
            "load",
            "--instance-name",
            instance_name,
            "--namespace",
            namespace,
            &fasta_path.to_string_lossy(),
        ])
        .output();

    let output = match output {
        Ok(o) if o.status.success() => o,
        _ => {
            // Fall back to plain seqrepo
            Command::new("seqrepo")
                .args([
                    "--root-directory",
                    &root_dir.to_string_lossy(),
                    "load",
                    "--instance-name",
                    instance_name,
                    "--namespace",
                    namespace,
                    &fasta_path.to_string_lossy(),
                ])
                .output()
                .map_err(|e| FerroError::Io {
                    msg: format!("Failed to run seqrepo load: {}", e),
                })?
        }
    };

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        // Some warnings are OK, check for fatal errors
        if stderr.contains("FATAL") || stderr.contains("RuntimeError") {
            return Err(FerroError::Io {
                msg: format!("seqrepo load failed: {}", stderr),
            });
        }
    }

    Ok(())
}

/// Set up SeqRepo by downloading sequence data.
///
/// This uses the `seqrepo` Python CLI to download sequence data from biocommons mirrors.
///
/// # Arguments
/// * `config` - Configuration specifying the SeqRepo directory and instance
/// * `force` - If true, re-download even if data exists
///
/// # Returns
/// A result with details about the setup operation.
pub fn setup_seqrepo(
    config: &BiocommonsLocalConfig,
    force: bool,
) -> Result<SeqRepoSetupResult, FerroError> {
    // Check if seqrepo CLI is available
    let seqrepo_available = Command::new("seqrepo")
        .args(["--version"])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);

    if !seqrepo_available {
        return Err(FerroError::Io {
            msg: "seqrepo CLI not found. Install with: pip install biocommons.seqrepo".to_string(),
        });
    }

    // Create the parent directory if it doesn't exist
    if let Some(parent) = config.seqrepo_dir.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory: {}", e),
        })?;
    }

    // Check if SeqRepo data already exists
    let data_exists = check_seqrepo(&config.seqrepo_dir);

    if data_exists && !force {
        eprintln!(
            "SeqRepo data already exists at: {}",
            config.seqrepo_dir.display()
        );
        return Ok(SeqRepoSetupResult {
            seqrepo_dir: config.seqrepo_dir.clone(),
            instance: config.seqrepo_instance.clone(),
            downloaded: false,
        });
    }

    // Download SeqRepo data using the seqrepo CLI
    eprintln!(
        "Downloading SeqRepo instance {} to {} (this will take a while, ~10GB)...",
        config.seqrepo_instance,
        config.seqrepo_dir.display()
    );

    // Find rsync executable - prefer one in PATH that's not openrsync
    let rsync_exe = find_gnu_rsync();

    let mut args = vec![
        "--root-directory".to_string(),
        config.seqrepo_dir.to_string_lossy().to_string(),
    ];

    if let Some(ref rsync_path) = rsync_exe {
        args.push("--rsync-exe".to_string());
        args.push(rsync_path.clone());
    }

    args.extend([
        "pull".to_string(),
        "-i".to_string(),
        config.seqrepo_instance.clone(),
    ]);

    let pull_output = Command::new("seqrepo")
        .args(&args)
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to run seqrepo pull: {}", e),
        })?;

    if !pull_output.status.success() {
        let stderr = String::from_utf8_lossy(&pull_output.stderr);
        return Err(FerroError::Io {
            msg: format!("Failed to download SeqRepo data: {}", stderr),
        });
    }

    eprintln!("SeqRepo data downloaded successfully.");

    Ok(SeqRepoSetupResult {
        seqrepo_dir: config.seqrepo_dir.clone(),
        instance: config.seqrepo_instance.clone(),
        downloaded: true,
    })
}

// ============================================================================
// Container Management
// ============================================================================

/// Start an existing UTA Docker container.
pub fn start_uta(container_name: &str) -> Result<(), FerroError> {
    let output = Command::new("docker")
        .args(["start", container_name])
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to start container: {}", e),
        })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(FerroError::Io {
            msg: format!("Failed to start container: {}", stderr),
        });
    }

    Ok(())
}

/// Stop a running UTA Docker container.
pub fn stop_uta(container_name: &str) -> Result<(), FerroError> {
    let output = Command::new("docker")
        .args(["stop", container_name])
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to stop container: {}", e),
        })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(FerroError::Io {
            msg: format!("Failed to stop container: {}", stderr),
        });
    }

    Ok(())
}

// ============================================================================
// Settings File Management
// ============================================================================

/// Write a biocommons settings file.
///
/// Creates a configuration file in key=value format that can be used
/// to configure biocommons/hgvs for local operation.
///
/// # Format
/// ```text
/// UTA_DB_URL = postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b
/// HGVS_SEQREPO_DIR = /path/to/seqrepo/2021-01-29
/// ```
pub fn write_biocommons_settings(
    settings_path: &Path,
    uta_db_url: &str,
    seqrepo_dir: &Path,
) -> Result<(), FerroError> {
    let mut file = std::fs::File::create(settings_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create settings file: {}", e),
    })?;

    writeln!(file, "# biocommons/hgvs local configuration").map_err(|e| FerroError::Io {
        msg: format!("Failed to write settings: {}", e),
    })?;
    writeln!(file, "# Generated by ferro-benchmark").map_err(|e| FerroError::Io {
        msg: format!("Failed to write settings: {}", e),
    })?;
    writeln!(file).map_err(|e| FerroError::Io {
        msg: format!("Failed to write settings: {}", e),
    })?;
    writeln!(file, "UTA_DB_URL = {}", uta_db_url).map_err(|e| FerroError::Io {
        msg: format!("Failed to write settings: {}", e),
    })?;
    writeln!(file, "HGVS_SEQREPO_DIR = {}", seqrepo_dir.display()).map_err(|e| FerroError::Io {
        msg: format!("Failed to write settings: {}", e),
    })?;

    Ok(())
}

/// Load biocommons settings from a configuration file.
///
/// Reads a key=value format configuration file and returns the settings.
pub fn load_biocommons_settings(settings_path: &Path) -> Result<BiocommonsSettings, FerroError> {
    let file = std::fs::File::open(settings_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open settings file: {}", e),
    })?;

    let reader = std::io::BufReader::new(file);
    let mut uta_db_url: Option<String> = None;
    let mut seqrepo_dir: Option<PathBuf> = None;

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read settings: {}", e),
        })?;

        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        if let Some((key, value)) = line.split_once('=') {
            let key = key.trim();
            let value = value.trim();

            match key {
                "UTA_DB_URL" => uta_db_url = Some(value.to_string()),
                "HGVS_SEQREPO_DIR" => seqrepo_dir = Some(PathBuf::from(value)),
                _ => {} // Ignore unknown keys
            }
        }
    }

    let uta_db_url = uta_db_url.ok_or_else(|| FerroError::Io {
        msg: "UTA_DB_URL not found in settings file".to_string(),
    })?;

    let seqrepo_dir = seqrepo_dir.ok_or_else(|| FerroError::Io {
        msg: "HGVS_SEQREPO_DIR not found in settings file".to_string(),
    })?;

    Ok(BiocommonsSettings {
        uta_db_url,
        seqrepo_dir,
    })
}

// ============================================================================
// Existing Functions (Remote Integration)
// ============================================================================

/// Check if the biocommons/hgvs Python package is available.
pub fn has_biocommons_normalizer() -> bool {
    Command::new("python3")
        .args([
            "-c",
            "import hgvs.parser; import hgvs.normalizer; import hgvs.dataproviders.uta",
        ])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Check if we can connect to the UTA database (remote or local).
///
/// This performs a lightweight connection test without loading any data.
pub fn check_uta_connection(uta_db_url: Option<&str>) -> bool {
    let check_script = if let Some(url) = uta_db_url {
        format!(
            r#"
import os
os.environ['UTA_DB_URL'] = '{}'
import hgvs.dataproviders.uta
hdp = hgvs.dataproviders.uta.connect()
print("OK")
"#,
            url
        )
    } else {
        r#"
import hgvs.dataproviders.uta
hdp = hgvs.dataproviders.uta.connect()
print("OK")
"#
        .to_string()
    };

    Command::new("python3")
        .args(["-c", &check_script])
        .output()
        .map(|o| o.status.success() && String::from_utf8_lossy(&o.stdout).contains("OK"))
        .unwrap_or(false)
}

/// Run biocommons/hgvs normalization via Python subprocess.
///
/// This calls the biocommons/hgvs package directly for normalization.
/// By default, it uses the remote UTA database at uta.biocommons.org.
/// Set `uta_db_url` to use a local UTA database for better performance.
/// Set `lrg_mapping_file` to enable LRG transcript pattern translation.
pub fn run_biocommons_normalizer_subprocess(
    input_file: &str,
    output_file: &str,
    uta_db_url: Option<&str>,
    seqrepo_dir: Option<&str>,
    lrg_mapping_file: Option<&str>,
) -> Result<(), FerroError> {
    // Python script for biocommons normalization with thread-based parallelism
    // Threads work well here since normalization is I/O bound (database queries)
    let python_code = r#"
import sys
import os
import json
import time
import re
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

# Optional: Set UTA database URL from environment or argument
if len(sys.argv) > 3 and sys.argv[3]:
    os.environ['UTA_DB_URL'] = sys.argv[3]

# Optional: Set SeqRepo directory from argument
if len(sys.argv) > 4 and sys.argv[4]:
    os.environ['HGVS_SEQREPO_DIR'] = sys.argv[4]

# Optional: Load LRG-to-RefSeq mapping file
# Maps LRG transcript accessions (e.g., LRG_199t1) to RefSeq (e.g., NM_004006.2)
LRG_MAPPING = {}
if len(sys.argv) > 5 and sys.argv[5]:
    try:
        with open(sys.argv[5], 'r') as f:
            for line in f:
                # Skip header and comment lines
                if line.startswith('#') or line.startswith('LRG\t'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    lrg, _, _, tx, refseq = parts[:5]
                    if refseq and refseq != '-':
                        # Map "LRG_199t1" -> "NM_004006.2"
                        LRG_MAPPING[f"{lrg}{tx}"] = refseq
    except Exception as e:
        print(f"Warning: Could not load LRG mapping file: {e}", file=sys.stderr)

def translate_lrg_pattern(pattern):
    """Translate LRG transcript pattern to RefSeq if mapping exists.

    Returns: (translated_pattern, original_lrg_accession or None)
    """
    # Match LRG transcript patterns like LRG_199t1:c.79del
    match = re.match(r'(LRG_\d+t\d+)(:.+)', pattern)
    if match:
        lrg_acc, rest = match.groups()
        if lrg_acc in LRG_MAPPING:
            return LRG_MAPPING[lrg_acc] + rest, lrg_acc
    return pattern, None

# Number of parallel workers - threads work well for I/O bound work
NUM_WORKERS = 8

try:
    import hgvs.parser
    import hgvs.normalizer
    import hgvs.dataproviders.uta
    import hgvs.variantmapper
except ImportError:
    print("ERROR: hgvs (biocommons) not installed. Run: pip install hgvs", file=sys.stderr)
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

def categorize_error(error_msg):
    """Categorize error message into a short category."""
    msg = str(error_msg)
    if 'HGVSParseError' in msg or 'parse' in msg.lower():
        return 'PARSE_ERROR'
    elif 'HGVSInvalidVariantError' in msg:
        return 'INVALID_VARIANT'
    elif 'HGVSDataNotAvailableError' in msg or 'not found' in msg.lower():
        return 'DATA_NOT_AVAILABLE'
    elif 'HGVSInvalidIntervalError' in msg:
        return 'INVALID_INTERVAL'
    elif 'HGVSNormalizationError' in msg:
        return 'NORMALIZATION_ERROR'
    elif 'HGVSUnsupportedOperationError' in msg:
        return 'UNSUPPORTED'
    elif 'connection' in msg.lower() or 'timeout' in msg.lower():
        return 'CONNECTION_ERROR'
    else:
        return 'OTHER'

# Thread-local storage for per-thread hgvs components
_thread_local = threading.local()

def get_normalizer():
    """Get or create thread-local hgvs normalizer and variant mapper."""
    if not hasattr(_thread_local, 'hn'):
        hdp = hgvs.dataproviders.uta.connect()
        _thread_local.hdp = hdp
        _thread_local.hp = hgvs.parser.Parser()
        _thread_local.hn = hgvs.normalizer.Normalizer(hdp)
        # Create variant mapper for reference replacement
        _thread_local.vm = hgvs.variantmapper.VariantMapper(hdp, replace_reference=True)
    return _thread_local.hp, _thread_local.hn, _thread_local.vm

def normalize_pattern(pattern):
    """Normalize a single pattern using thread-local components."""
    # Try to translate LRG transcript patterns to RefSeq
    translated_pattern, original_lrg = translate_lrg_pattern(pattern)

    try:
        hp, hn, vm = get_normalizer()
        var = hp.parse_hgvs_variant(translated_pattern)
        # Replace incorrect reference bases with actual sequence (skip for protein variants)
        if var.type != 'p':
            var = vm._replace_reference(var)
        normalized = hn.normalize(var)
        result = {
            "input": pattern,
            "success": True,
            "output": str(normalized)
        }
        # Include translation info if pattern was translated
        if original_lrg:
            result["translated_from"] = original_lrg
            result["translated_to"] = translated_pattern
        return result
    except Exception as e:
        error_msg = str(e)
        category = categorize_error(error_msg)
        result = {
            "input": pattern,
            "success": False,
            "error": error_msg[:200],
            "error_category": category
        }
        # Include translation info even on failure
        if original_lrg:
            result["translated_from"] = original_lrg
            result["translated_to"] = translated_pattern
        return result

with open(input_file, 'r') as f:
    patterns = [line.strip() for line in f if line.strip()]

start = time.perf_counter()

# Use thread pool for parallel I/O
if len(patterns) > 5:
    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
        results = list(executor.map(normalize_pattern, patterns))
else:
    # Sequential for tiny batches
    results = [normalize_pattern(p) for p in patterns]

elapsed = time.perf_counter() - start

successful = sum(1 for r in results if r['success'])
failed = len(results) - successful
error_counts = defaultdict(int)
for r in results:
    if not r['success'] and 'error_category' in r:
        error_counts[r['error_category']] += 1

output_data = {
    "tool": "biocommons-hgvs",
    "total_patterns": len(patterns),
    "successful": successful,
    "failed": failed,
    "elapsed_seconds": elapsed,
    "patterns_per_second": len(patterns) / elapsed if elapsed > 0 else 0,
    "error_counts": dict(error_counts),
    "results": results
}

with open(output_file, 'w') as f:
    json.dump(output_data, f, indent=2)
"#;

    // Try pixi run first to use pixi-managed Python environment,
    // then fall back to plain python3 if pixi is not available
    let (program, base_args): (&str, Vec<&str>) = if Command::new("pixi")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
    {
        (
            "pixi",
            vec!["run", "python3", "-c", python_code, input_file, output_file],
        )
    } else {
        ("python3", vec!["-c", python_code, input_file, output_file])
    };

    let mut cmd = Command::new(program);
    cmd.args(&base_args);

    if let Some(url) = uta_db_url {
        cmd.arg(url);
    } else {
        cmd.arg("");
    }

    if let Some(dir) = seqrepo_dir {
        cmd.arg(dir);
    } else {
        cmd.arg("");
    }

    if let Some(mapping) = lrg_mapping_file {
        cmd.arg(mapping);
    } else {
        cmd.arg("");
    }

    let output = cmd.output().map_err(|e| FerroError::Io {
        msg: format!("Failed to run Python: {}", e),
    })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(FerroError::Io {
            msg: format!("biocommons/hgvs normalizer failed: {}", stderr),
        });
    }

    Ok(())
}

/// Run biocommons/hgvs normalization in parallel using multiple subprocesses.
///
/// This shards the input patterns and runs N Python subprocesses in parallel,
/// which provides true parallelism since each subprocess has its own GIL.
pub fn run_biocommons_normalizer_parallel(
    input_file: &str,
    output_file: &str,
    uta_db_url: Option<&str>,
    seqrepo_dir: Option<&str>,
    lrg_mapping_file: Option<&str>,
    num_workers: usize,
) -> Result<(), FerroError> {
    use rayon::prelude::*;
    use std::io::{BufRead, BufReader, Write};
    use tempfile::NamedTempFile;

    // Read all patterns
    let file = std::fs::File::open(input_file).map_err(|e| FerroError::Io {
        msg: format!("Failed to open input file: {}", e),
    })?;
    let reader = BufReader::new(file);
    let patterns: Vec<String> = reader
        .lines()
        .map_while(Result::ok)
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();

    if patterns.is_empty() {
        // Write empty results
        std::fs::write(
            output_file,
            r#"{"tool":"biocommons-hgvs","total_patterns":0,"successful":0,"failed":0,"elapsed_seconds":0,"patterns_per_second":0,"error_counts":{},"results":[]}"#,
        )
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to write output: {}", e),
        })?;
        return Ok(());
    }

    let start = std::time::Instant::now();

    // Shard patterns across workers
    let chunk_size = patterns.len().div_ceil(num_workers);
    let chunks: Vec<Vec<String>> = patterns.chunks(chunk_size).map(|c| c.to_vec()).collect();

    eprintln!(
        "Running {} parallel biocommons workers on {} patterns...",
        chunks.len(),
        patterns.len()
    );

    // Run subprocesses in parallel
    let results: Vec<Result<Vec<ParseResult>, FerroError>> = chunks
        .par_iter()
        .map(|chunk| {
            // Create temp files for this shard
            let mut input_temp =
                NamedTempFile::new().map_err(|e| FerroError::Io { msg: e.to_string() })?;
            let output_temp =
                NamedTempFile::new().map_err(|e| FerroError::Io { msg: e.to_string() })?;

            // Write patterns to temp file
            for pattern in chunk {
                writeln!(input_temp, "{}", pattern)
                    .map_err(|e| FerroError::Io { msg: e.to_string() })?;
            }

            // Run normalizer
            run_biocommons_normalizer_subprocess(
                input_temp.path().to_str().unwrap(),
                output_temp.path().to_str().unwrap(),
                uta_db_url,
                seqrepo_dir,
                lrg_mapping_file,
            )?;

            // Parse results
            let output_data: serde_json::Value =
                serde_json::from_reader(std::fs::File::open(output_temp.path()).map_err(|e| {
                    FerroError::Io {
                        msg: format!("Failed to read output: {}", e),
                    }
                })?)
                .map_err(|e| FerroError::Io {
                    msg: format!("Failed to parse JSON: {}", e),
                })?;

            // Extract results
            let results_array =
                output_data["results"]
                    .as_array()
                    .ok_or_else(|| FerroError::Io {
                        msg: "Missing results array".to_string(),
                    })?;

            let parsed_results: Result<Vec<ParseResult>, FerroError> = results_array
                .iter()
                .enumerate()
                .map(|(idx, r)| {
                    // "input" is required - its absence indicates malformed JSON
                    let input = r["input"]
                        .as_str()
                        .ok_or_else(|| FerroError::Io {
                            msg: format!(
                                "Malformed subprocess output: missing 'input' field at index {}",
                                idx
                            ),
                        })?
                        .to_string();
                    // "success" defaults to false if missing (safe default)
                    let success = r["success"].as_bool().unwrap_or(false);
                    Ok(ParseResult {
                        input,
                        success,
                        output: r["output"].as_str().map(String::from),
                        error: r["error"].as_str().map(String::from),
                        error_category: r["error_category"].as_str().map(String::from),
                        ref_mismatch: None,
                        details: None,
                    })
                })
                .collect();
            let parsed_results = parsed_results?;

            Ok(parsed_results)
        })
        .collect();

    let elapsed = start.elapsed();

    // Merge all results
    let mut all_results: Vec<ParseResult> = Vec::new();
    for result in results {
        all_results.extend(result?);
    }

    let successful = all_results.iter().filter(|r| r.success).count();
    let failed = all_results.len() - successful;

    // Aggregate error counts
    let mut error_counts: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    for r in &all_results {
        if let Some(ref cat) = r.error_category {
            *error_counts.entry(cat.clone()).or_insert(0) += 1;
        }
    }

    // Write merged results
    let output_data = serde_json::json!({
        "tool": "biocommons-hgvs",
        "total_patterns": all_results.len(),
        "successful": successful,
        "failed": failed,
        "elapsed_seconds": elapsed.as_secs_f64(),
        "patterns_per_second": all_results.len() as f64 / elapsed.as_secs_f64(),
        "error_counts": error_counts,
        "results": all_results,
    });

    std::fs::write(
        output_file,
        serde_json::to_string_pretty(&output_data).unwrap(),
    )
    .map_err(|e| FerroError::Io {
        msg: format!("Failed to write output: {}", e),
    })?;

    Ok(())
}

/// Normalize a single HGVS expression using biocommons/hgvs.
///
/// This is a convenience function for normalizing individual patterns.
/// For bulk normalization, use `run_biocommons_normalizer_subprocess`.
pub fn normalize_single(
    hgvs: &str,
    uta_db_url: Option<&str>,
    seqrepo_dir: Option<&str>,
    lrg_mapping_file: Option<&str>,
) -> Result<ParseResult, FerroError> {
    use std::io::Write;
    use tempfile::NamedTempFile;

    // Create temp files for input and output
    let mut input_file = NamedTempFile::new().map_err(|e| FerroError::Io {
        msg: format!("Failed to create temp file: {}", e),
    })?;
    let output_file = NamedTempFile::new().map_err(|e| FerroError::Io {
        msg: format!("Failed to create temp file: {}", e),
    })?;

    // Write input pattern
    writeln!(input_file, "{}", hgvs).map_err(|e| FerroError::Io {
        msg: format!("Failed to write to temp file: {}", e),
    })?;

    // Run normalization
    run_biocommons_normalizer_subprocess(
        input_file.path().to_str().unwrap(),
        output_file.path().to_str().unwrap(),
        uta_db_url,
        seqrepo_dir,
        lrg_mapping_file,
    )?;

    // Read output
    let output_content =
        std::fs::read_to_string(output_file.path()).map_err(|e| FerroError::Io {
            msg: format!("Failed to read output file: {}", e),
        })?;

    // Parse JSON output
    let output_data: serde_json::Value =
        serde_json::from_str(&output_content).map_err(|e| FerroError::Io {
            msg: format!("Failed to parse JSON output: {}", e),
        })?;

    // Extract result
    let results = output_data["results"]
        .as_array()
        .ok_or_else(|| FerroError::Io {
            msg: "No results in output".to_string(),
        })?;

    if let Some(result) = results.first() {
        let success = result["success"].as_bool().unwrap_or(false);
        if success {
            Ok(ParseResult {
                input: hgvs.to_string(),
                success: true,
                output: result["output"].as_str().map(|s| s.to_string()),
                error: None,
                error_category: None,
                ref_mismatch: None,
                details: None,
            })
        } else {
            Ok(ParseResult {
                input: hgvs.to_string(),
                success: false,
                output: None,
                error: result["error"].as_str().map(|s| s.to_string()),
                error_category: result["error_category"].as_str().map(|s| s.to_string()),
                ref_mismatch: None,
                details: None,
            })
        }
    } else {
        Err(FerroError::Io {
            msg: "No results returned from biocommons/hgvs".to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_biocommons_available() {
        // This test just verifies the function doesn't panic
        let _ = has_biocommons_normalizer();
    }
}
