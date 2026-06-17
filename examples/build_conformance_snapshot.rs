//! Build the pinned, version-exact **transcript** reference snapshot for a
//! conformance corpus (#719, increment I2).
//!
//! For every transcript accession the corpus references (via the corpus
//! accession inventory, #707), harvest the transcript's bases at the **exact
//! pinned version** plus that version's **own** CDS bounds from a prepared
//! reference, and write the two committed snapshot files
//! (`transcripts.fna` + `transcripts.metadata.json`, see
//! [`conformance::reference_snapshot`]).
//!
//! Why this is correct, not a no-op recorder:
//! - **Bases** come only from an *exact-version* FASTA hit
//!   ([`MultiFastaProvider::contains_exact_sequence`]); the version-fuzzy
//!   `get_sequence` fallback can substitute a *sibling* version's bases, so we
//!   gate on the exact key first and skip (report) any accession whose exact
//!   bases are absent locally.
//! - **CDS** comes from [`CdotMapper::get_transcript_exact_any_build`] — the
//!   accession's OWN reading frame, sourced cross-build if needed (e.g.
//!   `NM_003002.2`'s CDS lives in the GRCh37 cdot), NEVER a fuzzy sibling
//!   version's frame (#714/#717/#720). CDS bounds in transcript coordinates are
//!   genome-build-independent, so a cross-build exact hit is the right frame.
//!
//! The snapshot is committed and reviewed (like the biocommons windows); it is
//! **not** regenerated in CI, because the build requires a prepared reference
//! manifest CI does not have. `--check` is a *local* guard: run against the same
//! manifest, it fails if the committed snapshot no longer matches what the
//! current builder + reference produce (i.e. the pinned bytes drifted).
//!
//! Run (write):
//!   cargo run --features dev --example build_conformance_snapshot -- \
//!     --manifest /path/to/prepared/reference/manifest.json
//! Check (local, against the same manifest):
//!   cargo run --features dev --example build_conformance_snapshot -- \
//!     --manifest /path/to/prepared/reference/manifest.json --check

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use sha2::{Digest, Sha256};

use ferro_hgvs::conformance::accession_inventory::{AccessionClass, Inventory};
use ferro_hgvs::conformance::mutalyzer::Fixture;
use ferro_hgvs::conformance::reference_snapshot::{
    render_fasta, Provenance, TranscriptEntry, TranscriptSnapshot, TRANSCRIPTS_FASTA,
    TRANSCRIPTS_METADATA,
};
use ferro_hgvs::reference::multi_fasta::MultiFastaProvider;
use ferro_hgvs::reference::ReferenceProvider;

#[derive(Parser, Debug)]
#[command(
    about = "Harvest the version-exact transcript reference snapshot for a conformance corpus"
)]
struct Cli {
    /// Prepared-reference manifest to harvest bases + CDS from.
    #[arg(long)]
    manifest: PathBuf,
    /// Corpus `cases.json` to inventory.
    #[arg(long, default_value = "tests/fixtures/mutalyzer-normalize/cases.json")]
    cases: PathBuf,
    /// Snapshot output directory (holds `transcripts.fna` +
    /// `transcripts.metadata.json`).
    #[arg(
        long,
        default_value = "tests/fixtures/mutalyzer-normalize/reference-snapshot"
    )]
    out_dir: PathBuf,
    /// Verify the committed snapshot matches a fresh harvest and exit non-zero
    /// on any difference, instead of writing. Local-only (needs `--manifest`).
    #[arg(long)]
    check: bool,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// The harvested snapshot: metadata + the FASTA bases (kept together so write
/// and `--check` operate on one consistent pair).
struct Harvest {
    snapshot: TranscriptSnapshot,
    sequences: BTreeMap<String, String>,
    /// Corpus transcript accessions whose exact-version bases are absent from
    /// the reference (reported, not fatal — these are the I2 residuals reserved
    /// for an authoritative fetch).
    absent: Vec<String>,
}

fn run(cli: &Cli) -> anyhow::Result<()> {
    let harvest = harvest(cli)?;

    let metadata_json = harvest.snapshot.to_json()?;
    let fasta = render_fasta(&harvest.sequences);
    let metadata_path = cli.out_dir.join(TRANSCRIPTS_METADATA);
    let fasta_path = cli.out_dir.join(TRANSCRIPTS_FASTA);

    println!(
        "Corpus:    {}\nReference: {}",
        cli.cases.display(),
        cli.manifest.display()
    );
    println!(
        "Harvested {} transcript(s) version-exact; {} absent from reference.",
        harvest.sequences.len(),
        harvest.absent.len()
    );
    if !harvest.absent.is_empty() {
        println!("  Absent (no exact-version bases locally; reserved for authoritative fetch):");
        for accession in &harvest.absent {
            println!("    {accession}");
        }
    }

    if cli.check {
        check_file(&metadata_path, &metadata_json)?;
        check_file(&fasta_path, &fasta)?;
        println!("OK: committed snapshot matches a fresh harvest.");
    } else {
        std::fs::create_dir_all(&cli.out_dir)
            .map_err(|e| anyhow::anyhow!("create {}: {e}", cli.out_dir.display()))?;
        std::fs::write(&metadata_path, &metadata_json)
            .map_err(|e| anyhow::anyhow!("write {}: {e}", metadata_path.display()))?;
        std::fs::write(&fasta_path, &fasta)
            .map_err(|e| anyhow::anyhow!("write {}: {e}", fasta_path.display()))?;
        println!(
            "Wrote {} and {}.",
            metadata_path.display(),
            fasta_path.display()
        );
    }
    Ok(())
}

/// Inventory the corpus, then harvest each version-exact transcript's bases +
/// own CDS frame from the prepared reference.
fn harvest(cli: &Cli) -> anyhow::Result<Harvest> {
    let content = std::fs::read_to_string(&cli.cases)
        .map_err(|e| anyhow::anyhow!("read {}: {e}", cli.cases.display()))?;
    let fixture: Fixture = serde_json::from_str(&content)
        .map_err(|e| anyhow::anyhow!("parse {}: {e}", cli.cases.display()))?;
    let inventory = Inventory::from_fixture(&fixture);

    let provider = MultiFastaProvider::from_manifest(&cli.manifest)
        .map_err(|e| anyhow::anyhow!("load manifest {}: {e}", cli.manifest.display()))?;
    let captured_from = read_prepared_at(&cli.manifest).unwrap_or_default();

    // A stable corpus label (the corpus directory name, e.g.
    // `mutalyzer-normalize`) — not the raw `--cases` path, which would bake a
    // machine-specific absolute path into the committed file.
    let corpus_label = cli
        .cases
        .parent()
        .and_then(|p| p.file_name())
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| cli.cases.display().to_string());

    let mut snapshot = TranscriptSnapshot {
        description: format!(
            "Pinned, version-exact transcript reference snapshot for the {corpus_label} corpus. \
             Generated by `cargo run --features dev --example build_conformance_snapshot`. \
             Bases + CDS harvested version-exact from a prepared reference (CDS via the \
             accession's own frame, cross-build if needed; never a fuzzy sibling version). \
             CDS coordinates are 1-based inclusive."
        ),
        captured_from,
        transcripts: BTreeMap::new(),
    };
    let mut sequences = BTreeMap::new();
    let mut absent = Vec::new();

    for entry in inventory.entries() {
        let acc_ref = &entry.acc_ref;
        // Only transcripts go in this (I2) snapshot; genomic/protein are I3.
        if acc_ref.class != AccessionClass::RefSeqTranscript {
            continue;
        }
        // A RefSeq transcript reference without a pinned version cannot be
        // harvested version-exact — skip and report.
        if acc_ref.version.is_none() {
            absent.push(acc_ref.versioned());
            continue;
        }
        let versioned = acc_ref.versioned();

        // Bases: exact-version only. The fuzzy `get_sequence` fallback would
        // substitute a sibling version's bases, so gate on the exact key.
        if !provider.contains_exact_sequence(&versioned) {
            absent.push(versioned);
            continue;
        }
        let length = provider
            .sequence_length(&versioned)
            .ok_or_else(|| anyhow::anyhow!("{versioned}: exact present but no length"))?;
        let bases = provider
            .get_sequence(&versioned, 0, length)
            .map_err(|e| anyhow::anyhow!("{versioned}: read bases: {e}"))?
            .to_ascii_uppercase();
        // The reported length and the bases we actually read must agree; a
        // mismatch is a broken reference, so fail fast rather than commit a
        // file the integrity test would later reject.
        if bases.len() as u64 != length {
            return Err(anyhow::anyhow!(
                "{versioned}: read {} bases but sequence_length is {length}",
                bases.len()
            ));
        }

        // CDS: the accession's OWN frame, cross-build exact if needed.
        let (cds_start, cds_end, gene_symbol, protein) = match provider
            .cdot_mapper()
            .and_then(|cdot| cdot.get_transcript_exact_any_build(&versioned))
        {
            Some(tx) => {
                // CDS only when BOTH bounds are present; a half-populated CDS
                // is a malformed reading frame, not a coding transcript. cdot
                // CDS is 0-based [start, end); store 1-based inclusive (start+1;
                // end is numerically the 1-based inclusive end).
                let (cds_start, cds_end) = match (tx.cds_start, tx.cds_end) {
                    (Some(start), Some(end)) => (Some(start + 1), Some(end)),
                    _ => (None, None),
                };
                (cds_start, cds_end, tx.gene_name.clone(), tx.protein.clone())
            }
            None => (None, None, None, None),
        };

        let sha256 = sha256_hex(bases.as_bytes());
        snapshot.transcripts.insert(
            versioned.clone(),
            TranscriptEntry {
                cds_start,
                cds_end,
                gene_symbol,
                protein,
                length,
                provenance: Provenance {
                    source: "manifest:transcript-fasta+cdot".to_string(),
                    sha256,
                },
            },
        );
        sequences.insert(versioned, bases);
    }

    absent.sort();
    absent.dedup();
    Ok(Harvest {
        snapshot,
        sequences,
        absent,
    })
}

/// Compare an on-disk file against freshly-rendered content; error with a hint
/// if it differs or is missing.
fn check_file(path: &std::path::Path, rendered: &str) -> anyhow::Result<()> {
    let on_disk = std::fs::read_to_string(path).map_err(|e| {
        anyhow::anyhow!(
            "read {}: {e}; rerun without --check to (re)generate the snapshot",
            path.display()
        )
    })?;
    if on_disk != *rendered {
        return Err(anyhow::anyhow!(
            "{} is out of date; rerun: cargo run --features dev --example \
             build_conformance_snapshot -- --manifest <manifest>",
            path.display()
        ));
    }
    Ok(())
}

/// Lowercase hex SHA-256 of `bytes`.
fn sha256_hex(bytes: &[u8]) -> String {
    let digest = Sha256::digest(bytes);
    let mut s = String::with_capacity(digest.len() * 2);
    for byte in digest {
        s.push_str(&format!("{byte:02x}"));
    }
    s
}

/// Read the manifest's `prepared_at` for snapshot provenance, if present.
fn read_prepared_at(manifest: &std::path::Path) -> Option<String> {
    let content = std::fs::read_to_string(manifest).ok()?;
    let value: serde_json::Value = serde_json::from_str(&content).ok()?;
    value
        .get("prepared_at")
        .and_then(|v| v.as_str())
        .map(|s| s.to_string())
}
