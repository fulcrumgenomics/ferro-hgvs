//! Produce the `derived_refseqgene_placements.json` artifact (#728): for each
//! pinned `NG_` version in an input list, fetch its GenBank record + sequence
//! from NCBI, derive its chromosomal placement from the record's exon
//! annotation + cdot, validate by full-sequence comparison against the genome,
//! and write the validated placements.
//!
//! This is the prepare-time **producer** for the runtime-consumption path
//! (`conformance`/`reference::derived_placement` + the manifest's
//! `derived_refseqgene_placements` field). It runs locally (it needs a prepared
//! reference for cdot + genome, and network for EFetch); the committed/manifest
//! artifact it writes is then consumed offline. Decline-rather-than-mis-anchor
//! is enforced by [`derive_ng_placement`] — an accession that cannot be proven
//! is simply omitted (reported), never guessed.
//!
//! Run:
//!   cargo run --features dev --example derive_ng_placements -- \
//!     --manifest <ref>/manifest.json --accessions ng_versions.txt \
//!     --out <ref>/derived_refseqgene_placements.json [--builds GRCh38,GRCh37]

use std::path::PathBuf;
use std::process::{Command, ExitCode};

use clap::Parser;

use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::hgvs::parser::accession::parse_accession;
use ferro_hgvs::hgvs::variant::Accession;
use ferro_hgvs::reference::derived_placement::{
    derive_ng_placement, DerivedPlacement, DerivedPlacements, GenomeSlice, NcExonSource,
};
use ferro_hgvs::reference::multi_fasta::MultiFastaProvider;
use ferro_hgvs::reference::{ReferenceProvider, Strand};

#[derive(Parser, Debug)]
#[command(about = "Derive version-independent NG_/LRG_ genomic placements (#728)")]
struct Cli {
    /// Prepared-reference manifest (for cdot + genome).
    #[arg(long)]
    manifest: PathBuf,
    /// Newline-delimited file of exact `NG_` versions to derive (e.g.
    /// `NG_012337.1`). Blank lines and `#` comments are ignored.
    #[arg(long)]
    accessions: PathBuf,
    /// Output JSON path (`derived_refseqgene_placements.json`).
    #[arg(long)]
    out: PathBuf,
    /// Comma-separated genome builds to derive for (default `GRCh38`).
    #[arg(long, default_value = "GRCh38")]
    builds: String,
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

/// `NcExonSource` over a cdot mapper for a fixed build. The loaded
/// `CdotTranscript.exons` entry `[genome_start, genome_end, ..]` is **0-based**
/// half-open (e.g. SDHD exon 1 = `[112086872, 112086959]`), so the 1-based
/// inclusive interval `derive_affine` expects is `(genome_start + 1, genome_end)`.
struct CdotNcExons<'a> {
    cdot: &'a CdotMapper,
    build: &'a str,
}

impl NcExonSource for CdotNcExons<'_> {
    fn nc_exons(&self, transcript_id: &str) -> Option<Vec<(u64, u64)>> {
        let tx = self
            .cdot
            .get_transcript_on_build_exact(transcript_id, self.build)?;
        Some(tx.exons.iter().map(|e| (e[0] + 1, e[1])).collect())
    }
    fn nc_contig(&self, transcript_id: &str) -> Option<Accession> {
        let tx = self
            .cdot
            .get_transcript_on_build_exact(transcript_id, self.build)?;
        match parse_accession(&tx.contig) {
            Ok(("", acc)) => Some(acc),
            _ => None,
        }
    }
}

/// `GenomeSlice` over a provider: 1-based inclusive `[start, end]` → the
/// provider's 0-based half-open `get_sequence`.
struct ProviderGenomeSlice<'a> {
    provider: &'a MultiFastaProvider,
}

impl GenomeSlice for ProviderGenomeSlice<'_> {
    fn slice(&self, nc: &Accession, start: u64, end: u64) -> Option<Vec<u8>> {
        if start == 0 || end < start {
            return None;
        }
        self.provider
            .get_sequence(&nc.full(), start - 1, end)
            .ok()
            .map(String::into_bytes)
    }
}

fn run(cli: &Cli) -> anyhow::Result<()> {
    let provider = MultiFastaProvider::from_manifest(&cli.manifest)
        .map_err(|e| anyhow::anyhow!("load manifest {}: {e}", cli.manifest.display()))?;
    let cdot = provider
        .cdot_mapper()
        .ok_or_else(|| anyhow::anyhow!("manifest has no cdot; cannot derive placements"))?;

    let accessions = read_accession_list(&cli.accessions)?;
    let builds: Vec<&str> = cli.builds.split(',').map(str::trim).collect();

    let mut placements = Vec::new();
    let mut declined = Vec::new();
    for ng in &accessions {
        let genbank = match efetch(ng, "gb") {
            Ok(t) => t,
            Err(e) => {
                declined.push(format!("{ng}: efetch gb failed: {e}"));
                continue;
            }
        };
        let ng_seq = match efetch(ng, "fasta") {
            Ok(fa) => fasta_bases(&fa),
            Err(e) => {
                declined.push(format!("{ng}: efetch fasta failed: {e}"));
                continue;
            }
        };
        let mut any_build = false;
        for build in &builds {
            let nc_source = CdotNcExons { cdot, build };
            let genome = ProviderGenomeSlice {
                provider: &provider,
            };
            if let Some(p) = derive_ng_placement(&genbank, &ng_seq, &nc_source, &genome) {
                placements.push(to_record(ng, &p));
                any_build = true;
            }
        }
        if !any_build {
            declined.push(format!(
                "{ng}: no validated placement on any requested build"
            ));
        }
    }

    let artifact = DerivedPlacements {
        description: format!(
            "GENERATED by `cargo run --features dev --example derive_ng_placements`. \
             Version-independent NG_/LRG_ placements (#728) derived from each NG_'s \
             GenBank exon annotation + cdot, validated by full-sequence comparison. \
             Source manifest: {}.",
            cli.manifest.display()
        ),
        placements,
    };
    let json = artifact
        .to_json()
        .map_err(|e| anyhow::anyhow!("serialize: {e}"))?;
    std::fs::write(&cli.out, &json)
        .map_err(|e| anyhow::anyhow!("write {}: {e}", cli.out.display()))?;

    println!(
        "Derived {} placement(s) for {} accession(s); {} declined. Wrote {}.",
        artifact.placements.len(),
        accessions.len(),
        declined.len(),
        cli.out.display()
    );
    for d in &declined {
        println!("  declined: {d}");
    }
    Ok(())
}

/// Build a serializable record from a derived placement. (Provenance: the
/// validation already passed, so the mismatch fraction is ~0; the producer does
/// not currently surface the exact value, so record 0.0.)
fn to_record(ng: &str, p: &ferro_hgvs::reference::GenomicPlacement) -> DerivedPlacement {
    DerivedPlacement {
        parent: ng.to_string(),
        nc: p.nc.full(),
        nc_start: p.nc_start,
        nc_end: p.nc_end,
        strand: match p.strand {
            Strand::Plus => "+",
            Strand::Minus => "-",
            _ => "?", // Strand is #[non_exhaustive]; a derived placement is only
                      // ever Plus/Minus (Unknown declines upstream).
        }
        .to_string(),
        anchored_by: String::new(),
        mismatch_fraction: 0.0,
    }
}

/// Read exact accessions from a newline-delimited file (skipping blanks/`#`).
fn read_accession_list(path: &std::path::Path) -> anyhow::Result<Vec<String>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| anyhow::anyhow!("read {}: {e}", path.display()))?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(str::to_string)
        .collect())
}

/// EFetch one accession from NCBI nuccore in the given `rettype` (`gb`/`fasta`).
fn efetch(accession: &str, rettype: &str) -> anyhow::Result<String> {
    let url = format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype={rettype}&retmode=text"
    );
    let out = Command::new("curl")
        .args(["-sS", "--fail", &url])
        .output()
        .map_err(|e| anyhow::anyhow!("spawn curl: {e}"))?;
    if !out.status.success() {
        anyhow::bail!("curl exit {:?}", out.status.code());
    }
    let body = String::from_utf8_lossy(&out.stdout).into_owned();
    if body.trim().is_empty() {
        anyhow::bail!("empty response");
    }
    Ok(body)
}

/// Extract uppercase bases from a FASTA string (strip headers + whitespace).
fn fasta_bases(fasta: &str) -> Vec<u8> {
    fasta
        .lines()
        .filter(|l| !l.starts_with('>'))
        .flat_map(|l| l.trim().bytes())
        .map(|b| b.to_ascii_uppercase())
        .collect()
}
