//! Produce `derived_transcript_placements.json` (#790): for each supplemental
//! old transcript version absent from GRCh38 cdot but with a sibling version
//! present, derive its exon→genome map by sibling-anchoring + genome readback,
//! and write the validated structures. Runs locally (needs a prepared
//! reference for cdot + genome + supplemental). Decline-rather-than-mis-anchor
//! is enforced by `derive_transcript_structure`.
//!
//! Usage:
//!   cargo run --features dev --example derive_tx_placements -- \
//!     --manifest <ref>/manifest.json \
//!     --out <ref>/derived_transcript_placements.json [--check] \
//!     [--accessions NM_003002.2,...]

use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;

use ferro_hgvs::reference::derived_tx_structure::{
    derive_transcript_structure, DerivedTxStructure, DerivedTxStructures, GenomeSlice,
};
use ferro_hgvs::reference::multi_fasta::MultiFastaProvider;
use ferro_hgvs::reference::ReferenceProvider;

/// Default candidate old version(s). A later refinement can enumerate
/// supplemental `NM_` versions absent from cdot with a sibling present; the
/// single proven case is the first cut (YAGNI).
const DEFAULT_ACCESSIONS: &str = "NM_003002.2";

/// How many sibling versions to probe in each direction for a version present
/// in cdot — `candidate_siblings` tries the next-higher versions first, then the
/// next-lower (e.g. for `NM_003002.2` this probes `.3`/`.4`/`.5`, then `.1`).
const MAX_SIBLING_PROBE: u32 = 3;

#[derive(Parser, Debug)]
#[command(about = "Derive cdot-absent old transcript exon→genome structures (#790)")]
struct Cli {
    /// Prepared-reference manifest (for cdot + genome + supplemental mRNA/CDS).
    #[arg(long)]
    manifest: PathBuf,
    /// Output JSON path (`derived_transcript_placements.json`).
    #[arg(long)]
    out: PathBuf,
    /// Comma-separated old `NM_` versions to derive (default `NM_003002.2`).
    #[arg(long, default_value = DEFAULT_ACCESSIONS)]
    accessions: String,
    /// Verify the on-disk `--out` is up to date; exit non-zero on drift instead
    /// of writing.
    #[arg(long)]
    check: bool,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(code) => code,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// `GenomeSlice` (Task-2 trait) over a provider. The trait contract is a
/// 1-based **inclusive** `start_1based` and 1-based **exclusive**
/// `end_exclusive` (the same convention as a loaded `CdotTranscript` exon's
/// `[genome_start, genome_end]`), which maps to the provider's 0-based
/// half-open `get_sequence` as `[start_1based - 1, end_exclusive - 1)`. This
/// mirrors the Task-2 `FakeGenome` fixture (`(start-1)..(end-1)`); a clean
/// genome readback of a sibling's exons reproduces its supplemental mRNA
/// exactly (verified: `NM_003002.3` → 1395 nt, equal).
struct ProviderGenomeSlice<'a> {
    provider: &'a MultiFastaProvider,
}

impl GenomeSlice for ProviderGenomeSlice<'_> {
    fn slice(&self, contig: &str, start_1based: u64, end_exclusive: u64) -> Option<String> {
        if start_1based == 0 || end_exclusive <= start_1based {
            return None;
        }
        self.provider
            .get_sequence(contig, start_1based - 1, end_exclusive - 1)
            .ok()
    }
}

/// Split an accession into `(base, version)` (`"NM_003002.2"` → `("NM_003002", 2)`).
fn split_versioned(accession: &str) -> Option<(&str, u32)> {
    let (base, ver) = accession.rsplit_once('.')?;
    Some((base, ver.parse().ok()?))
}

/// Candidate sibling versions of `old_id` present **exactly** in cdot, ordered
/// by preference: next-higher versions first (ascending), then next-lower
/// versions (descending). Each is `(sibling_id, &transcript)`. Because two
/// sibling versions can differ from the old version only at their poly-A tail
/// (which breaks the exact-overlap anchor for one but not the other), the
/// caller tries them in order and keeps the first that derives cleanly.
fn candidate_siblings<'a>(
    old_id: &str,
    cdot: &'a ferro_hgvs::data::cdot::CdotMapper,
) -> Vec<(String, &'a ferro_hgvs::data::cdot::CdotTranscript)> {
    let mut out = Vec::new();
    let Some((base, old_ver)) = split_versioned(old_id) else {
        return out;
    };
    // Next-higher first (ascending), then next-lower (descending). Accession
    // versions start at 1, so clamp the lower bound to >= 1 — a `.0` probe is
    // always a wasted `has_transcript_exact` miss.
    let higher = (old_ver + 1)..=(old_ver + MAX_SIBLING_PROBE);
    let lower = (old_ver.saturating_sub(MAX_SIBLING_PROBE).max(1)..old_ver).rev();
    for ver in higher.chain(lower) {
        let candidate = format!("{base}.{ver}");
        // `get_transcript` can version-fall-back, so only trust a sibling
        // confirmed present at the exact requested version.
        if cdot.has_transcript_exact(&candidate) {
            if let Some(tx) = cdot.get_transcript(&candidate) {
                out.push((candidate, tx));
            }
        }
    }
    out
}

fn run(cli: &Cli) -> anyhow::Result<ExitCode> {
    // Load WITHOUT injecting the manifest's `derived_transcript_placements`
    // artifact (#800): once that key is wired to this producer's own prior output,
    // injecting it would make the target accession present in cdot and the
    // derivation would decline (skip) it, silently producing an empty artifact on
    // a re-run or CI `--check`. The producer derives over real cdot records only.
    let provider = MultiFastaProvider::from_manifest_without_derived_tx(&cli.manifest)
        .map_err(|e| anyhow::anyhow!("load manifest {}: {e}", cli.manifest.display()))?;
    let cdot = provider
        .cdot_mapper()
        .ok_or_else(|| anyhow::anyhow!("manifest has no cdot; cannot derive tx placements"))?;
    let genome = ProviderGenomeSlice {
        provider: &provider,
    };

    let accessions: Vec<String> = cli
        .accessions
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();

    let mut structures: Vec<DerivedTxStructure> = Vec::new();
    let mut declined = 0usize;
    for old_id in &accessions {
        match derive_one(old_id, &provider, cdot, &genome) {
            Some(s) => structures.push(s),
            // `derive_one` logs the specific reason in each decline branch.
            None => declined += 1,
        }
    }
    // Deterministic order for a stable `--check` diff.
    structures.sort_by(|a, b| a.accession.cmp(&b.accession));

    // The description is intentionally path-free so `--check` compares only the
    // derived structures, not the absolute manifest path (which differs between
    // machines/CI and would otherwise trigger spurious staleness).
    let description = "GENERATED by `cargo run --features dev --example \
         derive_tx_placements`. Exon→genome structures for cdot-absent old \
         transcript versions (#790), each derived by anchoring to a sibling \
         cdot version and validated by genome readback."
        .to_string();
    let artifact = DerivedTxStructures {
        description,
        structures,
    };
    let json = artifact
        .to_json()
        .map_err(|e| anyhow::anyhow!("serialize: {e}"))?;

    if cli.check {
        let existing = std::fs::read_to_string(&cli.out)
            .map_err(|e| anyhow::anyhow!("read {} for --check: {e}", cli.out.display()))?;
        if existing == json {
            println!("{} is up to date.", cli.out.display());
            return Ok(ExitCode::SUCCESS);
        }
        eprintln!(
            "{} is STALE. Regenerate with:\n  cargo run --features dev --example \
             derive_tx_placements -- --manifest {} --out {} --accessions {}",
            cli.out.display(),
            cli.manifest.display(),
            cli.out.display(),
            cli.accessions,
        );
        return Ok(ExitCode::FAILURE);
    }

    std::fs::write(&cli.out, &json)
        .map_err(|e| anyhow::anyhow!("write {}: {e}", cli.out.display()))?;
    println!(
        "Derived {} structure(s) for {} accession(s), {} declined. Wrote {}.",
        artifact.structures.len(),
        accessions.len(),
        declined,
        cli.out.display()
    );
    Ok(ExitCode::SUCCESS)
}

/// Derive one old version's structure, or `None` if any input is missing or the
/// derivation declines.
fn derive_one(
    old_id: &str,
    provider: &MultiFastaProvider,
    cdot: &ferro_hgvs::data::cdot::CdotMapper,
    genome: &ProviderGenomeSlice,
) -> Option<DerivedTxStructure> {
    if cdot.has_transcript_exact(old_id) {
        eprintln!("  skipping {old_id}: present in cdot (not an old/absent version)");
        return None;
    }
    let siblings = candidate_siblings(old_id, cdot);
    if siblings.is_empty() {
        eprintln!("  declined: {old_id} (no sibling version present in cdot)");
        return None;
    }

    // Old mRNA bases from the supplemental FASTA (0-based half-open).
    let old_mrna = provider
        .get_sequence_length(old_id)
        .and_then(|len| provider.get_sequence(old_id, 0, len))
        .map_err(|e| eprintln!("  declined: {old_id} (no supplemental mRNA: {e})"))
        .ok()?;

    // Try each sibling in preference order; keep the first that derives cleanly.
    for (sibling_id, sibling) in &siblings {
        if let Some(s) = derive_transcript_structure(old_id, &old_mrna, sibling_id, sibling, genome)
        {
            return Some(s);
        }
        eprintln!("  {old_id}: sibling {sibling_id} did not derive cleanly; trying next");
    }
    eprintln!("  declined: {old_id} (no sibling anchored cleanly)");
    None
}
