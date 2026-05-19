//! Throughput / latency bench harness for `VariantProjector`.
//!
//! Two modes:
//!   gen-snp   <out.tsv>   — emit g.SNV strings across the CDS of fixed transcripts
//!   gen-indel <out.tsv>   — emit g.del/ins/dup/delins across the CDS of fixed transcripts
//!   bench     <in.tsv>    — run project_all over every line and report variants/sec
//!
//! Build:  cargo build --release --example projector-bench
//! Run:    target/release/examples/projector-bench <manifest.json> <subcmd> <path>
//!
//! Fixture transcripts (all on chr17 to keep the working set small):
//!   NM_000088.4  COL1A1
//!   NM_007294.4  BRCA1
//!   NM_000546.6  TP53

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use clap::{Parser, Subcommand};

use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{FerroError, MultiFastaProvider, ReferenceProvider, VariantProjector};

const TRANSCRIPTS: &[&str] = &["NM_000088.4", "NM_007294.4", "NM_000546.6"];

#[derive(Parser)]
#[command(about = "VariantProjector throughput bench")]
struct Cli {
    /// Path to a ferro-prepare manifest.json.
    manifest: PathBuf,

    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Generate a SNP-dense fixture (every CDS base, one alt per ref).
    GenSnp { out: PathBuf },
    /// Generate an indel-mixed fixture (every 3rd CDS base, random op).
    GenIndel { out: PathBuf },
    /// Run the projector bench over a fixture file.
    Bench {
        fixture: PathBuf,
        /// Stop after this many variants (for quick profiles). 0 = no limit.
        #[arg(long, default_value_t = 0)]
        limit: usize,
        /// Write a deterministic per-variant result line to this path
        /// (input \t n_projections \t sorted_tx_id:c|p|frameshift|intronic|utr;...).
        /// Used to byte-diff outputs across perf changes.
        #[arg(long)]
        dump: Option<PathBuf>,
        /// Run the fixture this many times in a row. Used to amortize the
        /// (large) manifest startup cost when capturing a profile, so the
        /// per-variant hot path dominates the sample distribution. 1 = single
        /// pass (default).
        #[arg(long, default_value_t = 1)]
        repeat: usize,
    },
}

#[derive(Clone)]
struct ArcProvider(Arc<MultiFastaProvider>);

impl ReferenceProvider for ArcProvider {
    fn get_transcript(&self, id: &str) -> Result<Transcript, FerroError> {
        self.0.get_transcript(id)
    }
    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        self.0.get_sequence(id, start, end)
    }
    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        self.0.get_genomic_sequence(contig, start, end)
    }
    fn has_genomic_data(&self) -> bool {
        self.0.has_genomic_data()
    }
    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        self.0.get_protein_sequence(accession, start, end)
    }
    fn has_protein_data(&self) -> bool {
        self.0.has_protein_data()
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    eprintln!("loading manifest: {}", cli.manifest.display());
    let t0 = Instant::now();
    let provider = Arc::new(MultiFastaProvider::from_manifest(&cli.manifest)?);
    eprintln!("manifest loaded in {:?}", t0.elapsed());

    match cli.cmd {
        Cmd::GenSnp { out } => gen_snp(&provider, &out),
        Cmd::GenIndel { out } => gen_indel(&provider, &out),
        Cmd::Bench {
            fixture,
            limit,
            dump,
            repeat,
        } => run_bench(provider, &fixture, limit, dump.as_deref(), repeat),
    }
}

fn write_dump_line(
    w: &mut BufWriter<File>,
    input: &str,
    projections: &[ferro_hgvs::VariantProjection],
) -> Result<(), Box<dyn std::error::Error>> {
    // Build a deterministic, sorted summary string for each variant.
    // Sort by transcript_id to be independent of projection ordering changes
    // that may happen during perf refactors (the projection set must be stable;
    // its ordering is asserted separately by the existing test suite).
    let mut entries: Vec<String> = projections
        .iter()
        .map(|p| {
            let c = p.coding.as_ref().map(|v| v.to_string()).unwrap_or_default();
            let pr = p
                .protein
                .as_ref()
                .map(|v| v.to_string())
                .unwrap_or_default();
            format!(
                "{}|c={}|p={}|fs={}|in={}|utr={}",
                p.transcript_id, c, pr, p.is_frameshift, p.is_intronic, p.is_utr
            )
        })
        .collect();
    entries.sort();
    writeln!(w, "{}\t{}\t{}", input, projections.len(), entries.join(";"))?;
    Ok(())
}

fn cds_segments(cdot: &CdotMapper, tx_id: &str) -> Result<Vec<(String, u64, u64)>, String> {
    let tx = cdot
        .get_transcript(tx_id)
        .ok_or_else(|| format!("missing transcript {tx_id} in cdot"))?;
    let cds_t0 = tx.cds_start.ok_or("transcript has no CDS start")?;
    let cds_t1 = tx.cds_end.ok_or("transcript has no CDS end")?;
    let contig = tx.contig.clone();
    let strand = tx.strand;
    let mut out = Vec::new();
    for exon in &tx.exons {
        let g0 = exon[0];
        let g1 = exon[1];
        let t0 = exon[2];
        let t1 = exon[3];
        // exon tx-range, clipped to CDS.
        let clip_t0 = t0.max(cds_t0);
        let clip_t1 = t1.min(cds_t1);
        if clip_t1 <= clip_t0 {
            continue;
        }
        // Map clipped tx range back to genome. cdot exons store
        // `[g_start, g_end)` in genomic order regardless of strand, and
        // `[t_start, t_end)` in transcript order; for minus-strand
        // transcripts the first tx base of the exon corresponds to the
        // last genomic base (`g_end - 1`), so the genomic range mirrors.
        let (seg_g0, seg_g1) = match strand {
            Strand::Plus => (g0 + (clip_t0 - t0), g1.min(g0 + (clip_t1 - t0))),
            Strand::Minus => (g1 - (clip_t1 - t0), g1 - (clip_t0 - t0)),
            other => return Err(format!("unsupported strand for tx {tx_id}: {other:?}")),
        };
        if seg_g1 > seg_g0 {
            out.push((contig.clone(), seg_g0, seg_g1));
        }
    }
    Ok(out)
}

fn gen_snp(
    provider: &Arc<MultiFastaProvider>,
    out: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let cdot = provider
        .cdot_mapper()
        .ok_or("cdot data not loaded from manifest")?;
    let mut w = BufWriter::new(File::create(out)?);
    let alt_map: HashMap<u8, u8> = [(b'A', b'G'), (b'C', b'T'), (b'G', b'A'), (b'T', b'C')]
        .into_iter()
        .collect();

    let mut total: u64 = 0;
    for &tx_id in TRANSCRIPTS {
        let segs = cds_segments(cdot, tx_id)?;
        for (contig, g0, g1) in segs {
            let seq = provider.get_genomic_sequence(&contig, g0, g1)?;
            for (i, ref_b) in seq.bytes().enumerate() {
                let upper = ref_b.to_ascii_uppercase();
                let Some(&alt) = alt_map.get(&upper) else {
                    continue;
                };
                let pos = g0 + i as u64 + 1; // 1-based for HGVS
                writeln!(w, "{}:g.{}{}>{}", contig, pos, upper as char, alt as char)?;
                total += 1;
            }
        }
    }
    eprintln!("snp fixture: {} variants -> {}", total, out.display());
    Ok(())
}

fn gen_indel(
    provider: &Arc<MultiFastaProvider>,
    out: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let cdot = provider
        .cdot_mapper()
        .ok_or("cdot data not loaded from manifest")?;
    let mut w = BufWriter::new(File::create(out)?);
    // Tiny deterministic LCG (so the fixture is reproducible).
    let mut rng: u64 = 0xc001_d00d_dead_beef;
    let mut next_u8 = || -> u8 {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (rng >> 33) as u8
    };
    let mut total: u64 = 0;
    for &tx_id in TRANSCRIPTS {
        let segs = cds_segments(cdot, tx_id)?;
        for (contig, g0, g1) in segs {
            // Step by 3 to keep the file from being absurdly large.
            for i in (0..(g1 - g0)).step_by(3) {
                let pos = g0 + i + 1; // 1-based
                let r = next_u8() % 4;
                let line = match r {
                    0 => format!("{}:g.{}del", contig, pos),
                    1 => format!("{}:g.{}_{}delinsACG", contig, pos, pos + 2),
                    2 => format!("{}:g.{}_{}insATG", contig, pos, pos + 1),
                    _ => format!("{}:g.{}dup", contig, pos),
                };
                writeln!(w, "{}", line)?;
                total += 1;
            }
        }
    }
    eprintln!("indel fixture: {} variants -> {}", total, out.display());
    Ok(())
}

fn run_bench(
    provider: Arc<MultiFastaProvider>,
    fixture: &Path,
    limit: usize,
    dump_path: Option<&std::path::Path>,
    repeat: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let cdot = provider
        .cdot_mapper()
        .ok_or("cdot data not loaded from manifest")?
        .clone();
    let projector = Projector::new(cdot);
    let vp = VariantProjector::new(projector, ArcProvider(provider.clone()));

    let mut lines: Vec<String> = BufReader::new(File::open(fixture)?)
        .lines()
        .collect::<Result<_, _>>()?;
    lines.retain(|s| !s.trim().is_empty());
    if limit > 0 && lines.len() > limit {
        lines.truncate(limit);
    }
    eprintln!("running over {} variants", lines.len());

    let t_warm = Instant::now();
    // small warmup: first 200 calls (caches FASTA pages, etc.)
    let warm_n = lines.len().min(200);
    for s in lines.iter().take(warm_n) {
        let _ = vp.project_all(s);
    }
    eprintln!("warmup ({} calls) in {:?}", warm_n, t_warm.elapsed());

    let mut latencies_us: Vec<u64> = Vec::with_capacity(lines.len());
    let mut ok = 0u64;
    let mut err = 0u64;
    let mut proj_total = 0u64;
    let mut dump_w = dump_path
        .map(|p| -> Result<BufWriter<File>, Box<dyn std::error::Error>> {
            Ok(BufWriter::new(File::create(p)?))
        })
        .transpose()?;

    let t0 = Instant::now();
    for pass in 0..repeat.max(1) {
        for s in &lines {
            let t = Instant::now();
            match vp.project_all(s) {
                Ok(v) => {
                    ok += 1;
                    proj_total += v.len() as u64;
                    // Only dump from the first pass — the result is
                    // deterministic so subsequent passes would just duplicate.
                    if pass == 0 {
                        if let Some(w) = dump_w.as_mut() {
                            write_dump_line(w, s, &v)?;
                        }
                    }
                }
                Err(e) => {
                    err += 1;
                    if pass == 0 {
                        if let Some(w) = dump_w.as_mut() {
                            writeln!(w, "{}\tERR\t{}", s, e)?;
                        }
                    }
                }
            }
            latencies_us.push(t.elapsed().as_micros() as u64);
        }
    }
    if let Some(mut w) = dump_w {
        w.flush()?;
    }
    let elapsed = t0.elapsed();
    let n = ok + err;
    let v_per_s = n as f64 / elapsed.as_secs_f64();
    eprintln!(
        "processed {} ({} ok, {} err) in {:?}\n  -> {:.1} variants/sec\n  -> {:.1} projections/sec (sum {} projections)",
        n, ok, err, elapsed, v_per_s, proj_total as f64 / elapsed.as_secs_f64(), proj_total
    );

    latencies_us.sort_unstable();
    let p = |q: f64| -> u64 {
        if latencies_us.is_empty() {
            0
        } else {
            let i = ((latencies_us.len() as f64) * q) as usize;
            latencies_us[i.min(latencies_us.len() - 1)]
        }
    };
    eprintln!(
        "latency µs  P50={} P90={} P99={} P999={} max={}",
        p(0.5),
        p(0.9),
        p(0.99),
        p(0.999),
        latencies_us.last().copied().unwrap_or(0)
    );

    // Print a single-line CSV summary that's easy to grep across runs.
    println!(
        "summary,{},{},{},{:.2},{},{},{},{}",
        fixture.display(),
        n,
        elapsed.as_micros(),
        v_per_s,
        p(0.5),
        p(0.9),
        p(0.99),
        latencies_us.last().copied().unwrap_or(0)
    );
    Ok(())
}
