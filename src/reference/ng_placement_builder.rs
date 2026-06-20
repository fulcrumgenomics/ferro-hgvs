//! Build `derived_refseqgene_placements` (#728): for each pinned `NG_` version,
//! fetch its GenBank record + sequence, derive its chromosomal placement per
//! genome build, and collect the validated placements. The pure core
//! (`build_placements`) takes `fetch`/`derive` closures so it is testable
//! without network or a concrete provider; the provider-bound wrapper and the
//! NCBI EFetch helper live alongside it (Tasks 2-3).

use std::process::Command;

use crate::data::cdot::CdotMapper;
use crate::hgvs::parser::accession::parse_accession;
use crate::hgvs::variant::Accession;
use crate::reference::derived_placement::{
    derive_ng_placement, parse_ng_gene_transcripts, DerivedPlacement, DerivedPlacements,
    GenomeSlice, NcExonSource,
};
use crate::reference::multi_fasta::MultiFastaProvider;
use crate::reference::{GenomicPlacement, ReferenceProvider, Strand};
use crate::FerroError;

/// Output of [`build_placements`]: the validated placement records, the
/// hosting map entries for `ng_hosted_transcripts.json` (#792), and a list of
/// human-readable decline messages.
pub struct BuildOutput {
    /// Validated `DerivedPlacement` records, one per validated (NG_, build) pair.
    pub placements: Vec<DerivedPlacement>,
    /// Per-NG_ hosting records: `(ng_acc_version, [(gene_upper, transcript_id)])`.
    /// Collected from the same GenBank fetch as the placements — no second fetch.
    pub hosted: Vec<(String, Vec<(String, String)>)>,
    /// Human-readable decline messages (fetch failures, no-build declines).
    pub declined: Vec<String>,
}

/// Derive placements for `accessions` over `builds`, fetching record/sequence
/// via `fetch(accession, rettype)` (`rettype` is `"gb"` or `"fasta"`) and
/// deriving via `derive(genbank, seq, build)`. Returns a [`BuildOutput`]
/// containing the validated records, per-NG_ hosting entries (#792), and
/// human-readable decline messages. Never panics on a single accession's
/// failure.
pub fn build_placements(
    accessions: &[String],
    builds: &[String],
    fetch: impl Fn(&str, &str) -> Result<String, FerroError>,
    derive: impl Fn(&str, &[u8], &str) -> Option<GenomicPlacement>,
) -> BuildOutput {
    let mut placements = Vec::new();
    let mut hosted = Vec::new();
    let mut declined = Vec::new();
    for ng in accessions {
        let genbank = match fetch(ng, "gb") {
            Ok(t) => t,
            Err(e) => {
                declined.push(format!("{ng}: efetch gb failed: {e}"));
                continue;
            }
        };
        // Collect hosted-transcript entries from this GenBank fetch (#792).
        let gene_txs = parse_ng_gene_transcripts(&genbank);
        if !gene_txs.is_empty() {
            hosted.push((ng.clone(), gene_txs));
        }
        let seq = match fetch(ng, "fasta") {
            Ok(fa) => fasta_bases(&fa),
            Err(e) => {
                declined.push(format!("{ng}: efetch fasta failed: {e}"));
                continue;
            }
        };
        let mut any_build = false;
        for build in builds {
            if let Some(p) = derive(&genbank, &seq, build) {
                placements.push(placement_to_record(ng, &p));
                any_build = true;
            }
        }
        if !any_build {
            declined.push(format!(
                "{ng}: no validated placement on any requested build"
            ));
        }
    }
    BuildOutput {
        placements,
        hosted,
        declined,
    }
}

/// Serializable record from a derived placement. Validation already passed, so
/// `mismatch_fraction` is recorded as 0.0 (the producer does not surface the
/// exact value).
fn placement_to_record(parent: &str, p: &GenomicPlacement) -> DerivedPlacement {
    DerivedPlacement {
        parent: parent.to_string(),
        nc: p.nc.full(),
        nc_start: p.nc_start,
        nc_end: p.nc_end,
        strand: match p.strand {
            Strand::Plus => "+",
            Strand::Minus => "-",
            _ => "?",
        }
        .to_string(),
        anchored_by: String::new(),
        mismatch_fraction: 0.0,
    }
}

/// NCBI nuccore EFetch URL for `accession` in `rettype` (`gb`/`fasta`), text mode.
fn efetch_url(accession: &str, rettype: &str) -> String {
    format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype={rettype}&retmode=text"
    )
}

/// EFetch one accession's record as text. Tries `curl` (present on the prepare
/// hosts), falling back to a blocking reqwest GET — matching `download_file`.
/// A non-empty body is required; an empty/failed fetch is an error.
pub fn efetch_text(accession: &str, rettype: &str) -> Result<String, FerroError> {
    let url = efetch_url(accession, rettype);
    // curl path
    if let Ok(out) = Command::new("curl")
        .args(["-sS", "--fail", "-L", &url])
        .output()
    {
        if out.status.success() {
            let body = String::from_utf8_lossy(&out.stdout).into_owned();
            if !body.trim().is_empty() {
                return Ok(body);
            }
        }
    }
    // reqwest fallback
    let body = reqwest::blocking::Client::builder()
        .build()
        .and_then(|c| c.get(&url).send())
        .and_then(|r| r.error_for_status())
        .and_then(|r| r.text())
        .map_err(|e| FerroError::Io {
            msg: format!("efetch {accession} ({rettype}): {e}"),
        })?;
    if body.trim().is_empty() {
        return Err(FerroError::Io {
            msg: format!("efetch {accession} ({rettype}): empty response"),
        });
    }
    Ok(body)
}

/// Uppercase bases from a FASTA string (strip `>` headers + whitespace).
fn fasta_bases(fasta: &str) -> Vec<u8> {
    fasta
        .lines()
        .filter(|l| !l.starts_with('>'))
        .flat_map(|l| l.trim().bytes())
        .map(|b| b.to_ascii_uppercase())
        .collect()
}

/// `NcExonSource` over a cdot mapper for one build. Since #742 the loaded
/// `CdotTranscript.exons` entry `[genome_start, genome_end, ..]` is in the HGVS
/// convention — `genome_start` 1-based inclusive, `genome_end` 1-based exclusive
/// — so the 1-based inclusive interval `derive_ng_placement` expects is
/// `(genome_start, genome_end - 1)`.
struct CdotNcExons<'a> {
    cdot: &'a CdotMapper,
    build: &'a str,
}

impl NcExonSource for CdotNcExons<'_> {
    fn nc_exons(&self, transcript_id: &str) -> Option<Vec<(u64, u64)>> {
        let tx = self
            .cdot
            .get_transcript_on_build_exact(transcript_id, self.build)?;
        Some(tx.exons.iter().map(|e| (e[0], e[1] - 1)).collect())
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

/// Output of [`derive_placements_for_accessions`]: the validated placement
/// artifact, the per-NG_ hosting records (#792), and the number of accessions
/// that declined on every requested build.
pub struct DeriveOutput {
    /// Validated `DerivedPlacements` artifact ready to write to disk.
    pub artifact: DerivedPlacements,
    /// Per-NG_ hosting records: `(ng_acc_version, [(gene_upper, transcript_id)])`.
    /// Collected from the same GenBank fetch as the placements — no second fetch.
    pub hosted: Vec<(String, Vec<(String, String)>)>,
    /// Number of accessions that declined on every requested build.
    pub declined_count: usize,
}

/// Derive validated placements for `accessions` over `builds` using `provider`
/// (cdot + genome) and live NCBI EFetch. Declines are logged via `warn!`.
/// Returns a [`DeriveOutput`] containing the validated artifact, per-NG_
/// hosting records (#792), and the number of declined accessions — so callers
/// can surface an all-declined run rather than reporting a bare success.
pub fn derive_placements_for_accessions(
    provider: &MultiFastaProvider,
    accessions: &[String],
    builds: &[String],
    description: String,
) -> Result<DeriveOutput, FerroError> {
    let cdot = provider.cdot_mapper().ok_or_else(|| FerroError::Io {
        msg: "manifest has no cdot; cannot derive NG_ placements".to_string(),
    })?;
    let genome = ProviderGenomeSlice { provider };
    let derive = |genbank: &str, seq: &[u8], build: &str| -> Option<GenomicPlacement> {
        let nc_source = CdotNcExons { cdot, build };
        derive_ng_placement(genbank, seq, &nc_source, &genome)
    };
    let out = build_placements(accessions, builds, efetch_text, derive);
    for d in &out.declined {
        log::warn!("derive NG_ placement declined: {d}");
    }
    let declined_count = out.declined.len();
    Ok(DeriveOutput {
        artifact: DerivedPlacements {
            description,
            placements: out.placements,
        },
        hosted: out.hosted,
        declined_count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::accession::parse_accession;
    use crate::reference::{GenomicPlacement, Strand};

    fn pl(nc: &str, start: u64, end: u64) -> GenomicPlacement {
        let (_, acc) = parse_accession(nc).unwrap();
        GenomicPlacement {
            nc: acc,
            parent_start: 1,
            nc_start: start,
            nc_end: end,
            strand: Strand::Plus,
        }
    }

    #[test]
    fn build_placements_aggregates_and_declines() {
        let accessions = vec!["NG_aaa.1".to_string(), "NG_bad.1".to_string()];
        let builds = vec!["GRCh38".to_string(), "GRCh37".to_string()];

        // fetch: succeeds for NG_aaa, fails for NG_bad on the gb call.
        let fetch = |acc: &str, rettype: &str| -> Result<String, FerroError> {
            if acc == "NG_bad.1" {
                return Err(FerroError::Io {
                    msg: "efetch failed".into(),
                });
            }
            // gb -> dummy genbank, fasta -> dummy fasta with bases
            Ok(if rettype == "fasta" {
                ">hdr\nACGT\n".to_string()
            } else {
                "LOCUS dummy".to_string()
            })
        };
        // derive: produces a placement only on GRCh38, declines GRCh37.
        let derive = |_gb: &str, _seq: &[u8], build: &str| -> Option<GenomicPlacement> {
            if build == "GRCh38" {
                Some(pl("NC_000011.10", 100, 200))
            } else {
                None
            }
        };

        let out = build_placements(&accessions, &builds, fetch, derive);

        assert_eq!(
            out.placements.len(),
            1,
            "one validated placement (NG_aaa on GRCh38)"
        );
        assert_eq!(out.placements[0].parent, "NG_aaa.1");
        assert_eq!(out.placements[0].nc_start, 100);
        assert_eq!(out.placements[0].strand, "+");
        // NG_bad declined (fetch failure); NG_aaa not fully declined (had a build).
        assert_eq!(out.declined.len(), 1);
        assert!(out.declined[0].contains("NG_bad.1"));
    }

    #[test]
    fn build_placements_declines_when_no_build_validates() {
        let accessions = vec!["NG_aaa.1".to_string()];
        let builds = vec!["GRCh38".to_string()];
        let fetch = |_a: &str, rettype: &str| -> Result<String, FerroError> {
            Ok(if rettype == "fasta" {
                ">h\nACGT\n".into()
            } else {
                "LOCUS".into()
            })
        };
        let derive = |_gb: &str, _seq: &[u8], _b: &str| -> Option<GenomicPlacement> { None };
        let out = build_placements(&accessions, &builds, fetch, derive);
        assert!(out.placements.is_empty());
        assert_eq!(out.declined.len(), 1);
        assert!(out.declined[0].contains("no validated placement"));
        assert!(out.hosted.is_empty());
    }

    #[test]
    fn fasta_bases_strips_headers_and_uppercases() {
        assert_eq!(fasta_bases(">hdr\nacgt\nTT\n"), b"ACGTTT".to_vec());
    }

    #[test]
    fn efetch_url_is_nuccore_text() {
        let url = efetch_url("NG_012337.3", "gb");
        assert!(url.contains("efetch.fcgi"));
        assert!(url.contains("db=nuccore"));
        assert!(url.contains("id=NG_012337.3"));
        assert!(url.contains("rettype=gb"));
        assert!(url.contains("retmode=text"));
    }

    #[test]
    fn derive_for_real_reference_when_available() {
        // Gated: only runs against a prepared reference (repo convention).
        let manifest = match std::env::var("FERRO_MANIFEST") {
            Ok(m) => std::path::PathBuf::from(m),
            Err(_) => return, // skip when no reference is configured
        };
        let provider = crate::reference::multi_fasta::MultiFastaProvider::from_manifest(&manifest)
            .expect("load manifest");
        let accs = vec!["NG_012337.3".to_string()];
        let builds = vec!["GRCh38".to_string()];
        let out = super::derive_placements_for_accessions(&provider, &accs, &builds, "test".into())
            .expect("derive");
        // Either a validated placement or a principled decline (no panic, valid artifact).
        assert!(out.artifact.placements.len() <= 1);
    }
}
