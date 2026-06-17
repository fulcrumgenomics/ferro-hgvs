//! Build `derived_refseqgene_placements` (#728): for each pinned `NG_` version,
//! fetch its GenBank record + sequence, derive its chromosomal placement per
//! genome build, and collect the validated placements. The pure core
//! (`build_placements`) takes `fetch`/`derive` closures so it is testable
//! without network or a concrete provider; the provider-bound wrapper and the
//! NCBI EFetch helper live alongside it (Tasks 2-3).

use std::process::Command;

use crate::reference::derived_placement::DerivedPlacement;
use crate::reference::{GenomicPlacement, Strand};
use crate::FerroError;

/// Derive placements for `accessions` over `builds`, fetching record/sequence
/// via `fetch(accession, rettype)` (`rettype` is `"gb"` or `"fasta"`) and
/// deriving via `derive(genbank, seq, build)`. Returns the validated records
/// and a list of human-readable decline messages. Never panics on a single
/// accession's failure.
pub fn build_placements(
    accessions: &[String],
    builds: &[String],
    fetch: impl Fn(&str, &str) -> Result<String, FerroError>,
    derive: impl Fn(&str, &[u8], &str) -> Option<GenomicPlacement>,
) -> (Vec<DerivedPlacement>, Vec<String>) {
    let mut records = Vec::new();
    let mut declined = Vec::new();
    for ng in accessions {
        let genbank = match fetch(ng, "gb") {
            Ok(t) => t,
            Err(e) => {
                declined.push(format!("{ng}: efetch gb failed: {e}"));
                continue;
            }
        };
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
                records.push(placement_to_record(ng, &p));
                any_build = true;
            }
        }
        if !any_build {
            declined.push(format!(
                "{ng}: no validated placement on any requested build"
            ));
        }
    }
    (records, declined)
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

        let (records, declined) = build_placements(&accessions, &builds, fetch, derive);

        assert_eq!(
            records.len(),
            1,
            "one validated placement (NG_aaa on GRCh38)"
        );
        assert_eq!(records[0].parent, "NG_aaa.1");
        assert_eq!(records[0].nc_start, 100);
        assert_eq!(records[0].strand, "+");
        // NG_bad declined (fetch failure); NG_aaa not fully declined (had a build).
        assert_eq!(declined.len(), 1);
        assert!(declined[0].contains("NG_bad.1"));
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
        let (records, declined) = build_placements(&accessions, &builds, fetch, derive);
        assert!(records.is_empty());
        assert_eq!(declined.len(), 1);
        assert!(declined[0].contains("no validated placement"));
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
}
