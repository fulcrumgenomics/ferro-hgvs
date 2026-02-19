//! Unified accession availability checking.
//!
//! This module provides utilities for checking what accessions are available
//! in a ferro reference directory, used by both mutalyzer and biocommons
//! preparation to avoid redundant NCBI fetches.

use crate::FerroError;
use std::collections::HashSet;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Categories of accession sources in ferro reference.
#[derive(Debug, Clone, Default)]
pub struct AccessionSources {
    /// Accessions in RefSeq transcript FASTAs (NM_, NR_, XM_, XR_)
    pub transcripts: HashSet<String>,
    /// Accessions in supplemental FASTAs (patterns, legacy, genbank)
    pub supplemental: HashSet<String>,
    /// Accessions in genome FASTAs (NC_)
    pub genomic: HashSet<String>,
    /// Accessions in RefSeqGene FASTAs (NG_)
    pub refseqgene: HashSet<String>,
    /// Accessions in LRG FASTAs
    pub lrg: HashSet<String>,
}

impl AccessionSources {
    /// Load all accession sources from ferro reference directory.
    pub fn from_ferro_reference(ferro_ref: &Path) -> Result<Self, FerroError> {
        Ok(Self {
            transcripts: load_accessions_from_dir(&ferro_ref.join("transcripts"))?,
            supplemental: load_accessions_from_dir(&ferro_ref.join("supplemental"))?,
            genomic: load_accessions_from_dir(&ferro_ref.join("genome"))?,
            refseqgene: load_accessions_from_dir(&ferro_ref.join("refseqgene"))?,
            lrg: load_accessions_from_dir(&ferro_ref.join("lrg"))?,
        })
    }

    /// Check if an accession is available in any source.
    pub fn contains(&self, accession: &str) -> bool {
        self.transcripts.contains(accession)
            || self.supplemental.contains(accession)
            || self.genomic.contains(accession)
            || self.refseqgene.contains(accession)
            || self.lrg.contains(accession)
    }

    /// Get total count of unique available accessions.
    pub fn total_count(&self) -> usize {
        // Use a HashSet to deduplicate across sources
        let mut all = HashSet::new();
        all.extend(self.transcripts.iter().cloned());
        all.extend(self.supplemental.iter().cloned());
        all.extend(self.genomic.iter().cloned());
        all.extend(self.refseqgene.iter().cloned());
        all.extend(self.lrg.iter().cloned());
        all.len()
    }

    /// Find which accessions from a list are missing from all sources.
    pub fn find_missing<'a>(&self, accessions: &'a [String]) -> Vec<&'a String> {
        accessions
            .iter()
            .filter(|acc| !self.contains(acc))
            .collect()
    }

    /// Check which source contains a specific accession.
    pub fn source_for(&self, accession: &str) -> Option<&'static str> {
        if self.transcripts.contains(accession) {
            Some("transcripts")
        } else if self.supplemental.contains(accession) {
            Some("supplemental")
        } else if self.genomic.contains(accession) {
            Some("genomic")
        } else if self.refseqgene.contains(accession) {
            Some("refseqgene")
        } else if self.lrg.contains(accession) {
            Some("lrg")
        } else {
            None
        }
    }

    /// Print summary of available accessions.
    pub fn print_summary(&self) {
        eprintln!("  Ferro reference contains:");
        eprintln!(
            "    {} transcripts (NM_, NR_, XM_, XR_)",
            self.transcripts.len()
        );
        eprintln!("    {} supplemental sequences", self.supplemental.len());
        eprintln!("    {} genomic sequences (NC_)", self.genomic.len());
        eprintln!("    {} RefSeqGene sequences (NG_)", self.refseqgene.len());
        eprintln!("    {} LRG sequences", self.lrg.len());
        eprintln!("    {} total unique accessions", self.total_count());
    }
}

/// Load accessions from all .fai files in a directory.
fn load_accessions_from_dir(dir: &Path) -> Result<HashSet<String>, FerroError> {
    let mut accessions = HashSet::new();

    if !dir.exists() {
        return Ok(accessions);
    }

    for entry in std::fs::read_dir(dir).map_err(|e| FerroError::Io {
        msg: format!("Failed to read directory {}: {}", dir.display(), e),
    })? {
        let path = entry
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to read entry: {}", e),
            })?
            .path();

        if path.extension().is_some_and(|e| e == "fai") {
            let file = std::fs::File::open(&path).map_err(|e| FerroError::Io {
                msg: format!("Failed to open {}: {}", path.display(), e),
            })?;

            for line in BufReader::new(file).lines() {
                let line = line.map_err(|e| FerroError::Io {
                    msg: format!("Failed to read line: {}", e),
                })?;
                if let Some(acc) = line.split('\t').next() {
                    accessions.insert(acc.to_string());
                }
            }
        }
    }

    Ok(accessions)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    #[test]
    fn test_load_accessions_from_dir() {
        let dir = TempDir::new().unwrap();

        // Create a mock FAI file
        let fai_path = dir.path().join("test.fna.fai");
        let mut f = std::fs::File::create(&fai_path).unwrap();
        writeln!(f, "NM_000001.1\t1000\t0\t80\t81").unwrap();
        writeln!(f, "NM_000002.2\t2000\t1000\t80\t81").unwrap();

        let accessions = load_accessions_from_dir(dir.path()).unwrap();
        assert_eq!(accessions.len(), 2);
        assert!(accessions.contains("NM_000001.1"));
        assert!(accessions.contains("NM_000002.2"));
    }

    #[test]
    fn test_accession_sources_contains() {
        let mut sources = AccessionSources::default();
        sources.transcripts.insert("NM_000001.1".to_string());
        sources.supplemental.insert("U31929.1".to_string());
        sources.genomic.insert("NC_000001.11".to_string());

        assert!(sources.contains("NM_000001.1"));
        assert!(sources.contains("U31929.1"));
        assert!(sources.contains("NC_000001.11"));
        assert!(!sources.contains("NM_999999.1"));
    }

    #[test]
    fn test_find_missing() {
        let mut sources = AccessionSources::default();
        sources.transcripts.insert("NM_000001.1".to_string());
        sources.transcripts.insert("NM_000002.2".to_string());

        let to_check = vec![
            "NM_000001.1".to_string(),
            "NM_000002.2".to_string(),
            "NM_000003.3".to_string(),
        ];
        let missing = sources.find_missing(&to_check);
        assert_eq!(missing.len(), 1);
        assert_eq!(missing[0], "NM_000003.3");
    }

    #[test]
    fn test_total_count_deduplicates() {
        let mut sources = AccessionSources::default();
        sources.transcripts.insert("NM_000001.1".to_string());
        sources.supplemental.insert("NM_000001.1".to_string()); // Same as transcripts
        sources.supplemental.insert("U31929.1".to_string());

        // Should deduplicate NM_000001.1
        assert_eq!(sources.total_count(), 2);
    }
}
