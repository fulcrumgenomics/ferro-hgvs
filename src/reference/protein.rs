//! Protein sequence provider
//!
//! Fetches protein sequences from NCBI E-utilities for reference validation.
//!
//! This module requires the `protein-fetch` feature to be enabled for
//! fetching sequences from NCBI.

use crate::error::FerroError;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::RwLock;

/// A protein sequence cache that stores sequences locally
///
/// This provider fetches protein sequences from NCBI E-utilities and caches
/// them locally for faster subsequent access.
pub struct ProteinCache {
    /// Cache directory for storing sequences
    cache_dir: PathBuf,
    /// In-memory cache for frequently accessed sequences
    memory_cache: RwLock<HashMap<String, String>>,
    /// Maximum size of in-memory cache
    max_memory_entries: usize,
}

impl ProteinCache {
    /// Create a new protein cache with the specified cache directory
    pub fn new(cache_dir: impl AsRef<Path>) -> Result<Self, FerroError> {
        let cache_dir = cache_dir.as_ref().to_path_buf();
        fs::create_dir_all(&cache_dir)?;

        Ok(Self {
            cache_dir,
            memory_cache: RwLock::new(HashMap::new()),
            max_memory_entries: 10000,
        })
    }

    /// Get a protein sequence by accession
    ///
    /// First checks the in-memory cache, then disk cache, then fetches from NCBI.
    pub fn get_sequence(&self, accession: &str) -> Result<String, FerroError> {
        // Check memory cache first
        {
            let cache = self.memory_cache.read().unwrap();
            if let Some(seq) = cache.get(accession) {
                return Ok(seq.clone());
            }
        }

        // Check disk cache
        if let Ok(seq) = self.load_from_disk(accession) {
            // Add to memory cache
            self.add_to_memory_cache(accession, &seq);
            return Ok(seq);
        }

        // Fetch from NCBI
        let seq = self.fetch_from_ncbi(accession)?;

        // Cache to disk and memory
        self.save_to_disk(accession, &seq)?;
        self.add_to_memory_cache(accession, &seq);

        Ok(seq)
    }

    /// Get a subsequence of a protein
    pub fn get_subsequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        let seq = self.get_sequence(accession)?;
        let start = start as usize;
        let end = end as usize;

        if start >= seq.len() || end > seq.len() || start > end {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "Position {}:{}-{} out of bounds for protein {} (length {})",
                    accession,
                    start,
                    end,
                    accession,
                    seq.len()
                ),
            });
        }

        Ok(seq[start..end].to_string())
    }

    /// Add a sequence to the memory cache
    fn add_to_memory_cache(&self, accession: &str, sequence: &str) {
        let mut cache = self.memory_cache.write().unwrap();
        if cache.len() >= self.max_memory_entries {
            // Simple eviction: clear half the cache
            let to_remove: Vec<_> = cache
                .keys()
                .take(self.max_memory_entries / 2)
                .cloned()
                .collect();
            for key in to_remove {
                cache.remove(&key);
            }
        }
        cache.insert(accession.to_string(), sequence.to_string());
    }

    /// Load sequence from disk cache
    fn load_from_disk(&self, accession: &str) -> Result<String, FerroError> {
        let path = self.cache_path(accession);
        let file = File::open(&path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open cache file {:?}: {}", path, e),
        })?;
        let reader = BufReader::new(file);

        let mut sequence = String::new();
        for line in reader.lines() {
            let line = line?;
            if !line.starts_with('>') {
                sequence.push_str(line.trim());
            }
        }

        if sequence.is_empty() {
            return Err(FerroError::ProteinReferenceNotAvailable {
                accession: accession.to_string(),
                start: 0,
                end: 0,
            });
        }

        Ok(sequence)
    }

    /// Save sequence to disk cache
    fn save_to_disk(&self, accession: &str, sequence: &str) -> Result<(), FerroError> {
        let path = self.cache_path(accession);
        let file = File::create(&path)?;
        let mut writer = BufWriter::new(file);

        // Write in FASTA format
        writeln!(writer, ">{}", accession)?;
        // Write sequence in 60-char lines
        for chunk in sequence.as_bytes().chunks(60) {
            writeln!(writer, "{}", std::str::from_utf8(chunk).unwrap_or(""))?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Get the cache file path for an accession
    fn cache_path(&self, accession: &str) -> PathBuf {
        // Sanitize accession for filesystem
        let safe_name = accession.replace(['/', '\\', ':', '*', '?', '"', '<', '>', '|'], "_");
        self.cache_dir.join(format!("{}.fa", safe_name))
    }

    /// Fetch protein sequence from NCBI E-utilities
    #[cfg(feature = "protein-fetch")]
    fn fetch_from_ncbi(&self, accession: &str) -> Result<String, FerroError> {
        // Use efetch to get FASTA format
        let url = format!(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={}&rettype=fasta&retmode=text",
            accession
        );

        // Use ureq for HTTP requests (synchronous, simpler than reqwest)
        let response = ureq::get(&url).call().map_err(|e| FerroError::Io {
            msg: format!("Failed to fetch protein {} from NCBI: {}", accession, e),
        })?;

        let body = response.into_string().map_err(|e| FerroError::Io {
            msg: format!("Failed to read response for {}: {}", accession, e),
        })?;

        // Parse FASTA
        let mut sequence = String::new();
        for line in body.lines() {
            if !line.starts_with('>') && !line.is_empty() {
                sequence.push_str(line.trim());
            }
        }

        if sequence.is_empty() {
            return Err(FerroError::ProteinReferenceNotAvailable {
                accession: accession.to_string(),
                start: 0,
                end: 0,
            });
        }

        Ok(sequence)
    }

    /// Fetch protein sequence from NCBI E-utilities (stub when feature not enabled)
    #[cfg(not(feature = "protein-fetch"))]
    fn fetch_from_ncbi(&self, accession: &str) -> Result<String, FerroError> {
        Err(FerroError::Io {
            msg: format!(
                "Cannot fetch protein {} from NCBI: protein-fetch feature not enabled",
                accession
            ),
        })
    }

    /// Pre-fetch multiple protein sequences in parallel
    ///
    /// Returns the number of sequences successfully fetched.
    pub fn prefetch(&self, accessions: &[&str]) -> usize {
        let mut success_count = 0;
        for accession in accessions {
            if self.get_sequence(accession).is_ok() {
                success_count += 1;
            }
        }
        success_count
    }

    /// Check if a sequence is cached
    pub fn is_cached(&self, accession: &str) -> bool {
        self.cache_path(accession).exists()
    }

    /// Get the number of cached sequences
    pub fn cached_count(&self) -> usize {
        fs::read_dir(&self.cache_dir)
            .map(|entries| entries.filter_map(Result::ok).count())
            .unwrap_or(0)
    }

    /// Clear the in-memory cache
    pub fn clear_memory_cache(&self) {
        let mut cache = self.memory_cache.write().unwrap();
        cache.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_protein_cache_new() {
        let temp_dir = TempDir::new().unwrap();
        let cache = ProteinCache::new(temp_dir.path()).unwrap();
        assert_eq!(cache.cached_count(), 0);
    }

    #[test]
    fn test_cache_path_sanitization() {
        let temp_dir = TempDir::new().unwrap();
        let cache = ProteinCache::new(temp_dir.path()).unwrap();

        // Accession with version
        let path = cache.cache_path("NP_000079.2");
        assert!(path.to_string_lossy().contains("NP_000079.2.fa"));

        // Accession with special chars (shouldn't happen but test sanitization)
        let path = cache.cache_path("test/acc:1");
        // The filename should not contain special chars - check file_name() only
        let filename = path.file_name().unwrap().to_string_lossy();
        assert!(!filename.contains('/'));
        assert!(!filename.contains(':'));
        assert!(filename.contains("test_acc_1.fa"));
    }

    #[test]
    fn test_save_and_load_sequence() {
        let temp_dir = TempDir::new().unwrap();
        let cache = ProteinCache::new(temp_dir.path()).unwrap();

        let accession = "TEST_PROTEIN.1";
        let sequence = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH";

        // Save
        cache.save_to_disk(accession, sequence).unwrap();
        assert!(cache.is_cached(accession));

        // Load
        let loaded = cache.load_from_disk(accession).unwrap();
        assert_eq!(loaded, sequence);
    }

    #[test]
    fn test_get_subsequence() {
        let temp_dir = TempDir::new().unwrap();
        let cache = ProteinCache::new(temp_dir.path()).unwrap();

        let accession = "TEST_PROTEIN.1";
        let sequence = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH";
        cache.save_to_disk(accession, sequence).unwrap();

        // Get subsequence
        let subseq = cache.get_subsequence(accession, 0, 5).unwrap();
        assert_eq!(subseq, "MVLSP");

        let subseq = cache.get_subsequence(accession, 10, 20).unwrap();
        assert_eq!(subseq, "VKAAWGKVGA");
    }

    #[test]
    fn test_get_subsequence_out_of_bounds() {
        let temp_dir = TempDir::new().unwrap();
        let cache = ProteinCache::new(temp_dir.path()).unwrap();

        let accession = "TEST_PROTEIN.1";
        let sequence = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH";
        cache.save_to_disk(accession, sequence).unwrap();

        // Out of bounds
        let result = cache.get_subsequence(accession, 0, 1000);
        assert!(result.is_err());
    }

    #[test]
    fn test_memory_cache() {
        let temp_dir = TempDir::new().unwrap();
        let cache = ProteinCache::new(temp_dir.path()).unwrap();

        let accession = "TEST_PROTEIN.1";
        let sequence = "MVLSPADKTN";
        cache.save_to_disk(accession, sequence).unwrap();

        // First access loads from disk and caches in memory
        let _ = cache.get_sequence(accession).unwrap();

        // Verify in memory cache
        {
            let mem_cache = cache.memory_cache.read().unwrap();
            assert!(mem_cache.contains_key(accession));
        }

        // Clear memory cache
        cache.clear_memory_cache();

        {
            let mem_cache = cache.memory_cache.read().unwrap();
            assert!(mem_cache.is_empty());
        }
    }
}
