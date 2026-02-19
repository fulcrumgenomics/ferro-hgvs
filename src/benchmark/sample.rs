//! Stratified sampling for benchmark datasets.

use crate::benchmark::types::PatternCategory;
use crate::FerroError;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Create a stratified sample that maintains the distribution of pattern types.
///
/// This is important for normalization benchmarks where we want to test
/// a representative subset without processing millions of patterns through
/// a slow API.
///
/// If `exclude_protein` is true, protein patterns (NP_*, :p.) are filtered out
/// before sampling. This is recommended for normalization benchmarks since
/// protein normalization is handled inconsistently across tools.
pub fn stratified_sample<P: AsRef<Path>>(
    input: P,
    output: P,
    sample_size: usize,
    seed: u64,
    exclude_protein: bool,
) -> Result<SampleStats, FerroError> {
    let input = input.as_ref();
    let output = output.as_ref();

    // First pass: categorize all patterns
    let file = File::open(input).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", input.display(), e),
    })?;
    let reader = BufReader::new(file);

    let mut by_category: HashMap<PatternCategory, Vec<String>> = HashMap::new();
    let mut total = 0usize;
    let mut protein_excluded = 0usize;

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Error reading: {}", e),
        })?;

        let pattern = line.trim();
        if pattern.is_empty() {
            continue;
        }

        let category = PatternCategory::categorize(pattern);

        // Skip protein patterns if requested
        if exclude_protein && category.is_protein() {
            protein_excluded += 1;
            continue;
        }

        by_category
            .entry(category)
            .or_default()
            .push(pattern.to_string());
        total += 1;
    }

    if exclude_protein && protein_excluded > 0 {
        eprintln!(
            "  Excluded {} protein patterns (--exclude-protein)",
            protein_excluded
        );
    }

    if total == 0 {
        return Err(FerroError::Io {
            msg: "Input file is empty".to_string(),
        });
    }

    // Calculate proportional samples for each category
    let mut sampled = Vec::new();
    let mut stats = SampleStats {
        total_patterns: total,
        sample_size: 0,
        by_category: HashMap::new(),
    };

    // Use a simple PRNG for reproducibility
    let mut rng = SimpleRng::new(seed);

    // Sort categories for deterministic iteration order
    let mut categories: Vec<_> = by_category.into_iter().collect();
    categories.sort_by_key(|(cat, _)| format!("{:?}", cat));

    for (category, mut patterns) in categories {
        let proportion = patterns.len() as f64 / total as f64;
        let cat_sample_size = ((sample_size as f64 * proportion).ceil() as usize).max(1);
        let actual_sample = cat_sample_size.min(patterns.len());

        // Sort patterns within category for deterministic shuffle input
        patterns.sort();

        // Shuffle and take sample
        rng.shuffle(&mut patterns);
        sampled.extend(patterns.into_iter().take(actual_sample));

        stats.by_category.insert(
            category,
            CategoryStats {
                total,
                sampled: actual_sample,
                proportion,
            },
        );
    }

    // Shuffle final sample and trim to exact size
    rng.shuffle(&mut sampled);
    sampled.truncate(sample_size);
    stats.sample_size = sampled.len();

    // Create output directory if needed
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory {}: {}", parent.display(), e),
        })?;
    }

    // Write output
    let out_file = File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output.display(), e),
    })?;
    let mut writer = BufWriter::new(out_file);

    for pattern in &sampled {
        writeln!(writer, "{}", pattern).map_err(|e| FerroError::Io {
            msg: format!("Error writing: {}", e),
        })?;
    }

    eprintln!(
        "Created stratified sample: {} patterns from {} total",
        stats.sample_size, stats.total_patterns
    );

    Ok(stats)
}

/// Create a stratified sample from a vector of patterns.
///
/// This is useful when patterns are already loaded in memory.
///
/// If `exclude_protein` is true, protein patterns (NP_*, :p.) are filtered out
/// before sampling.
pub fn stratified_sample_vec(
    patterns: &[String],
    sample_size: usize,
    seed: u64,
    exclude_protein: bool,
) -> Result<Vec<String>, FerroError> {
    if patterns.is_empty() {
        return Err(FerroError::Io {
            msg: "Input is empty".to_string(),
        });
    }

    // Categorize all patterns, optionally excluding proteins
    let mut by_category: HashMap<PatternCategory, Vec<String>> = HashMap::new();
    let mut protein_excluded = 0usize;

    for pattern in patterns {
        let category = PatternCategory::categorize(pattern);

        if exclude_protein && category.is_protein() {
            protein_excluded += 1;
            continue;
        }

        by_category
            .entry(category)
            .or_default()
            .push(pattern.clone());
    }

    if exclude_protein && protein_excluded > 0 {
        eprintln!("  Excluded {} protein patterns", protein_excluded);
    }

    let total: usize = by_category.values().map(|v| v.len()).sum();

    if total == 0 {
        return Err(FerroError::Io {
            msg: "No patterns remaining after filtering".to_string(),
        });
    }

    // Calculate proportional samples for each category
    let mut sampled = Vec::new();
    let mut rng = SimpleRng::new(seed);

    // Sort categories for deterministic iteration order
    let mut categories: Vec<_> = by_category.into_iter().collect();
    categories.sort_by_key(|(cat, _)| format!("{:?}", cat));

    for (_category, mut cat_patterns) in categories {
        let proportion = cat_patterns.len() as f64 / total as f64;
        let cat_sample_size = ((sample_size as f64 * proportion).ceil() as usize).max(1);
        let actual_sample = cat_sample_size.min(cat_patterns.len());

        // Sort patterns within category for deterministic shuffle input
        cat_patterns.sort();

        // Shuffle and take sample
        rng.shuffle(&mut cat_patterns);
        sampled.extend(cat_patterns.into_iter().take(actual_sample));
    }

    // Shuffle final sample and trim to exact size
    rng.shuffle(&mut sampled);
    sampled.truncate(sample_size);

    Ok(sampled)
}

/// Statistics about a stratified sample.
#[derive(Debug, Clone)]
pub struct SampleStats {
    pub total_patterns: usize,
    pub sample_size: usize,
    pub by_category: HashMap<PatternCategory, CategoryStats>,
}

/// Statistics for a single category.
#[derive(Debug, Clone)]
pub struct CategoryStats {
    pub total: usize,
    pub sampled: usize,
    pub proportion: f64,
}

/// Simple PRNG for reproducible sampling (xorshift64).
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    fn next(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn shuffle<T>(&mut self, slice: &mut [T]) {
        let len = slice.len();
        for i in (1..len).rev() {
            let j = (self.next() as usize) % (i + 1);
            slice.swap(i, j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_stratified_sample() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("input.txt");
        let output = dir.path().join("sample.txt");

        // Create input with mixed pattern types
        let mut f = File::create(&input).unwrap();
        // 5 genomic SNVs
        for i in 0..5 {
            writeln!(f, "NC_000001.11:g.{}A>G", 1000 + i).unwrap();
        }
        // 3 coding SNVs
        for i in 0..3 {
            writeln!(f, "NM_000088.3:c.{}G>T", 100 + i).unwrap();
        }
        // 2 protein subs
        writeln!(f, "NP_000079.2:p.Val200Met").unwrap();
        writeln!(f, "NP_000079.2:p.Arg300Gln").unwrap();

        // Test with proteins included
        let stats = stratified_sample(&input, &output, 5, 42, false).unwrap();

        assert_eq!(stats.total_patterns, 10);
        assert!(stats.sample_size <= 5);

        // Verify output exists and has correct number of lines
        let content = std::fs::read_to_string(&output).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert!(lines.len() <= 5);

        // Test with proteins excluded
        let output_no_protein = dir.path().join("sample_no_protein.txt");
        let stats_no_protein = stratified_sample(&input, &output_no_protein, 5, 42, true).unwrap();

        // Should have 8 patterns (10 - 2 proteins)
        assert_eq!(stats_no_protein.total_patterns, 8);

        // Verify no protein patterns in output
        let content = std::fs::read_to_string(&output_no_protein).unwrap();
        for line in content.lines() {
            assert!(
                !line.contains(":p."),
                "Found protein pattern in exclude_protein output: {}",
                line
            );
        }
    }

    #[test]
    fn test_pattern_categorization() {
        assert_eq!(
            PatternCategory::categorize("NC_000001.11:g.12345A>G"),
            PatternCategory::GenomicSnv
        );
        assert_eq!(
            PatternCategory::categorize("NC_000001.11:g.12345del"),
            PatternCategory::GenomicDel
        );
        assert_eq!(
            PatternCategory::categorize("NM_000088.3:c.589G>T"),
            PatternCategory::CodingSnv
        );
        assert_eq!(
            PatternCategory::categorize("NM_000088.3:c.589+5G>T"),
            PatternCategory::CodingIntronic
        );
        assert_eq!(
            PatternCategory::categorize("NP_000079.2:p.Val200Met"),
            PatternCategory::ProteinSub
        );
        assert_eq!(
            PatternCategory::categorize("NC_012920.1:m.8993T>G"),
            PatternCategory::Mitochondrial
        );
    }
}
