//! Comprehensive external pattern tests from cloned repositories
//!
//! Extracts and tests patterns from:
//! - biocommons/hgvs-eval
//! - openvar/vv_hgvs
//! - openvar/variantValidator

use ferro_hgvs::parse_hgvs;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

/// Extract HGVS patterns from a TSV file with HGVSg, HGVSc, HGVSp columns
fn extract_hgvs_from_gcp_tsv(path: &Path) -> Vec<String> {
    let mut patterns = Vec::new();

    if let Ok(content) = fs::read_to_string(path) {
        for line in content.lines().skip(1) {
            // Skip header
            let parts: Vec<&str> = line.split('\t').collect();
            // Look for HGVS columns (typically columns 1, 2, 3 for HGVSg, HGVSc, HGVSp)
            for part in parts.iter().skip(1).take(3) {
                let part = part.trim();
                if !part.is_empty()
                    && (part.contains(':')
                        && (part.starts_with("NC_")
                            || part.starts_with("NM_")
                            || part.starts_with("NP_")
                            || part.starts_with("NG_")
                            || part.starts_with("NR_")
                            || part.starts_with("LRG_")
                            || part.starts_with("ENST")))
                {
                    // Filter out error messages and malformed patterns
                    if !is_error_message(part) && !is_malformed_pattern(part) {
                        patterns.push(part.to_string());
                    }
                }
            }
        }
    }

    patterns
}

/// Extract HGVS patterns from hgvs-eval TSV (input and output columns)
fn extract_hgvs_eval_patterns(path: &Path) -> Vec<String> {
    let mut patterns = Vec::new();

    if let Ok(content) = fs::read_to_string(path) {
        for line in content.lines().skip(1) {
            // Skip header
            let parts: Vec<&str> = line.split('\t').collect();

            // Column 3 is input, columns 4-5 are output_accepted/output_preferred
            for &col in &[3, 4, 5] {
                if let Some(part) = parts.get(col) {
                    // Split on | for multiple accepted outputs
                    for variant in part.split('|') {
                        let variant = variant.trim();
                        if !variant.is_empty()
                            && variant.contains(':')
                            && (variant.starts_with("NC_")
                                || variant.starts_with("NM_")
                                || variant.starts_with("NP_")
                                || variant.starts_with("NG_")
                                || variant.starts_with("NR_")
                                || variant.starts_with("LRG_")
                                || variant.starts_with("ENST"))
                        {
                            // Filter out error messages and malformed patterns
                            if !is_error_message(variant) && !is_malformed_pattern(variant) {
                                patterns.push(variant.to_string());
                            }
                        }
                    }
                }
            }
        }
    }

    patterns
}

/// Check if a string looks like an error message rather than a valid HGVS pattern
fn is_error_message(s: &str) -> bool {
    // Common phrases found in error messages extracted from Python test files
    let error_phrases = [
        "may also be written as",
        "updated to",
        "does not agree with reference sequence",
        "cannot be mapped directly",
        "automapped to equivalent",
        "auto-mapped to",
        "is not associated with",
        "is not part of genome build",
        "ExonBoundaryError",
        "HGVSInvalidVariantError",
        "HGVSParseError",
        "HGVSDataNotAvailableError",
        "HGVSInvalidIntervalError",
        "HGVSError",
        "Normalization of",
        "is not in the expected",
        "boundary error",
        "does not exist on",
        "Unable to",
        "not supported",
        "Invalid",
        "where int requires",            // Description text
        "wrong",                         // Test case with appended "wrong"
        "automapped to genome position", // Mapping messages
    ];

    error_phrases.iter().any(|phrase| s.contains(phrase))
}

/// Check if a pattern is malformed and should be excluded from valid pattern tests
fn is_malformed_pattern(s: &str) -> bool {
    // Control characters (except normal whitespace)
    if s.chars().any(|c| c.is_control() && c != ' ' && c != '\t') {
        return true;
    }

    // Trailing quotes or other punctuation that shouldn't be there
    if s.ends_with('"') || s.ends_with('\'') {
        return true;
    }

    // Double colons (malformed accession separator)
    if s.contains("::") {
        return true;
    }

    // Whitespace in accession (e.g., "NM_015120  .4")
    if s.contains("  ") || s.contains("\t") {
        return true;
    }

    // Parentheses around REMOVE or other text in accession
    if s.contains("(REMOVE)") {
        return true;
    }

    // Invalid bases in sequence (X, I are not valid IUPAC codes for DNA)
    // But be careful - X is valid in protein (Xaa) and some contexts
    // Check for patterns like GGI[, CXXX[ which are clearly invalid
    if s.contains("GGI[") || s.contains("XXX[") {
        return true;
    }

    // Unclosed brackets
    let open_brackets = s.chars().filter(|&c| c == '[').count();
    let close_brackets = s.chars().filter(|&c| c == ']').count();
    if open_brackets != close_brackets {
        return true;
    }

    // Wrong case for edit types (DUp, INSC, inSC should be dup, ins)
    if s.contains("DUp") || s.contains("INSC") || s.contains("inSC") {
        return true;
    }

    // "con" is not a valid edit type
    if s.ends_with("con") || s.contains("con;") || s.contains("con]") {
        return true;
    }

    // Protein position without amino acid (p.34= is invalid)
    if s.contains(":p.") {
        // Check if p. is followed by a digit (invalid - needs amino acid)
        if let Some(idx) = s.find(":p.") {
            let after_p = &s[idx + 3..];
            if after_p.chars().next().is_some_and(|c| c.is_ascii_digit()) {
                return true;
            }
        }
    }

    // Empty delins (c.34_34delins with nothing after)
    if s.ends_with("delins") {
        return true;
    }

    // Incomplete variant (ends with offset but no edit, e.g., c.22+6)
    if let Some(idx) = s.rfind('+') {
        let after_plus = &s[idx + 1..];
        // If it's just digits after the +, it's incomplete
        if !after_plus.is_empty() && after_plus.chars().all(|c| c.is_ascii_digit()) {
            return true;
        }
    }

    false
}

/// Extract HGVS patterns from Python test files using regex-like pattern matching
fn extract_patterns_from_python(path: &Path) -> Vec<String> {
    let mut patterns = Vec::new();

    if let Ok(content) = fs::read_to_string(path) {
        // Look for patterns like 'NM_...', "NC_...", etc.
        let accession_prefixes = ["NC_", "NM_", "NP_", "NG_", "NR_", "LRG_", "ENST"];

        for line in content.lines() {
            // Skip lines that are clearly error/exception handling
            if line.contains("raise")
                || line.contains("Exception")
                || line.contains("assert") && line.contains("Error")
            {
                continue;
            }

            // Find quoted strings that look like HGVS variants
            let mut in_quote = false;
            let mut quote_char = ' ';
            let mut current = String::new();

            for c in line.chars() {
                if (c == '\'' || c == '"') && !in_quote {
                    in_quote = true;
                    quote_char = c;
                    current.clear();
                } else if c == quote_char && in_quote {
                    in_quote = false;
                    // Check if this looks like an HGVS variant
                    if current.contains(':')
                        && accession_prefixes.iter().any(|p| current.starts_with(p))
                    {
                        // Basic validation - should have coordinate type
                        if current.contains(":c.")
                            || current.contains(":g.")
                            || current.contains(":p.")
                            || current.contains(":n.")
                            || current.contains(":m.")
                            || current.contains(":r.")
                        {
                            // Filter out error messages and malformed patterns
                            if !is_error_message(&current) && !is_malformed_pattern(&current) {
                                patterns.push(current.clone());
                            }
                        }
                    }
                } else if in_quote {
                    current.push(c);
                }
            }
        }
    }

    patterns
}

/// Recursively find all TSV files in a directory
fn find_tsv_files(dir: &Path) -> Vec<std::path::PathBuf> {
    let mut files = Vec::new();

    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                files.extend(find_tsv_files(&path));
            } else if path.extension().is_some_and(|e| e == "tsv") {
                files.push(path);
            }
        }
    }

    files
}

/// Recursively find all Python test files
fn find_python_test_files(dir: &Path) -> Vec<std::path::PathBuf> {
    let mut files = Vec::new();

    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                files.extend(find_python_test_files(&path));
            } else if path.extension().is_some_and(|e| e == "py") {
                if let Some(name) = path.file_name() {
                    if name.to_string_lossy().starts_with("test_") {
                        files.push(path);
                    }
                }
            }
        }
    }

    files
}

#[test]
fn test_vv_hgvs_patterns() {
    let base_path = Path::new("external-repos/vv_hgvs/tests/data");

    if !base_path.exists() {
        println!("Skipping vv_hgvs tests - repo not cloned");
        return;
    }

    println!("\n=== vv_hgvs Pattern Extraction ===\n");

    let tsv_files = find_tsv_files(base_path);
    let mut all_patterns: HashSet<String> = HashSet::new();
    let mut by_file: HashMap<String, usize> = HashMap::new();

    for file in &tsv_files {
        let patterns = extract_hgvs_from_gcp_tsv(file);
        let count = patterns.len();
        if count > 0 {
            by_file.insert(
                file.file_name().unwrap().to_string_lossy().to_string(),
                count,
            );
        }
        all_patterns.extend(patterns);
    }

    println!(
        "Found {} unique patterns from {} TSV files:",
        all_patterns.len(),
        tsv_files.len()
    );
    for (file, count) in by_file.iter() {
        println!("  {}: {} patterns", file, count);
    }

    // Test all patterns
    let mut supported = 0;
    let mut unsupported = Vec::new();

    for pattern in &all_patterns {
        if parse_hgvs(pattern).is_ok() {
            supported += 1;
        } else {
            unsupported.push(pattern.clone());
        }
    }

    let total = all_patterns.len();
    let rate = if total > 0 {
        supported as f64 / total as f64 * 100.0
    } else {
        0.0
    };

    println!("\nResults:");
    println!("  Supported: {}/{} ({:.1}%)", supported, total, rate);
    println!("  Unsupported: {}", unsupported.len());

    if !unsupported.is_empty() && unsupported.len() <= 50 {
        println!("\nUnsupported patterns:");
        for p in &unsupported {
            println!("  {}", p);
        }
    } else if unsupported.len() > 50 {
        println!("\nFirst 50 unsupported patterns:");
        for p in unsupported.iter().take(50) {
            println!("  {}", p);
        }
    }
}

#[test]
fn test_hgvs_eval_patterns() {
    let eval_file =
        Path::new("external-repos/hgvs-eval/hgvseval/tests/data/hgvs_eval_criteria_tests.tsv");

    if !eval_file.exists() {
        println!("Skipping hgvs-eval tests - repo not cloned");
        return;
    }

    println!("\n=== hgvs-eval Pattern Extraction ===\n");

    let patterns: HashSet<String> = extract_hgvs_eval_patterns(eval_file).into_iter().collect();

    println!("Found {} unique patterns from hgvs-eval", patterns.len());

    let mut supported = 0;
    let mut unsupported = Vec::new();

    for pattern in &patterns {
        if parse_hgvs(pattern).is_ok() {
            supported += 1;
        } else {
            unsupported.push(pattern.clone());
        }
    }

    let total = patterns.len();
    let rate = if total > 0 {
        supported as f64 / total as f64 * 100.0
    } else {
        0.0
    };

    println!("\nResults:");
    println!("  Supported: {}/{} ({:.1}%)", supported, total, rate);

    if !unsupported.is_empty() {
        println!("\nUnsupported patterns:");
        for p in &unsupported {
            println!("  {}", p);
        }
    }
}

#[test]
fn test_variant_validator_patterns() {
    let tests_path = Path::new("external-repos/variantValidator/tests");

    if !tests_path.exists() {
        println!("Skipping variantValidator tests - repo not cloned");
        return;
    }

    println!("\n=== variantValidator Pattern Extraction ===\n");

    let py_files = find_python_test_files(tests_path);
    let mut all_patterns: HashSet<String> = HashSet::new();

    for file in &py_files {
        let patterns = extract_patterns_from_python(file);
        all_patterns.extend(patterns);
    }

    println!(
        "Found {} unique patterns from {} Python test files",
        all_patterns.len(),
        py_files.len()
    );

    let mut supported = 0;
    let mut unsupported = Vec::new();

    for pattern in &all_patterns {
        if parse_hgvs(pattern).is_ok() {
            supported += 1;
        } else {
            unsupported.push(pattern.clone());
        }
    }

    let total = all_patterns.len();
    let rate = if total > 0 {
        supported as f64 / total as f64 * 100.0
    } else {
        0.0
    };

    println!("\nResults:");
    println!("  Supported: {}/{} ({:.1}%)", supported, total, rate);

    if !unsupported.is_empty() && unsupported.len() <= 30 {
        println!("\nUnsupported patterns:");
        for p in &unsupported {
            println!("  {}", p);
        }
    }
}

#[test]
fn test_comprehensive_external_summary() {
    println!("\n============================================================");
    println!("     COMPREHENSIVE EXTERNAL PATTERN TEST SUMMARY");
    println!("============================================================\n");

    let mut total_patterns = 0;
    let mut total_supported = 0;
    let mut all_unsupported: Vec<String> = Vec::new();

    // vv_hgvs
    let vv_path = Path::new("external-repos/vv_hgvs/tests/data");
    if vv_path.exists() {
        let tsv_files = find_tsv_files(vv_path);
        let patterns: HashSet<String> = tsv_files
            .iter()
            .flat_map(|f| extract_hgvs_from_gcp_tsv(f))
            .collect();

        let supported = patterns.iter().filter(|p| parse_hgvs(p).is_ok()).count();
        let unsupported: Vec<_> = patterns
            .iter()
            .filter(|p| parse_hgvs(p).is_err())
            .cloned()
            .collect();

        println!("vv_hgvs TSV files:");
        println!("  Patterns: {}", patterns.len());
        println!(
            "  Supported: {} ({:.1}%)",
            supported,
            supported as f64 / patterns.len() as f64 * 100.0
        );

        total_patterns += patterns.len();
        total_supported += supported;
        all_unsupported.extend(unsupported);
    }

    // hgvs-eval
    let eval_file =
        Path::new("external-repos/hgvs-eval/hgvseval/tests/data/hgvs_eval_criteria_tests.tsv");
    if eval_file.exists() {
        let patterns: HashSet<String> = extract_hgvs_eval_patterns(eval_file).into_iter().collect();
        let supported = patterns.iter().filter(|p| parse_hgvs(p).is_ok()).count();
        let unsupported: Vec<_> = patterns
            .iter()
            .filter(|p| parse_hgvs(p).is_err())
            .cloned()
            .collect();

        println!("\nhgvs-eval:");
        println!("  Patterns: {}", patterns.len());
        println!(
            "  Supported: {} ({:.1}%)",
            supported,
            supported as f64 / patterns.len() as f64 * 100.0
        );

        total_patterns += patterns.len();
        total_supported += supported;
        all_unsupported.extend(unsupported);
    }

    // variantValidator
    let vv_tests = Path::new("external-repos/variantValidator/tests");
    if vv_tests.exists() {
        let py_files = find_python_test_files(vv_tests);
        let patterns: HashSet<String> = py_files
            .iter()
            .flat_map(|f| extract_patterns_from_python(f))
            .collect();

        let supported = patterns.iter().filter(|p| parse_hgvs(p).is_ok()).count();
        let unsupported: Vec<_> = patterns
            .iter()
            .filter(|p| parse_hgvs(p).is_err())
            .cloned()
            .collect();

        println!("\nvariantValidator Python tests:");
        println!("  Patterns: {}", patterns.len());
        println!(
            "  Supported: {} ({:.1}%)",
            supported,
            if patterns.is_empty() {
                0.0
            } else {
                supported as f64 / patterns.len() as f64 * 100.0
            }
        );

        total_patterns += patterns.len();
        total_supported += supported;
        all_unsupported.extend(unsupported);
    }

    // Deduplicate unsupported patterns
    let unique_unsupported: HashSet<String> = all_unsupported.into_iter().collect();

    println!("\n------------------------------------------------------------");
    println!("TOTALS (new sources only):");
    println!("  Total patterns: {}", total_patterns);
    println!("  Supported: {}", total_supported);
    println!("  Unsupported: {}", unique_unsupported.len());
    if total_patterns > 0 {
        println!(
            "  Support rate: {:.1}%",
            total_supported as f64 / total_patterns as f64 * 100.0
        );
    }

    if !unique_unsupported.is_empty() {
        println!("\n------------------------------------------------------------");
        println!(
            "UNIQUE UNSUPPORTED PATTERNS ({}):",
            unique_unsupported.len()
        );
        let mut sorted: Vec<_> = unique_unsupported.into_iter().collect();
        sorted.sort();
        for p in sorted.iter().take(100) {
            println!("  {}", p);
        }
        if sorted.len() > 100 {
            println!("  ... and {} more", sorted.len() - 100);
        }
    }

    println!("\n============================================================");
}
