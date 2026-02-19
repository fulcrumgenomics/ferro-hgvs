//! Extract HGVS patterns from various source formats.

use crate::benchmark::types::DatasetFormat;
use crate::FerroError;
use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Extract HGVS patterns from a ClinVar TSV file.
///
/// The ClinVar hgvs4variation.txt.gz file has the following columns:
/// - Column 6 (0-indexed): NucleotideExpression
/// - Column 8 (0-indexed): ProteinExpression
///
/// Both columns may contain HGVS expressions separated by '|'.
pub fn extract_clinvar<P: AsRef<Path>>(input: P, output: P) -> Result<usize, FerroError> {
    let input = input.as_ref();
    let output = output.as_ref();

    // Open input (handle gzip)
    let file = File::open(input).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", input.display(), e),
    })?;

    let reader: Box<dyn BufRead> = if input.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    // Create output directory if needed
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory {}: {}", parent.display(), e),
        })?;
    }

    let out_file = File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output.display(), e),
    })?;
    let mut writer = BufWriter::new(out_file);

    let mut patterns_seen = HashSet::new();
    let mut pattern_count = 0usize;
    let mut line_count = 0usize;

    // Columns to extract (0-indexed)
    let nucleotide_column = 6;
    let protein_column = 8;

    // Progress bar
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Error reading line: {}", e),
        })?;

        // Skip header lines
        if line.starts_with('#') {
            continue;
        }

        line_count += 1;

        let parts: Vec<&str> = line.split('\t').collect();

        // Extract from both nucleotide and protein columns
        for &column in &[nucleotide_column, protein_column] {
            if let Some(field) = parts.get(column) {
                for pattern in field.split('|') {
                    let pattern = pattern.trim();
                    if pattern.is_empty() || pattern == "-" {
                        continue;
                    }

                    // Basic validation: should contain : and .
                    if pattern.contains(':')
                        && pattern.contains('.')
                        && !patterns_seen.contains(pattern)
                    {
                        patterns_seen.insert(pattern.to_string());
                        writeln!(writer, "{}", pattern).map_err(|e| FerroError::Io {
                            msg: format!("Error writing: {}", e),
                        })?;
                        pattern_count += 1;
                    }
                }
            }
        }

        if line_count.is_multiple_of(100_000) {
            pb.set_message(format!(
                "Processed {:>10} lines, found {:>10} unique patterns",
                line_count, pattern_count
            ));
        }
    }

    pb.finish_with_message(format!(
        "Extracted {} unique patterns from {} lines",
        pattern_count, line_count
    ));

    Ok(pattern_count)
}

/// Extract HGVS patterns from a JSON file.
///
/// Supports multiple formats:
/// - test_cases_json: `{"test_cases": [{"input": "..."}, ...]}`
/// - json_array: `["pattern1", "pattern2", ...]` or `[{"hgvs": "..."}, ...]`
pub fn extract_json<P: AsRef<Path>>(
    input: P,
    output: P,
    format: DatasetFormat,
    field_name: &str,
) -> Result<usize, FerroError> {
    let input = input.as_ref();
    let output = output.as_ref();

    // Open and read JSON (handle gzip)
    let file = File::open(input).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", input.display(), e),
    })?;

    let data: serde_json::Value = if input.extension().is_some_and(|e| e == "gz") {
        let reader = BufReader::new(GzDecoder::new(file));
        serde_json::from_reader(reader).map_err(|e| FerroError::Json {
            msg: format!("Failed to parse JSON: {}", e),
        })?
    } else {
        let reader = BufReader::new(file);
        serde_json::from_reader(reader).map_err(|e| FerroError::Json {
            msg: format!("Failed to parse JSON: {}", e),
        })?
    };

    // Create output directory if needed
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory {}: {}", parent.display(), e),
        })?;
    }

    let out_file = File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output.display(), e),
    })?;
    let mut writer = BufWriter::new(out_file);

    let mut patterns = Vec::new();

    match format {
        DatasetFormat::TestCasesJson => {
            // Format: {"test_cases": [{"input": "...", ...}, ...]}
            if let Some(test_cases) = data.get("test_cases").and_then(|v| v.as_array()) {
                for item in test_cases {
                    if let Some(pattern) = item.get(field_name).and_then(|v| v.as_str()) {
                        patterns.push(pattern.to_string());
                    }
                }
            }
        }
        DatasetFormat::JsonArray => {
            // Format: ["pattern1", ...] or [{"field": "pattern1"}, ...]
            if let Some(arr) = data.as_array() {
                for item in arr {
                    if let Some(s) = item.as_str() {
                        patterns.push(s.to_string());
                    } else if let Some(pattern) = item.get(field_name).and_then(|v| v.as_str()) {
                        patterns.push(pattern.to_string());
                    }
                }
            }
        }
        _ => {
            return Err(FerroError::Io {
                msg: format!("Unsupported format for JSON extraction: {:?}", format),
            });
        }
    }

    for pattern in &patterns {
        writeln!(writer, "{}", pattern.trim()).map_err(|e| FerroError::Io {
            msg: format!("Error writing: {}", e),
        })?;
    }

    eprintln!(
        "Extracted {} patterns from {}",
        patterns.len(),
        input.display()
    );

    Ok(patterns.len())
}

/// Copy a plain text file (one pattern per line).
pub fn copy_text<P: AsRef<Path>>(input: P, output: P) -> Result<usize, FerroError> {
    let input = input.as_ref();
    let output = output.as_ref();

    // Create output directory if needed
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory {}: {}", parent.display(), e),
        })?;
    }

    let in_file = File::open(input).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", input.display(), e),
    })?;
    let out_file = File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output.display(), e),
    })?;

    let reader = BufReader::new(in_file);
    let mut writer = BufWriter::new(out_file);

    let mut count = 0usize;
    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Error reading: {}", e),
        })?;
        let trimmed = line.trim();
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            writeln!(writer, "{}", trimmed).map_err(|e| FerroError::Io {
                msg: format!("Error writing: {}", e),
            })?;
            count += 1;
        }
    }

    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    #[test]
    fn test_extract_text() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("input.txt");
        let output = dir.path().join("output.txt");

        // Create input file
        let mut f = File::create(&input).unwrap();
        writeln!(f, "NC_000001.11:g.12345A>G").unwrap();
        writeln!(f, "# comment").unwrap();
        writeln!(f, "NM_000088.3:c.589G>T").unwrap();
        writeln!(f).unwrap(); // Empty line
        writeln!(f, "NP_000079.2:p.Val200Met").unwrap();

        let count = copy_text(&input, &output).unwrap();
        assert_eq!(count, 3);

        // Verify output
        let content = std::fs::read_to_string(&output).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0], "NC_000001.11:g.12345A>G");
        assert_eq!(lines[1], "NM_000088.3:c.589G>T");
        assert_eq!(lines[2], "NP_000079.2:p.Val200Met");
    }
}
