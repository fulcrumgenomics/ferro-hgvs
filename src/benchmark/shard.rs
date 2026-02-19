//! Dataset sharding for parallel processing.

use crate::FerroError;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Shard a dataset into N files using round-robin distribution.
///
/// This ensures roughly equal distribution of patterns across shards,
/// regardless of the order in the input file.
pub fn shard_dataset<P: AsRef<Path>>(
    input: P,
    output_dir: P,
    num_shards: usize,
) -> Result<Vec<usize>, FerroError> {
    let input = input.as_ref();
    let output_dir = output_dir.as_ref();

    if num_shards == 0 {
        return Err(FerroError::Io {
            msg: "Number of shards must be > 0".to_string(),
        });
    }

    // Create output directory
    std::fs::create_dir_all(output_dir).map_err(|e| FerroError::Io {
        msg: format!("Failed to create directory {}: {}", output_dir.display(), e),
    })?;

    // Open input file
    let file = File::open(input).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", input.display(), e),
    })?;
    let reader = BufReader::new(file);

    // Create output files
    let mut writers: Vec<BufWriter<File>> = Vec::with_capacity(num_shards);
    for i in 0..num_shards {
        let shard_path = output_dir.join(format!("shard_{}.txt", i));
        let file = File::create(&shard_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create {}: {}", shard_path.display(), e),
        })?;
        writers.push(BufWriter::new(file));
    }

    // Distribute lines round-robin
    let mut counts = vec![0usize; num_shards];
    let mut line_index = 0usize;

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Error reading: {}", e),
        })?;

        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        let shard = line_index % num_shards;
        writeln!(writers[shard], "{}", trimmed).map_err(|e| FerroError::Io {
            msg: format!("Error writing to shard {}: {}", shard, e),
        })?;
        counts[shard] += 1;
        line_index += 1;
    }

    // Flush all writers
    for writer in &mut writers {
        writer.flush().map_err(|e| FerroError::Io {
            msg: format!("Error flushing: {}", e),
        })?;
    }

    let total: usize = counts.iter().sum();
    eprintln!(
        "Sharded {} patterns into {} files in {}",
        total,
        num_shards,
        output_dir.display()
    );

    Ok(counts)
}

/// Get the path to a specific shard file.
pub fn shard_path<P: AsRef<Path>>(output_dir: P, shard_index: usize) -> std::path::PathBuf {
    output_dir
        .as_ref()
        .join(format!("shard_{}.txt", shard_index))
}

/// Count lines in a file (for verification).
pub fn count_lines<P: AsRef<Path>>(path: P) -> Result<usize, FerroError> {
    let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", path.as_ref().display(), e),
    })?;
    let reader = BufReader::new(file);
    Ok(reader.lines().count())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    #[test]
    fn test_shard_dataset() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("input.txt");
        let output_dir = dir.path().join("shards");

        // Create input with 10 lines
        let mut f = File::create(&input).unwrap();
        for i in 0..10 {
            writeln!(f, "pattern_{}", i).unwrap();
        }

        // Shard into 3 files
        let counts = shard_dataset(&input, &output_dir, 3).unwrap();

        // Round-robin: 0,3,6,9 -> shard 0 (4 items)
        //              1,4,7   -> shard 1 (3 items)
        //              2,5,8   -> shard 2 (3 items)
        assert_eq!(counts, vec![4, 3, 3]);

        // Verify files exist and have correct line counts
        assert_eq!(count_lines(shard_path(&output_dir, 0)).unwrap(), 4);
        assert_eq!(count_lines(shard_path(&output_dir, 1)).unwrap(), 3);
        assert_eq!(count_lines(shard_path(&output_dir, 2)).unwrap(), 3);
    }
}
