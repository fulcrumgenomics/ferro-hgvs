//! UCSC chain file parser.
//!
//! This module parses UCSC chain files that define alignments between
//! genome assemblies for coordinate liftover.

use crate::error::FerroError;
use crate::reference::Strand;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// A single chain (alignment between two sequences).
#[derive(Debug, Clone)]
pub struct Chain {
    /// Chain ID.
    pub id: u64,
    /// Alignment score.
    pub score: u64,
    /// Target (reference) sequence name.
    pub target_name: String,
    /// Target sequence size.
    pub target_size: u64,
    /// Target strand.
    pub target_strand: Strand,
    /// Start position in target (0-based).
    pub target_start: u64,
    /// End position in target (0-based, exclusive).
    pub target_end: u64,
    /// Query sequence name.
    pub query_name: String,
    /// Query sequence size.
    pub query_size: u64,
    /// Query strand.
    pub query_strand: Strand,
    /// Start position in query (0-based).
    pub query_start: u64,
    /// End position in query (0-based, exclusive).
    pub query_end: u64,
    /// Alignment blocks.
    pub blocks: Vec<ChainBlock>,
}

/// An alignment block within a chain.
#[derive(Debug, Clone, Copy)]
pub struct ChainBlock {
    /// Size of the aligned block.
    pub size: u64,
    /// Gap in target after this block.
    pub target_gap: u64,
    /// Gap in query after this block.
    pub query_gap: u64,
}

impl Chain {
    /// Check if a target position falls within this chain's range.
    pub fn contains_target_pos(&self, pos: u64) -> bool {
        pos >= self.target_start && pos < self.target_end
    }

    /// Lift a position from target to query coordinates.
    ///
    /// Returns None if the position falls in a gap.
    pub fn lift_position(&self, target_pos: u64) -> Option<u64> {
        if !self.contains_target_pos(target_pos) {
            return None;
        }

        let mut t_pos = self.target_start;
        let mut q_pos = self.query_start;

        for block in &self.blocks {
            let block_end = t_pos + block.size;

            if target_pos < block_end {
                // Position is within this block
                let offset = target_pos - t_pos;
                return Some(match self.query_strand {
                    Strand::Plus => q_pos + offset,
                    Strand::Minus => self.query_size - (q_pos + offset) - 1,
                });
            }

            // Move to next block
            t_pos = block_end + block.target_gap;
            q_pos += block.size + block.query_gap;

            // Check if position is in a gap
            if target_pos < t_pos {
                return None; // In target gap
            }
        }

        None
    }

    /// Check if a position falls in a gap.
    pub fn is_in_gap(&self, target_pos: u64) -> bool {
        if !self.contains_target_pos(target_pos) {
            return false;
        }
        self.lift_position(target_pos).is_none()
    }
}

/// Chain file containing multiple chains indexed by target contig.
#[derive(Debug, Clone, Default)]
pub struct ChainFile {
    /// Chains indexed by target contig name.
    chains: HashMap<String, Vec<Chain>>,
}

impl ChainFile {
    /// Create a new empty chain file.
    pub fn new() -> Self {
        Self::default()
    }

    /// Load chain file from a path (supports .chain and .chain.gz).
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let path = path.as_ref();
        let file = File::open(path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open chain file: {}", e),
        })?;

        let path_str = path.to_string_lossy();
        if path_str.ends_with(".gz") {
            let decoder = flate2::read::GzDecoder::new(file);
            let reader = BufReader::new(decoder);
            Self::parse(reader)
        } else {
            let reader = BufReader::new(file);
            Self::parse(reader)
        }
    }

    /// Parse chain file from a reader.
    pub fn parse<R: Read>(reader: R) -> Result<Self, FerroError> {
        let buf_reader = BufReader::new(reader);
        let mut chains = ChainFile::new();
        let mut current_chain: Option<Chain> = None;
        let mut line_num = 0;

        for line_result in buf_reader.lines() {
            line_num += 1;
            let line = line_result.map_err(|e| FerroError::Io {
                msg: format!("Failed to read line {}: {}", line_num, e),
            })?;

            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            if line.starts_with("chain") {
                // Save previous chain if any
                if let Some(chain) = current_chain.take() {
                    chains.add_chain(chain);
                }

                // Parse chain header
                current_chain = Some(Self::parse_chain_header(line, line_num)?);
            } else if let Some(ref mut chain) = current_chain {
                // Parse alignment block
                if let Some(block) = Self::parse_block(line, line_num)? {
                    chain.blocks.push(block);
                }
            }
        }

        // Don't forget the last chain
        if let Some(chain) = current_chain {
            chains.add_chain(chain);
        }

        Ok(chains)
    }

    /// Parse a chain header line.
    fn parse_chain_header(line: &str, line_num: usize) -> Result<Chain, FerroError> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 12 {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "Invalid chain header at line {}: expected 12+ fields, got {}",
                    line_num,
                    parts.len()
                ),
            });
        }

        let score = parts[1]
            .parse::<u64>()
            .map_err(|_| FerroError::InvalidCoordinates {
                msg: format!("Invalid score at line {}", line_num),
            })?;

        let target_size = parts[3]
            .parse::<u64>()
            .map_err(|_| FerroError::InvalidCoordinates {
                msg: format!("Invalid target size at line {}", line_num),
            })?;

        let target_strand = match parts[4] {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            _ => {
                return Err(FerroError::InvalidCoordinates {
                    msg: format!("Invalid target strand '{}' at line {}", parts[4], line_num),
                })
            }
        };

        let target_start = parts[5]
            .parse::<u64>()
            .map_err(|_| FerroError::InvalidCoordinates {
                msg: format!("Invalid target start at line {}", line_num),
            })?;

        let target_end = parts[6]
            .parse::<u64>()
            .map_err(|_| FerroError::InvalidCoordinates {
                msg: format!("Invalid target end at line {}", line_num),
            })?;

        let query_size = parts[8]
            .parse::<u64>()
            .map_err(|_| FerroError::InvalidCoordinates {
                msg: format!("Invalid query size at line {}", line_num),
            })?;

        let query_strand = match parts[9] {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            _ => {
                return Err(FerroError::InvalidCoordinates {
                    msg: format!("Invalid query strand '{}' at line {}", parts[9], line_num),
                })
            }
        };

        let query_start = parts[10]
            .parse::<u64>()
            .map_err(|_| FerroError::InvalidCoordinates {
                msg: format!("Invalid query start at line {}", line_num),
            })?;

        let query_end = parts[11]
            .parse::<u64>()
            .map_err(|_| FerroError::InvalidCoordinates {
                msg: format!("Invalid query end at line {}", line_num),
            })?;

        let id = if parts.len() > 12 {
            parts[12].parse::<u64>().unwrap_or(0)
        } else {
            0
        };

        Ok(Chain {
            id,
            score,
            target_name: parts[2].to_string(),
            target_size,
            target_strand,
            target_start,
            target_end,
            query_name: parts[7].to_string(),
            query_size,
            query_strand,
            query_start,
            query_end,
            blocks: Vec::new(),
        })
    }

    /// Parse an alignment block line.
    fn parse_block(line: &str, _line_num: usize) -> Result<Option<ChainBlock>, FerroError> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            return Ok(None);
        }

        let size = parts[0]
            .parse::<u64>()
            .map_err(|_| FerroError::InvalidCoordinates {
                msg: format!("Invalid block size: {}", parts[0]),
            })?;

        // Last block only has size
        let (target_gap, query_gap) = if parts.len() >= 3 {
            let tg = parts[1].parse::<u64>().unwrap_or(0);
            let qg = parts[2].parse::<u64>().unwrap_or(0);
            (tg, qg)
        } else {
            (0, 0)
        };

        Ok(Some(ChainBlock {
            size,
            target_gap,
            query_gap,
        }))
    }

    /// Add a chain to the file.
    pub fn add_chain(&mut self, chain: Chain) {
        self.chains
            .entry(chain.target_name.clone())
            .or_default()
            .push(chain);
    }

    /// Find all chains covering a position on a contig.
    pub fn find_chains(&self, contig: &str, pos: u64) -> Vec<&Chain> {
        self.chains
            .get(contig)
            .map(|chains| {
                chains
                    .iter()
                    .filter(|c| c.contains_target_pos(pos))
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Get the best chain (highest score) covering a position.
    pub fn best_chain(&self, contig: &str, pos: u64) -> Option<&Chain> {
        self.find_chains(contig, pos)
            .into_iter()
            .max_by_key(|c| c.score)
    }

    /// Get all chains for a contig.
    pub fn chains_for_contig(&self, contig: &str) -> Option<&Vec<Chain>> {
        self.chains.get(contig)
    }

    /// Get all contig names.
    pub fn contig_names(&self) -> impl Iterator<Item = &str> {
        self.chains.keys().map(|s| s.as_str())
    }

    /// Get total number of chains.
    pub fn chain_count(&self) -> usize {
        self.chains.values().map(|v| v.len()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_chain_data() -> &'static str {
        r#"chain 1000 chr1 1000 + 0 1000 chr1 1100 + 0 1100 1
100	10	20
200	5	5
500

"#
    }

    #[test]
    fn test_parse_chain_file() {
        let chain_file = ChainFile::parse(simple_chain_data().as_bytes()).unwrap();
        assert_eq!(chain_file.chain_count(), 1);

        let chains = chain_file.chains_for_contig("chr1").unwrap();
        assert_eq!(chains.len(), 1);

        let chain = &chains[0];
        assert_eq!(chain.id, 1);
        assert_eq!(chain.score, 1000);
        assert_eq!(chain.target_name, "chr1");
        assert_eq!(chain.target_size, 1000);
        assert_eq!(chain.target_start, 0);
        assert_eq!(chain.target_end, 1000);
        assert_eq!(chain.query_name, "chr1");
        assert_eq!(chain.query_size, 1100);
        assert_eq!(chain.query_start, 0);
        assert_eq!(chain.query_end, 1100);
        assert_eq!(chain.blocks.len(), 3);
    }

    #[test]
    fn test_chain_blocks() {
        let chain_file = ChainFile::parse(simple_chain_data().as_bytes()).unwrap();
        let chain = &chain_file.chains_for_contig("chr1").unwrap()[0];

        assert_eq!(chain.blocks[0].size, 100);
        assert_eq!(chain.blocks[0].target_gap, 10);
        assert_eq!(chain.blocks[0].query_gap, 20);

        assert_eq!(chain.blocks[1].size, 200);
        assert_eq!(chain.blocks[1].target_gap, 5);
        assert_eq!(chain.blocks[1].query_gap, 5);

        assert_eq!(chain.blocks[2].size, 500);
        assert_eq!(chain.blocks[2].target_gap, 0);
        assert_eq!(chain.blocks[2].query_gap, 0);
    }

    #[test]
    fn test_lift_position_in_first_block() {
        let chain_file = ChainFile::parse(simple_chain_data().as_bytes()).unwrap();
        let chain = &chain_file.chains_for_contig("chr1").unwrap()[0];

        // Position 50 is in first block (0-99)
        let lifted = chain.lift_position(50).unwrap();
        assert_eq!(lifted, 50);
    }

    #[test]
    fn test_lift_position_in_gap() {
        let chain_file = ChainFile::parse(simple_chain_data().as_bytes()).unwrap();
        let chain = &chain_file.chains_for_contig("chr1").unwrap()[0];

        // Position 105 is in the gap after first block (100-109)
        assert!(chain.is_in_gap(105));
        assert!(chain.lift_position(105).is_none());
    }

    #[test]
    fn test_lift_position_in_second_block() {
        let chain_file = ChainFile::parse(simple_chain_data().as_bytes()).unwrap();
        let chain = &chain_file.chains_for_contig("chr1").unwrap()[0];

        // Position 110 is start of second block
        // Target: 0-99 (block1), 100-109 (gap), 110-309 (block2)
        // Query:  0-99 (block1), 100-119 (gap), 120-319 (block2)
        let lifted = chain.lift_position(110).unwrap();
        assert_eq!(lifted, 120);
    }

    #[test]
    fn test_find_chains() {
        let chain_file = ChainFile::parse(simple_chain_data().as_bytes()).unwrap();

        // Position in chain
        let chains = chain_file.find_chains("chr1", 50);
        assert_eq!(chains.len(), 1);

        // Position outside chain
        let chains = chain_file.find_chains("chr1", 2000);
        assert_eq!(chains.len(), 0);

        // Unknown contig
        let chains = chain_file.find_chains("chr2", 50);
        assert_eq!(chains.len(), 0);
    }

    #[test]
    fn test_best_chain() {
        let chain_data = r#"chain 1000 chr1 1000 + 0 500 chr1 1000 + 0 500 1
500

chain 2000 chr1 1000 + 0 800 chr1 1000 + 0 800 2
800

"#;
        let chain_file = ChainFile::parse(chain_data.as_bytes()).unwrap();

        // Position covered by both chains - should return higher score
        let best = chain_file.best_chain("chr1", 100).unwrap();
        assert_eq!(best.score, 2000);

        // Position only in second chain
        let best = chain_file.best_chain("chr1", 600).unwrap();
        assert_eq!(best.score, 2000);
    }

    #[test]
    fn test_chain_contains_pos() {
        let chain_file = ChainFile::parse(simple_chain_data().as_bytes()).unwrap();
        let chain = &chain_file.chains_for_contig("chr1").unwrap()[0];

        assert!(chain.contains_target_pos(0));
        assert!(chain.contains_target_pos(500));
        assert!(chain.contains_target_pos(999));
        assert!(!chain.contains_target_pos(1000)); // exclusive end
        assert!(!chain.contains_target_pos(1500));
    }
}
