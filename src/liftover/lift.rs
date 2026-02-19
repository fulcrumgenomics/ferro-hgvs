//! Coordinate liftover between genome builds.
//!
//! This module provides the main liftover engine for converting coordinates
//! between GRCh37 and GRCh38.
//!
//! # Coordinate System
//!
//! | Context | Basis | Notes |
//! |---------|-------|-------|
//! | Public API (`pos`, `start`, `end`) | 1-based | HGVS convention |
//! | Chain file operations | 0-based | Converted internally via `saturating_sub(1)` |
//! | Lifted result | 1-based | Converted back via `+ 1` |
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use super::chain::ChainFile;
use crate::coords::{hgvs_pos_to_index, index_to_hgvs_pos};
use crate::error::FerroError;
use crate::reference::transcript::GenomeBuild;
use std::path::Path;

/// Result of a liftover operation.
#[derive(Debug, Clone)]
pub struct LiftoverResult {
    /// Source genome build.
    pub source_build: GenomeBuild,
    /// Source contig/chromosome name.
    pub source_contig: String,
    /// Source position (1-based).
    pub source_pos: u64,
    /// Target genome build.
    pub target_build: GenomeBuild,
    /// Target contig/chromosome name.
    pub target_contig: String,
    /// Target position (1-based).
    pub target_pos: u64,
    /// Chain ID used for the mapping.
    pub chain_id: u64,
    /// Chain score (higher = more reliable).
    pub chain_score: u64,
    /// True if position falls in a gap (unmappable but interpolated).
    pub in_gap: bool,
    /// True if multiple chains cover this position.
    pub multiple_mappings: bool,
}

/// Result of an interval liftover.
#[derive(Debug, Clone)]
pub struct IntervalLiftoverResult {
    /// Source interval (start, end) in 1-based coordinates.
    pub source_interval: (u64, u64),
    /// Target interval if successfully lifted.
    pub target_interval: Option<(u64, u64)>,
    /// Status of the liftover.
    pub status: LiftoverStatus,
    /// Chain used for the mapping.
    pub chain_id: Option<u64>,
}

/// Status of an interval liftover.
#[derive(Debug, Clone, PartialEq)]
pub enum LiftoverStatus {
    /// Successfully lifted.
    Success,
    /// Partially mapped (some bases fell in gaps).
    PartiallyMapped {
        /// Fraction of interval that was mapped.
        mapped_fraction: f64,
    },
    /// Interval split across multiple chains.
    SplitAcrossChains,
    /// Position is unmappable.
    Unmappable,
    /// Contig not found.
    ContigNotFound,
}

/// Liftover engine for coordinate conversion between genome builds.
#[derive(Debug)]
pub struct Liftover {
    /// GRCh37 -> GRCh38 chain file.
    forward: ChainFile,
    /// GRCh38 -> GRCh37 chain file.
    reverse: ChainFile,
}

impl Liftover {
    /// Create a new liftover engine from chain files.
    pub fn new(forward: ChainFile, reverse: ChainFile) -> Self {
        Self { forward, reverse }
    }

    /// Load liftover from chain file paths.
    ///
    /// # Arguments
    ///
    /// * `grch37_to_38` - Path to hg19ToHg38.over.chain.gz
    /// * `grch38_to_37` - Path to hg38ToHg19.over.chain.gz
    pub fn from_files<P: AsRef<Path>>(
        grch37_to_38: P,
        grch38_to_37: P,
    ) -> Result<Self, FerroError> {
        let forward = ChainFile::from_file(grch37_to_38)?;
        let reverse = ChainFile::from_file(grch38_to_37)?;
        Ok(Self::new(forward, reverse))
    }

    /// Create a one-way liftover (forward direction only).
    pub fn one_way(chain_file: ChainFile) -> Self {
        Self {
            forward: chain_file,
            reverse: ChainFile::new(),
        }
    }

    /// Get the appropriate chain file for a liftover direction.
    fn get_chain_file(&self, from: GenomeBuild, to: GenomeBuild) -> Result<&ChainFile, FerroError> {
        match (from, to) {
            (GenomeBuild::GRCh37, GenomeBuild::GRCh38) => Ok(&self.forward),
            (GenomeBuild::GRCh38, GenomeBuild::GRCh37) => Ok(&self.reverse),
            _ => Err(FerroError::InvalidCoordinates {
                msg: format!("Unsupported liftover direction: {:?} -> {:?}", from, to),
            }),
        }
    }

    /// Lift a single position from one build to another.
    ///
    /// Positions are 1-based (HGVS convention).
    pub fn lift(
        &self,
        from: GenomeBuild,
        to: GenomeBuild,
        contig: &str,
        pos: u64,
    ) -> Result<LiftoverResult, FerroError> {
        let chain_file = self.get_chain_file(from, to)?;

        // Convert from 1-based HGVS to 0-based for chain operations
        let pos_0based = hgvs_pos_to_index(pos) as u64;

        // Find chains covering this position
        let chains = chain_file.find_chains(contig, pos_0based);

        if chains.is_empty() {
            // Suggest alternatives when position cannot be lifted
            let suggestion = if from == GenomeBuild::GRCh37 {
                "Position may not exist in GRCh38 (deleted region) or contig name may differ."
            } else {
                "Position may not exist in GRCh37 (new sequence in GRCh38) or contig name may differ."
            };
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "No chain found for {}:{} in {:?} -> {:?} liftover. {}",
                    contig, pos, from, to, suggestion
                ),
            });
        }

        let multiple_mappings = chains.len() > 1;

        // Use best chain (highest score)
        let chain = chains.into_iter().max_by_key(|c| c.score).unwrap();

        // Try to lift the position
        match chain.lift_position(pos_0based) {
            Some(lifted_0based) => Ok(LiftoverResult {
                source_build: from,
                source_contig: contig.to_string(),
                source_pos: pos,
                target_build: to,
                target_contig: chain.query_name.clone(),
                target_pos: index_to_hgvs_pos(lifted_0based as usize), // Convert 0-based to 1-based
                chain_id: chain.id,
                chain_score: chain.score,
                in_gap: false,
                multiple_mappings,
            }),
            None => {
                // Position is in a gap - try to interpolate or return error
                Err(FerroError::InvalidCoordinates {
                    msg: format!("Position {}:{} falls in an alignment gap", contig, pos),
                })
            }
        }
    }

    /// Lift from GRCh37 to GRCh38.
    pub fn lift_37_to_38(&self, contig: &str, pos: u64) -> Result<LiftoverResult, FerroError> {
        self.lift(GenomeBuild::GRCh37, GenomeBuild::GRCh38, contig, pos)
    }

    /// Lift from GRCh38 to GRCh37.
    pub fn lift_38_to_37(&self, contig: &str, pos: u64) -> Result<LiftoverResult, FerroError> {
        self.lift(GenomeBuild::GRCh38, GenomeBuild::GRCh37, contig, pos)
    }

    /// Lift an interval from one build to another.
    pub fn lift_interval(
        &self,
        from: GenomeBuild,
        to: GenomeBuild,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<IntervalLiftoverResult, FerroError> {
        let chain_file = self.get_chain_file(from, to)?;

        // Convert 1-based HGVS positions to 0-based for chain operations
        let start_0based = hgvs_pos_to_index(start) as u64;
        let end_0based = hgvs_pos_to_index(end) as u64;

        // Find chains covering the start and end
        let start_chains = chain_file.find_chains(contig, start_0based);
        let end_chains = chain_file.find_chains(contig, end_0based);

        if start_chains.is_empty() || end_chains.is_empty() {
            return Ok(IntervalLiftoverResult {
                source_interval: (start, end),
                target_interval: None,
                status: if start_chains.is_empty() && end_chains.is_empty() {
                    LiftoverStatus::ContigNotFound
                } else {
                    LiftoverStatus::Unmappable
                },
                chain_id: None,
            });
        }

        // Check if both ends map through the same chain
        let common_chain = start_chains
            .iter()
            .find(|sc| end_chains.iter().any(|ec| ec.id == sc.id));

        match common_chain {
            Some(chain) => {
                let lifted_start = chain.lift_position(start_0based);
                let lifted_end = chain.lift_position(end_0based);

                match (lifted_start, lifted_end) {
                    (Some(ls), Some(le)) => {
                        // Both ends mapped successfully
                        let (target_start, target_end) = if ls <= le {
                            (ls + 1, le + 1)
                        } else {
                            (le + 1, ls + 1) // Minus strand reversal
                        };

                        Ok(IntervalLiftoverResult {
                            source_interval: (start, end),
                            target_interval: Some((target_start, target_end)),
                            status: LiftoverStatus::Success,
                            chain_id: Some(chain.id),
                        })
                    }
                    (Some(ls), None) | (None, Some(ls)) => {
                        // One end in a gap
                        Ok(IntervalLiftoverResult {
                            source_interval: (start, end),
                            target_interval: Some((ls + 1, ls + 1)),
                            status: LiftoverStatus::PartiallyMapped {
                                mapped_fraction: 0.5,
                            },
                            chain_id: Some(chain.id),
                        })
                    }
                    (None, None) => {
                        // Both ends in gaps
                        Ok(IntervalLiftoverResult {
                            source_interval: (start, end),
                            target_interval: None,
                            status: LiftoverStatus::Unmappable,
                            chain_id: Some(chain.id),
                        })
                    }
                }
            }
            None => {
                // Interval spans multiple chains
                Ok(IntervalLiftoverResult {
                    source_interval: (start, end),
                    target_interval: None,
                    status: LiftoverStatus::SplitAcrossChains,
                    chain_id: None,
                })
            }
        }
    }

    /// Get the forward chain file.
    pub fn forward_chains(&self) -> &ChainFile {
        &self.forward
    }

    /// Get the reverse chain file.
    pub fn reverse_chains(&self) -> &ChainFile {
        &self.reverse
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_chain_file() -> ChainFile {
        let chain_data = r#"chain 1000 chr1 1000000 + 10000 20000 chr1 1000100 + 10050 20150 1
5000
5000

"#;
        ChainFile::parse(chain_data.as_bytes()).unwrap()
    }

    #[test]
    fn test_lift_position() {
        let chain_file = test_chain_file();
        let liftover = Liftover::one_way(chain_file);

        // Position 10001 (1-based) should map to query
        let result = liftover
            .lift(GenomeBuild::GRCh37, GenomeBuild::GRCh38, "chr1", 10001)
            .unwrap();

        assert_eq!(result.source_pos, 10001);
        assert_eq!(result.target_contig, "chr1");
        assert!(!result.in_gap);
    }

    #[test]
    fn test_lift_position_not_found() {
        let chain_file = test_chain_file();
        let liftover = Liftover::one_way(chain_file);

        // Position outside chain range
        let result = liftover.lift(GenomeBuild::GRCh37, GenomeBuild::GRCh38, "chr1", 1);

        assert!(result.is_err());
    }

    #[test]
    fn test_lift_interval() {
        let chain_file = test_chain_file();
        let liftover = Liftover::one_way(chain_file);

        let result = liftover
            .lift_interval(
                GenomeBuild::GRCh37,
                GenomeBuild::GRCh38,
                "chr1",
                10001,
                10100,
            )
            .unwrap();

        assert_eq!(result.status, LiftoverStatus::Success);
        assert!(result.target_interval.is_some());
    }

    #[test]
    fn test_lift_unknown_contig() {
        let chain_file = test_chain_file();
        let liftover = Liftover::one_way(chain_file);

        let result = liftover.lift(GenomeBuild::GRCh37, GenomeBuild::GRCh38, "chrZ", 1000);

        assert!(result.is_err());
    }

    // Additional liftover tests

    #[test]
    fn test_liftover_result_fields() {
        let result = LiftoverResult {
            source_build: GenomeBuild::GRCh37,
            source_contig: "chr1".to_string(),
            source_pos: 12345,
            target_build: GenomeBuild::GRCh38,
            target_contig: "chr1".to_string(),
            target_pos: 12445,
            chain_id: 1,
            chain_score: 1000,
            in_gap: false,
            multiple_mappings: false,
        };

        assert_eq!(result.source_pos, 12345);
        assert_eq!(result.target_pos, 12445);
        assert_eq!(result.chain_score, 1000);
        assert!(!result.in_gap);
        assert!(!result.multiple_mappings);
    }

    #[test]
    fn test_liftover_result_clone() {
        let result = LiftoverResult {
            source_build: GenomeBuild::GRCh37,
            source_contig: "chr1".to_string(),
            source_pos: 100,
            target_build: GenomeBuild::GRCh38,
            target_contig: "chr1".to_string(),
            target_pos: 200,
            chain_id: 1,
            chain_score: 500,
            in_gap: true,
            multiple_mappings: true,
        };

        let cloned = result.clone();
        assert_eq!(cloned.source_pos, result.source_pos);
        assert_eq!(cloned.target_pos, result.target_pos);
        assert_eq!(cloned.in_gap, result.in_gap);
    }

    #[test]
    fn test_interval_liftover_result_success() {
        let result = IntervalLiftoverResult {
            source_interval: (1000, 2000),
            target_interval: Some((1100, 2100)),
            status: LiftoverStatus::Success,
            chain_id: Some(1),
        };

        assert_eq!(result.source_interval, (1000, 2000));
        assert_eq!(result.target_interval, Some((1100, 2100)));
        assert_eq!(result.status, LiftoverStatus::Success);
        assert_eq!(result.chain_id, Some(1));
    }

    #[test]
    fn test_interval_liftover_result_unmappable() {
        let result = IntervalLiftoverResult {
            source_interval: (1000, 2000),
            target_interval: None,
            status: LiftoverStatus::Unmappable,
            chain_id: None,
        };

        assert!(result.target_interval.is_none());
        assert_eq!(result.status, LiftoverStatus::Unmappable);
    }

    #[test]
    fn test_liftover_status_partially_mapped() {
        let status = LiftoverStatus::PartiallyMapped {
            mapped_fraction: 0.75,
        };

        if let LiftoverStatus::PartiallyMapped { mapped_fraction } = status {
            assert!((mapped_fraction - 0.75).abs() < 0.01);
        } else {
            panic!("Expected PartiallyMapped status");
        }
    }

    #[test]
    fn test_liftover_status_equality() {
        assert_eq!(LiftoverStatus::Success, LiftoverStatus::Success);
        assert_eq!(LiftoverStatus::Unmappable, LiftoverStatus::Unmappable);
        assert_eq!(
            LiftoverStatus::ContigNotFound,
            LiftoverStatus::ContigNotFound
        );
        assert_ne!(LiftoverStatus::Success, LiftoverStatus::Unmappable);
    }

    #[test]
    fn test_liftover_status_split_across_chains() {
        let status = LiftoverStatus::SplitAcrossChains;
        assert_eq!(status, LiftoverStatus::SplitAcrossChains);
    }

    #[test]
    fn test_liftover_new() {
        let forward = ChainFile::new();
        let reverse = ChainFile::new();
        let liftover = Liftover::new(forward, reverse);

        // Just verify we can access the chains
        let _ = liftover.forward_chains();
        let _ = liftover.reverse_chains();
    }

    #[test]
    fn test_liftover_one_way() {
        let chain_file = test_chain_file();
        let liftover = Liftover::one_way(chain_file);

        // Forward should have chains (can find positions)
        let result = liftover.lift(GenomeBuild::GRCh37, GenomeBuild::GRCh38, "chr1", 10001);
        assert!(result.is_ok());

        // Reverse should be empty for one-way (lift fails)
        let reverse_result = liftover.lift(GenomeBuild::GRCh38, GenomeBuild::GRCh37, "chr1", 10001);
        assert!(reverse_result.is_err());
    }

    #[test]
    fn test_lift_37_to_38_convenience() {
        let chain_file = test_chain_file();
        let liftover = Liftover::one_way(chain_file);

        let result = liftover.lift_37_to_38("chr1", 10001);
        assert!(result.is_ok());

        let r = result.unwrap();
        assert_eq!(r.source_build, GenomeBuild::GRCh37);
        assert_eq!(r.target_build, GenomeBuild::GRCh38);
    }

    #[test]
    fn test_lift_38_to_37_no_reverse_chain() {
        let chain_file = test_chain_file();
        let liftover = Liftover::one_way(chain_file);

        // One-way liftover has no reverse chain
        let result = liftover.lift_38_to_37("chr1", 10001);
        assert!(result.is_err());
    }

    #[test]
    fn test_lift_interval_contig_not_found() {
        let chain_file = test_chain_file();
        let liftover = Liftover::one_way(chain_file);

        let result = liftover
            .lift_interval(
                GenomeBuild::GRCh37,
                GenomeBuild::GRCh38,
                "chrUnknown",
                1000,
                2000,
            )
            .unwrap();

        assert_eq!(result.status, LiftoverStatus::ContigNotFound);
        assert!(result.target_interval.is_none());
    }

    #[test]
    fn test_lift_interval_partial_unmappable() {
        let chain_file = test_chain_file();
        let liftover = Liftover::one_way(chain_file);

        // Start in range, end out of range
        let result = liftover
            .lift_interval(
                GenomeBuild::GRCh37,
                GenomeBuild::GRCh38,
                "chr1",
                10001,
                99999, // Way outside chain range
            )
            .unwrap();

        // Should be unmappable or partially mapped
        assert_ne!(result.status, LiftoverStatus::Success);
    }

    #[test]
    fn test_lift_position_boundary() {
        let chain_file = test_chain_file();
        let liftover = Liftover::one_way(chain_file);

        // Test at exact chain start boundary
        let result = liftover.lift(GenomeBuild::GRCh37, GenomeBuild::GRCh38, "chr1", 10001);
        assert!(result.is_ok());
    }

    #[test]
    fn test_get_chain_file_unsupported_direction() {
        let liftover = Liftover::new(ChainFile::new(), ChainFile::new());

        // Same build to same build should fail
        let result = liftover.lift(GenomeBuild::GRCh37, GenomeBuild::GRCh37, "chr1", 1000);
        assert!(result.is_err());
    }

    #[test]
    fn test_liftover_debug() {
        let result = LiftoverResult {
            source_build: GenomeBuild::GRCh37,
            source_contig: "chr1".to_string(),
            source_pos: 100,
            target_build: GenomeBuild::GRCh38,
            target_contig: "chr1".to_string(),
            target_pos: 200,
            chain_id: 1,
            chain_score: 500,
            in_gap: false,
            multiple_mappings: false,
        };

        let debug_str = format!("{:?}", result);
        assert!(debug_str.contains("LiftoverResult"));
        assert!(debug_str.contains("source_pos"));
    }

    #[test]
    fn test_interval_result_debug() {
        let result = IntervalLiftoverResult {
            source_interval: (100, 200),
            target_interval: Some((150, 250)),
            status: LiftoverStatus::Success,
            chain_id: Some(1),
        };

        let debug_str = format!("{:?}", result);
        assert!(debug_str.contains("IntervalLiftoverResult"));
    }

    #[test]
    fn test_liftover_status_debug() {
        let status = LiftoverStatus::PartiallyMapped {
            mapped_fraction: 0.5,
        };
        let debug_str = format!("{:?}", status);
        assert!(debug_str.contains("PartiallyMapped"));
        assert!(debug_str.contains("0.5"));
    }
}
