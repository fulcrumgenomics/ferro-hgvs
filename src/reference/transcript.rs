//! Transcript and exon models
//!
//! # Coordinate System
//!
//! All coordinates in this module are **1-based inclusive**:
//!
//! | Field | Basis | Notes |
//! |-------|-------|-------|
//! | `Exon.start`, `Exon.end` | 1-based | Transcript coordinates (inclusive) |
//! | `Exon.genomic_start`, `Exon.genomic_end` | 1-based | Genomic coordinates (inclusive) |
//! | `Intron.genomic_start`, `Intron.genomic_end` | 1-based | First/last intronic base |
//! | `Transcript.cds_start`, `Transcript.cds_end` | 1-based | CDS boundaries in transcript space |
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::sync::OnceLock;

/// Genome build/assembly version
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum GenomeBuild {
    /// GRCh37 / hg19
    GRCh37,
    /// GRCh38 / hg38
    #[default]
    GRCh38,
    /// Unknown build
    Unknown,
}

impl std::fmt::Display for GenomeBuild {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GenomeBuild::GRCh37 => write!(f, "GRCh37"),
            GenomeBuild::GRCh38 => write!(f, "GRCh38"),
            GenomeBuild::Unknown => write!(f, "Unknown"),
        }
    }
}

/// Strand orientation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum Strand {
    #[serde(rename = "+")]
    #[default]
    Plus,
    #[serde(rename = "-")]
    Minus,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
        }
    }
}

/// An exon in a transcript
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Exon {
    /// Exon number (1-based)
    pub number: u32,
    /// Start position in transcript coordinates (1-based, inclusive)
    pub start: u64,
    /// End position in transcript coordinates (1-based, inclusive)
    pub end: u64,
    /// Genomic start position (1-based, inclusive)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genomic_start: Option<u64>,
    /// Genomic end position (1-based, inclusive)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genomic_end: Option<u64>,
}

/// An intron between two exons
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Intron {
    /// Intron number (1-based, intron 1 is between exon 1 and exon 2)
    pub number: u32,
    /// 5' exon number (upstream exon)
    pub upstream_exon: u32,
    /// 3' exon number (downstream exon)
    pub downstream_exon: u32,
    /// Genomic start position (1-based, first intronic base)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genomic_start: Option<u64>,
    /// Genomic end position (1-based, last intronic base)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genomic_end: Option<u64>,
    /// Transcript position of last exonic base before intron (5' boundary)
    pub tx_5prime_boundary: u64,
    /// Transcript position of first exonic base after intron (3' boundary)
    pub tx_3prime_boundary: u64,
}

impl Intron {
    /// Create a new intron
    pub fn new(
        number: u32,
        upstream_exon: u32,
        downstream_exon: u32,
        tx_5prime_boundary: u64,
        tx_3prime_boundary: u64,
    ) -> Self {
        Self {
            number,
            upstream_exon,
            downstream_exon,
            genomic_start: None,
            genomic_end: None,
            tx_5prime_boundary,
            tx_3prime_boundary,
        }
    }

    /// Create a new intron with genomic coordinates
    pub fn with_genomic(
        number: u32,
        upstream_exon: u32,
        downstream_exon: u32,
        tx_5prime_boundary: u64,
        tx_3prime_boundary: u64,
        genomic_start: u64,
        genomic_end: u64,
    ) -> Self {
        Self {
            number,
            upstream_exon,
            downstream_exon,
            genomic_start: Some(genomic_start),
            genomic_end: Some(genomic_end),
            tx_5prime_boundary,
            tx_3prime_boundary,
        }
    }

    /// Get the genomic length of the intron
    ///
    /// Uses saturating arithmetic to prevent overflow in edge cases.
    pub fn genomic_length(&self) -> Option<u64> {
        match (self.genomic_start, self.genomic_end) {
            (Some(start), Some(end)) if end >= start => Some((end - start).saturating_add(1)),
            _ => None,
        }
    }

    /// Check if a genomic position is within this intron
    pub fn contains_genomic(&self, pos: u64) -> bool {
        match (self.genomic_start, self.genomic_end) {
            (Some(start), Some(end)) => pos >= start && pos <= end,
            _ => false,
        }
    }
}

impl Exon {
    /// Create a new exon with transcript coordinates only
    pub fn new(number: u32, start: u64, end: u64) -> Self {
        Self {
            number,
            start,
            end,
            genomic_start: None,
            genomic_end: None,
        }
    }

    /// Create a new exon with both transcript and genomic coordinates
    pub fn with_genomic(
        number: u32,
        start: u64,
        end: u64,
        genomic_start: u64,
        genomic_end: u64,
    ) -> Self {
        Self {
            number,
            start,
            end,
            genomic_start: Some(genomic_start),
            genomic_end: Some(genomic_end),
        }
    }

    /// Length of the exon
    pub fn len(&self) -> u64 {
        if self.end >= self.start {
            self.end - self.start + 1
        } else {
            0
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Check if a position is within this exon
    pub fn contains(&self, pos: u64) -> bool {
        pos >= self.start && pos <= self.end
    }
}

/// MANE (Matched Annotation from NCBI and EBI) transcript status
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum ManeStatus {
    /// Not a MANE transcript
    #[default]
    None,
    /// MANE Select - single representative transcript per gene
    Select,
    /// MANE Plus Clinical - additional clinically relevant transcripts
    PlusClinical,
}

impl ManeStatus {
    /// Check if this is a MANE transcript of any type
    pub fn is_mane(&self) -> bool {
        !matches!(self, ManeStatus::None)
    }

    /// Check if this is MANE Select
    pub fn is_select(&self) -> bool {
        matches!(self, ManeStatus::Select)
    }

    /// Check if this is MANE Plus Clinical
    pub fn is_plus_clinical(&self) -> bool {
        matches!(self, ManeStatus::PlusClinical)
    }

    /// Get priority score for sorting (lower is better)
    pub fn priority(&self) -> u8 {
        match self {
            ManeStatus::Select => 0,
            ManeStatus::PlusClinical => 1,
            ManeStatus::None => 2,
        }
    }
}

impl std::fmt::Display for ManeStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ManeStatus::None => write!(f, ""),
            ManeStatus::Select => write!(f, "MANE Select"),
            ManeStatus::PlusClinical => write!(f, "MANE Plus Clinical"),
        }
    }
}

/// A transcript with its exon structure and sequence
#[derive(Debug, Serialize, Deserialize)]
pub struct Transcript {
    /// Transcript accession (e.g., "NM_000088.3")
    pub id: String,

    /// Gene symbol (e.g., "COL1A1")
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gene_symbol: Option<String>,

    /// Strand orientation
    pub strand: Strand,

    /// Full transcript sequence
    pub sequence: String,

    /// CDS start position (1-based, in transcript coordinates)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cds_start: Option<u64>,

    /// CDS end position (1-based, in transcript coordinates)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cds_end: Option<u64>,

    /// List of exons
    pub exons: Vec<Exon>,

    /// Chromosome name (e.g., "chr1", "1", "X")
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chromosome: Option<String>,

    /// Genomic start position of transcript (1-based, inclusive)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genomic_start: Option<u64>,

    /// Genomic end position of transcript (1-based, inclusive)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genomic_end: Option<u64>,

    /// Genome build/assembly version
    #[serde(default)]
    pub genome_build: GenomeBuild,

    /// MANE (Matched Annotation from NCBI and EBI) status
    #[serde(default)]
    pub mane_status: ManeStatus,

    /// RefSeq accession matched to this transcript (for Ensembl transcripts)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub refseq_match: Option<String>,

    /// Ensembl accession matched to this transcript (for RefSeq transcripts)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ensembl_match: Option<String>,

    /// Cached introns (computed lazily from exons)
    /// This field is for internal use and will be initialized automatically.
    #[serde(skip)]
    pub(crate) cached_introns: OnceLock<Vec<Intron>>,
}

impl Clone for Transcript {
    fn clone(&self) -> Self {
        Self {
            id: self.id.clone(),
            gene_symbol: self.gene_symbol.clone(),
            strand: self.strand,
            sequence: self.sequence.clone(),
            cds_start: self.cds_start,
            cds_end: self.cds_end,
            exons: self.exons.clone(),
            chromosome: self.chromosome.clone(),
            genomic_start: self.genomic_start,
            genomic_end: self.genomic_end,
            genome_build: self.genome_build,
            mane_status: self.mane_status,
            refseq_match: self.refseq_match.clone(),
            ensembl_match: self.ensembl_match.clone(),
            // Cache is reset on clone - will be lazily re-initialized
            cached_introns: OnceLock::new(),
        }
    }
}

impl PartialEq for Transcript {
    fn eq(&self, other: &Self) -> bool {
        // Compare all fields except the cache
        self.id == other.id
            && self.gene_symbol == other.gene_symbol
            && self.strand == other.strand
            && self.sequence == other.sequence
            && self.cds_start == other.cds_start
            && self.cds_end == other.cds_end
            && self.exons == other.exons
            && self.chromosome == other.chromosome
            && self.genomic_start == other.genomic_start
            && self.genomic_end == other.genomic_end
            && self.genome_build == other.genome_build
            && self.mane_status == other.mane_status
            && self.refseq_match == other.refseq_match
            && self.ensembl_match == other.ensembl_match
    }
}

impl Eq for Transcript {}

impl Transcript {
    /// Create a new Transcript with the given fields
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: String,
        gene_symbol: Option<String>,
        strand: Strand,
        sequence: String,
        cds_start: Option<u64>,
        cds_end: Option<u64>,
        exons: Vec<Exon>,
        chromosome: Option<String>,
        genomic_start: Option<u64>,
        genomic_end: Option<u64>,
        genome_build: GenomeBuild,
        mane_status: ManeStatus,
        refseq_match: Option<String>,
        ensembl_match: Option<String>,
    ) -> Self {
        Self {
            id,
            gene_symbol,
            strand,
            sequence,
            cds_start,
            cds_end,
            exons,
            chromosome,
            genomic_start,
            genomic_end,
            genome_build,
            mane_status,
            refseq_match,
            ensembl_match,
            cached_introns: OnceLock::new(),
        }
    }

    /// Get the length of the transcript sequence
    pub fn sequence_length(&self) -> u64 {
        self.sequence.len() as u64
    }

    /// Check if this is a coding transcript
    pub fn is_coding(&self) -> bool {
        self.cds_start.is_some() && self.cds_end.is_some()
    }

    /// Get the CDS length
    pub fn cds_length(&self) -> Option<u64> {
        match (self.cds_start, self.cds_end) {
            (Some(start), Some(end)) if end >= start => Some(end - start + 1),
            _ => None,
        }
    }

    /// Get sequence at a position range (0-based)
    pub fn get_sequence(&self, start: u64, end: u64) -> Option<&str> {
        let start = start as usize;
        let end = end as usize;
        if end <= self.sequence.len() && start < end {
            Some(&self.sequence[start..end])
        } else {
            None
        }
    }

    /// Find which exon contains a position using binary search
    ///
    /// Assumes exons are sorted by start position (which they typically are).
    /// Falls back to linear search if exons are not sorted.
    pub fn exon_at(&self, pos: u64) -> Option<&Exon> {
        // Use binary search for O(log n) lookup
        self.exons
            .binary_search_by(|e| {
                if pos < e.start {
                    Ordering::Greater
                } else if pos > e.end {
                    Ordering::Less
                } else {
                    Ordering::Equal
                }
            })
            .ok()
            .map(|i| &self.exons[i])
    }

    /// Get the 5' UTR length
    pub fn utr5_length(&self) -> Option<u64> {
        self.cds_start.map(|s| s.saturating_sub(1))
    }

    /// Get the 3' UTR length
    pub fn utr3_length(&self) -> Option<u64> {
        self.cds_end
            .map(|e| self.sequence_length().saturating_sub(e))
    }

    /// Check if this transcript has genomic coordinates
    pub fn has_genomic_coords(&self) -> bool {
        self.chromosome.is_some() && self.genomic_start.is_some() && self.genomic_end.is_some()
    }

    /// Get the genomic length (span) of the transcript
    pub fn genomic_length(&self) -> Option<u64> {
        match (self.genomic_start, self.genomic_end) {
            (Some(start), Some(end)) if end >= start => Some(end - start + 1),
            _ => None,
        }
    }

    /// Check if a genomic position falls within this transcript's span
    pub fn contains_genomic_pos(&self, pos: u64) -> bool {
        match (self.genomic_start, self.genomic_end) {
            (Some(start), Some(end)) => pos >= start && pos <= end,
            _ => false,
        }
    }

    /// Check if this is a MANE Select transcript
    pub fn is_mane_select(&self) -> bool {
        self.mane_status.is_select()
    }

    /// Check if this is a MANE Plus Clinical transcript
    pub fn is_mane_plus_clinical(&self) -> bool {
        self.mane_status.is_plus_clinical()
    }

    /// Check if this is any type of MANE transcript
    pub fn is_mane(&self) -> bool {
        self.mane_status.is_mane()
    }

    /// Get cached introns, computing them lazily if not already cached
    ///
    /// This method provides O(1) access to introns after the first call,
    /// avoiding recalculation on every lookup.
    pub fn introns(&self) -> &[Intron] {
        self.cached_introns.get_or_init(|| self.compute_introns())
    }

    /// Calculate introns from exon boundaries
    ///
    /// Returns a vector of Intron structs derived from adjacent exon pairs.
    /// Exons should be sorted by transcript position.
    ///
    /// Note: This method uses caching internally. For repeated access,
    /// prefer `introns()` which returns a slice reference.
    pub fn calculate_introns(&self) -> Vec<Intron> {
        self.introns().to_vec()
    }

    /// Compute introns from exon boundaries (internal, uncached)
    fn compute_introns(&self) -> Vec<Intron> {
        let mut introns = Vec::new();

        // Sort exons by transcript position
        let mut sorted_exons: Vec<_> = self.exons.iter().collect();
        sorted_exons.sort_by_key(|e| e.start);

        for (i, window) in sorted_exons.windows(2).enumerate() {
            let upstream = window[0];
            let downstream = window[1];

            let intron_number = (i + 1) as u32;

            // Calculate genomic coordinates for the intron if available
            let (genomic_start, genomic_end) = match self.strand {
                Strand::Plus => {
                    // Plus strand: intron starts after upstream exon ends genomically
                    let g_start = upstream.genomic_end.map(|e| e + 1);
                    let g_end = downstream.genomic_start.map(|s| s - 1);
                    (g_start, g_end)
                }
                Strand::Minus => {
                    // Minus strand: genomic coordinates are reversed
                    // Intron is between downstream's genomic_end+1 and upstream's genomic_start-1
                    let g_start = downstream.genomic_end.map(|e| e + 1);
                    let g_end = upstream.genomic_start.map(|s| s - 1);
                    (g_start, g_end)
                }
            };

            let mut intron = Intron::new(
                intron_number,
                upstream.number,
                downstream.number,
                upstream.end,     // tx position of last base in upstream exon
                downstream.start, // tx position of first base in downstream exon
            );

            if let (Some(gs), Some(ge)) = (genomic_start, genomic_end) {
                if ge >= gs {
                    intron.genomic_start = Some(gs);
                    intron.genomic_end = Some(ge);
                }
            }

            introns.push(intron);
        }

        introns
    }

    /// Find which intron contains a genomic position
    ///
    /// Returns the intron and the offset from the nearest exon boundary.
    /// Positive offset means downstream from 5' boundary (c.N+offset notation).
    /// Negative offset means upstream from 3' boundary (c.N-offset notation).
    pub fn find_intron_at_genomic(&self, genomic_pos: u64) -> Option<(Intron, IntronPosition)> {
        // Use cached introns for O(1) access after first call
        for intron in self.introns() {
            if intron.contains_genomic(genomic_pos) {
                let (g_start, g_end) = (intron.genomic_start?, intron.genomic_end?);
                // Use saturating_add to prevent overflow at u64::MAX
                let intron_length = (g_end - g_start).saturating_add(1);

                // Calculate distance from each boundary
                // Use saturating arithmetic to prevent overflow
                let (dist_to_5prime, dist_to_3prime) = match self.strand {
                    Strand::Plus => {
                        // Plus strand: genomic start is 5' boundary
                        let from_5prime = (genomic_pos - g_start).saturating_add(1); // 1-based offset
                        let from_3prime = (g_end - genomic_pos).saturating_add(1);
                        (from_5prime, from_3prime)
                    }
                    Strand::Minus => {
                        // Minus strand: genomic end is 5' boundary
                        let from_5prime = (g_end - genomic_pos).saturating_add(1);
                        let from_3prime = (genomic_pos - g_start).saturating_add(1);
                        (from_5prime, from_3prime)
                    }
                };

                // Determine which boundary is closer and create position
                let position = if dist_to_5prime <= dist_to_3prime {
                    // Closer to 5' boundary (or equal): use +offset notation
                    IntronPosition {
                        intron_number: intron.number,
                        boundary: IntronBoundary::FivePrime,
                        offset: dist_to_5prime as i64,
                        tx_boundary_pos: intron.tx_5prime_boundary,
                        intron_length,
                    }
                } else {
                    // Closer to 3' boundary: use -offset notation
                    IntronPosition {
                        intron_number: intron.number,
                        boundary: IntronBoundary::ThreePrime,
                        offset: -(dist_to_3prime as i64),
                        tx_boundary_pos: intron.tx_3prime_boundary,
                        intron_length,
                    }
                };

                return Some((intron.clone(), position));
            }
        }

        None
    }

    /// Get the number of introns in this transcript
    pub fn intron_count(&self) -> usize {
        if self.exons.len() > 1 {
            self.exons.len() - 1
        } else {
            0
        }
    }

    /// Find an intron given a transcript boundary position and offset
    ///
    /// This is used to convert intronic positions like c.100+5 or c.200-10
    /// to find which intron they're in.
    ///
    /// # Arguments
    /// * `tx_boundary` - The transcript position of the exon boundary
    /// * `offset` - The offset into the intron (positive = after exon, negative = before exon)
    ///
    /// # Returns
    /// The intron if found, along with whether this is a 5' or 3' boundary reference
    pub fn find_intron_at_tx_boundary(&self, tx_boundary: u64, offset: i64) -> Option<&Intron> {
        for intron in self.introns() {
            if offset > 0 && intron.tx_5prime_boundary == tx_boundary {
                // Positive offset: c.N+offset means we're after the 5' boundary (end of upstream exon)
                return Some(intron);
            } else if offset < 0 && intron.tx_3prime_boundary == tx_boundary {
                // Negative offset: c.N-offset means we're before the 3' boundary (start of downstream exon)
                return Some(intron);
            }
        }
        None
    }

    /// Convert an intronic position to a genomic coordinate
    ///
    /// # Arguments
    /// * `tx_boundary` - The transcript position of the exon boundary (the base part of c.N+offset)
    /// * `offset` - The offset into the intron
    ///
    /// # Returns
    /// The genomic position if the transcript has genomic coordinates
    pub fn intronic_to_genomic(&self, tx_boundary: u64, offset: i64) -> Option<u64> {
        let intron = self.find_intron_at_tx_boundary(tx_boundary, offset)?;

        // Get the genomic coordinates of the intron
        let (g_start, g_end) = (intron.genomic_start?, intron.genomic_end?);

        match self.strand {
            Strand::Plus => {
                if offset > 0 {
                    // c.N+offset: offset bases after the 5' exon boundary
                    // genomic = intron start + offset - 1 (offset is 1-based)
                    Some(g_start + offset as u64 - 1)
                } else {
                    // c.N-offset: |offset| bases before the 3' exon boundary
                    // genomic = intron end - |offset| + 1
                    Some(g_end - (-offset) as u64 + 1)
                }
            }
            Strand::Minus => {
                // On minus strand, genomic coordinates are reversed relative to transcript
                if offset > 0 {
                    // c.N+offset: offset bases after the 5' exon boundary
                    // For minus strand: intron end - offset + 1
                    Some(g_end - offset as u64 + 1)
                } else {
                    // c.N-offset: |offset| bases before the 3' exon boundary
                    // For minus strand: intron start + |offset| - 1
                    Some(g_start + (-offset) as u64 - 1)
                }
            }
        }
    }

    /// Convert a genomic position to intronic transcript notation
    ///
    /// # Arguments
    /// * `genomic_pos` - The genomic position
    ///
    /// # Returns
    /// A tuple of (tx_boundary_position, offset) where:
    /// - tx_boundary_position is the CDS/transcript position of the nearest exon boundary
    /// - offset is positive (c.N+offset) or negative (c.N-offset)
    pub fn genomic_to_intronic(&self, genomic_pos: u64) -> Option<(u64, i64)> {
        let (intron, intron_pos) = self.find_intron_at_genomic(genomic_pos)?;

        match intron_pos.boundary {
            IntronBoundary::FivePrime => {
                // Reference the 5' exon boundary with positive offset
                Some((intron.tx_5prime_boundary, intron_pos.offset))
            }
            IntronBoundary::ThreePrime => {
                // Reference the 3' exon boundary with negative offset
                Some((intron.tx_3prime_boundary, intron_pos.offset))
            }
        }
    }
}

/// Position within an intron
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IntronPosition {
    /// Intron number (1-based)
    pub intron_number: u32,
    /// Which exon boundary this position is relative to
    pub boundary: IntronBoundary,
    /// Offset from the boundary (positive for 5' boundary, negative for 3')
    pub offset: i64,
    /// Transcript position of the exon boundary
    pub tx_boundary_pos: u64,
    /// Total length of the intron in bases
    pub intron_length: u64,
}

impl IntronPosition {
    /// Check if this is a deep intronic position (>50bp from nearest exon)
    pub fn is_deep_intronic(&self) -> bool {
        self.offset.abs() > 50
    }

    /// Check if this is at a canonical splice site (within 2bp of exon)
    pub fn is_canonical_splice_site(&self) -> bool {
        self.offset.abs() <= 2
    }

    /// Check if this is near a splice site (within 10bp of exon)
    pub fn is_near_splice_site(&self) -> bool {
        self.offset.abs() <= 10
    }

    /// Check if this is in the extended splice region (within 20bp of exon)
    pub fn is_extended_splice_region(&self) -> bool {
        self.offset.abs() <= 20
    }

    /// Get the splice site type based on position
    pub fn splice_site_type(&self) -> SpliceSiteType {
        let abs_offset = self.offset.abs();
        match self.boundary {
            IntronBoundary::FivePrime => {
                // 5' end of intron = splice donor site
                if abs_offset <= 2 {
                    SpliceSiteType::DonorCanonical
                } else if abs_offset <= 6 {
                    SpliceSiteType::DonorExtended
                } else if abs_offset <= 20 {
                    SpliceSiteType::DonorRegion
                } else if abs_offset <= 50 {
                    SpliceSiteType::NearSplice
                } else {
                    SpliceSiteType::DeepIntronic
                }
            }
            IntronBoundary::ThreePrime => {
                // 3' end of intron = splice acceptor site
                if abs_offset <= 2 {
                    SpliceSiteType::AcceptorCanonical
                } else if abs_offset <= 12 {
                    // Branch point region is typically -18 to -35
                    SpliceSiteType::AcceptorExtended
                } else if abs_offset <= 20 {
                    SpliceSiteType::AcceptorRegion
                } else if abs_offset <= 50 {
                    SpliceSiteType::NearSplice
                } else {
                    SpliceSiteType::DeepIntronic
                }
            }
        }
    }

    /// Get distance from splice donor (5' end of intron)
    pub fn distance_from_donor(&self) -> Option<u64> {
        match self.boundary {
            IntronBoundary::FivePrime => Some(self.offset.unsigned_abs()),
            IntronBoundary::ThreePrime => {
                // Distance = intron_length - distance_from_3prime
                Some(self.intron_length - self.offset.unsigned_abs())
            }
        }
    }

    /// Get distance from splice acceptor (3' end of intron)
    pub fn distance_from_acceptor(&self) -> Option<u64> {
        match self.boundary {
            IntronBoundary::ThreePrime => Some(self.offset.unsigned_abs()),
            IntronBoundary::FivePrime => {
                // Distance = intron_length - distance_from_5prime
                Some(self.intron_length - self.offset.unsigned_abs())
            }
        }
    }
}

/// Which boundary of an intron a position is relative to
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntronBoundary {
    /// 5' boundary (end of upstream exon) - uses + offset notation
    FivePrime,
    /// 3' boundary (start of downstream exon) - uses - offset notation
    ThreePrime,
}

/// Type of splice site position
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpliceSiteType {
    /// Canonical splice donor (GT, positions +1 to +2)
    DonorCanonical,
    /// Extended splice donor consensus (positions +3 to +6)
    DonorExtended,
    /// Splice donor region (positions +7 to +20)
    DonorRegion,
    /// Canonical splice acceptor (AG, positions -1 to -2)
    AcceptorCanonical,
    /// Extended splice acceptor including polypyrimidine tract (positions -3 to -12)
    AcceptorExtended,
    /// Splice acceptor region (positions -13 to -20)
    AcceptorRegion,
    /// Near splice site but not in critical region (21-50bp from exon)
    NearSplice,
    /// Deep intronic (>50bp from nearest exon)
    DeepIntronic,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_transcript() -> Transcript {
        Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGCATGCATGCATGCATGC".to_string(), // 20 bases
            cds_start: Some(5),
            cds_end: Some(15),
            exons: vec![
                Exon::new(1, 1, 7),
                Exon::new(2, 8, 14),
                Exon::new(3, 15, 20),
            ],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        }
    }

    fn make_test_transcript_with_genomic() -> Transcript {
        Transcript {
            id: "NM_000088.4".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGCATGCATGCATGCATGC".to_string(),
            cds_start: Some(5),
            cds_end: Some(15),
            exons: vec![
                Exon::with_genomic(1, 1, 7, 50189542, 50189548),
                Exon::with_genomic(2, 8, 14, 50190100, 50190106),
                Exon::with_genomic(3, 15, 20, 50190500, 50190505),
            ],
            chromosome: Some("chr17".to_string()),
            genomic_start: Some(50189542),
            genomic_end: Some(50190505),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: Some("ENST00000123456.5".to_string()),
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn test_transcript_sequence_length() {
        let tx = make_test_transcript();
        assert_eq!(tx.sequence_length(), 20);
    }

    #[test]
    fn test_transcript_is_coding() {
        let tx = make_test_transcript();
        assert!(tx.is_coding());

        let noncoding = Transcript {
            cds_start: None,
            cds_end: None,
            ..tx
        };
        assert!(!noncoding.is_coding());
    }

    #[test]
    fn test_transcript_cds_length() {
        let tx = make_test_transcript();
        assert_eq!(tx.cds_length(), Some(11));
    }

    #[test]
    fn test_exon_contains() {
        let exon = Exon::new(1, 10, 20);
        assert!(exon.contains(10));
        assert!(exon.contains(15));
        assert!(exon.contains(20));
        assert!(!exon.contains(9));
        assert!(!exon.contains(21));
    }

    #[test]
    fn test_exon_at() {
        let tx = make_test_transcript();
        let exon = tx.exon_at(10).unwrap();
        assert_eq!(exon.number, 2);
    }

    #[test]
    fn test_get_sequence() {
        let tx = make_test_transcript();
        assert_eq!(tx.get_sequence(0, 3), Some("ATG"));
    }

    #[test]
    fn test_strand_display() {
        assert_eq!(format!("{}", Strand::Plus), "+");
        assert_eq!(format!("{}", Strand::Minus), "-");
    }

    #[test]
    fn test_genome_build_display() {
        assert_eq!(format!("{}", GenomeBuild::GRCh37), "GRCh37");
        assert_eq!(format!("{}", GenomeBuild::GRCh38), "GRCh38");
        assert_eq!(format!("{}", GenomeBuild::Unknown), "Unknown");
    }

    #[test]
    fn test_genome_build_default() {
        assert_eq!(GenomeBuild::default(), GenomeBuild::GRCh38);
    }

    #[test]
    fn test_exon_with_genomic() {
        let exon = Exon::with_genomic(1, 1, 100, 50000, 50099);
        assert_eq!(exon.number, 1);
        assert_eq!(exon.start, 1);
        assert_eq!(exon.end, 100);
        assert_eq!(exon.genomic_start, Some(50000));
        assert_eq!(exon.genomic_end, Some(50099));
    }

    #[test]
    fn test_transcript_has_genomic_coords() {
        let tx = make_test_transcript();
        assert!(!tx.has_genomic_coords());

        let tx_genomic = make_test_transcript_with_genomic();
        assert!(tx_genomic.has_genomic_coords());
    }

    #[test]
    fn test_transcript_genomic_length() {
        let tx = make_test_transcript();
        assert_eq!(tx.genomic_length(), None);

        let tx_genomic = make_test_transcript_with_genomic();
        // 50190505 - 50189542 + 1 = 964
        assert_eq!(tx_genomic.genomic_length(), Some(964));
    }

    #[test]
    fn test_transcript_contains_genomic_pos() {
        let tx = make_test_transcript();
        assert!(!tx.contains_genomic_pos(50189542));

        let tx_genomic = make_test_transcript_with_genomic();
        assert!(tx_genomic.contains_genomic_pos(50189542)); // start
        assert!(tx_genomic.contains_genomic_pos(50190000)); // middle
        assert!(tx_genomic.contains_genomic_pos(50190505)); // end
        assert!(!tx_genomic.contains_genomic_pos(50189541)); // before
        assert!(!tx_genomic.contains_genomic_pos(50190506)); // after
    }

    #[test]
    fn test_mane_status_default() {
        assert_eq!(ManeStatus::default(), ManeStatus::None);
    }

    #[test]
    fn test_mane_status_methods() {
        assert!(!ManeStatus::None.is_mane());
        assert!(!ManeStatus::None.is_select());
        assert!(!ManeStatus::None.is_plus_clinical());

        assert!(ManeStatus::Select.is_mane());
        assert!(ManeStatus::Select.is_select());
        assert!(!ManeStatus::Select.is_plus_clinical());

        assert!(ManeStatus::PlusClinical.is_mane());
        assert!(!ManeStatus::PlusClinical.is_select());
        assert!(ManeStatus::PlusClinical.is_plus_clinical());
    }

    #[test]
    fn test_mane_status_priority() {
        assert!(ManeStatus::Select.priority() < ManeStatus::PlusClinical.priority());
        assert!(ManeStatus::PlusClinical.priority() < ManeStatus::None.priority());
    }

    #[test]
    fn test_mane_status_display() {
        assert_eq!(format!("{}", ManeStatus::None), "");
        assert_eq!(format!("{}", ManeStatus::Select), "MANE Select");
        assert_eq!(
            format!("{}", ManeStatus::PlusClinical),
            "MANE Plus Clinical"
        );
    }

    #[test]
    fn test_transcript_mane_methods() {
        let tx = make_test_transcript();
        assert!(!tx.is_mane());
        assert!(!tx.is_mane_select());
        assert!(!tx.is_mane_plus_clinical());

        let tx_mane = make_test_transcript_with_genomic();
        assert!(tx_mane.is_mane());
        assert!(tx_mane.is_mane_select());
        assert!(!tx_mane.is_mane_plus_clinical());
    }

    #[test]
    fn test_transcript_matched_accessions() {
        let tx = make_test_transcript_with_genomic();
        assert_eq!(tx.refseq_match, None);
        assert_eq!(tx.ensembl_match, Some("ENST00000123456.5".to_string()));
    }

    #[test]
    fn test_calculate_introns() {
        let tx = make_test_transcript_with_genomic();
        let introns = tx.calculate_introns();

        // Should have 2 introns for 3 exons
        assert_eq!(introns.len(), 2);

        // Intron 1 is between exon 1 and exon 2
        let intron1 = &introns[0];
        assert_eq!(intron1.number, 1);
        assert_eq!(intron1.upstream_exon, 1);
        assert_eq!(intron1.downstream_exon, 2);
        assert_eq!(intron1.tx_5prime_boundary, 7); // end of exon 1
        assert_eq!(intron1.tx_3prime_boundary, 8); // start of exon 2

        // Intron 2 is between exon 2 and exon 3
        let intron2 = &introns[1];
        assert_eq!(intron2.number, 2);
        assert_eq!(intron2.upstream_exon, 2);
        assert_eq!(intron2.downstream_exon, 3);
    }

    #[test]
    fn test_intron_genomic_coords() {
        let tx = make_test_transcript_with_genomic();
        let introns = tx.calculate_introns();

        // Intron 1 genomic: 50189549 to 50190099 (between exon 1 end 50189548 and exon 2 start 50190100)
        let intron1 = &introns[0];
        assert_eq!(intron1.genomic_start, Some(50189549));
        assert_eq!(intron1.genomic_end, Some(50190099));
        assert!(intron1.genomic_length().is_some());
        assert_eq!(intron1.genomic_length().unwrap(), 551);
    }

    #[test]
    fn test_intron_contains_genomic() {
        let tx = make_test_transcript_with_genomic();
        let introns = tx.calculate_introns();
        let intron1 = &introns[0];

        // Position inside intron 1
        assert!(intron1.contains_genomic(50189600));
        assert!(intron1.contains_genomic(50190000));

        // Position outside intron 1
        assert!(!intron1.contains_genomic(50189548)); // in exon 1
        assert!(!intron1.contains_genomic(50190100)); // in exon 2
    }

    #[test]
    fn test_find_intron_at_genomic() {
        let tx = make_test_transcript_with_genomic();

        // Position in intron 1 (close to 5' end)
        let result = tx.find_intron_at_genomic(50189550);
        assert!(result.is_some());
        let (intron, pos) = result.unwrap();
        assert_eq!(intron.number, 1);
        assert_eq!(pos.intron_number, 1);
        assert_eq!(pos.boundary, IntronBoundary::FivePrime);
        assert_eq!(pos.offset, 2); // 2 bases into intron

        // Position in intron 1 (close to 3' end)
        let result = tx.find_intron_at_genomic(50190098);
        assert!(result.is_some());
        let (_, pos) = result.unwrap();
        assert_eq!(pos.boundary, IntronBoundary::ThreePrime);
        assert!(pos.offset < 0);
    }

    #[test]
    fn test_intron_position_splice_site_type() {
        // Test canonical donor (+1, +2)
        let pos = IntronPosition {
            intron_number: 1,
            boundary: IntronBoundary::FivePrime,
            offset: 1,
            tx_boundary_pos: 100,
            intron_length: 500,
        };
        assert_eq!(pos.splice_site_type(), SpliceSiteType::DonorCanonical);
        assert!(pos.is_canonical_splice_site());
        assert!(!pos.is_deep_intronic());

        // Test canonical acceptor (-1, -2)
        let pos = IntronPosition {
            intron_number: 1,
            boundary: IntronBoundary::ThreePrime,
            offset: -2,
            tx_boundary_pos: 200,
            intron_length: 500,
        };
        assert_eq!(pos.splice_site_type(), SpliceSiteType::AcceptorCanonical);

        // Test deep intronic (>50bp)
        let pos = IntronPosition {
            intron_number: 1,
            boundary: IntronBoundary::FivePrime,
            offset: 100,
            tx_boundary_pos: 100,
            intron_length: 500,
        };
        assert_eq!(pos.splice_site_type(), SpliceSiteType::DeepIntronic);
        assert!(pos.is_deep_intronic());
    }

    #[test]
    fn test_intron_position_distances() {
        let pos = IntronPosition {
            intron_number: 1,
            boundary: IntronBoundary::FivePrime,
            offset: 10,
            tx_boundary_pos: 100,
            intron_length: 500,
        };

        // Distance from donor should be the offset
        assert_eq!(pos.distance_from_donor(), Some(10));
        // Distance from acceptor should be intron_length - offset
        assert_eq!(pos.distance_from_acceptor(), Some(490));
    }

    #[test]
    fn test_intron_count() {
        let tx = make_test_transcript_with_genomic();
        assert_eq!(tx.intron_count(), 2);

        // Single exon transcript should have 0 introns
        let single_exon = Transcript {
            id: "NR_TEST.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: "A".repeat(100),
            cds_start: None,
            cds_end: None,
            exons: vec![Exon::new(1, 1, 100)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            cached_introns: std::sync::OnceLock::new(),
        };
        assert_eq!(single_exon.intron_count(), 0);
    }
}
