//! Protein effect prediction from DNA/RNA variants.
//!
//! This module provides functionality to predict the protein-level effects
//! of coding sequence variants.
//!
//! # Example
//!
//! ```
//! use ferro_hgvs::effect::{Consequence, Impact, EffectPredictor};
//! use ferro_hgvs::hgvs::location::AminoAcid;
//!
//! let predictor = EffectPredictor::new();
//!
//! // Classify a missense change
//! let effect = predictor.classify_amino_acid_change(
//!     &AminoAcid::Val,
//!     &AminoAcid::Glu,
//!     600,
//! );
//!
//! assert_eq!(effect.consequences[0], Consequence::MissenseVariant);
//! assert_eq!(effect.impact, Impact::Moderate);
//! ```

use crate::backtranslate::{Codon, CodonTable};
use crate::hgvs::location::AminoAcid;

/// Sequence Ontology consequence term.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Consequence {
    /// Complete transcript deletion.
    TranscriptAblation,
    /// Within 2bp of splice acceptor site (AG).
    SpliceAcceptorVariant,
    /// Within 2bp of splice donor site (GT).
    SpliceDonorVariant,
    /// Introduces a premature stop codon.
    StopGained,
    /// Insertion/deletion causing frameshift.
    FrameshiftVariant,
    /// Stop codon changed to amino acid.
    StopLost,
    /// Start codon changed to other.
    StartLost,
    /// Amino acid substitution.
    MissenseVariant,
    /// In-frame insertion of amino acids.
    InframeInsertion,
    /// In-frame deletion of amino acids.
    InframeDeletion,
    /// Generic protein-altering variant.
    ProteinAlteringVariant,
    /// Within 3-8bp of splice site.
    SpliceRegionVariant,
    /// Silent change (codon change, same amino acid).
    SynonymousVariant,
    /// Start codon unchanged.
    StartRetainedVariant,
    /// Stop codon unchanged.
    StopRetainedVariant,
    /// Variant in 5' UTR.
    FivePrimeUtrVariant,
    /// Variant in 3' UTR.
    ThreePrimeUtrVariant,
    /// Variant in intron.
    IntronVariant,
    /// Coding sequence variant (general).
    CodingSequenceVariant,
}

impl Consequence {
    /// Get the Sequence Ontology term.
    pub fn so_term(&self) -> &'static str {
        match self {
            Consequence::TranscriptAblation => "transcript_ablation",
            Consequence::SpliceAcceptorVariant => "splice_acceptor_variant",
            Consequence::SpliceDonorVariant => "splice_donor_variant",
            Consequence::StopGained => "stop_gained",
            Consequence::FrameshiftVariant => "frameshift_variant",
            Consequence::StopLost => "stop_lost",
            Consequence::StartLost => "start_lost",
            Consequence::MissenseVariant => "missense_variant",
            Consequence::InframeInsertion => "inframe_insertion",
            Consequence::InframeDeletion => "inframe_deletion",
            Consequence::ProteinAlteringVariant => "protein_altering_variant",
            Consequence::SpliceRegionVariant => "splice_region_variant",
            Consequence::SynonymousVariant => "synonymous_variant",
            Consequence::StartRetainedVariant => "start_retained_variant",
            Consequence::StopRetainedVariant => "stop_retained_variant",
            Consequence::FivePrimeUtrVariant => "5_prime_UTR_variant",
            Consequence::ThreePrimeUtrVariant => "3_prime_UTR_variant",
            Consequence::IntronVariant => "intron_variant",
            Consequence::CodingSequenceVariant => "coding_sequence_variant",
        }
    }

    /// Get the Sequence Ontology ID.
    pub fn so_id(&self) -> &'static str {
        match self {
            Consequence::TranscriptAblation => "SO:0001893",
            Consequence::SpliceAcceptorVariant => "SO:0001574",
            Consequence::SpliceDonorVariant => "SO:0001575",
            Consequence::StopGained => "SO:0001587",
            Consequence::FrameshiftVariant => "SO:0001589",
            Consequence::StopLost => "SO:0001578",
            Consequence::StartLost => "SO:0002012",
            Consequence::MissenseVariant => "SO:0001583",
            Consequence::InframeInsertion => "SO:0001821",
            Consequence::InframeDeletion => "SO:0001822",
            Consequence::ProteinAlteringVariant => "SO:0001818",
            Consequence::SpliceRegionVariant => "SO:0001630",
            Consequence::SynonymousVariant => "SO:0001819",
            Consequence::StartRetainedVariant => "SO:0002019",
            Consequence::StopRetainedVariant => "SO:0001567",
            Consequence::FivePrimeUtrVariant => "SO:0001623",
            Consequence::ThreePrimeUtrVariant => "SO:0001624",
            Consequence::IntronVariant => "SO:0001627",
            Consequence::CodingSequenceVariant => "SO:0001580",
        }
    }

    /// Get the impact level.
    pub fn impact(&self) -> Impact {
        match self {
            Consequence::TranscriptAblation
            | Consequence::SpliceAcceptorVariant
            | Consequence::SpliceDonorVariant
            | Consequence::StopGained
            | Consequence::FrameshiftVariant
            | Consequence::StopLost
            | Consequence::StartLost => Impact::High,

            Consequence::MissenseVariant
            | Consequence::InframeInsertion
            | Consequence::InframeDeletion
            | Consequence::ProteinAlteringVariant => Impact::Moderate,

            Consequence::SpliceRegionVariant
            | Consequence::SynonymousVariant
            | Consequence::StartRetainedVariant
            | Consequence::StopRetainedVariant => Impact::Low,

            Consequence::FivePrimeUtrVariant
            | Consequence::ThreePrimeUtrVariant
            | Consequence::IntronVariant
            | Consequence::CodingSequenceVariant => Impact::Modifier,
        }
    }

    /// Get a human-readable description.
    pub fn description(&self) -> &'static str {
        match self {
            Consequence::TranscriptAblation => "Complete transcript deletion",
            Consequence::SpliceAcceptorVariant => "Variant in splice acceptor site (2bp)",
            Consequence::SpliceDonorVariant => "Variant in splice donor site (2bp)",
            Consequence::StopGained => "Premature stop codon introduced",
            Consequence::FrameshiftVariant => "Frameshift causing protein truncation",
            Consequence::StopLost => "Stop codon changed to amino acid",
            Consequence::StartLost => "Start codon changed",
            Consequence::MissenseVariant => "Amino acid substitution",
            Consequence::InframeInsertion => "In-frame amino acid insertion",
            Consequence::InframeDeletion => "In-frame amino acid deletion",
            Consequence::ProteinAlteringVariant => "Protein-altering variant",
            Consequence::SpliceRegionVariant => "Variant in splice region (3-8bp)",
            Consequence::SynonymousVariant => "Silent change (same amino acid)",
            Consequence::StartRetainedVariant => "Start codon preserved",
            Consequence::StopRetainedVariant => "Stop codon preserved",
            Consequence::FivePrimeUtrVariant => "Variant in 5' UTR",
            Consequence::ThreePrimeUtrVariant => "Variant in 3' UTR",
            Consequence::IntronVariant => "Variant in intron",
            Consequence::CodingSequenceVariant => "Variant in coding sequence",
        }
    }
}

impl std::fmt::Display for Consequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.so_term())
    }
}

/// Variant impact level (VEP-style).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Impact {
    /// Modifier - minimal predicted impact.
    Modifier,
    /// Low impact.
    Low,
    /// Moderate impact.
    Moderate,
    /// High impact (likely deleterious).
    High,
}

impl Impact {
    /// Get the impact as a string.
    pub fn as_str(&self) -> &'static str {
        match self {
            Impact::High => "HIGH",
            Impact::Moderate => "MODERATE",
            Impact::Low => "LOW",
            Impact::Modifier => "MODIFIER",
        }
    }
}

impl std::fmt::Display for Impact {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Amino acid change details.
#[derive(Debug, Clone)]
pub struct AminoAcidChange {
    /// Position in protein (1-based).
    pub position: u64,
    /// Reference amino acid.
    pub ref_aa: AminoAcid,
    /// Alternate amino acid.
    pub alt_aa: AminoAcid,
    /// Reference codon (if known).
    pub ref_codon: Option<Codon>,
    /// Alternate codon (if known).
    pub alt_codon: Option<Codon>,
}

/// Protein effect prediction result.
#[derive(Debug, Clone)]
pub struct ProteinEffect {
    /// All applicable consequences.
    pub consequences: Vec<Consequence>,
    /// Highest impact level.
    pub impact: Impact,
    /// Amino acid change details (for coding variants).
    pub amino_acid_change: Option<AminoAcidChange>,
    /// Intronic offset (for splice variants).
    pub intronic_offset: Option<i64>,
}

impl ProteinEffect {
    /// Get the most severe consequence.
    pub fn most_severe(&self) -> Option<&Consequence> {
        self.consequences.iter().max_by_key(|c| c.impact())
    }

    /// Check if this is a high-impact variant.
    pub fn is_high_impact(&self) -> bool {
        self.impact == Impact::High
    }

    /// Check if this is a protein-altering variant.
    pub fn is_protein_altering(&self) -> bool {
        self.consequences.iter().any(|c| {
            matches!(
                c,
                Consequence::MissenseVariant
                    | Consequence::StopGained
                    | Consequence::StopLost
                    | Consequence::StartLost
                    | Consequence::FrameshiftVariant
                    | Consequence::InframeInsertion
                    | Consequence::InframeDeletion
                    | Consequence::ProteinAlteringVariant
            )
        })
    }
}

/// Protein effect predictor.
#[derive(Debug, Clone)]
pub struct EffectPredictor {
    codon_table: CodonTable,
}

impl EffectPredictor {
    /// Create a new effect predictor.
    pub fn new() -> Self {
        Self {
            codon_table: CodonTable::standard(),
        }
    }

    /// Classify an amino acid change.
    pub fn classify_amino_acid_change(
        &self,
        ref_aa: &AminoAcid,
        alt_aa: &AminoAcid,
        position: u64,
    ) -> ProteinEffect {
        let consequence = if ref_aa == alt_aa {
            Consequence::SynonymousVariant
        } else if *alt_aa == AminoAcid::Ter {
            Consequence::StopGained
        } else if *ref_aa == AminoAcid::Ter {
            Consequence::StopLost
        } else if *ref_aa == AminoAcid::Met && position == 1 {
            Consequence::StartLost
        } else {
            Consequence::MissenseVariant
        };

        ProteinEffect {
            consequences: vec![consequence],
            impact: consequence.impact(),
            amino_acid_change: Some(AminoAcidChange {
                position,
                ref_aa: *ref_aa,
                alt_aa: *alt_aa,
                ref_codon: None,
                alt_codon: None,
            }),
            intronic_offset: None,
        }
    }

    /// Classify an indel by frame effect.
    pub fn classify_indel(&self, ref_len: usize, alt_len: usize) -> ProteinEffect {
        let net_change = alt_len as i64 - ref_len as i64;

        let consequence = if net_change % 3 == 0 {
            if net_change > 0 {
                Consequence::InframeInsertion
            } else if net_change < 0 {
                Consequence::InframeDeletion
            } else {
                Consequence::CodingSequenceVariant
            }
        } else {
            Consequence::FrameshiftVariant
        };

        ProteinEffect {
            consequences: vec![consequence],
            impact: consequence.impact(),
            amino_acid_change: None,
            intronic_offset: None,
        }
    }

    /// Classify a splice site variant by distance.
    pub fn classify_splice_variant(&self, offset: i64) -> ProteinEffect {
        let abs_offset = offset.abs();

        let consequence = if abs_offset <= 2 {
            if offset > 0 {
                Consequence::SpliceDonorVariant
            } else {
                Consequence::SpliceAcceptorVariant
            }
        } else if abs_offset <= 8 {
            Consequence::SpliceRegionVariant
        } else {
            Consequence::IntronVariant
        };

        ProteinEffect {
            consequences: vec![consequence],
            impact: consequence.impact(),
            amino_acid_change: None,
            intronic_offset: Some(offset),
        }
    }

    /// Classify a UTR variant.
    pub fn classify_utr_variant(&self, is_5_prime: bool) -> ProteinEffect {
        let consequence = if is_5_prime {
            Consequence::FivePrimeUtrVariant
        } else {
            Consequence::ThreePrimeUtrVariant
        };

        ProteinEffect {
            consequences: vec![consequence],
            impact: consequence.impact(),
            amino_acid_change: None,
            intronic_offset: None,
        }
    }

    /// Get the codon table.
    pub fn codon_table(&self) -> &CodonTable {
        &self.codon_table
    }
}

impl Default for EffectPredictor {
    fn default() -> Self {
        Self::new()
    }
}

// =============================================================================
// NMD (Nonsense-Mediated Decay) Prediction
// =============================================================================

/// NMD prediction result.
///
/// Nonsense-mediated decay is triggered when a premature termination codon (PTC)
/// is located more than 50-55 nucleotides upstream of the last exon-exon junction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NmdPrediction {
    /// NMD is predicted to be triggered.
    Triggered {
        /// Distance from PTC to last exon-exon junction (in nucleotides).
        distance_to_junction: u64,
        /// Confidence level of prediction.
        confidence: NmdConfidence,
    },
    /// NMD is not predicted to be triggered.
    NotTriggered {
        /// Reason NMD is not triggered.
        reason: NmdNotTriggeredReason,
    },
    /// NMD prediction is uncertain.
    Uncertain {
        /// Reason for uncertainty.
        reason: String,
    },
}

/// Confidence level for NMD prediction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NmdConfidence {
    /// High confidence - PTC is >55nt upstream of last junction.
    High,
    /// Medium confidence - PTC is 50-55nt upstream of last junction.
    Medium,
    /// Low confidence - borderline case or limited transcript information.
    Low,
}

/// Reason why NMD is not triggered.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NmdNotTriggeredReason {
    /// PTC is in the last exon (no downstream exon-exon junction).
    PtcInLastExon,
    /// PTC is within 50nt of the last exon-exon junction.
    PtcNearLastJunction {
        /// Distance to last junction.
        distance: u64,
    },
    /// Transcript has only one exon (no exon-exon junctions).
    SingleExonTranscript,
    /// No premature termination codon present.
    NoPtc,
    /// PTC position is at or after the normal stop codon.
    PtcAtOrAfterNormalStop,
}

/// Transcript information needed for NMD prediction.
#[derive(Debug, Clone)]
pub struct TranscriptForNmd {
    /// Exon boundaries (list of (start, end) positions in CDS coordinates).
    /// These are 1-based, inclusive.
    pub exon_boundaries: Vec<(u64, u64)>,
    /// Position of normal stop codon in CDS coordinates (1-based).
    pub normal_stop_position: u64,
    /// CDS length in nucleotides.
    pub cds_length: u64,
}

/// NMD predictor.
#[derive(Debug, Clone)]
pub struct NmdPredictor {
    /// NMD threshold distance (default: 50-55 nucleotides).
    /// PTCs more than this distance upstream of the last exon-exon junction trigger NMD.
    pub threshold_low: u64,
    pub threshold_high: u64,
}

impl NmdPredictor {
    /// Create a new NMD predictor with default thresholds.
    pub fn new() -> Self {
        Self {
            threshold_low: 50,
            threshold_high: 55,
        }
    }

    /// Create an NMD predictor with custom thresholds.
    pub fn with_thresholds(threshold_low: u64, threshold_high: u64) -> Self {
        Self {
            threshold_low,
            threshold_high,
        }
    }

    /// Predict whether NMD will be triggered for a premature stop codon.
    ///
    /// # Arguments
    /// * `ptc_position` - Position of the premature termination codon in CDS coordinates (1-based, first codon position).
    /// * `transcript` - Transcript information including exon boundaries.
    ///
    /// # Returns
    /// NMD prediction result.
    pub fn predict(&self, ptc_position: u64, transcript: &TranscriptForNmd) -> NmdPrediction {
        // Check if transcript has exon junctions
        if transcript.exon_boundaries.len() <= 1 {
            return NmdPrediction::NotTriggered {
                reason: NmdNotTriggeredReason::SingleExonTranscript,
            };
        }

        // Check if PTC is at or after normal stop
        if ptc_position >= transcript.normal_stop_position {
            return NmdPrediction::NotTriggered {
                reason: NmdNotTriggeredReason::PtcAtOrAfterNormalStop,
            };
        }

        // Find the last exon-exon junction position
        // The last junction is at the start of the last exon
        let last_junction = transcript.exon_boundaries.last().map(|(start, _)| *start);

        match last_junction {
            Some(junction_pos) => {
                // Check if PTC is in the last exon
                if ptc_position >= junction_pos {
                    return NmdPrediction::NotTriggered {
                        reason: NmdNotTriggeredReason::PtcInLastExon,
                    };
                }

                // Calculate distance from PTC to last exon-exon junction
                let distance = junction_pos.saturating_sub(ptc_position);

                if distance > self.threshold_high {
                    NmdPrediction::Triggered {
                        distance_to_junction: distance,
                        confidence: NmdConfidence::High,
                    }
                } else if distance >= self.threshold_low {
                    NmdPrediction::Triggered {
                        distance_to_junction: distance,
                        confidence: NmdConfidence::Medium,
                    }
                } else {
                    NmdPrediction::NotTriggered {
                        reason: NmdNotTriggeredReason::PtcNearLastJunction { distance },
                    }
                }
            }
            None => NmdPrediction::Uncertain {
                reason: "Could not determine last exon-exon junction position".to_string(),
            },
        }
    }

    /// Predict NMD for a frameshift variant.
    ///
    /// # Arguments
    /// * `frameshift_position` - Position where frameshift occurs (CDS coordinate, 1-based).
    /// * `new_stop_position` - Position of the new stop codon after frameshift (if known).
    /// * `transcript` - Transcript information.
    pub fn predict_for_frameshift(
        &self,
        frameshift_position: u64,
        new_stop_position: Option<u64>,
        transcript: &TranscriptForNmd,
    ) -> NmdPrediction {
        match new_stop_position {
            Some(stop_pos) => self.predict(stop_pos, transcript),
            None => {
                // If we don't know where the new stop is, make a conservative prediction
                // based on the frameshift position
                NmdPrediction::Uncertain {
                    reason: format!(
                        "Frameshift at position {}; new stop codon position unknown",
                        frameshift_position
                    ),
                }
            }
        }
    }
}

impl Default for NmdPredictor {
    fn default() -> Self {
        Self::new()
    }
}

// =============================================================================
// Kozak Sequence Analysis
// =============================================================================

/// Kozak sequence strength.
///
/// The Kozak consensus sequence is (gcc)gccRccAUGG where:
/// - R is a purine (A or G) at position -3
/// - G at position +4 (after AUG)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KozakStrength {
    /// Strong Kozak context: A/G at -3 AND G at +4.
    Strong,
    /// Adequate/Moderate Kozak context: A/G at -3 OR G at +4 (but not both).
    Adequate,
    /// Weak Kozak context: neither A/G at -3 nor G at +4.
    Weak,
    /// Unknown - insufficient sequence information.
    Unknown,
}

/// Result of Kozak sequence analysis.
#[derive(Debug, Clone)]
pub struct KozakAnalysis {
    /// Kozak strength of the original sequence.
    pub original_strength: KozakStrength,
    /// Kozak strength after variant (if applicable).
    pub variant_strength: Option<KozakStrength>,
    /// Whether the variant affects the Kozak sequence.
    pub affects_kozak: bool,
    /// Whether the variant damages translation initiation.
    pub damages_initiation: bool,
    /// The Kozak sequence positions affected (-6 to +4 relative to ATG).
    pub positions_affected: Vec<i8>,
    /// Clinical interpretation.
    pub interpretation: KozakInterpretation,
}

/// Clinical interpretation of Kozak sequence variant effect.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KozakInterpretation {
    /// Variant has no effect on translation initiation.
    NoEffect,
    /// Variant may reduce translation efficiency.
    ReducedEfficiency,
    /// Variant likely abolishes translation from this start codon.
    LikelyAbolished,
    /// Effect is uncertain.
    Uncertain,
}

/// Kozak sequence analyzer.
#[derive(Debug, Clone)]
pub struct KozakAnalyzer {
    /// Kozak consensus positions to check (relative to A of ATG = position 1).
    /// Position -3 (R = A/G) and position +4 (G) are most critical.
    _consensus: Vec<(i8, Vec<char>)>,
}

impl KozakAnalyzer {
    /// Create a new Kozak analyzer.
    pub fn new() -> Self {
        Self {
            // Kozak consensus: (gcc)gccRccAUGG
            // Positions relative to A of ATG (A=1, T=2, G=3)
            // -6 to -1 is upstream context, +4 is downstream
            _consensus: vec![
                (-3, vec!['A', 'G']), // R (purine) at -3 is critical
                (4, vec!['G']),       // G at +4 is critical
            ],
        }
    }

    /// Check if a variant affects the Kozak sequence region.
    ///
    /// The Kozak sequence spans positions -6 to +4 relative to the ATG start codon.
    ///
    /// # Arguments
    /// * `variant_start` - Start position of variant in CDS coordinates (negative for 5' UTR).
    /// * `variant_end` - End position of variant.
    ///
    /// # Returns
    /// True if the variant overlaps the Kozak region.
    pub fn affects_kozak_region(&self, variant_start: i64, variant_end: i64) -> bool {
        // Kozak region is -6 to +4 relative to ATG (positions -6 to -1 in 5' UTR, 1-4 in CDS)
        let kozak_start: i64 = -6;
        let kozak_end: i64 = 4;

        variant_start <= kozak_end && variant_end >= kozak_start
    }

    /// Analyze the Kozak sequence context.
    ///
    /// # Arguments
    /// * `sequence` - Sequence around the start codon. Should include positions -6 to +4.
    ///   Index 0 should correspond to position -6.
    ///
    /// # Returns
    /// Kozak strength assessment.
    pub fn analyze_strength(&self, sequence: &str) -> KozakStrength {
        if sequence.len() < 10 {
            return KozakStrength::Unknown;
        }

        let chars: Vec<char> = sequence.chars().collect();

        // Position -3 is at index 3 (if sequence starts at -6)
        // Position +4 is at index 9
        let pos_minus3 = chars.get(3).map(|c| c.to_ascii_uppercase());
        let pos_plus4 = chars.get(9).map(|c| c.to_ascii_uppercase());

        let has_purine_minus3 = pos_minus3.is_some_and(|c| c == 'A' || c == 'G');
        let has_g_plus4 = pos_plus4 == Some('G');

        match (has_purine_minus3, has_g_plus4) {
            (true, true) => KozakStrength::Strong,
            (true, false) | (false, true) => KozakStrength::Adequate,
            (false, false) => KozakStrength::Weak,
        }
    }

    /// Analyze the effect of a variant on Kozak sequence.
    ///
    /// # Arguments
    /// * `variant_start` - Start position of variant relative to ATG (A=1).
    /// * `variant_end` - End position of variant.
    /// * `original_sequence` - Original sequence around start codon (-6 to +4).
    /// * `variant_sequence` - Sequence after variant (if known).
    pub fn analyze_variant(
        &self,
        variant_start: i64,
        variant_end: i64,
        original_sequence: Option<&str>,
        variant_sequence: Option<&str>,
    ) -> KozakAnalysis {
        let affects_kozak = self.affects_kozak_region(variant_start, variant_end);

        let original_strength = original_sequence
            .map(|s| self.analyze_strength(s))
            .unwrap_or(KozakStrength::Unknown);

        let variant_strength = variant_sequence.map(|s| self.analyze_strength(s));

        // Determine which positions are affected
        let mut positions_affected = Vec::new();
        for pos in -6..=4i8 {
            if (pos as i64) >= variant_start && (pos as i64) <= variant_end {
                positions_affected.push(pos);
            }
        }

        // Determine if initiation is damaged
        let damages_initiation = if affects_kozak {
            // Check if critical positions are affected
            let affects_critical =
                positions_affected.contains(&-3) || positions_affected.contains(&4);

            // Check if strength decreases
            let strength_decreased = matches!(
                (&original_strength, &variant_strength),
                (KozakStrength::Strong, Some(KozakStrength::Adequate))
                    | (KozakStrength::Strong, Some(KozakStrength::Weak))
                    | (KozakStrength::Adequate, Some(KozakStrength::Weak))
            );

            affects_critical || strength_decreased
        } else {
            false
        };

        let interpretation = if !affects_kozak {
            KozakInterpretation::NoEffect
        } else if damages_initiation {
            match variant_strength {
                Some(KozakStrength::Weak) => KozakInterpretation::LikelyAbolished,
                Some(KozakStrength::Adequate) => KozakInterpretation::ReducedEfficiency,
                Some(KozakStrength::Strong) => KozakInterpretation::NoEffect,
                _ => KozakInterpretation::Uncertain,
            }
        } else {
            KozakInterpretation::Uncertain
        };

        KozakAnalysis {
            original_strength,
            variant_strength,
            affects_kozak,
            damages_initiation,
            positions_affected,
            interpretation,
        }
    }
}

impl Default for KozakAnalyzer {
    fn default() -> Self {
        Self::new()
    }
}

// =============================================================================
// Inversion Sequence Validation
// =============================================================================

/// Inversion validation result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InversionValidation {
    /// Inversion is valid - sequence matches expected reverse complement.
    Valid,
    /// Inversion sequence doesn't match reference.
    SequenceMismatch {
        /// Expected sequence (reverse complement of reference).
        expected: String,
        /// Provided sequence in HGVS.
        provided: String,
    },
    /// Length mismatch between specified range and sequence.
    LengthMismatch {
        /// Expected length from range.
        range_length: usize,
        /// Provided sequence length.
        sequence_length: usize,
    },
    /// Cannot validate - reference sequence not available.
    CannotValidate {
        /// Reason.
        reason: String,
    },
}

/// Validate that an inversion sequence matches the reference.
///
/// # Arguments
/// * `reference_sequence` - The reference sequence at the inversion position.
/// * `provided_sequence` - The sequence provided in the HGVS notation (optional).
/// * `start` - Start position (1-based).
/// * `end` - End position (1-based, inclusive).
///
/// # Returns
/// Validation result.
pub fn validate_inversion(
    reference_sequence: Option<&str>,
    provided_sequence: Option<&str>,
    start: u64,
    end: u64,
) -> InversionValidation {
    let expected_length = (end - start + 1) as usize;

    // If a sequence is provided, check its length
    if let Some(provided) = provided_sequence {
        if provided.len() != expected_length {
            return InversionValidation::LengthMismatch {
                range_length: expected_length,
                sequence_length: provided.len(),
            };
        }
    }

    // If we have the reference, validate the provided sequence
    match (reference_sequence, provided_sequence) {
        (Some(reference), Some(provided)) => {
            // The provided sequence should be the reverse complement of reference
            let expected_revcomp = reverse_complement(reference);

            if provided.to_uppercase() == expected_revcomp.to_uppercase() {
                InversionValidation::Valid
            } else {
                InversionValidation::SequenceMismatch {
                    expected: expected_revcomp,
                    provided: provided.to_string(),
                }
            }
        }
        (Some(_reference), None) => {
            // No sequence provided, but we have reference - this is fine for validation
            InversionValidation::Valid
        }
        (None, Some(_)) => InversionValidation::CannotValidate {
            reason: "Reference sequence not available to validate inversion".to_string(),
        },
        (None, None) => InversionValidation::Valid, // Nothing to validate
    }
}

/// Compute reverse complement of a DNA sequence.
fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            other => other,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // =========================================================================
    // NMD Predictor Tests
    // =========================================================================

    #[test]
    fn test_nmd_triggered_high_confidence() {
        let predictor = NmdPredictor::new();
        let transcript = TranscriptForNmd {
            exon_boundaries: vec![(1, 100), (101, 200), (201, 300)],
            normal_stop_position: 300,
            cds_length: 300,
        };

        // PTC at position 100, last junction at 201 -> distance = 101 > 55
        let result = predictor.predict(100, &transcript);
        match result {
            NmdPrediction::Triggered {
                distance_to_junction,
                confidence,
            } => {
                assert_eq!(distance_to_junction, 101);
                assert_eq!(confidence, NmdConfidence::High);
            }
            _ => panic!("Expected NMD to be triggered with high confidence"),
        }
    }

    #[test]
    fn test_nmd_triggered_medium_confidence() {
        let predictor = NmdPredictor::new();
        let transcript = TranscriptForNmd {
            exon_boundaries: vec![(1, 100), (101, 200), (201, 300)],
            normal_stop_position: 300,
            cds_length: 300,
        };

        // PTC at position 148, last junction at 201 -> distance = 53 (between 50 and 55)
        let result = predictor.predict(148, &transcript);
        match result {
            NmdPrediction::Triggered {
                distance_to_junction,
                confidence,
            } => {
                assert_eq!(distance_to_junction, 53);
                assert_eq!(confidence, NmdConfidence::Medium);
            }
            _ => panic!("Expected NMD to be triggered with medium confidence"),
        }
    }

    #[test]
    fn test_nmd_not_triggered_near_junction() {
        let predictor = NmdPredictor::new();
        let transcript = TranscriptForNmd {
            exon_boundaries: vec![(1, 100), (101, 200), (201, 300)],
            normal_stop_position: 300,
            cds_length: 300,
        };

        // PTC at position 180, last junction at 201 -> distance = 21 < 50
        let result = predictor.predict(180, &transcript);
        match result {
            NmdPrediction::NotTriggered {
                reason: NmdNotTriggeredReason::PtcNearLastJunction { distance },
            } => {
                assert_eq!(distance, 21);
            }
            _ => panic!("Expected NMD not to be triggered due to proximity to junction"),
        }
    }

    #[test]
    fn test_nmd_not_triggered_last_exon() {
        let predictor = NmdPredictor::new();
        let transcript = TranscriptForNmd {
            exon_boundaries: vec![(1, 100), (101, 200), (201, 300)],
            normal_stop_position: 300,
            cds_length: 300,
        };

        // PTC at position 250, which is in the last exon (201-300)
        let result = predictor.predict(250, &transcript);
        assert_eq!(
            result,
            NmdPrediction::NotTriggered {
                reason: NmdNotTriggeredReason::PtcInLastExon
            }
        );
    }

    #[test]
    fn test_nmd_single_exon_transcript() {
        let predictor = NmdPredictor::new();
        let transcript = TranscriptForNmd {
            exon_boundaries: vec![(1, 300)],
            normal_stop_position: 300,
            cds_length: 300,
        };

        let result = predictor.predict(100, &transcript);
        assert_eq!(
            result,
            NmdPrediction::NotTriggered {
                reason: NmdNotTriggeredReason::SingleExonTranscript
            }
        );
    }

    #[test]
    fn test_nmd_ptc_after_normal_stop() {
        let predictor = NmdPredictor::new();
        let transcript = TranscriptForNmd {
            exon_boundaries: vec![(1, 100), (101, 200), (201, 300)],
            normal_stop_position: 300,
            cds_length: 300,
        };

        let result = predictor.predict(350, &transcript);
        assert_eq!(
            result,
            NmdPrediction::NotTriggered {
                reason: NmdNotTriggeredReason::PtcAtOrAfterNormalStop
            }
        );
    }

    // =========================================================================
    // Kozak Analyzer Tests
    // =========================================================================

    #[test]
    fn test_kozak_strong() {
        let analyzer = KozakAnalyzer::new();
        // Strong Kozak: A at -3 and G at +4
        // Sequence: -6 -5 -4 -3 -2 -1  A  T  G +4
        //            G  C  C  A  C  C  A  T  G  G
        let strength = analyzer.analyze_strength("GCCACCATGG");
        assert_eq!(strength, KozakStrength::Strong);
    }

    #[test]
    fn test_kozak_adequate_purine_only() {
        let analyzer = KozakAnalyzer::new();
        // Adequate: A at -3 but not G at +4
        let strength = analyzer.analyze_strength("GCCACCATGT");
        assert_eq!(strength, KozakStrength::Adequate);
    }

    #[test]
    fn test_kozak_adequate_g_plus4_only() {
        let analyzer = KozakAnalyzer::new();
        // Adequate: G at +4 but not purine at -3
        let strength = analyzer.analyze_strength("GCCTCCATGG");
        assert_eq!(strength, KozakStrength::Adequate);
    }

    #[test]
    fn test_kozak_weak() {
        let analyzer = KozakAnalyzer::new();
        // Weak: neither purine at -3 nor G at +4
        let strength = analyzer.analyze_strength("GCCTCCATGT");
        assert_eq!(strength, KozakStrength::Weak);
    }

    #[test]
    fn test_kozak_affects_region() {
        let analyzer = KozakAnalyzer::new();

        // Variant at position -5 (within Kozak region)
        assert!(analyzer.affects_kozak_region(-5, -5));

        // Variant at position 3 (within Kozak region)
        assert!(analyzer.affects_kozak_region(3, 3));

        // Variant spanning -3 to +4 (within Kozak region)
        assert!(analyzer.affects_kozak_region(-3, 4));

        // Variant at position -10 (outside Kozak region)
        assert!(!analyzer.affects_kozak_region(-10, -10));

        // Variant at position 10 (outside Kozak region)
        assert!(!analyzer.affects_kozak_region(10, 10));
    }

    #[test]
    fn test_kozak_variant_analysis() {
        let analyzer = KozakAnalyzer::new();

        // Variant in Kozak region (-5 to -1 deletion)
        let analysis = analyzer.analyze_variant(
            -5,
            -1,
            Some("GCCACCATGG"), // Strong original
            Some("GCCATATGG"),  // Weak after variant (hypothetical)
        );

        assert!(analysis.affects_kozak);
        assert_eq!(analysis.original_strength, KozakStrength::Strong);
    }

    // =========================================================================
    // Inversion Validation Tests
    // =========================================================================

    #[test]
    fn test_inversion_valid() {
        // Reference: ATGC, reverse complement: GCAT
        let result = validate_inversion(Some("ATGC"), Some("GCAT"), 1, 4);
        assert_eq!(result, InversionValidation::Valid);
    }

    #[test]
    fn test_inversion_sequence_mismatch() {
        // Reference: ATGC, reverse complement should be GCAT, but provided AAAA
        let result = validate_inversion(Some("ATGC"), Some("AAAA"), 1, 4);
        match result {
            InversionValidation::SequenceMismatch { expected, provided } => {
                assert_eq!(expected, "GCAT");
                assert_eq!(provided, "AAAA");
            }
            _ => panic!("Expected sequence mismatch"),
        }
    }

    #[test]
    fn test_inversion_length_mismatch() {
        let result = validate_inversion(Some("ATGC"), Some("GC"), 1, 4);
        match result {
            InversionValidation::LengthMismatch {
                range_length,
                sequence_length,
            } => {
                assert_eq!(range_length, 4);
                assert_eq!(sequence_length, 2);
            }
            _ => panic!("Expected length mismatch"),
        }
    }

    #[test]
    fn test_inversion_no_sequence_provided() {
        let result = validate_inversion(Some("ATGC"), None, 1, 4);
        assert_eq!(result, InversionValidation::Valid);
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATGC"), "GCAT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("GCGC"), "GCGC");
        assert_eq!(reverse_complement("A"), "T");
    }

    // =========================================================================
    // Original Tests
    // =========================================================================

    #[test]
    fn test_consequence_impact() {
        assert_eq!(Consequence::StopGained.impact(), Impact::High);
        assert_eq!(Consequence::MissenseVariant.impact(), Impact::Moderate);
        assert_eq!(Consequence::SynonymousVariant.impact(), Impact::Low);
        assert_eq!(Consequence::IntronVariant.impact(), Impact::Modifier);
    }

    #[test]
    fn test_consequence_so_terms() {
        assert_eq!(Consequence::MissenseVariant.so_term(), "missense_variant");
        assert_eq!(Consequence::MissenseVariant.so_id(), "SO:0001583");
    }

    #[test]
    fn test_classify_missense() {
        let predictor = EffectPredictor::new();
        let effect = predictor.classify_amino_acid_change(&AminoAcid::Val, &AminoAcid::Glu, 600);

        assert_eq!(effect.consequences[0], Consequence::MissenseVariant);
        assert_eq!(effect.impact, Impact::Moderate);
    }

    #[test]
    fn test_classify_synonymous() {
        let predictor = EffectPredictor::new();
        let effect = predictor.classify_amino_acid_change(&AminoAcid::Val, &AminoAcid::Val, 100);

        assert_eq!(effect.consequences[0], Consequence::SynonymousVariant);
        assert_eq!(effect.impact, Impact::Low);
    }

    #[test]
    fn test_classify_stop_gained() {
        let predictor = EffectPredictor::new();
        let effect = predictor.classify_amino_acid_change(&AminoAcid::Gln, &AminoAcid::Ter, 100);

        assert_eq!(effect.consequences[0], Consequence::StopGained);
        assert_eq!(effect.impact, Impact::High);
    }

    #[test]
    fn test_classify_stop_lost() {
        let predictor = EffectPredictor::new();
        let effect = predictor.classify_amino_acid_change(&AminoAcid::Ter, &AminoAcid::Gln, 500);

        assert_eq!(effect.consequences[0], Consequence::StopLost);
        assert_eq!(effect.impact, Impact::High);
    }

    #[test]
    fn test_classify_start_lost() {
        let predictor = EffectPredictor::new();
        let effect = predictor.classify_amino_acid_change(&AminoAcid::Met, &AminoAcid::Val, 1);

        assert_eq!(effect.consequences[0], Consequence::StartLost);
        assert_eq!(effect.impact, Impact::High);
    }

    #[test]
    fn test_classify_frameshift() {
        let predictor = EffectPredictor::new();

        // 1bp deletion = frameshift
        let effect = predictor.classify_indel(1, 0);
        assert_eq!(effect.consequences[0], Consequence::FrameshiftVariant);
        assert_eq!(effect.impact, Impact::High);

        // 2bp insertion = frameshift
        let effect = predictor.classify_indel(0, 2);
        assert_eq!(effect.consequences[0], Consequence::FrameshiftVariant);
    }

    #[test]
    fn test_classify_inframe() {
        let predictor = EffectPredictor::new();

        // 3bp deletion = in-frame
        let effect = predictor.classify_indel(3, 0);
        assert_eq!(effect.consequences[0], Consequence::InframeDeletion);
        assert_eq!(effect.impact, Impact::Moderate);

        // 6bp insertion = in-frame
        let effect = predictor.classify_indel(0, 6);
        assert_eq!(effect.consequences[0], Consequence::InframeInsertion);
    }

    #[test]
    fn test_classify_splice_donor() {
        let predictor = EffectPredictor::new();
        let effect = predictor.classify_splice_variant(1);

        assert_eq!(effect.consequences[0], Consequence::SpliceDonorVariant);
        assert_eq!(effect.impact, Impact::High);
    }

    #[test]
    fn test_classify_splice_acceptor() {
        let predictor = EffectPredictor::new();
        let effect = predictor.classify_splice_variant(-2);

        assert_eq!(effect.consequences[0], Consequence::SpliceAcceptorVariant);
        assert_eq!(effect.impact, Impact::High);
    }

    #[test]
    fn test_classify_splice_region() {
        let predictor = EffectPredictor::new();
        let effect = predictor.classify_splice_variant(5);

        assert_eq!(effect.consequences[0], Consequence::SpliceRegionVariant);
        assert_eq!(effect.impact, Impact::Low);
    }

    #[test]
    fn test_classify_intron() {
        let predictor = EffectPredictor::new();
        let effect = predictor.classify_splice_variant(50);

        assert_eq!(effect.consequences[0], Consequence::IntronVariant);
        assert_eq!(effect.impact, Impact::Modifier);
    }

    #[test]
    fn test_classify_utr() {
        let predictor = EffectPredictor::new();

        let effect = predictor.classify_utr_variant(true);
        assert_eq!(effect.consequences[0], Consequence::FivePrimeUtrVariant);

        let effect = predictor.classify_utr_variant(false);
        assert_eq!(effect.consequences[0], Consequence::ThreePrimeUtrVariant);
    }

    #[test]
    fn test_protein_effect_methods() {
        let effect = ProteinEffect {
            consequences: vec![Consequence::StopGained, Consequence::SpliceRegionVariant],
            impact: Impact::High,
            amino_acid_change: None,
            intronic_offset: None,
        };

        assert_eq!(effect.most_severe(), Some(&Consequence::StopGained));
        assert!(effect.is_high_impact());
        assert!(effect.is_protein_altering());
    }

    #[test]
    fn test_impact_ordering() {
        assert!(Impact::High > Impact::Moderate);
        assert!(Impact::Moderate > Impact::Low);
        assert!(Impact::Low > Impact::Modifier);
    }
}
