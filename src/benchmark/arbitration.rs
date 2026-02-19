//! Arbitration system for comparing normalization results between tools.
//!
//! This module provides a framework for:
//! - Recording differences between ferro and mutalyzer
//! - Categorizing each difference (correct, equivalent, incorrect, unknown)
//! - Storing arbitration decisions persistently
//! - Generating reports on normalization quality

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use crate::error::FerroError;

/// Category of an arbitration decision
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ArbitrationCategory {
    /// Ferro output is correct, mutalyzer is wrong
    FerroCorrect,
    /// Mutalyzer output is correct, ferro is wrong
    MutalyzerCorrect,
    /// Both outputs are equivalent (valid alternative representations)
    Equivalent,
    /// Both outputs are incorrect
    BothIncorrect,
    /// Unable to determine which is correct
    Unknown,
}

impl std::fmt::Display for ArbitrationCategory {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ArbitrationCategory::FerroCorrect => write!(f, "ferro_correct"),
            ArbitrationCategory::MutalyzerCorrect => write!(f, "mutalyzer_correct"),
            ArbitrationCategory::Equivalent => write!(f, "equivalent"),
            ArbitrationCategory::BothIncorrect => write!(f, "both_incorrect"),
            ArbitrationCategory::Unknown => write!(f, "unknown"),
        }
    }
}

/// An arbitration decision for a single variant
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ArbitrationDecision {
    /// The input HGVS variant
    pub input: String,
    /// Ferro's normalized output
    pub ferro_output: String,
    /// Mutalyzer's normalized output
    pub mutalyzer_output: String,
    /// The arbitration category
    pub category: ArbitrationCategory,
    /// Explanation for the decision
    pub reason: String,
    /// Optional reference to HGVS spec section
    pub spec_reference: Option<String>,
    /// When this decision was made
    pub timestamp: String,
}

/// Collection of arbitration decisions
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ArbitrationDatabase {
    /// Version of the database format
    pub version: String,
    /// Decisions keyed by input variant
    pub decisions: HashMap<String, ArbitrationDecision>,
    /// Summary statistics
    #[serde(skip_serializing_if = "Option::is_none")]
    pub stats: Option<ArbitrationStats>,
}

/// Summary statistics for arbitration decisions
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ArbitrationStats {
    pub total: usize,
    pub ferro_correct: usize,
    pub mutalyzer_correct: usize,
    pub equivalent: usize,
    pub both_incorrect: usize,
    pub unknown: usize,
}

impl ArbitrationDatabase {
    /// Create a new empty database
    pub fn new() -> Self {
        Self {
            version: "1.0".to_string(),
            decisions: HashMap::new(),
            stats: None,
        }
    }

    /// Load database from JSON file
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let path = path.as_ref();
        if !path.exists() {
            return Ok(Self::new());
        }
        let file = File::open(path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open {}: {}", path.display(), e),
        })?;
        let reader = BufReader::new(file);
        serde_json::from_reader(reader).map_err(|e| FerroError::Io {
            msg: format!("Failed to parse {}: {}", path.display(), e),
        })
    }

    /// Save database to JSON file
    pub fn save<P: AsRef<Path>>(&mut self, path: P) -> Result<(), FerroError> {
        self.update_stats();
        let path = path.as_ref();
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
                msg: format!("Failed to create directory: {}", e),
            })?;
        }
        let file = File::create(path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create {}: {}", path.display(), e),
        })?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, self).map_err(|e| FerroError::Io {
            msg: format!("Failed to write {}: {}", path.display(), e),
        })
    }

    /// Add or update an arbitration decision
    pub fn add_decision(&mut self, decision: ArbitrationDecision) {
        self.decisions.insert(decision.input.clone(), decision);
    }

    /// Get a decision by input variant
    pub fn get_decision(&self, input: &str) -> Option<&ArbitrationDecision> {
        self.decisions.get(input)
    }

    /// Check if we have a decision for an input
    pub fn has_decision(&self, input: &str) -> bool {
        self.decisions.contains_key(input)
    }

    /// Update statistics
    pub fn update_stats(&mut self) {
        let mut stats = ArbitrationStats {
            total: self.decisions.len(),
            ..Default::default()
        };
        for decision in self.decisions.values() {
            match decision.category {
                ArbitrationCategory::FerroCorrect => stats.ferro_correct += 1,
                ArbitrationCategory::MutalyzerCorrect => stats.mutalyzer_correct += 1,
                ArbitrationCategory::Equivalent => stats.equivalent += 1,
                ArbitrationCategory::BothIncorrect => stats.both_incorrect += 1,
                ArbitrationCategory::Unknown => stats.unknown += 1,
            }
        }
        self.stats = Some(stats);
    }

    /// Get summary statistics
    pub fn stats(&self) -> ArbitrationStats {
        self.stats.clone().unwrap_or_else(|| {
            let mut stats = ArbitrationStats {
                total: self.decisions.len(),
                ..Default::default()
            };
            for decision in self.decisions.values() {
                match decision.category {
                    ArbitrationCategory::FerroCorrect => stats.ferro_correct += 1,
                    ArbitrationCategory::MutalyzerCorrect => stats.mutalyzer_correct += 1,
                    ArbitrationCategory::Equivalent => stats.equivalent += 1,
                    ArbitrationCategory::BothIncorrect => stats.both_incorrect += 1,
                    ArbitrationCategory::Unknown => stats.unknown += 1,
                }
            }
            stats
        })
    }

    /// Get decisions by category
    pub fn by_category(&self, category: ArbitrationCategory) -> Vec<&ArbitrationDecision> {
        self.decisions
            .values()
            .filter(|d| d.category == category)
            .collect()
    }
}

/// Builder for creating arbitration decisions
pub struct ArbitrationBuilder {
    input: String,
    ferro_output: String,
    mutalyzer_output: String,
    category: Option<ArbitrationCategory>,
    reason: String,
    spec_reference: Option<String>,
}

impl ArbitrationBuilder {
    pub fn new(input: &str, ferro_output: &str, mutalyzer_output: &str) -> Self {
        Self {
            input: input.to_string(),
            ferro_output: ferro_output.to_string(),
            mutalyzer_output: mutalyzer_output.to_string(),
            category: None,
            reason: String::new(),
            spec_reference: None,
        }
    }

    pub fn ferro_correct(mut self, reason: &str) -> Self {
        self.category = Some(ArbitrationCategory::FerroCorrect);
        self.reason = reason.to_string();
        self
    }

    pub fn mutalyzer_correct(mut self, reason: &str) -> Self {
        self.category = Some(ArbitrationCategory::MutalyzerCorrect);
        self.reason = reason.to_string();
        self
    }

    pub fn equivalent(mut self, reason: &str) -> Self {
        self.category = Some(ArbitrationCategory::Equivalent);
        self.reason = reason.to_string();
        self
    }

    pub fn both_incorrect(mut self, reason: &str) -> Self {
        self.category = Some(ArbitrationCategory::BothIncorrect);
        self.reason = reason.to_string();
        self
    }

    pub fn unknown(mut self, reason: &str) -> Self {
        self.category = Some(ArbitrationCategory::Unknown);
        self.reason = reason.to_string();
        self
    }

    pub fn spec_ref(mut self, reference: &str) -> Self {
        self.spec_reference = Some(reference.to_string());
        self
    }

    pub fn build(self) -> ArbitrationDecision {
        ArbitrationDecision {
            input: self.input,
            ferro_output: self.ferro_output,
            mutalyzer_output: self.mutalyzer_output,
            category: self.category.unwrap_or(ArbitrationCategory::Unknown),
            reason: self.reason,
            spec_reference: self.spec_reference,
            timestamp: chrono::Utc::now().to_rfc3339(),
        }
    }
}

/// Apply known arbitrations to a comparison result
pub fn apply_arbitrations(
    db: &ArbitrationDatabase,
    differences: &[(String, String, String)], // (input, ferro, mutalyzer)
) -> ArbitrationSummary {
    let mut summary = ArbitrationSummary::default();

    for (input, ferro, mutalyzer) in differences {
        if let Some(decision) = db.get_decision(input) {
            match decision.category {
                ArbitrationCategory::FerroCorrect => summary.ferro_correct.push(input.clone()),
                ArbitrationCategory::MutalyzerCorrect => {
                    summary.mutalyzer_correct.push(input.clone())
                }
                ArbitrationCategory::Equivalent => summary.equivalent.push(input.clone()),
                ArbitrationCategory::BothIncorrect => summary.both_incorrect.push(input.clone()),
                ArbitrationCategory::Unknown => summary.unknown.push(input.clone()),
            }
        } else {
            summary
                .unarbitrated
                .push((input.clone(), ferro.clone(), mutalyzer.clone()));
        }
    }

    summary
}

/// Summary of applying arbitrations to differences
#[derive(Debug, Default)]
pub struct ArbitrationSummary {
    pub ferro_correct: Vec<String>,
    pub mutalyzer_correct: Vec<String>,
    pub equivalent: Vec<String>,
    pub both_incorrect: Vec<String>,
    pub unknown: Vec<String>,
    pub unarbitrated: Vec<(String, String, String)>,
}

impl ArbitrationSummary {
    pub fn print_report(&self) {
        println!("\n=== Arbitration Summary ===");
        println!("Ferro correct:     {}", self.ferro_correct.len());
        println!("Mutalyzer correct: {}", self.mutalyzer_correct.len());
        println!("Equivalent:        {}", self.equivalent.len());
        println!("Both incorrect:    {}", self.both_incorrect.len());
        println!("Unknown:           {}", self.unknown.len());
        println!("Unarbitrated:      {}", self.unarbitrated.len());

        if !self.unarbitrated.is_empty() {
            println!("\nUnarbitrated differences:");
            for (i, (input, ferro, mutalyzer)) in self.unarbitrated.iter().take(10).enumerate() {
                println!("  {}. {}", i + 1, input);
                println!("     ferro:     {}", ferro);
                println!("     mutalyzer: {}", mutalyzer);
            }
            if self.unarbitrated.len() > 10 {
                println!("  ... and {} more", self.unarbitrated.len() - 10);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arbitration_database() {
        let mut db = ArbitrationDatabase::new();

        let decision = ArbitrationBuilder::new(
            "NM_000255.4:c.445_446insA",
            "NM_000255.4:c.446dup",
            "NM_000255.4:c.445dup",
        )
        .ferro_correct("Insertion at c.445_446 of A should shift 3' to c.446 since c.446=A")
        .spec_ref("HGVS 3.1: 3' rule for insertions")
        .build();

        db.add_decision(decision);

        assert!(db.has_decision("NM_000255.4:c.445_446insA"));
        assert!(!db.has_decision("NM_000255.4:c.100A>G"));

        let stats = db.stats();
        assert_eq!(stats.total, 1);
        assert_eq!(stats.ferro_correct, 1);
    }
}
