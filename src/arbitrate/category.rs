use serde::{Deserialize, Serialize};

/// Category of an arbitration decision (who is correct).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ArbitrationCategory {
    /// Ferro output is correct, the other tool is wrong.
    FerroCorrect,
    /// The other tool's output is correct, ferro is wrong.
    MutalyzerCorrect,
    /// Both outputs are equivalent (valid alternative representations).
    Equivalent,
    /// Both outputs are incorrect.
    BothIncorrect,
    /// Unable to determine which is correct.
    Unknown,
}

impl std::fmt::Display for ArbitrationCategory {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            ArbitrationCategory::FerroCorrect => "ferro_correct",
            ArbitrationCategory::MutalyzerCorrect => "mutalyzer_correct",
            ArbitrationCategory::Equivalent => "equivalent",
            ArbitrationCategory::BothIncorrect => "both_incorrect",
            ArbitrationCategory::Unknown => "unknown",
        };
        write!(f, "{s}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_matches_serde_snake_case() {
        assert_eq!(
            ArbitrationCategory::FerroCorrect.to_string(),
            "ferro_correct"
        );
        let json = serde_json::to_string(&ArbitrationCategory::MutalyzerCorrect).unwrap();
        assert_eq!(json, "\"mutalyzer_correct\"");
    }
}
