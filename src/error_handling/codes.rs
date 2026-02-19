//! Error and warning code definitions with rich documentation.
//!
//! This module provides comprehensive information about each error (E-codes)
//! and warning (W-codes) in ferro-hgvs, including examples, explanations,
//! and links to relevant documentation.

use std::fmt;

/// The action taken for an error/warning in a given mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModeAction {
    /// Reject the input with an error.
    Reject,
    /// Auto-correct and emit a warning.
    WarnAndCorrect,
    /// Auto-correct silently.
    SilentCorrect,
    /// Warn but accept as-is (no correction).
    WarnAndAccept,
    /// Accept as-is without warning or correction.
    Accept,
}

impl fmt::Display for ModeAction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ModeAction::Reject => write!(f, "Rejected"),
            ModeAction::WarnAndCorrect => write!(f, "Warn + Correct"),
            ModeAction::SilentCorrect => write!(f, "Silent Correct"),
            ModeAction::WarnAndAccept => write!(f, "Warn + Accept"),
            ModeAction::Accept => write!(f, "Accept"),
        }
    }
}

/// Behavior of an error/warning in each mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ModeBehavior {
    /// Behavior in Strict mode.
    pub strict: ModeAction,
    /// Behavior in Lenient mode.
    pub lenient: ModeAction,
    /// Behavior in Silent mode.
    pub silent: ModeAction,
}

impl ModeBehavior {
    /// Create a new mode behavior specification.
    pub const fn new(strict: ModeAction, lenient: ModeAction, silent: ModeAction) -> Self {
        Self {
            strict,
            lenient,
            silent,
        }
    }

    /// Standard behavior: Reject in strict, warn+correct in lenient, silent correct in silent.
    pub const fn standard_correctable() -> Self {
        Self::new(
            ModeAction::Reject,
            ModeAction::WarnAndCorrect,
            ModeAction::SilentCorrect,
        )
    }

    /// Always reject regardless of mode.
    pub const fn always_reject() -> Self {
        Self::new(ModeAction::Reject, ModeAction::Reject, ModeAction::Reject)
    }

    /// Reject in strict, warn+accept in lenient/silent (for things that can't be auto-corrected).
    pub const fn warn_accept() -> Self {
        Self::new(
            ModeAction::Reject,
            ModeAction::WarnAndAccept,
            ModeAction::Accept,
        )
    }

    /// Reject in strict, warn in lenient, warn in silent (always warn if not rejected).
    pub const fn always_warn_if_not_rejected() -> Self {
        Self::new(
            ModeAction::Reject,
            ModeAction::WarnAndAccept,
            ModeAction::WarnAndAccept,
        )
    }
}

/// Category of error or warning code.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CodeCategory {
    /// Parse errors (E1xxx)
    Parse,
    /// Reference/data errors (E2xxx)
    Reference,
    /// Validation errors (E3xxx)
    Validation,
    /// Normalization errors (E4xxx)
    Normalization,
    /// Conversion errors (E5xxx)
    Conversion,
    /// IO errors (E9xxx)
    Io,
    /// Case/capitalization warnings (W1xxx)
    Case,
    /// Character warnings (W2xxx)
    Character,
    /// Format/syntax warnings (W3xxx)
    Format,
    /// Position/range warnings (W4xxx)
    Position,
    /// Semantic warnings (W5xxx)
    Semantic,
}

impl CodeCategory {
    /// Returns the display name for this category.
    pub fn name(&self) -> &'static str {
        match self {
            CodeCategory::Parse => "Parse Errors",
            CodeCategory::Reference => "Reference Errors",
            CodeCategory::Validation => "Validation Errors",
            CodeCategory::Normalization => "Normalization Errors",
            CodeCategory::Conversion => "Conversion Errors",
            CodeCategory::Io => "IO Errors",
            CodeCategory::Case => "Case/Capitalization",
            CodeCategory::Character => "Character Issues",
            CodeCategory::Format => "Format/Syntax",
            CodeCategory::Position => "Position/Range",
            CodeCategory::Semantic => "Semantic Issues",
        }
    }
}

impl fmt::Display for CodeCategory {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}

/// Complete information about an error or warning code.
#[derive(Debug, Clone)]
pub struct CodeInfo {
    /// Code identifier (e.g., "E1001", "W2001").
    pub code: &'static str,
    /// Short name (e.g., "InvalidAccession", "WrongDashCharacter").
    pub name: &'static str,
    /// One-line summary.
    pub summary: &'static str,
    /// Detailed multi-line explanation.
    pub explanation: &'static str,
    /// Category of this code.
    pub category: CodeCategory,
    /// Examples of incorrect input.
    pub bad_examples: &'static [&'static str],
    /// Examples of correct input (or corrected output).
    pub good_examples: &'static [&'static str],
    /// Behavior in each error mode (None for pure errors that aren't configurable).
    pub mode_behavior: Option<ModeBehavior>,
    /// HGVS specification URL (if applicable).
    pub hgvs_spec_url: Option<&'static str>,
    /// Related error/warning codes.
    pub related_codes: &'static [&'static str],
}

impl CodeInfo {
    /// Returns true if this is an error code (E-prefix).
    pub fn is_error(&self) -> bool {
        self.code.starts_with('E')
    }

    /// Returns true if this is a warning code (W-prefix).
    pub fn is_warning(&self) -> bool {
        self.code.starts_with('W')
    }

    /// Returns the documentation URL for this code.
    pub fn docs_url(&self) -> String {
        format!(
            "https://docs.rs/ferro-hgvs/latest/ferro_hgvs/error_handling/#{}",
            self.code
        )
    }

    /// Format this code info for terminal display.
    pub fn format_terminal(&self, use_color: bool) -> String {
        let mut output = String::new();

        // Header line
        let separator = "━".repeat(78);
        if use_color {
            output.push_str(&format!("\x1b[1m{}\x1b[0m\n", separator));
            output.push_str(&format!(
                "\x1b[1;33m{}\x1b[0m: \x1b[1m{}\x1b[0m\n",
                self.code, self.name
            ));
            output.push_str(&format!("\x1b[1m{}\x1b[0m\n\n", separator));
        } else {
            output.push_str(&format!("{}\n", separator));
            output.push_str(&format!("{}: {}\n", self.code, self.name));
            output.push_str(&format!("{}\n\n", separator));
        }

        // Summary and explanation
        output.push_str(self.summary);
        output.push_str("\n\n");
        output.push_str(self.explanation);
        output.push_str("\n\n");

        // Mode behavior (for warnings)
        if let Some(behavior) = &self.mode_behavior {
            output.push_str("Mode Behavior:\n");
            output.push_str(&format!("  • Strict:  {}\n", behavior.strict));
            output.push_str(&format!("  • Lenient: {}\n", behavior.lenient));
            output.push_str(&format!("  • Silent:  {}\n", behavior.silent));
            output.push('\n');
        }

        // Examples
        if !self.bad_examples.is_empty() {
            output.push_str("Examples:\n");
            for (bad, good) in self.bad_examples.iter().zip(self.good_examples.iter()) {
                if use_color {
                    output.push_str(&format!("  \x1b[31m✗\x1b[0m Incorrect: {}\n", bad));
                    output.push_str(&format!("  \x1b[32m✓\x1b[0m Correct:   {}\n\n", good));
                } else {
                    output.push_str(&format!("  ✗ Incorrect: {}\n", bad));
                    output.push_str(&format!("  ✓ Correct:   {}\n\n", good));
                }
            }
        }

        // Ignore instructions (for warnings)
        if self.is_warning() {
            output.push_str("To ignore this warning:\n");
            output.push_str(&format!(
                "  ferro parse --ignore {} \"<input>\"\n\n",
                self.code
            ));
            output.push_str("Or in .ferro.toml:\n");
            output.push_str("  [error-handling]\n");
            output.push_str(&format!("  ignore = [\"{}\"]\n\n", self.code));
        }

        // Related codes
        if !self.related_codes.is_empty() {
            output.push_str("Related Codes:\n");
            for related in self.related_codes {
                output.push_str(&format!("  {}\n", related));
            }
            output.push('\n');
        }

        // URLs
        if let Some(hgvs_url) = self.hgvs_spec_url {
            output.push_str("HGVS Specification:\n");
            output.push_str(&format!("  {}\n\n", hgvs_url));
        }

        output.push_str("Documentation:\n");
        output.push_str(&format!("  {}\n", self.docs_url()));

        if use_color {
            output.push_str(&format!("\x1b[1m{}\x1b[0m\n", separator));
        } else {
            output.push_str(&format!("{}\n", separator));
        }

        output
    }

    /// Format this code info as JSON.
    pub fn format_json(&self) -> String {
        let mode_behavior = self.mode_behavior.map(|b| {
            format!(
                r#"{{"strict":"{}","lenient":"{}","silent":"{}"}}"#,
                b.strict, b.lenient, b.silent
            )
        });

        let bad_examples: Vec<String> = self
            .bad_examples
            .iter()
            .map(|e| format!("\"{}\"", e.replace('\"', "\\\"")))
            .collect();
        let good_examples: Vec<String> = self
            .good_examples
            .iter()
            .map(|e| format!("\"{}\"", e.replace('\"', "\\\"")))
            .collect();
        let related: Vec<String> = self
            .related_codes
            .iter()
            .map(|c| format!("\"{}\"", c))
            .collect();

        format!(
            r#"{{"code":"{}","name":"{}","summary":"{}","category":"{}","bad_examples":[{}],"good_examples":[{}],"mode_behavior":{},"hgvs_spec_url":{},"related_codes":[{}],"docs_url":"{}"}}"#,
            self.code,
            self.name,
            self.summary.replace('\"', "\\\""),
            self.category,
            bad_examples.join(","),
            good_examples.join(","),
            mode_behavior.unwrap_or_else(|| "null".to_string()),
            self.hgvs_spec_url
                .map(|u| format!("\"{}\"", u))
                .unwrap_or_else(|| "null".to_string()),
            related.join(","),
            self.docs_url()
        )
    }

    /// Format this code info as Markdown.
    pub fn format_markdown(&self) -> String {
        let mut output = String::new();

        output.push_str(&format!("## {}: {}\n\n", self.code, self.name));
        output.push_str(&format!("**Category:** {}\n\n", self.category));
        output.push_str(&format!("{}\n\n", self.summary));
        output.push_str(&format!("{}\n\n", self.explanation));

        if let Some(behavior) = &self.mode_behavior {
            output.push_str("### Mode Behavior\n\n");
            output.push_str("| Mode | Action |\n");
            output.push_str("|------|--------|\n");
            output.push_str(&format!("| Strict | {} |\n", behavior.strict));
            output.push_str(&format!("| Lenient | {} |\n", behavior.lenient));
            output.push_str(&format!("| Silent | {} |\n\n", behavior.silent));
        }

        if !self.bad_examples.is_empty() {
            output.push_str("### Examples\n\n");
            for (bad, good) in self.bad_examples.iter().zip(self.good_examples.iter()) {
                output.push_str(&format!("- ❌ `{}`\n", bad));
                output.push_str(&format!("- ✅ `{}`\n\n", good));
            }
        }

        if !self.related_codes.is_empty() {
            output.push_str("### Related Codes\n\n");
            for related in self.related_codes {
                output.push_str(&format!("- {}\n", related));
            }
            output.push('\n');
        }

        if let Some(url) = self.hgvs_spec_url {
            output.push_str(&format!("**HGVS Spec:** {}\n\n", url));
        }
        output.push_str(&format!("**Documentation:** {}\n", self.docs_url()));

        output
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mode_action_display() {
        assert_eq!(format!("{}", ModeAction::Reject), "Rejected");
        assert_eq!(format!("{}", ModeAction::WarnAndCorrect), "Warn + Correct");
        assert_eq!(format!("{}", ModeAction::SilentCorrect), "Silent Correct");
    }

    #[test]
    fn test_mode_behavior_standard() {
        let behavior = ModeBehavior::standard_correctable();
        assert_eq!(behavior.strict, ModeAction::Reject);
        assert_eq!(behavior.lenient, ModeAction::WarnAndCorrect);
        assert_eq!(behavior.silent, ModeAction::SilentCorrect);
    }

    #[test]
    fn test_mode_behavior_always_reject() {
        let behavior = ModeBehavior::always_reject();
        assert_eq!(behavior.strict, ModeAction::Reject);
        assert_eq!(behavior.lenient, ModeAction::Reject);
        assert_eq!(behavior.silent, ModeAction::Reject);
    }

    #[test]
    fn test_code_category_display() {
        assert_eq!(format!("{}", CodeCategory::Parse), "Parse Errors");
        assert_eq!(format!("{}", CodeCategory::Case), "Case/Capitalization");
    }

    #[test]
    fn test_code_info_is_error() {
        let code = CodeInfo {
            code: "E1001",
            name: "Test",
            summary: "Test",
            explanation: "Test",
            category: CodeCategory::Parse,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &[],
        };
        assert!(code.is_error());
        assert!(!code.is_warning());
    }

    #[test]
    fn test_code_info_is_warning() {
        let code = CodeInfo {
            code: "W1001",
            name: "Test",
            summary: "Test",
            explanation: "Test",
            category: CodeCategory::Case,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: None,
            related_codes: &[],
        };
        assert!(!code.is_error());
        assert!(code.is_warning());
    }

    #[test]
    fn test_code_info_docs_url() {
        let code = CodeInfo {
            code: "W1001",
            name: "Test",
            summary: "Test",
            explanation: "Test",
            category: CodeCategory::Case,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &[],
        };
        assert!(code.docs_url().contains("W1001"));
    }
}
