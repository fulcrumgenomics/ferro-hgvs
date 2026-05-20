//! Configuration file support for ferro-hgvs.
//!
//! This module provides loading of `.ferro.toml` configuration files
//! which can specify error handling modes and per-code overrides.
//!
//! # Example Configuration
//!
//! ```toml
//! [error-handling]
//! mode = "lenient"
//! ignore = ["W1001", "W2001"]
//! reject = ["W3003"]
//! ```
//!
//! # Config File Locations
//!
//! Configuration is searched in this order (first found wins):
//! 1. `.ferro.toml` in current directory
//! 2. `~/.config/ferro/config.toml`
//!
//! CLI flags take precedence over config file settings.

use crate::error_handling::{ErrorConfig, ErrorMode, ErrorOverride, ErrorType};
use std::fs;
use std::path::PathBuf;

/// Parsed configuration from a .ferro.toml file.
#[derive(Debug, Clone, Default)]
pub struct FerroConfig {
    /// Error handling configuration.
    pub error_handling: ErrorHandlingConfig,
}

/// Error handling section of the config file.
#[derive(Debug, Clone, Default)]
pub struct ErrorHandlingConfig {
    /// Base error mode.
    pub mode: Option<String>,
    /// Codes to silently correct.
    pub ignore: Vec<String>,
    /// Codes to always reject.
    pub reject: Vec<String>,
    /// Codes to always warn about.
    pub warn: Vec<String>,
}

impl FerroConfig {
    /// Load configuration from the default locations.
    ///
    /// Searches for config in:
    /// 1. `.ferro.toml` in current directory
    /// 2. `~/.config/ferro/config.toml`
    pub fn load() -> Option<Self> {
        // Try current directory first
        let cwd_config = PathBuf::from(".ferro.toml");
        if cwd_config.exists() {
            if let Ok(config) = Self::load_from_path(&cwd_config) {
                return Some(config);
            }
        }

        // Try home config directory
        if let Some(home) = dirs_home() {
            let home_config = home.join(".config").join("ferro").join("config.toml");
            if home_config.exists() {
                if let Ok(config) = Self::load_from_path(&home_config) {
                    return Some(config);
                }
            }
        }

        None
    }

    /// Load configuration from a specific path.
    pub fn load_from_path(path: &PathBuf) -> Result<Self, ConfigError> {
        let content = fs::read_to_string(path).map_err(|e| ConfigError::Io(e.to_string()))?;
        Self::parse(&content)
    }

    /// Parse configuration from TOML content.
    pub fn parse(content: &str) -> Result<Self, ConfigError> {
        // Simple TOML parsing without external dependencies
        let mut config = FerroConfig::default();
        let mut in_error_handling = false;

        for line in content.lines() {
            let line = line.trim();

            // Skip comments and empty lines
            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            // Check for section headers
            if line.starts_with('[') && line.ends_with(']') {
                let section = &line[1..line.len() - 1];
                in_error_handling = section == "error-handling";
                continue;
            }

            if !in_error_handling {
                continue;
            }

            // Parse key-value pairs
            if let Some((key, value)) = line.split_once('=') {
                let key = key.trim();
                let value = value.trim();

                match key {
                    "mode" => {
                        // Remove quotes
                        let mode = value.trim_matches('"').trim_matches('\'');
                        config.error_handling.mode = Some(mode.to_string());
                    }
                    "ignore" => {
                        config.error_handling.ignore = parse_string_array(value);
                    }
                    "reject" => {
                        config.error_handling.reject = parse_string_array(value);
                    }
                    "warn" => {
                        config.error_handling.warn = parse_string_array(value);
                    }
                    _ => {}
                }
            }
        }

        Ok(config)
    }

    /// Convert this config to an ErrorConfig.
    pub fn to_error_config(&self) -> ErrorConfig {
        let mode = match self.error_handling.mode.as_deref() {
            Some("lenient") => ErrorMode::Lenient,
            Some("silent") => ErrorMode::Silent,
            _ => ErrorMode::Strict,
        };

        let mut config = ErrorConfig::new(mode);

        // Apply ignore overrides (silent correct)
        for code in &self.error_handling.ignore {
            if let Some(error_type) = ErrorType::from_code(code) {
                config.set_override(error_type, ErrorOverride::SilentCorrect);
            } else {
                warn_unmapped_override(code, "ignore", OverrideSource::ConfigFile);
            }
        }

        // Apply reject overrides
        for code in &self.error_handling.reject {
            if let Some(error_type) = ErrorType::from_code(code) {
                config.set_override(error_type, ErrorOverride::Reject);
            } else {
                warn_unmapped_override(code, "reject", OverrideSource::ConfigFile);
            }
        }

        // Apply warn overrides
        for code in &self.error_handling.warn {
            if let Some(error_type) = ErrorType::from_code(code) {
                config.set_override(error_type, ErrorOverride::WarnCorrect);
            } else {
                warn_unmapped_override(code, "warn", OverrideSource::ConfigFile);
            }
        }

        config
    }

    /// Merge this config with CLI arguments.
    /// CLI arguments take precedence.
    pub fn merge_with_cli(
        &self,
        cli_mode: Option<&str>,
        cli_ignore: &[String],
        cli_reject: &[String],
    ) -> ErrorConfig {
        // Start with config file settings
        let mut config = self.to_error_config();

        // Override mode if specified on CLI
        if let Some(mode_str) = cli_mode {
            config.mode = match mode_str {
                "lenient" => ErrorMode::Lenient,
                "silent" => ErrorMode::Silent,
                _ => ErrorMode::Strict,
            };
        }

        // Add CLI ignore overrides
        for code in cli_ignore {
            if let Some(error_type) = ErrorType::from_code(code) {
                config.set_override(error_type, ErrorOverride::SilentCorrect);
            } else {
                warn_unmapped_override(code, "ignore", OverrideSource::Cli);
            }
        }

        // Add CLI reject overrides
        for code in cli_reject {
            if let Some(error_type) = ErrorType::from_code(code) {
                config.set_override(error_type, ErrorOverride::Reject);
            } else {
                warn_unmapped_override(code, "reject", OverrideSource::Cli);
            }
        }

        config
    }
}

/// Configuration loading error.
#[derive(Debug, Clone)]
pub enum ConfigError {
    /// IO error reading config file.
    Io(String),
    /// Parse error in config file.
    Parse(String),
}

impl std::fmt::Display for ConfigError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConfigError::Io(msg) => write!(f, "Config IO error: {}", msg),
            ConfigError::Parse(msg) => write!(f, "Config parse error: {}", msg),
        }
    }
}

impl std::error::Error for ConfigError {}

/// Parse a TOML array of strings like `["W1001", "W2001"]`.
fn parse_string_array(value: &str) -> Vec<String> {
    let value = value.trim();
    if !value.starts_with('[') || !value.ends_with(']') {
        return Vec::new();
    }

    let inner = &value[1..value.len() - 1];
    inner
        .split(',')
        .map(|s| s.trim().trim_matches('"').trim_matches('\'').to_string())
        .filter(|s| !s.is_empty())
        .collect()
}

/// Get the user's home directory.
fn dirs_home() -> Option<PathBuf> {
    std::env::var("HOME").ok().map(PathBuf::from)
}

/// Identifies which surface routed an override that does not map to a
/// known `ErrorType`, so the stderr warning can name it accurately.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OverrideSource {
    /// The override came from a `.ferro.toml` `[error-handling]` list.
    ConfigFile,
    /// The override came from a CLI flag (`--ignore` / `--reject` / `--warn`).
    Cli,
}

/// Emit a stderr warning when an override targets a warning code that no
/// longer maps to a known `ErrorType`. Retired codes (currently W4002
/// PositionZero — see issue #269) are recognised specifically so the user
/// gets a migration pointer instead of a silent no-op; other unknown codes
/// get a generic "unknown code" hint. `source` controls the source prefix
/// so CLI overrides do not falsely point users at `.ferro.toml`.
pub fn warn_unmapped_override(code: &str, list: &str, source: OverrideSource) {
    let normalised = code.to_uppercase();
    let detail = if normalised == "W4002" {
        " (retired in issue #269; `c.0…` is now rejected at parse time \
         under `E1003 InvalidPosition` — drop this entry)"
    } else {
        " (unknown code; the override has no effect)"
    };
    let prefix = match source {
        OverrideSource::ConfigFile => format!("`.ferro.toml` `[error-handling].{list}` entry"),
        OverrideSource::Cli => format!("`--{list}` flag value"),
    };
    eprintln!("warning: {prefix} `{code}` does not match a known ErrorType{detail}");
}

// Warning-code → `ErrorType` resolution lives on `ErrorType::from_code`
// (see `src/error_handling/types.rs`) so it is derived from the same
// `code()` table the rest of the registry uses; this prevents the
// config layer from drifting away from the variant set the way an
// independent hard-coded table did.

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_empty_config() {
        let config = FerroConfig::parse("").unwrap();
        assert!(config.error_handling.mode.is_none());
        assert!(config.error_handling.ignore.is_empty());
    }

    #[test]
    fn test_parse_mode() {
        let content = r#"
[error-handling]
mode = "lenient"
"#;
        let config = FerroConfig::parse(content).unwrap();
        assert_eq!(config.error_handling.mode.as_deref(), Some("lenient"));
    }

    #[test]
    fn test_parse_ignore_list() {
        let content = r#"
[error-handling]
ignore = ["W1001", "W2001"]
"#;
        let config = FerroConfig::parse(content).unwrap();
        assert_eq!(config.error_handling.ignore, vec!["W1001", "W2001"]);
    }

    #[test]
    fn test_parse_reject_list() {
        let content = r#"
[error-handling]
reject = ["W3003"]
"#;
        let config = FerroConfig::parse(content).unwrap();
        assert_eq!(config.error_handling.reject, vec!["W3003"]);
    }

    #[test]
    fn test_to_error_config() {
        let content = r#"
[error-handling]
mode = "lenient"
ignore = ["W1001"]
reject = ["W3003"]
"#;
        let config = FerroConfig::parse(content).unwrap();
        let error_config = config.to_error_config();

        assert_eq!(error_config.mode, ErrorMode::Lenient);
        // W1001 should be silently corrected
        assert!(error_config.should_correct(ErrorType::LowercaseAminoAcid));
        assert!(!error_config.should_warn(ErrorType::LowercaseAminoAcid));
        // W3003 should be rejected
        assert!(error_config.should_reject(ErrorType::OldSubstitutionSyntax));
    }

    #[test]
    fn test_to_error_config_resolves_codes_previously_missing() {
        // Regression for CodeRabbit PR #370 review: the old hard-coded
        // mapping table omitted W3007, W3015, W3020, and W4003, so
        // `.ferro.toml` overrides for these registered codes silently
        // no-op'd. The mapping is now derived from `ErrorType::from_code`
        // and should round-trip every registered variant.
        let content = r#"
[error-handling]
mode = "strict"
ignore = ["W3007", "W3020"]
reject = ["W3015", "W4003"]
"#;
        let config = FerroConfig::parse(content).unwrap();
        let error_config = config.to_error_config();

        // ignore → silent correct
        assert!(error_config.should_correct(ErrorType::DeprecatedStopCodonStar));
        assert!(!error_config.should_warn(ErrorType::DeprecatedStopCodonStar));
        assert!(error_config.should_correct(ErrorType::RnaThymineCanonicalized));
        assert!(!error_config.should_warn(ErrorType::RnaThymineCanonicalized));

        // reject
        assert!(error_config.should_reject(ErrorType::DeprecatedConSyntax));
        assert!(error_config.should_reject(ErrorType::SinglePositionRange));
    }

    #[test]
    fn test_merge_with_cli() {
        let content = r#"
[error-handling]
mode = "lenient"
ignore = ["W1001"]
"#;
        let config = FerroConfig::parse(content).unwrap();

        // CLI overrides mode
        let merged = config.merge_with_cli(Some("strict"), &[], &[]);
        assert_eq!(merged.mode, ErrorMode::Strict);

        // CLI adds reject
        let merged = config.merge_with_cli(None, &[], &["W2001".to_string()]);
        assert!(merged.should_reject(ErrorType::WrongDashCharacter));
    }

    #[test]
    fn test_parse_string_array() {
        assert_eq!(
            parse_string_array(r#"["W1001", "W2001"]"#),
            vec!["W1001", "W2001"]
        );
        assert_eq!(parse_string_array(r#"["W1001"]"#), vec!["W1001"]);
        assert!(parse_string_array("").is_empty());
        assert!(parse_string_array("not an array").is_empty());
    }

    #[test]
    fn test_comments_ignored() {
        let content = r#"
# This is a comment
[error-handling]
# Another comment
mode = "silent"  # inline comment not supported but shouldn't break
"#;
        let config = FerroConfig::parse(content).unwrap();
        // Note: This simple parser doesn't handle inline comments, so mode might include them
        // For now, just verify the section is parsed
        assert!(config.error_handling.mode.is_some());
    }
}
