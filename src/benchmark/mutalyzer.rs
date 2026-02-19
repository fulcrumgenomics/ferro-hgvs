//! Mutalyzer API client for normalization benchmarks.

use crate::benchmark::types::ParseResult;
use crate::FerroError;
use reqwest::blocking::Client;
use serde::Deserialize;
use std::time::Duration;

/// HTTP client for the Mutalyzer API.
pub struct MutalyzerClient {
    client: Client,
    base_url: String,
}

/// Response from Mutalyzer normalize endpoint.
#[derive(Debug, Deserialize)]
pub struct NormalizeResponse {
    /// The normalized HGVS description
    pub normalized_description: Option<String>,

    /// Original input
    pub input_description: Option<String>,

    /// Errors if any
    #[serde(default)]
    pub errors: Vec<MutalyzerError>,

    /// Warnings if any
    #[serde(default)]
    pub warnings: Vec<MutalyzerWarning>,
}

#[derive(Debug, Deserialize)]
pub struct MutalyzerError {
    pub code: Option<String>,
    pub message: Option<String>,
}

#[derive(Debug, Deserialize)]
pub struct MutalyzerWarning {
    pub code: Option<String>,
    pub message: Option<String>,
}

impl MutalyzerClient {
    /// Create a new Mutalyzer client.
    pub fn new(base_url: &str) -> Result<Self, FerroError> {
        let client = Client::builder()
            .timeout(Duration::from_secs(30))
            .build()
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to create HTTP client: {}", e),
            })?;

        Ok(Self {
            client,
            base_url: base_url.trim_end_matches('/').to_string(),
        })
    }

    /// Check if the Mutalyzer API is reachable.
    pub fn health_check(&self) -> Result<bool, FerroError> {
        let url = format!("{}/", self.base_url);
        match self.client.get(&url).send() {
            Ok(resp) => Ok(resp.status().is_success()),
            Err(_) => Ok(false),
        }
    }

    /// Wait for the Mutalyzer API to become available.
    pub fn wait_for_ready(&self, timeout: Duration) -> Result<(), FerroError> {
        let start = std::time::Instant::now();
        let check_interval = Duration::from_secs(5);

        while start.elapsed() < timeout {
            if self.health_check()? {
                eprintln!("Mutalyzer API is ready at {}", self.base_url);
                return Ok(());
            }
            eprintln!("Waiting for Mutalyzer API...");
            std::thread::sleep(check_interval);
        }

        Err(FerroError::Io {
            msg: format!(
                "Mutalyzer API not ready after {:?} at {}",
                timeout, self.base_url
            ),
        })
    }

    /// Normalize a single HGVS expression.
    pub fn normalize(&self, hgvs: &str) -> Result<ParseResult, FerroError> {
        let encoded = urlencoding::encode(hgvs);
        let url = format!("{}/normalize/{}", self.base_url, encoded);

        let response = self.client.get(&url).send().map_err(|e| FerroError::Io {
            msg: format!("HTTP request failed: {}", e),
        })?;

        if !response.status().is_success() {
            return Ok(ParseResult {
                input: hgvs.to_string(),
                success: false,
                output: None,
                error: Some(format!("HTTP {}", response.status())),
                error_category: Some("http_error".to_string()),
                ref_mismatch: None,
                details: None,
            });
        }

        let body: NormalizeResponse = response.json().map_err(|e| FerroError::Io {
            msg: format!("Failed to parse response: {}", e),
        })?;

        if !body.errors.is_empty() {
            let error_msg = body
                .errors
                .iter()
                .filter_map(|e| e.message.as_ref())
                .cloned()
                .collect::<Vec<_>>()
                .join("; ");
            return Ok(ParseResult {
                input: hgvs.to_string(),
                success: false,
                output: None,
                error: Some(error_msg),
                error_category: Some("mutalyzer_error".to_string()),
                ref_mismatch: None,
                details: None,
            });
        }

        match body.normalized_description {
            Some(normalized) => Ok(ParseResult {
                input: hgvs.to_string(),
                success: true,
                output: Some(normalized),
                error: None,
                error_category: None,
                ref_mismatch: None,
                details: None,
            }),
            None => Ok(ParseResult {
                input: hgvs.to_string(),
                success: false,
                output: None,
                error: Some("No normalized description returned".to_string()),
                error_category: Some("empty_response".to_string()),
                ref_mismatch: None,
                details: None,
            }),
        }
    }

    /// Normalize a batch of HGVS expressions with rate limiting.
    pub fn normalize_batch(
        &self,
        patterns: &[String],
        rate_limit_ms: Option<u64>,
    ) -> Vec<ParseResult> {
        let delay = rate_limit_ms.map(Duration::from_millis);

        patterns
            .iter()
            .map(|pattern| {
                let result = self.normalize(pattern).unwrap_or_else(|e| ParseResult {
                    input: pattern.clone(),
                    success: false,
                    output: None,
                    error: Some(format!("{}", e)),
                    error_category: Some("client_error".to_string()),
                    ref_mismatch: None,
                    details: None,
                });

                if let Some(d) = delay {
                    std::thread::sleep(d);
                }

                result
            })
            .collect()
    }
}

/// Run Mutalyzer parsing via Python subprocess.
///
/// This calls the mutalyzer-hgvs-parser Python package which is separate
/// from the web API.
pub fn run_mutalyzer_parser_subprocess(
    input_file: &str,
    output_file: &str,
) -> Result<(), FerroError> {
    use std::process::Command;

    // Python script to run - outputs unified ToolParseOutput format
    let python_code = r#"
import sys
import json
import time

try:
    from mutalyzer_hgvs_parser import to_model
except ImportError:
    print("ERROR: mutalyzer-hgvs-parser not installed", file=sys.stderr)
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

results = []
successful = 0

start = time.perf_counter()

with open(input_file, 'r') as f:
    patterns = [line.strip() for line in f if line.strip()]

for pattern in patterns:
    try:
        model = to_model(pattern)
        # Output the original pattern on success - mutalyzer-hgvs-parser doesn't have
        # a built-in serializer, so we return the input if it parses successfully
        results.append({
            "input": pattern,
            "success": True,
            "output": pattern,
            "error": None
        })
        successful += 1
    except Exception as e:
        results.append({
            "input": pattern,
            "success": False,
            "output": None,
            "error": str(e)[:200]
        })

elapsed = time.perf_counter() - start

# Unified output format matching ToolParseOutput
output = {
    "tool": "mutalyzer",
    "total_patterns": len(patterns),
    "successful": successful,
    "failed": len(patterns) - successful,
    "elapsed_seconds": elapsed,
    "throughput": len(patterns) / elapsed if elapsed > 0 else 0,
    "results": results
}

with open(output_file, 'w') as f:
    json.dump(output, f, indent=2)
"#;

    let output = Command::new("python3")
        .args(["-c", python_code, input_file, output_file])
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to run Python: {}", e),
        })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(FerroError::Io {
            msg: format!("Mutalyzer parser failed: {}", stderr),
        });
    }

    Ok(())
}

/// Check if the mutalyzer-hgvs-parser Python package is available.
pub fn has_mutalyzer_parser() -> bool {
    use std::process::Command;

    Command::new("python3")
        .args(["-c", "from mutalyzer_hgvs_parser import to_model"])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Check if the mutalyzer package (for local normalization) is available.
pub fn has_mutalyzer_normalizer() -> bool {
    use std::process::Command;

    Command::new("python3")
        .args(["-c", "from mutalyzer.description import Description"])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Run local mutalyzer normalization via Python subprocess.
///
/// This calls the mutalyzer package directly for local normalization,
/// avoiding HTTP API calls for fair performance comparison.
pub fn run_mutalyzer_normalizer_subprocess(
    input_file: &str,
    output_file: &str,
    settings_file: Option<&str>,
) -> Result<(), FerroError> {
    use std::process::Command;

    // Python script for local normalization
    let python_code = r#"
import sys
import json
import time
import os

# Network call tracking - must be set up before any imports that use urllib
network_calls = []
_original_urlopen = None

def _tracking_urlopen(url, *args, **kwargs):
    """Track network calls made by urllib."""
    url_str = url if isinstance(url, str) else url.full_url if hasattr(url, 'full_url') else str(url)
    start = time.perf_counter()
    try:
        result = _original_urlopen(url, *args, **kwargs)
        elapsed = time.perf_counter() - start
        network_calls.append({"url": url_str[:200], "elapsed_seconds": elapsed, "success": True})
        return result
    except Exception as e:
        elapsed = time.perf_counter() - start
        network_calls.append({"url": url_str[:200], "elapsed_seconds": elapsed, "success": False, "error": str(e)[:100]})
        raise

# Install network tracking
import urllib.request
_original_urlopen = urllib.request.urlopen
urllib.request.urlopen = _tracking_urlopen

# Set up mutalyzer settings if provided
# We must configure the cache BEFORE importing mutalyzer, because
# mutalyzer_retriever.configuration loads settings at import time.
cache_dir = None
cache_add = False
if len(sys.argv) > 3 and sys.argv[3]:
    settings_file = sys.argv[3]
    try:
        with open(settings_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('MUTALYZER_CACHE_DIR'):
                    parts = line.split('=', 1)
                    if len(parts) == 2:
                        cache_dir = parts[1].strip()
                elif line.startswith('MUTALYZER_FILE_CACHE_ADD'):
                    parts = line.split('=', 1)
                    if len(parts) == 2:
                        cache_add = parts[1].strip().lower() in ('true', '1', 'yes')
    except Exception as e:
        print(f"WARNING: Failed to parse settings file: {e}", file=sys.stderr)

# Import mutalyzer_retriever first and configure the cache directory
# directly in its settings dict (env vars don't work after import)
from mutalyzer_retriever import configuration
if cache_dir:
    configuration.settings['MUTALYZER_CACHE_DIR'] = cache_dir
if cache_add:
    configuration.settings['MUTALYZER_FILE_CACHE_ADD'] = True

try:
    from mutalyzer.description import Description
except ImportError:
    print("ERROR: mutalyzer not installed", file=sys.stderr)
    print("Install with: pip install mutalyzer", file=sys.stderr)
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

results = []
successful = 0
failed = 0

with open(input_file, 'r') as f:
    patterns = [line.strip() for line in f if line.strip()]

start = time.perf_counter()

for pattern in patterns:
    try:
        d = Description(description=pattern)
        d.normalize()
        if not d.errors:
            output = d.output()
            results.append({
                "input": pattern,
                "success": True,
                "output": output.get("normalized_description", pattern)
            })
            successful += 1
        else:
            error_msg = str(d.errors[0]) if d.errors else "Unknown error"
            results.append({
                "input": pattern,
                "success": False,
                "error": error_msg[:200]
            })
            failed += 1
    except Exception as e:
        results.append({
            "input": pattern,
            "success": False,
            "error": str(e)[:200]
        })
        failed += 1

elapsed = time.perf_counter() - start

output_data = {
    "tool": "mutalyzer",
    "total_patterns": len(patterns),
    "successful": successful,
    "failed": failed,
    "elapsed_seconds": elapsed,
    "patterns_per_second": len(patterns) / elapsed if elapsed > 0 else 0,
    "network_calls": len(network_calls),
    "network_stats": network_calls[:50] if network_calls else [],  # First 50 for debugging
    "results": results
}

with open(output_file, 'w') as f:
    json.dump(output_data, f, indent=2)

# Report network statistics
if network_calls:
    total_network_time = sum(c.get('elapsed_seconds', 0) for c in network_calls)
    print(f"\n=== Network Statistics ===", file=sys.stderr)
    print(f"Total network calls: {len(network_calls)}", file=sys.stderr)
    print(f"Network calls per pattern: {len(network_calls) / len(patterns):.2f}", file=sys.stderr)
    print(f"Total network time: {total_network_time:.2f}s", file=sys.stderr)
    print(f"Average time per call: {total_network_time / len(network_calls):.2f}s", file=sys.stderr)
    # Show unique URLs called
    unique_urls = set()
    for c in network_calls:
        url = c.get('url', '')
        # Extract domain/path pattern
        if 'ncbi.nlm.nih.gov' in url:
            unique_urls.add('NCBI')
        elif 'ebi.ac.uk' in url:
            unique_urls.add('EBI')
        else:
            unique_urls.add(url[:60])
    print(f"Endpoints hit: {', '.join(sorted(unique_urls))}", file=sys.stderr)
else:
    print(f"\n=== Network Statistics ===", file=sys.stderr)
    print(f"Total network calls: 0 (all cached)", file=sys.stderr)
"#;

    let mut cmd = Command::new("python3");
    cmd.args(["-c", python_code, input_file, output_file]);

    if let Some(settings) = settings_file {
        cmd.arg(settings);
    } else {
        cmd.arg("");
    }

    let output = cmd.output().map_err(|e| FerroError::Io {
        msg: format!("Failed to run Python: {}", e),
    })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(FerroError::Io {
            msg: format!("Mutalyzer normalizer failed: {}", stderr),
        });
    }

    Ok(())
}

/// Normalize a single HGVS variant using local mutalyzer subprocess.
///
/// This is a convenience function for normalizing a single variant,
/// following the pattern used by biocommons::normalize_single().
pub fn normalize_single(
    hgvs: &str,
    settings_file: Option<&str>,
    allow_network: bool,
) -> Result<ParseResult, FerroError> {
    use std::io::Write;
    use tempfile::NamedTempFile;

    // Create temp files for input and output
    let mut input_file = NamedTempFile::new().map_err(|e| FerroError::Io {
        msg: format!("Failed to create temp file: {}", e),
    })?;
    let output_file = NamedTempFile::new().map_err(|e| FerroError::Io {
        msg: format!("Failed to create temp file: {}", e),
    })?;

    // Write input pattern
    writeln!(input_file, "{}", hgvs).map_err(|e| FerroError::Io {
        msg: format!("Failed to write to temp file: {}", e),
    })?;

    // Run normalization
    if allow_network {
        run_mutalyzer_normalizer_subprocess(
            input_file.path().to_str().unwrap(),
            output_file.path().to_str().unwrap(),
            settings_file,
        )?;
    } else {
        run_mutalyzer_normalizer_subprocess_offline(
            input_file.path().to_str().unwrap(),
            output_file.path().to_str().unwrap(),
            settings_file,
        )?;
    }

    // Read output
    let output_content =
        std::fs::read_to_string(output_file.path()).map_err(|e| FerroError::Io {
            msg: format!("Failed to read output file: {}", e),
        })?;

    // Parse JSON output
    let output_data: serde_json::Value =
        serde_json::from_str(&output_content).map_err(|e| FerroError::Io {
            msg: format!("Failed to parse JSON output: {}", e),
        })?;

    // Extract result
    let results = output_data["results"]
        .as_array()
        .ok_or_else(|| FerroError::Io {
            msg: "No results in output".to_string(),
        })?;

    if let Some(result) = results.first() {
        let success = result["success"].as_bool().unwrap_or(false);
        if success {
            Ok(ParseResult {
                input: hgvs.to_string(),
                success: true,
                output: result["output"].as_str().map(|s| s.to_string()),
                error: None,
                error_category: None,
                ref_mismatch: None,
                details: None,
            })
        } else {
            Ok(ParseResult {
                input: hgvs.to_string(),
                success: false,
                output: None,
                error: result["error"].as_str().map(|s| s.to_string()),
                error_category: Some("mutalyzer_error".to_string()),
                ref_mismatch: None,
                details: None,
            })
        }
    } else {
        Err(FerroError::Io {
            msg: "No results returned from mutalyzer".to_string(),
        })
    }
}

/// Run local mutalyzer normalization via Python subprocess (offline mode).
///
/// This version disables network access to force cache-only operation.
pub fn run_mutalyzer_normalizer_subprocess_offline(
    input_file: &str,
    output_file: &str,
    settings_file: Option<&str>,
) -> Result<(), FerroError> {
    use std::process::Command;

    // Python script for local normalization with network disabled
    let python_code = r#"
import sys
import json
import time
import os

# Block network access
def _blocked_urlopen(*args, **kwargs):
    raise ConnectionError("Network access disabled in offline mode")

import urllib.request
urllib.request.urlopen = _blocked_urlopen

# Set up mutalyzer settings if provided
cache_dir = None
cache_add = False
if len(sys.argv) > 3 and sys.argv[3]:
    settings_file = sys.argv[3]
    try:
        with open(settings_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('MUTALYZER_CACHE_DIR'):
                    parts = line.split('=', 1)
                    if len(parts) == 2:
                        cache_dir = parts[1].strip()
                elif line.startswith('MUTALYZER_FILE_CACHE_ADD'):
                    parts = line.split('=', 1)
                    if len(parts) == 2:
                        cache_add = parts[1].strip().lower() in ('true', '1', 'yes')
    except Exception as e:
        print(f"WARNING: Failed to parse settings file: {e}", file=sys.stderr)

# Import mutalyzer_retriever first and configure the cache directory
from mutalyzer_retriever import configuration
if cache_dir:
    configuration.settings['MUTALYZER_CACHE_DIR'] = cache_dir
if cache_add:
    configuration.settings['MUTALYZER_FILE_CACHE_ADD'] = True

try:
    from mutalyzer.description import Description
except ImportError:
    print("ERROR: mutalyzer not installed", file=sys.stderr)
    print("Install with: pip install mutalyzer", file=sys.stderr)
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

results = []
successful = 0
failed = 0

with open(input_file, 'r') as f:
    patterns = [line.strip() for line in f if line.strip()]

start = time.perf_counter()

for pattern in patterns:
    try:
        d = Description(description=pattern)
        d.normalize()
        if not d.errors:
            output = d.output()
            results.append({
                "input": pattern,
                "success": True,
                "output": output.get("normalized_description", pattern)
            })
            successful += 1
        else:
            error_msg = str(d.errors[0]) if d.errors else "Unknown error"
            results.append({
                "input": pattern,
                "success": False,
                "error": error_msg[:200]
            })
            failed += 1
    except ConnectionError as e:
        results.append({
            "input": pattern,
            "success": False,
            "error": f"Network access required but disabled: {str(e)[:150]}"
        })
        failed += 1
    except Exception as e:
        results.append({
            "input": pattern,
            "success": False,
            "error": str(e)[:200]
        })
        failed += 1

elapsed = time.perf_counter() - start

output_data = {
    "tool": "mutalyzer",
    "total_patterns": len(patterns),
    "successful": successful,
    "failed": failed,
    "elapsed_seconds": elapsed,
    "patterns_per_second": len(patterns) / elapsed if elapsed > 0 else 0,
    "network_calls": 0,
    "results": results
}

with open(output_file, 'w') as f:
    json.dump(output_data, f, indent=2)
"#;

    let mut cmd = Command::new("python3");
    cmd.args(["-c", python_code, input_file, output_file]);

    if let Some(settings) = settings_file {
        cmd.arg(settings);
    } else {
        cmd.arg("");
    }

    let output = cmd.output().map_err(|e| FerroError::Io {
        msg: format!("Failed to run Python: {}", e),
    })?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(FerroError::Io {
            msg: format!("Mutalyzer normalizer (offline) failed: {}", stderr),
        });
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_mutalyzer_available() {
        // This test just verifies the function doesn't panic
        let _ = has_mutalyzer_parser();
    }
}
