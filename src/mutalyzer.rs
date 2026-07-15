//! Minimal blocking Mutalyzer client for `ferro arbitrate`.
//!
//! Deliberately un-gated and blocking (uses the already-enabled
//! `reqwest::blocking`), so it ships in the default `ferro` binary without
//! dragging in the async/web-service stack. Parsing is split into a pure
//! `from_response_json` for testability (no network in CI).

use crate::error::FerroError;
use std::time::Duration;

/// Outcome category of a Mutalyzer normalize call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MutalyzerStatus {
    /// A normalized description was returned on the input's basis.
    Ok,
    /// Intronic-on-transcript: Mutalyzer redirected to genomic options.
    GenomicRedirect,
    /// Mutalyzer rejected the input (syntax/reference error).
    ParseError,
    /// The endpoint could not be reached (network/timeout).
    Unavailable,
}

/// A genomic-framed alternative Mutalyzer offered (EINTRONIC redirect).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicOption {
    pub assembly_id: String,
    pub description: String,
}

/// Parsed result of a Mutalyzer normalize call.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MutalyzerResult {
    pub status: MutalyzerStatus,
    pub normalized: Option<String>,
    pub genomic_options: Vec<GenomicOption>,
    pub errors: Vec<String>,
}

impl MutalyzerResult {
    fn unavailable(reason: impl Into<String>) -> Self {
        MutalyzerResult {
            status: MutalyzerStatus::Unavailable,
            normalized: None,
            genomic_options: Vec::new(),
            errors: vec![reason.into()],
        }
    }

    /// Pure parser over a Mutalyzer JSON body. See §5.1.4a of the design spec.
    pub fn from_response_json(_input: &str, body: &serde_json::Value) -> Self {
        // 1. Clean case: top-level normalized_description.
        if let Some(norm) = body.get("normalized_description").and_then(|v| v.as_str()) {
            return MutalyzerResult {
                status: MutalyzerStatus::Ok,
                normalized: Some(norm.to_string()),
                genomic_options: Vec::new(),
                errors: Vec::new(),
            };
        }
        // 2. Genomic-redirect case: errors[].options[].
        let mut genomic_options = Vec::new();
        let mut errors = Vec::new();
        if let Some(errs) = body.get("errors").and_then(|v| v.as_array()) {
            for e in errs {
                if let (Some(code), details) = (
                    e.get("code").and_then(|c| c.as_str()),
                    e.get("details").and_then(|d| d.as_str()),
                ) {
                    errors.push(match details {
                        Some(d) => format!("{code}: {d}"),
                        None => code.to_string(),
                    });
                }
                if let Some(opts) = e.get("options").and_then(|o| o.as_array()) {
                    for o in opts {
                        if let (Some(a), Some(d)) = (
                            o.get("assembly_id").and_then(|x| x.as_str()),
                            o.get("description").and_then(|x| x.as_str()),
                        ) {
                            genomic_options.push(GenomicOption {
                                assembly_id: a.to_string(),
                                description: d.to_string(),
                            });
                        }
                    }
                }
            }
        }
        // 3. custom.errors as an additional error source.
        if let Some(cerrs) = body
            .get("custom")
            .and_then(|c| c.get("errors"))
            .and_then(|v| v.as_array())
        {
            for e in cerrs {
                if let Some(code) = e.get("code").and_then(|c| c.as_str()) {
                    let msg = match e.get("details").and_then(|d| d.as_str()) {
                        Some(d) => format!("{code}: {d}"),
                        None => code.to_string(),
                    };
                    if !errors.contains(&msg) {
                        errors.push(msg);
                    }
                }
            }
        }
        if !genomic_options.is_empty() {
            MutalyzerResult {
                status: MutalyzerStatus::GenomicRedirect,
                normalized: None,
                genomic_options,
                errors,
            }
        } else {
            MutalyzerResult {
                status: MutalyzerStatus::ParseError,
                normalized: None,
                genomic_options: Vec::new(),
                errors: if errors.is_empty() {
                    vec!["Mutalyzer returned no normalized description".to_string()]
                } else {
                    errors
                },
            }
        }
    }
}

/// Blocking HTTP client for the Mutalyzer normalize endpoint.
pub struct MutalyzerClient {
    base_url: String,
    client: reqwest::blocking::Client,
}

impl MutalyzerClient {
    /// Create a client. `base_url` defaults to the public API when empty.
    pub fn new(base_url: impl Into<String>) -> Result<Self, FerroError> {
        let mut base_url = base_url.into();
        if base_url.is_empty() {
            base_url = "https://mutalyzer.nl".to_string();
        }
        let client = reqwest::blocking::Client::builder()
            .timeout(Duration::from_secs(30))
            .build()
            .map_err(|e| FerroError::Io {
                msg: format!("HTTP client build failed: {e}"),
            })?;
        Ok(MutalyzerClient {
            base_url: base_url.trim_end_matches('/').to_string(),
            client,
        })
    }

    /// Normalize one HGVS expression. Network failures return
    /// `MutalyzerStatus::Unavailable` rather than erroring, so the caller can
    /// fall back to a pasted result.
    pub fn normalize(&self, hgvs: &str) -> Result<MutalyzerResult, FerroError> {
        let encoded = urlencode(hgvs);
        let url = if self.base_url.contains("mutalyzer.nl") {
            format!("{}/api/normalize/{}", self.base_url, encoded)
        } else {
            format!("{}/normalize/{}", self.base_url, encoded)
        };
        let resp = match self.client.get(&url).send() {
            Ok(r) => r,
            Err(e) => return Ok(MutalyzerResult::unavailable(format!("request failed: {e}"))),
        };
        let body: serde_json::Value = match resp.json() {
            Ok(v) => v,
            Err(e) => {
                return Ok(MutalyzerResult::unavailable(format!(
                    "bad response body: {e}"
                )))
            }
        };
        Ok(MutalyzerResult::from_response_json(hgvs, &body))
    }
}

/// Percent-encode a path segment (avoids adding the `urlencoding` optional dep,
/// which is gated behind the benchmark/web-service features).
///
/// Intentionally hand-rolled rather than reusing the already-present `url`
/// crate: `url`'s path-segment API (`Url::path_segments_mut`) leaves `:` and
/// `+` unescaped, and its `form_urlencoded` query-pair encoder leaves `*`
/// unescaped (verified against `NM_000059.3:c.68-7_316+7del` and HGVS `*`
/// UTR/stop-codon notation, e.g. `c.*4A>G`) — neither is byte-identical to
/// this encoding, and Mutalyzer's URL path segment needs `:`/`+`/`*` escaped.
fn urlencode(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for b in s.bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_' | b'.' | b'~' => {
                out.push(b as char)
            }
            _ => out.push_str(&format!("%{b:02X}")),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_clean_normalized_description() {
        let body: serde_json::Value = serde_json::json!({
            "normalized_description": "NM_004006.2:c.20del",
            "input_description": "NM_004006.2:c.20del"
        });
        let r = MutalyzerResult::from_response_json("NM_004006.2:c.20del", &body);
        assert_eq!(r.status, MutalyzerStatus::Ok);
        assert_eq!(r.normalized.as_deref(), Some("NM_004006.2:c.20del"));
        assert!(r.genomic_options.is_empty());
    }

    #[test]
    fn parses_eintronic_genomic_redirect() {
        // Shape verified live 2026-07-09 for NM_000059.3:c.68-7_316+7del.
        let body: serde_json::Value = serde_json::json!({
            "errors": [{
                "code": "EINTRONIC",
                "details": "Intronic positions were used with a non intronic reference.",
                "options": [
                    {"assembly_id": "GRCh38", "description": "NC_000013.11(NM_000059.3):c.68-7_316+7del"},
                    {"assembly_id": "GRCh37", "description": "NC_000013.10(NM_000059.3):c.68-7_316+7del"}
                ]
            }]
        });
        let r = MutalyzerResult::from_response_json("NM_000059.3:c.68-7_316+7del", &body);
        assert_eq!(r.status, MutalyzerStatus::GenomicRedirect);
        assert_eq!(r.genomic_options.len(), 2);
        assert_eq!(r.genomic_options[0].assembly_id, "GRCh38");
        assert_eq!(
            r.genomic_options[0].description,
            "NC_000013.11(NM_000059.3):c.68-7_316+7del"
        );
    }

    #[test]
    fn parses_plain_error() {
        let body: serde_json::Value = serde_json::json!({
            "errors": [{"code": "ESYNTAX", "details": "bad syntax"}],
            "custom": {"errors": [{"code": "ESYNTAX", "details": "bad syntax"}]}
        });
        let r = MutalyzerResult::from_response_json("garbage", &body);
        assert_eq!(r.status, MutalyzerStatus::ParseError);
        // Exact-match: `custom.errors` repeats the same `ESYNTAX: bad
        // syntax` message already present in top-level `errors[]`, so the
        // `custom.errors` dedup (`!errors.contains(&msg)`) must collapse
        // this to a single entry. A loose `.any(...)` assertion would not
        // catch a regression that let the dedup double this entry.
        assert_eq!(r.errors, vec!["ESYNTAX: bad syntax".to_string()]);
    }
}
