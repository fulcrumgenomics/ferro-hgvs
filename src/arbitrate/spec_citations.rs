//! Vendored HGVS spec citations, embedded so external users (who lack the
//! `assets/hgvs-nomenclature` submodule) still get quoted rules offline.
//! Generated/verified by `examples/generate_spec_citations.rs` (Task 4).

use serde::Deserialize;
use std::collections::HashMap;
use std::sync::OnceLock;

/// A quoted, versioned spec passage.
#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
pub struct SpecCitation {
    /// The HGVS spec release the excerpt was verified against. Filled from the
    /// map's top-level `spec_version` when a record does not carry its own.
    #[serde(default)]
    pub spec_version: String,
    /// The spec source file, relative to the spec repository root.
    pub file: String,
    /// The section heading the excerpt sits under.
    pub heading: String,
    /// The quoted passage, checked against `file` by
    /// `examples/generate_spec_citations.rs`.
    pub excerpt: String,
}

/// The normalization/projection operation a disagreement is about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operation {
    /// The 3' rule: which of several equivalent positions is reported.
    ThreePrimeShift,
    /// Whether an insertion that copies its 5' flank must be a `dup`.
    DupVsIns,
    /// Whether nearby changes are merged into one `delins`. No citation is
    /// mapped yet, so [`lookup`] returns an empty list.
    DelinsCoalesce,
    /// Re-expressing a variant on another axis or reference. No citation is
    /// mapped yet, so [`lookup`] returns an empty list.
    Projection,
}

impl Operation {
    fn key(self) -> &'static str {
        match self {
            Operation::ThreePrimeShift => "three_prime_shift",
            Operation::DupVsIns => "dup_vs_ins",
            Operation::DelinsCoalesce => "delins_coalesce",
            Operation::Projection => "projection",
        }
    }
}

#[derive(Deserialize)]
struct RawMap {
    spec_version: String,
    entries: HashMap<String, Vec<SpecCitation>>,
}

const EMBEDDED: &str = include_str!("spec_citations.json");

fn map() -> &'static RawMap {
    static M: OnceLock<RawMap> = OnceLock::new();
    M.get_or_init(|| serde_json::from_str(EMBEDDED).expect("embedded spec_citations.json is valid"))
}

/// Look up the governing citations for an operation. Empty if unmapped.
pub fn lookup(op: Operation) -> Vec<SpecCitation> {
    let m = map();
    m.entries
        .get(op.key())
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .map(|mut c| {
            if c.spec_version.is_empty() {
                c.spec_version = m.spec_version.clone();
            }
            c
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_launch_operation_resolves_a_nonempty_excerpt() {
        for op in [Operation::ThreePrimeShift, Operation::DupVsIns] {
            let cites = lookup(op);
            assert!(!cites.is_empty(), "no citation for {op:?}");
            assert!(cites.iter().all(|c| !c.excerpt.trim().is_empty()));
            assert!(cites.iter().all(|c| !c.spec_version.is_empty()));
        }
    }
}
