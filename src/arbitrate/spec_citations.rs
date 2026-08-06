//! Vendored HGVS spec citations, embedded so external users (who lack the
//! `assets/hgvs-nomenclature` submodule) still get quoted rules offline.
//! Generated/verified by `examples/generate_spec_citations.rs` (Task 4).

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::OnceLock;

/// A quoted, versioned spec passage.
#[derive(Debug, Clone, Deserialize, Serialize, PartialEq, Eq)]
pub struct SpecCitation {
    #[serde(default)]
    pub spec_version: String,
    pub file: String,
    pub heading: String,
    pub excerpt: String,
}

/// The normalization/projection operation a disagreement is about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operation {
    ThreePrimeShift,
    DupVsIns,
    DelinsCoalesce,
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
