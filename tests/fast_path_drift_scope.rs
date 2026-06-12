//! Scope-enumeration (not an assertion gate): bucket every way the fast path
//! diverges from the standard parser, so we can size the parity gap before
//! choosing a fix. Two diagnostics:
//!
//!   * `enumerate_fast_path_drift` — a widened accession x coord-system x edit
//!     grid (mito / ncRNA / genomic-region / protein refs, UTR + intronic +
//!     repeat edits).
//!   * `replay_corpus_drift` — replays a real HGVS corpus (newline-delimited
//!     text at `FERRO_DIFF_CORPUS_TXT`, e.g. the 500k ClinVar set) and
//!     histograms the divergences by class and by the standard parser's
//!     rejection reason.
//!
//! Run with:
//!   cargo test --release --test fast_path_drift_scope --features dev -- \
//!     --nocapture --ignored

use ferro_hgvs::{parse_hgvs, parse_hgvs_fast};
use std::collections::BTreeMap;

/// Divergence class for a single input.
enum Class {
    Agree,
    FastAcceptsStdRejects(String), // carries the standard parser's reason
    StdAcceptsFastRejects,
    BothOkDiffer,
}

fn classify(input: &str) -> Class {
    match (parse_hgvs(input), parse_hgvs_fast(input)) {
        (Ok(a), Ok(b)) if a == b => Class::Agree,
        (Err(_), Err(_)) => Class::Agree,
        (Err(e), Ok(_)) => {
            // Normalize the reason to a stable bucket key: take the text up to
            // the first reference-specific token so NR_/XR_/accession numbers
            // collapse together.
            let msg = e.to_string();
            let key = msg
                .split(" (")
                .next()
                .unwrap_or(&msg)
                .replace("position 0: ", "")
                .trim()
                .to_string();
            Class::FastAcceptsStdRejects(key)
        }
        (Ok(_), Err(_)) => Class::StdAcceptsFastRejects,
        (Ok(_), Ok(_)) => Class::BothOkDiffer,
    }
}

const ACCESSIONS: &[&str] = &[
    "NC_000001.11", // genomic
    "NC_012920.1",  // mitochondrial genome
    "NG_012232.1",  // genomic region (RefSeqGene)
    "NM_000088.3",  // mRNA
    "NR_003051.3",  // non-coding RNA
    "XR_001234.1",  // predicted non-coding RNA
    "NP_000079.2",  // protein
    "ENST00000357033.8",
    "ENSG00000139618.15",
    "ENSP00000369497.3", // Ensembl protein
    "LRG_1",
    "LRG_1t1",
    "GRCh38(chr1)",
];
const SYSTEMS: &[&str] = &["g", "c", "n", "r", "m", "p"];
const EDITS: &[&str] = &[
    "100A>G",       // substitution
    "100del",       // deletion
    "100_101insAC", // insertion
    "100dup",       // duplication
    "100_101inv",   // inversion
    "*100A>G",      // 3' UTR substitution
    "-50A>G",       // 5' UTR substitution
    "100+5A>G",     // intronic substitution
    "100-5A>G",     // intronic substitution (acceptor side)
    "100GT[3]",     // repeat
    "Arg100Gly",    // protein substitution (3-letter)
];

#[test]
#[ignore = "diagnostic enumeration, run explicitly with --ignored --nocapture"]
fn enumerate_fast_path_drift() {
    let mut total = 0usize;
    let mut by_reason: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut std_accepts_fast_rejects = Vec::new();
    let mut both_ok_differ = Vec::new();

    for acc in ACCESSIONS {
        for sys in SYSTEMS {
            for edit in EDITS {
                let input = format!("{acc}:{sys}.{edit}");
                total += 1;
                match classify(&input) {
                    Class::Agree => {}
                    Class::FastAcceptsStdRejects(reason) => {
                        by_reason.entry(reason).or_default().push(input);
                    }
                    Class::StdAcceptsFastRejects => std_accepts_fast_rejects.push(input),
                    Class::BothOkDiffer => both_ok_differ.push(input),
                }
            }
        }
    }

    let diverged: usize = by_reason.values().map(Vec::len).sum::<usize>()
        + std_accepts_fast_rejects.len()
        + both_ok_differ.len();

    println!("\n=== GRID drift: {diverged}/{total} inputs diverge ===\n");
    println!("[A] fast ACCEPTS, standard REJECTS — grouped by standard's reason:");
    for (reason, inputs) in &by_reason {
        println!("  ({}) {reason}", inputs.len());
        for input in inputs {
            println!("        {input}");
        }
    }
    println!(
        "\n[B] standard ACCEPTS, fast REJECTS ({}):",
        std_accepts_fast_rejects.len()
    );
    for input in &std_accepts_fast_rejects {
        println!("        {input}");
    }
    println!("\n[C] both accept, ASTs DIFFER ({}):", both_ok_differ.len());
    for input in &both_ok_differ {
        println!("        {input}");
    }
    println!();
}

#[test]
#[ignore = "corpus replay, run explicitly with --ignored --nocapture; needs FERRO_DIFF_CORPUS_TXT"]
fn replay_corpus_drift() {
    let Ok(path) = std::env::var("FERRO_DIFF_CORPUS_TXT") else {
        println!(
            "\nSKIP replay_corpus_drift: set FERRO_DIFF_CORPUS_TXT to a \
             newline-delimited HGVS file (e.g. the 500k ClinVar inputs)\n"
        );
        return;
    };
    let text = std::fs::read_to_string(&path).expect("read corpus txt");

    let mut total = 0usize;
    let mut by_reason: BTreeMap<String, (usize, Vec<String>)> = BTreeMap::new();
    let mut std_accepts_fast_rejects: (usize, Vec<String>) = (0, Vec::new());
    let mut both_ok_differ: (usize, Vec<String>) = (0, Vec::new());

    // Keep at most this many example inputs per bucket (count all of them).
    const MAX_EXAMPLES: usize = 8;

    for line in text.lines() {
        let input = line.trim();
        if input.is_empty() {
            continue;
        }
        total += 1;
        match classify(input) {
            Class::Agree => {}
            Class::FastAcceptsStdRejects(reason) => {
                let entry = by_reason.entry(reason).or_insert((0, Vec::new()));
                entry.0 += 1;
                if entry.1.len() < MAX_EXAMPLES {
                    entry.1.push(input.to_string());
                }
            }
            Class::StdAcceptsFastRejects => {
                std_accepts_fast_rejects.0 += 1;
                if std_accepts_fast_rejects.1.len() < MAX_EXAMPLES {
                    std_accepts_fast_rejects.1.push(input.to_string());
                }
            }
            Class::BothOkDiffer => {
                both_ok_differ.0 += 1;
                if both_ok_differ.1.len() < MAX_EXAMPLES {
                    both_ok_differ.1.push(input.to_string());
                }
            }
        }
    }

    let diverged: usize = by_reason.values().map(|(n, _)| *n).sum::<usize>()
        + std_accepts_fast_rejects.0
        + both_ok_differ.0;

    println!("\n=== CORPUS replay: {diverged}/{total} inputs diverge ===\n");
    println!("[A] fast ACCEPTS, standard REJECTS — grouped by standard's reason:");
    for (reason, (count, examples)) in &by_reason {
        println!("  ({count}) {reason}");
        for ex in examples {
            println!("        e.g. {ex}");
        }
    }
    println!(
        "\n[B] standard ACCEPTS, fast REJECTS ({}):",
        std_accepts_fast_rejects.0
    );
    for ex in &std_accepts_fast_rejects.1 {
        println!("        e.g. {ex}");
    }
    println!("\n[C] both accept, ASTs DIFFER ({}):", both_ok_differ.0);
    for ex in &both_ok_differ.1 {
        println!("        e.g. {ex}");
    }
    println!();
}
