//! Report the reference accessions a conformance corpus needs that a prepared
//! reference does **not** carry — the operational form of the
//! `reference_unavailable` triage (#706) and a discovery step toward the
//! version-exact reference snapshot (#719).
//!
//! It builds the corpus accession inventory (`conformance::accession_inventory`)
//! and diffs it against the set of accessions present in the reference's FASTA
//! `.fai` indexes (`Inventory::missing_against`). Each missing entry is a corpus
//! row that no-ops for lack of a reference — not a ferro defect — and is exactly
//! what a version-complete snapshot must provision.
//!
//! Run:
//!   cargo run --features dev --example report_conformance_reference_gaps -- \
//!     --reference /path/to/prepared/reference
//!   (defaults to the mutalyzer corpus; pass --cases for another corpus)

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use ferro_hgvs::conformance::accession_inventory::{AccessionClass, Inventory};
use ferro_hgvs::conformance::mutalyzer::Fixture;

#[derive(Parser, Debug)]
#[command(about = "Report corpus-referenced accessions absent from a prepared reference")]
struct Cli {
    /// Corpus `cases.json` to inventory.
    #[arg(long, default_value = "tests/fixtures/mutalyzer-normalize/cases.json")]
    cases: PathBuf,
    /// Prepared-reference directory to scan for FASTA `.fai` indexes
    /// (searched recursively).
    #[arg(long)]
    reference: PathBuf,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: &Cli) -> Result<(), String> {
    let content = std::fs::read_to_string(&cli.cases)
        .map_err(|e| format!("read {}: {e}", cli.cases.display()))?;
    let fixture: Fixture = serde_json::from_str(&content)
        .map_err(|e| format!("parse {}: {e}", cli.cases.display()))?;
    let inventory = Inventory::from_fixture(&fixture);

    let available = collect_fai_accessions(&cli.reference)?;
    if available.is_empty() {
        eprintln!(
            "warning: no `.fai` indexes found under {} — every accession will report missing",
            cli.reference.display()
        );
    }
    let missing = inventory.missing_against(&available);

    println!("Corpus:    {}", cli.cases.display());
    println!("Reference: {}", cli.reference.display());
    println!(
        "Referenced accessions: {} ({})",
        inventory.len(),
        format_class_counts(&inventory.counts_by_class())
    );
    println!("Available in reference: {}", available.len());
    println!("MISSING (corpus needs, reference lacks): {}", missing.len());

    if missing.is_empty() {
        println!("\nNo gaps — the reference is version-complete for this corpus.");
        return Ok(());
    }

    println!(
        "\nNote: matching is exact on the versioned key. Entries marked \
         [unversioned reference] (bare `LRG_<n>`, or an accession the corpus \
         cites without a version) won't exact-match a versioned/suffixed index \
         (e.g. `LRG_199` vs the indexed `LRG_199g`/`t1`) and may be present in \
         another form — verify these manually. Versioned RefSeq/Ensembl entries \
         are genuine gaps."
    );

    // Group missing entries by class for a readable report.
    let mut by_class: std::collections::BTreeMap<AccessionClass, Vec<&_>> =
        std::collections::BTreeMap::new();
    for entry in &missing {
        by_class.entry(entry.acc_ref.class).or_default().push(entry);
    }
    println!();
    for (class, entries) in &by_class {
        println!("  {} ({}):", class.as_str(), entries.len());
        for entry in entries {
            let versioned = entry.acc_ref.versioned();
            let unversioned = if entry.acc_ref.version.is_none() {
                "  [unversioned reference]"
            } else {
                ""
            };
            println!(
                "    {:<28} {} input(s){}",
                versioned,
                entry.referencing_inputs.len(),
                unversioned
            );
        }
    }
    Ok(())
}

/// Stable `class=count` summary, e.g. `refseq_transcript=30, refseq_genomic=15`.
fn format_class_counts(counts: &std::collections::BTreeMap<AccessionClass, usize>) -> String {
    counts
        .iter()
        .map(|(c, n)| format!("{}={}", c.as_str(), n))
        .collect::<Vec<_>>()
        .join(", ")
}

/// Recursively collect the accessions indexed by every `*.fai` file under `dir`.
/// A FASTA index is tab-separated; the first column is the sequence name
/// (the versioned accession as it keys lookups).
fn collect_fai_accessions(dir: &Path) -> Result<HashSet<String>, String> {
    let mut acc = HashSet::new();
    collect_fai_accessions_into(dir, &mut acc)?;
    Ok(acc)
}

fn collect_fai_accessions_into(dir: &Path, acc: &mut HashSet<String>) -> Result<(), String> {
    let entries = std::fs::read_dir(dir).map_err(|e| format!("read_dir {}: {e}", dir.display()))?;
    for entry in entries {
        let entry = entry.map_err(|e| format!("dir entry under {}: {e}", dir.display()))?;
        let path = entry.path();
        let file_type = entry
            .file_type()
            .map_err(|e| format!("file_type {}: {e}", path.display()))?;
        if file_type.is_dir() {
            collect_fai_accessions_into(&path, acc)?;
        } else if path.extension().and_then(|e| e.to_str()) == Some("fai") {
            let text = std::fs::read_to_string(&path)
                .map_err(|e| format!("read {}: {e}", path.display()))?;
            acc.extend(fai_names(&text));
        }
    }
    Ok(())
}

/// The sequence names (first column) of a FASTA `.fai` index.
fn fai_names(fai_text: &str) -> impl Iterator<Item = String> + '_ {
    fai_text
        .lines()
        .filter_map(|line| line.split('\t').next())
        .filter(|name| !name.is_empty())
        .map(|name| name.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fai_names_takes_first_tab_column() {
        // `name\tlength\toffset\tlinebases\tlinewidth` per `.fai` format.
        let fai = "NM_003002.2\t1382\t12\t70\t71\nNG_012337.3\t100\t0\t60\t61\n\n";
        let names: Vec<String> = fai_names(fai).collect();
        assert_eq!(
            names,
            vec!["NM_003002.2".to_string(), "NG_012337.3".to_string()]
        );
    }
}
