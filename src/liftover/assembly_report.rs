//! Parser for NCBI `*_assembly_report.txt` files.
//!
//! NCBI ships an `assembly_report.txt` alongside each genome assembly (in the
//! same directory as the `*_genomic.fna.gz` that `ferro prepare` downloads).
//! It is the authoritative description of every sequence in the assembly:
//! which RefSeq accession maps to which molecule, and — from the header — the
//! assembly's name (e.g. `GRCh38.p14`).
//!
//! This parser turns that file into a data-driven `accession → assembly` map,
//! so build inference can consult the prepared reference's own report instead
//! of a hardcoded per-chromosome version table (#716).
//!
//! # Format
//!
//! Header lines begin with `#`. The assembly name is carried on a
//! `# Assembly name:  <name>` line. Data rows are tab-separated with the
//! columns (NCBI's stable layout):
//!
//! | idx | column |
//! |-----|--------|
//! | 0 | Sequence-Name (`1`, `X`, `MT`) |
//! | 1 | Sequence-Role (`assembled-molecule`, `alt-scaffold`, …) |
//! | 2 | Assigned-Molecule (`1`, `X`, `MT`) |
//! | 3 | Assigned-Molecule-Location/Type (`Chromosome`, `Mitochondrion`) |
//! | 4 | GenBank-Accn (`CM000663.2`) |
//! | 5 | Relationship (`=`) |
//! | 6 | RefSeq-Accn (`NC_000001.11`) |
//! | 7 | Assembly-Unit (`Primary Assembly`, `non-nuclear`) |
//! | 8 | Sequence-Length |
//! | 9 | UCSC-style-name (`chr1`, or `na`) |

/// A single sequence row from an NCBI assembly report.
///
/// Only the columns relevant to contig/build inference are retained. A `na`
/// in any source column is preserved verbatim (callers decide how to treat
/// missing UCSC names, etc.).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssemblyReportEntry {
    /// Sequence-Name (column 0), e.g. `1`, `X`, `MT`, `HSCHR1_CTG1_UNLOCALIZED`.
    pub sequence_name: String,
    /// Sequence-Role (column 1), e.g. `assembled-molecule`, `alt-scaffold`.
    pub sequence_role: String,
    /// Assigned-Molecule (column 2), e.g. `1`, `X`, `MT`, `na`.
    pub assigned_molecule: String,
    /// GenBank-Accn (column 4), e.g. `CM000663.2`.
    pub genbank_accession: String,
    /// RefSeq-Accn (column 6), e.g. `NC_000001.11` (may be `na`).
    pub refseq_accession: String,
    /// UCSC-style-name (column 9), e.g. `chr1` (may be `na`).
    pub ucsc_name: String,
}

impl AssemblyReportEntry {
    /// `true` for rows whose Sequence-Role is `assembled-molecule` — the
    /// primary chromosomes (1–22, X, Y) plus the mitochondrion. These are the
    /// rows that carry build-distinguishing `NC_` accessions; alt/patch/
    /// unplaced scaffolds are excluded from build inference.
    pub fn is_assembled_molecule(&self) -> bool {
        self.sequence_role == "assembled-molecule"
    }

    /// `true` when the RefSeq accession column carries a real accession (not
    /// the `na` placeholder NCBI uses for GenBank-only sequences).
    pub fn has_refseq_accession(&self) -> bool {
        !self.refseq_accession.is_empty() && self.refseq_accession != "na"
    }
}

/// A parsed NCBI assembly report: the assembly name from the header plus every
/// data row.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AssemblyReport {
    /// Assembly name from the `# Assembly name:` header line, e.g.
    /// `GRCh38.p14`. Empty if the header line was absent.
    pub assembly_name: String,
    /// All data (non-comment) rows that parsed into at least the columns
    /// needed for inference.
    pub entries: Vec<AssemblyReportEntry>,
}

impl AssemblyReport {
    /// Iterator over the assembled-molecule entries that carry a RefSeq
    /// accession — the `(NC_ accession, …)` rows build inference keys on.
    pub fn assembled_molecules_with_refseq(&self) -> impl Iterator<Item = &AssemblyReportEntry> {
        self.entries
            .iter()
            .filter(|e| e.is_assembled_molecule() && e.has_refseq_accession())
    }
}

/// The `#` prefix that marks a header / comment line.
const COMMENT_PREFIX: &str = "#";
/// Header key carrying the assembly name.
const ASSEMBLY_NAME_KEY: &str = "# Assembly name:";
/// Number of tab-separated columns NCBI emits; rows with fewer are skipped.
const MIN_COLUMNS: usize = 10;

/// Parse the text of an NCBI `assembly_report.txt`.
///
/// Comment/header lines (`#`-prefixed) are scanned for the assembly name;
/// everything else is parsed as a tab-separated data row. Rows with fewer than
/// [`MIN_COLUMNS`] columns are skipped (defensive against truncated or
/// reformatted files). The parse never fails — a malformed file simply yields
/// fewer entries — because a missing/garbled report should degrade to the
/// hardcoded fallback, not abort reference loading.
pub fn parse_assembly_report(text: &str) -> AssemblyReport {
    let mut assembly_name = String::new();
    let mut entries = Vec::new();

    for line in text.lines() {
        if line.starts_with(COMMENT_PREFIX) {
            // Header line. Capture the assembly name; ignore everything else
            // (including the `##`-prefixed Assembly-Units block, which still
            // starts with `#`).
            if let Some(name) = line.strip_prefix(ASSEMBLY_NAME_KEY) {
                let name = name.trim();
                if !name.is_empty() {
                    assembly_name = name.to_string();
                }
            }
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < MIN_COLUMNS {
            continue;
        }
        entries.push(AssemblyReportEntry {
            sequence_name: cols[0].to_string(),
            sequence_role: cols[1].to_string(),
            assigned_molecule: cols[2].to_string(),
            genbank_accession: cols[4].to_string(),
            refseq_accession: cols[6].to_string(),
            ucsc_name: cols[9].to_string(),
        });
    }

    AssemblyReport {
        assembly_name,
        entries,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A faithful slice of an NCBI GRCh38 assembly_report.txt: the header
    /// block (including the `##` Assembly-Units sub-block), a few
    /// assembled-molecule chromosomes, the mitochondrion, and one
    /// alt-scaffold row that must be excluded from build inference.
    const GRCH38_SAMPLE: &str = "\
# Assembly name:  GRCh38.p14
# Organism name:  Homo sapiens (human)
# Assembly level: Chromosome
# RefSeq assembly accession: GCF_000001405.40
#
## Assembly-Units:
## GenBank Unit Accession\tRefSeq Unit Accession\tAssembly-Unit name
## GCA_000001305.2\tGCF_000001305.15\tPrimary Assembly
# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name
1\tassembled-molecule\t1\tChromosome\tCM000663.2\t=\tNC_000001.11\tPrimary Assembly\t248956422\tchr1
17\tassembled-molecule\t17\tChromosome\tCM000679.2\t=\tNC_000017.11\tPrimary Assembly\t83257441\tchr17
X\tassembled-molecule\tX\tChromosome\tCM000685.2\t=\tNC_000023.11\tPrimary Assembly\t156040895\tchrX
MT\tassembled-molecule\tMT\tMitochondrion\tJ01415.2\t=\tNC_012920.1\tnon-nuclear\t16569\tchrM
HSCHR1_CTG1\talt-scaffold\t1\tChromosome\tKI270762.1\t=\tNT_187515.1\tALT_REF_LOCI_1\t354444\tchr1_KI270762v1_alt
";

    #[test]
    fn parses_assembly_name_from_header() {
        let report = parse_assembly_report(GRCH38_SAMPLE);
        assert_eq!(report.assembly_name, "GRCh38.p14");
    }

    #[test]
    fn parses_assembled_molecule_rows() {
        let report = parse_assembly_report(GRCH38_SAMPLE);
        // 5 data rows total (4 assembled-molecule + 1 alt-scaffold).
        assert_eq!(report.entries.len(), 5);

        let chr1 = &report.entries[0];
        assert_eq!(chr1.sequence_name, "1");
        assert_eq!(chr1.sequence_role, "assembled-molecule");
        assert_eq!(chr1.refseq_accession, "NC_000001.11");
        assert_eq!(chr1.genbank_accession, "CM000663.2");
        assert_eq!(chr1.ucsc_name, "chr1");
        assert!(chr1.is_assembled_molecule());
        assert!(chr1.has_refseq_accession());
    }

    #[test]
    fn assembled_molecules_with_refseq_excludes_alt_scaffolds() {
        let report = parse_assembly_report(GRCH38_SAMPLE);
        let refseqs: Vec<&str> = report
            .assembled_molecules_with_refseq()
            .map(|e| e.refseq_accession.as_str())
            .collect();
        assert_eq!(
            refseqs,
            vec![
                "NC_000001.11",
                "NC_000017.11",
                "NC_000023.11",
                "NC_012920.1"
            ],
            "alt-scaffold NT_187515.1 must be excluded from inference rows"
        );
    }

    #[test]
    fn skips_comment_and_short_rows() {
        // A header without an Assembly-name line, a too-short row, and a blank
        // line: the parser must tolerate all three and still find the one
        // valid data row.
        let text = "\
# some banner
1\tassembled-molecule\t1\tChromosome\tCM000663.2\t=\tNC_000001.11\tPrimary Assembly\t248956422\tchr1

short\trow\twith\tfew\tcols
";
        let report = parse_assembly_report(text);
        assert_eq!(report.assembly_name, "");
        assert_eq!(report.entries.len(), 1);
        assert_eq!(report.entries[0].refseq_accession, "NC_000001.11");
    }

    #[test]
    fn treats_na_refseq_as_absent() {
        let text = "\
# Assembly name:  GRCh38.p14
HSCHR1_RANDOM\tunplaced-scaffold\tna\tna\tGL000008.2\t<>\tna\tPrimary Assembly\t209709\tchr1_GL000008v2_random
";
        let report = parse_assembly_report(text);
        assert_eq!(report.entries.len(), 1);
        assert!(!report.entries[0].has_refseq_accession());
        assert_eq!(report.assembled_molecules_with_refseq().count(), 0);
    }

    #[test]
    fn empty_input_yields_empty_report() {
        let report = parse_assembly_report("");
        assert_eq!(report, AssemblyReport::default());
    }
}
