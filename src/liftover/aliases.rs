//! Contig name aliases for different naming conventions.
//!
//! Handles mapping between different contig naming conventions:
//! - UCSC: chr1, chr2, ..., chrX, chrY, chrM
//! - RefSeq: NC_000001.10 (GRCh37), NC_000001.11 (GRCh38)
//! - Ensembl: 1, 2, ..., X, Y, MT

use crate::liftover::assembly_report::AssemblyReport;
use crate::reference::transcript::GenomeBuild;
use std::collections::HashMap;

/// Contig alias mappings.
#[derive(Debug, Clone, Default)]
pub struct ContigAliases {
    /// Maps (contig_name, build) -> canonical RefSeq name
    to_refseq: HashMap<(String, GenomeBuild), String>,
    /// Maps RefSeq name -> UCSC name
    refseq_to_ucsc: HashMap<String, String>,
    /// Maps RefSeq name -> Ensembl name
    refseq_to_ensembl: HashMap<String, String>,
}

impl ContigAliases {
    /// Create empty alias mapping.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create default human genome aliases.
    pub fn default_human() -> Self {
        let mut aliases = Self::new();

        // Autosomes
        for i in 1..=22 {
            let ucsc = format!("chr{}", i);
            let ensembl = format!("{}", i);
            let refseq_37 = format!("NC_{:06}.{}", i, Self::grch37_version(i));
            let refseq_38 = format!("NC_{:06}.{}", i, Self::grch38_version(i));

            // GRCh37 mappings
            aliases.add_alias(&ucsc, GenomeBuild::GRCh37, &refseq_37);
            aliases.add_alias(&ensembl, GenomeBuild::GRCh37, &refseq_37);
            aliases.add_alias(&refseq_37, GenomeBuild::GRCh37, &refseq_37);

            // GRCh38 mappings
            aliases.add_alias(&ucsc, GenomeBuild::GRCh38, &refseq_38);
            aliases.add_alias(&ensembl, GenomeBuild::GRCh38, &refseq_38);
            aliases.add_alias(&refseq_38, GenomeBuild::GRCh38, &refseq_38);

            // Reverse lookups
            aliases
                .refseq_to_ucsc
                .insert(refseq_37.clone(), ucsc.clone());
            aliases
                .refseq_to_ucsc
                .insert(refseq_38.clone(), ucsc.clone());
            aliases.refseq_to_ensembl.insert(refseq_37, ensembl.clone());
            aliases.refseq_to_ensembl.insert(refseq_38, ensembl);
        }

        // X chromosome
        aliases.add_alias("chrX", GenomeBuild::GRCh37, "NC_000023.10");
        aliases.add_alias("X", GenomeBuild::GRCh37, "NC_000023.10");
        aliases.add_alias("NC_000023.10", GenomeBuild::GRCh37, "NC_000023.10");
        aliases.add_alias("chrX", GenomeBuild::GRCh38, "NC_000023.11");
        aliases.add_alias("X", GenomeBuild::GRCh38, "NC_000023.11");
        aliases.add_alias("NC_000023.11", GenomeBuild::GRCh38, "NC_000023.11");
        aliases
            .refseq_to_ucsc
            .insert("NC_000023.10".to_string(), "chrX".to_string());
        aliases
            .refseq_to_ucsc
            .insert("NC_000023.11".to_string(), "chrX".to_string());
        aliases
            .refseq_to_ensembl
            .insert("NC_000023.10".to_string(), "X".to_string());
        aliases
            .refseq_to_ensembl
            .insert("NC_000023.11".to_string(), "X".to_string());

        // Y chromosome
        aliases.add_alias("chrY", GenomeBuild::GRCh37, "NC_000024.9");
        aliases.add_alias("Y", GenomeBuild::GRCh37, "NC_000024.9");
        aliases.add_alias("NC_000024.9", GenomeBuild::GRCh37, "NC_000024.9");
        aliases.add_alias("chrY", GenomeBuild::GRCh38, "NC_000024.10");
        aliases.add_alias("Y", GenomeBuild::GRCh38, "NC_000024.10");
        aliases.add_alias("NC_000024.10", GenomeBuild::GRCh38, "NC_000024.10");
        aliases
            .refseq_to_ucsc
            .insert("NC_000024.9".to_string(), "chrY".to_string());
        aliases
            .refseq_to_ucsc
            .insert("NC_000024.10".to_string(), "chrY".to_string());
        aliases
            .refseq_to_ensembl
            .insert("NC_000024.9".to_string(), "Y".to_string());
        aliases
            .refseq_to_ensembl
            .insert("NC_000024.10".to_string(), "Y".to_string());

        // Mitochondrial — NC_012920.1 is the canonical rCRS accession on
        // both GRCh37 and GRCh38 (mtDNA was not re-assembled between
        // builds), so register a self-alias under each build. Without
        // the GRCh38 self-alias `infer_genome_build_from_accession` would
        // misclassify NC_012920.1 inputs as GRCh37-only.
        aliases.add_alias("chrM", GenomeBuild::GRCh37, "NC_012920.1");
        aliases.add_alias("MT", GenomeBuild::GRCh37, "NC_012920.1");
        aliases.add_alias("NC_012920.1", GenomeBuild::GRCh37, "NC_012920.1");
        aliases.add_alias("chrM", GenomeBuild::GRCh38, "NC_012920.1");
        aliases.add_alias("MT", GenomeBuild::GRCh38, "NC_012920.1");
        aliases.add_alias("NC_012920.1", GenomeBuild::GRCh38, "NC_012920.1");
        aliases
            .refseq_to_ucsc
            .insert("NC_012920.1".to_string(), "chrM".to_string());
        aliases
            .refseq_to_ensembl
            .insert("NC_012920.1".to_string(), "MT".to_string());

        aliases
    }

    /// Build a contig-alias table from parsed NCBI assembly reports, one per
    /// genome build.
    ///
    /// Data-driven counterpart to [`default_human`](Self::default_human): the
    /// `(name, build) → RefSeq` and reverse maps are derived from the report's
    /// own `assembled-molecule` rows rather than a hardcoded per-chromosome
    /// version table (#716). For each such row the RefSeq accession is
    /// registered as a self-alias (so [`infer_genome_build_from_accession_with`]
    /// can classify it), alongside its UCSC-style name (column 9) and Ensembl
    /// name (the Assigned-Molecule, column 2). Rows whose UCSC name is the `na`
    /// placeholder contribute no UCSC alias.
    ///
    /// Multiple reports are merged, so passing both the GRCh37 and GRCh38
    /// reports yields one table covering both builds — an accession present in
    /// only one build resolves under that build alone, which is exactly the
    /// signal build inference needs.
    pub fn from_assembly_reports(reports: &[(GenomeBuild, &AssemblyReport)]) -> Self {
        let mut aliases = Self::new();
        for (build, report) in reports {
            for entry in report.assembled_molecules_with_refseq() {
                let refseq = entry.refseq_accession.as_str();
                let ensembl = entry.assigned_molecule.as_str();
                let ucsc = entry.ucsc_name.as_str();

                aliases.add_alias(refseq, *build, refseq);
                if !ensembl.is_empty() && ensembl != "na" {
                    aliases.add_alias(ensembl, *build, refseq);
                    aliases
                        .refseq_to_ensembl
                        .insert(refseq.to_string(), ensembl.to_string());
                }
                if !ucsc.is_empty() && ucsc != "na" {
                    aliases.add_alias(ucsc, *build, refseq);
                    aliases
                        .refseq_to_ucsc
                        .insert(refseq.to_string(), ucsc.to_string());
                }
            }
        }
        aliases
    }

    /// Get GRCh37 version number for a chromosome.
    fn grch37_version(chr: i32) -> i32 {
        match chr {
            1 => 10,
            2 => 11,
            3 => 11,
            4 => 11,
            5 => 9,
            6 => 11,
            7 => 13,
            8 => 10,
            9 => 11,
            10 => 10,
            11 => 9,
            12 => 11,
            13 => 10,
            14 => 8,
            15 => 9,
            16 => 9,
            17 => 10,
            18 => 9,
            19 => 9,
            20 => 10,
            21 => 8,
            22 => 10,
            _ => 1,
        }
    }

    /// Get GRCh38 version number for a chromosome.
    fn grch38_version(chr: i32) -> i32 {
        match chr {
            1 => 11,
            2 => 12,
            3 => 12,
            4 => 12,
            5 => 10,
            6 => 12,
            7 => 14,
            8 => 11,
            9 => 12,
            10 => 11,
            11 => 10,
            12 => 12,
            13 => 11,
            14 => 9,
            15 => 10,
            16 => 10,
            17 => 11,
            18 => 10,
            19 => 10,
            20 => 11,
            21 => 9,
            22 => 11,
            _ => 1,
        }
    }

    /// Add an alias mapping.
    pub fn add_alias(&mut self, name: &str, build: GenomeBuild, refseq: &str) {
        self.to_refseq
            .insert((name.to_string(), build), refseq.to_string());
    }

    /// Resolve a contig name to its canonical RefSeq name.
    pub fn resolve_to_refseq(&self, name: &str, build: GenomeBuild) -> Option<&str> {
        self.to_refseq
            .get(&(name.to_string(), build))
            .map(|s| s.as_str())
    }

    /// Resolve a RefSeq name to UCSC format.
    pub fn refseq_to_ucsc(&self, refseq: &str) -> Option<&str> {
        self.refseq_to_ucsc.get(refseq).map(|s| s.as_str())
    }

    /// Resolve a RefSeq name to Ensembl format.
    pub fn refseq_to_ensembl(&self, refseq: &str) -> Option<&str> {
        self.refseq_to_ensembl.get(refseq).map(|s| s.as_str())
    }

    /// Check if two contig names are equivalent for a given build.
    pub fn are_equivalent(&self, name1: &str, name2: &str, build: GenomeBuild) -> bool {
        let refseq1 = self.resolve_to_refseq(name1, build);
        let refseq2 = self.resolve_to_refseq(name2, build);
        match (refseq1, refseq2) {
            (Some(r1), Some(r2)) => r1 == r2,
            _ => name1 == name2,
        }
    }

    /// Normalize a contig name to the preferred format for a build.
    ///
    /// Returns the RefSeq accession if found, otherwise returns the input unchanged.
    pub fn normalize(&self, name: &str, build: GenomeBuild) -> String {
        self.resolve_to_refseq(name, build)
            .map(|s| s.to_string())
            .unwrap_or_else(|| name.to_string())
    }
}

/// Lazily-constructed singleton of [`ContigAliases::default_human`] so
/// repeated build inference (per projection) does not rebuild the
/// ~150-entry alias table every call.
fn default_human_aliases() -> &'static ContigAliases {
    use std::sync::OnceLock;
    static ALIASES: OnceLock<ContigAliases> = OnceLock::new();
    ALIASES.get_or_init(ContigAliases::default_human)
}

/// Infer the canonical genome build name for a chromosome accession.
///
/// Sources of build information consulted, in order:
///   1. `accession.assembly` — assembly-style references like
///      `GRCh37(chr1):g.…` carry the build verbatim.
///   2. `NC_*` accessions are looked up in [`ContigAliases::default_human`].
///      `NC_*` carries a build-distinguishing version (e.g.
///      `NC_000017.10` is GRCh37 chr17, `NC_000017.11` is GRCh38). NC
///      accessions present under both builds (e.g. mtDNA `NC_012920.1`)
///      resolve to GRCh38, our default-when-ambiguous.
///   3. Everything else (`NG_*`, `LRG_*`, `NW_*`, unknown prefixes, NC
///      accessions outside the human alias table) returns `None` so the
///      caller falls back to the build-agnostic probe order.
///
/// Centralized here so call sites that need to disambiguate cdot lookups
/// against multi-build data (e.g. `MultiFastaProvider`, `VariantProjector`)
/// share one implementation rather than duplicating the GRCh37/GRCh38
/// inference logic.
pub fn infer_genome_build_from_accession(
    accession: &crate::hgvs::variant::Accession,
) -> Option<&'static str> {
    infer_genome_build_from_accession_with(default_human_aliases(), accession)
}

/// Like [`infer_genome_build_from_accession`] but consulting an explicitly
/// supplied [`ContigAliases`] table rather than the bundled
/// [`default_human`](ContigAliases::default_human) heuristic.
///
/// This is the data-driven entry point (#716): pass a table built from the
/// prepared reference's own `assembly_report.txt`
/// ([`ContigAliases::from_assembly_reports`]) and build inference becomes
/// authoritative — it classifies any `NC_` accession the report describes,
/// including versions absent from the hardcoded table. The `accession.assembly`
/// explicit-build field and the `NC_`-prefix gate are honored identically to
/// the bundled path, so the only behavioral difference is the source of the
/// accession→build mapping.
pub fn infer_genome_build_from_accession_with(
    aliases: &ContigAliases,
    accession: &crate::hgvs::variant::Accession,
) -> Option<&'static str> {
    // Assembly-style references (`GRCh37(chr1)`, `GRCh38(chrX)`) carry the
    // build name explicitly; honor it before the prefix check.
    if let Some(asm) = accession.assembly.as_deref() {
        match asm {
            "GRCh37" => return Some("GRCh37"),
            "GRCh38" => return Some("GRCh38"),
            _ => {}
        }
    }
    if &*accession.prefix != "NC" {
        return None;
    }
    let acc_str = accession.full();
    // GRCh38 is checked first so an accession present under both builds (e.g.
    // mtDNA `NC_012920.1`) resolves to GRCh38, ferro's default-when-ambiguous.
    if aliases
        .resolve_to_refseq(&acc_str, GenomeBuild::GRCh38)
        .is_some()
    {
        return Some("GRCh38");
    }
    if aliases
        .resolve_to_refseq(&acc_str, GenomeBuild::GRCh37)
        .is_some()
    {
        return Some("GRCh37");
    }
    None
}

/// Infer the genome build for an accession, consulting a data-driven
/// [`ContigAliases`] (from the prepared reference's own assembly report, #716)
/// before falling back to the bundled hardcoded table.
///
/// `table` is `None` when the prepared reference carries no assembly report
/// (old manifest, bare library use), in which case this is exactly
/// [`infer_genome_build_from_accession`]. When `Some`, the data-driven table is
/// authoritative for any accession it describes; accessions it omits still
/// classify via the hardcoded fallback. Returns the first build found, so the
/// data-driven layer wins on disagreement — safe because the two only overlap
/// on accessions identical across builds (e.g. mtDNA), which the GRCh38-first
/// check order in [`infer_genome_build_from_accession_with`] resolves the same
/// way either source would.
pub fn infer_genome_build_layered(
    table: Option<&ContigAliases>,
    accession: &crate::hgvs::variant::Accession,
) -> Option<&'static str> {
    if let Some(table) = table {
        if let Some(build) = infer_genome_build_from_accession_with(table, accession) {
            return Some(build);
        }
    }
    infer_genome_build_from_accession(accession)
}

/// Normalize a user-supplied assembly name to ferro's canonical build string
/// (`"GRCh37"` / `"GRCh38"`), or `None` if unrecognized (#715).
///
/// Accepts the common aliases (`hg19`/`b37`/`37` → `GRCh37`,
/// `hg38`/`b38`/`38` → `GRCh38`), case-insensitively. The canonical strings
/// returned here are the exact `&'static str` values
/// [`infer_genome_build_from_accession`] yields and that
/// `select_placement_for_build` / cdot's `get_transcript_on_build` compare
/// against, so callers should validate at the boundary (CLI flag, Python kwarg)
/// and only ever hand the canonical value to
/// [`crate::project::VariantProjector::with_assembly`].
pub fn normalize_assembly_name(name: &str) -> Option<&'static str> {
    match name.trim().to_ascii_lowercase().as_str() {
        "grch37" | "hg19" | "b37" | "37" => Some("GRCh37"),
        "grch38" | "hg38" | "b38" | "38" => Some("GRCh38"),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_assembly_name() {
        for name in ["GRCh37", "grch37", "hg19", "HG19", "b37", "37", " grch37 "] {
            assert_eq!(normalize_assembly_name(name), Some("GRCh37"), "{name:?}");
        }
        for name in ["GRCh38", "grch38", "hg38", "b38", "38"] {
            assert_eq!(normalize_assembly_name(name), Some("GRCh38"), "{name:?}");
        }
        for name in ["", "hg20", "GRCh39", "chm13", "t2t"] {
            assert_eq!(normalize_assembly_name(name), None, "{name:?}");
        }
    }

    #[test]
    fn test_default_human_autosomes() {
        let aliases = ContigAliases::default_human();

        // Chr1 GRCh37
        assert_eq!(
            aliases.resolve_to_refseq("chr1", GenomeBuild::GRCh37),
            Some("NC_000001.10")
        );
        assert_eq!(
            aliases.resolve_to_refseq("1", GenomeBuild::GRCh37),
            Some("NC_000001.10")
        );

        // Chr1 GRCh38
        assert_eq!(
            aliases.resolve_to_refseq("chr1", GenomeBuild::GRCh38),
            Some("NC_000001.11")
        );
        assert_eq!(
            aliases.resolve_to_refseq("1", GenomeBuild::GRCh38),
            Some("NC_000001.11")
        );
    }

    #[test]
    fn test_sex_chromosomes() {
        let aliases = ContigAliases::default_human();

        assert_eq!(
            aliases.resolve_to_refseq("chrX", GenomeBuild::GRCh37),
            Some("NC_000023.10")
        );
        assert_eq!(
            aliases.resolve_to_refseq("X", GenomeBuild::GRCh38),
            Some("NC_000023.11")
        );
        assert_eq!(
            aliases.resolve_to_refseq("chrY", GenomeBuild::GRCh37),
            Some("NC_000024.9")
        );
    }

    #[test]
    fn test_mitochondrial() {
        let aliases = ContigAliases::default_human();

        assert_eq!(
            aliases.resolve_to_refseq("chrM", GenomeBuild::GRCh37),
            Some("NC_012920.1")
        );
        assert_eq!(
            aliases.resolve_to_refseq("MT", GenomeBuild::GRCh38),
            Some("NC_012920.1")
        );
    }

    #[test]
    fn test_reverse_lookup() {
        let aliases = ContigAliases::default_human();

        assert_eq!(aliases.refseq_to_ucsc("NC_000001.10"), Some("chr1"));
        assert_eq!(aliases.refseq_to_ensembl("NC_000001.10"), Some("1"));
    }

    #[test]
    fn test_equivalence() {
        let aliases = ContigAliases::default_human();

        assert!(aliases.are_equivalent("chr1", "1", GenomeBuild::GRCh37));
        assert!(aliases.are_equivalent("chr1", "NC_000001.10", GenomeBuild::GRCh37));
        assert!(!aliases.are_equivalent("chr1", "chr2", GenomeBuild::GRCh37));
    }

    #[test]
    fn test_normalize() {
        let aliases = ContigAliases::default_human();

        assert_eq!(
            aliases.normalize("chr1", GenomeBuild::GRCh37),
            "NC_000001.10"
        );
        assert_eq!(aliases.normalize("unknown", GenomeBuild::GRCh37), "unknown");
    }

    // ----- #716: data-driven build inference from assembly reports -----

    use crate::hgvs::variant::Accession;
    use crate::liftover::assembly_report::parse_assembly_report;

    /// A minimal but format-faithful assembly_report for one or more
    /// chromosomes: `(ucsc, ensembl, refseq)` rows under the given name.
    fn synthetic_report(name: &str, rows: &[(&str, &str, &str)]) -> String {
        let mut s = format!("# Assembly name:  {name}\n");
        for (ucsc, ensembl, refseq) in rows {
            // cols: name role assigned type genbank rel refseq unit length ucsc
            s.push_str(&format!(
                "{ensembl}\tassembled-molecule\t{ensembl}\tChromosome\tCM000000.1\t=\t{refseq}\tPrimary Assembly\t1000\t{ucsc}\n"
            ));
        }
        s
    }

    #[test]
    fn from_assembly_reports_builds_resolvable_aliases() {
        let report = parse_assembly_report(&synthetic_report(
            "GRCh38.p14",
            &[("chr1", "1", "NC_000001.11"), ("chrX", "X", "NC_000023.11")],
        ));
        let aliases = ContigAliases::from_assembly_reports(&[(GenomeBuild::GRCh38, &report)]);

        assert_eq!(
            aliases.resolve_to_refseq("chr1", GenomeBuild::GRCh38),
            Some("NC_000001.11")
        );
        assert_eq!(
            aliases.resolve_to_refseq("1", GenomeBuild::GRCh38),
            Some("NC_000001.11")
        );
        // Self-alias powers build inference.
        assert_eq!(
            aliases.resolve_to_refseq("NC_000023.11", GenomeBuild::GRCh38),
            Some("NC_000023.11")
        );
        assert_eq!(aliases.refseq_to_ucsc("NC_000001.11"), Some("chr1"));
        assert_eq!(aliases.refseq_to_ensembl("NC_000023.11"), Some("X"));
    }

    #[test]
    fn infer_with_classifies_from_merged_reports() {
        let g38 = parse_assembly_report(&synthetic_report(
            "GRCh38.p14",
            &[("chr1", "1", "NC_000001.11")],
        ));
        let g37 = parse_assembly_report(&synthetic_report(
            "GRCh37.p13",
            &[("chr1", "1", "NC_000001.10")],
        ));
        let aliases = ContigAliases::from_assembly_reports(&[
            (GenomeBuild::GRCh38, &g38),
            (GenomeBuild::GRCh37, &g37),
        ]);

        assert_eq!(
            infer_genome_build_from_accession_with(
                &aliases,
                &Accession::new("NC", "000001", Some(11))
            ),
            Some("GRCh38")
        );
        assert_eq!(
            infer_genome_build_from_accession_with(
                &aliases,
                &Accession::new("NC", "000001", Some(10))
            ),
            Some("GRCh37")
        );
    }

    #[test]
    fn infer_with_hardens_beyond_the_hardcoded_table() {
        // A version the bundled table does not know (`NC_000001.12` is a
        // hypothetical future GRCh38 chr1 build). The hardcoded default
        // declines; a report-built table classifies it authoritatively.
        let unknown_to_default = Accession::new("NC", "000001", Some(12));
        assert_eq!(
            infer_genome_build_from_accession(&unknown_to_default),
            None,
            "the bundled version table must not know NC_000001.12"
        );

        let report = parse_assembly_report(&synthetic_report(
            "GRCh38.p15",
            &[("chr1", "1", "NC_000001.12")],
        ));
        let aliases = ContigAliases::from_assembly_reports(&[(GenomeBuild::GRCh38, &report)]);
        assert_eq!(
            infer_genome_build_from_accession_with(&aliases, &unknown_to_default),
            Some("GRCh38"),
            "report-derived aliases must classify accessions absent from the bundled table"
        );
    }

    #[test]
    fn infer_default_path_is_unchanged() {
        // The public no-arg entry point delegates to the bundled table; its
        // behavior must be identical to before the #716 refactor.
        assert_eq!(
            infer_genome_build_from_accession(&Accession::new("NC", "000017", Some(11))),
            Some("GRCh38")
        );
        assert_eq!(
            infer_genome_build_from_accession(&Accession::new("NC", "000017", Some(10))),
            Some("GRCh37")
        );
        // mtDNA is shared across builds and resolves to the GRCh38 default.
        assert_eq!(
            infer_genome_build_from_accession(&Accession::new("NC", "012920", Some(1))),
            Some("GRCh38")
        );
        // Non-NC prefixes still decline.
        assert_eq!(
            infer_genome_build_from_accession(&Accession::new("NG", "012337", Some(1))),
            None
        );
    }

    /// A GRCh38-only report whose chr17 row carries a version (`.99`) the
    /// hardcoded table does not know — proves the data-driven layer adds reach.
    const GRCH38_FUTURE_SAMPLE: &str = "\
# Assembly name:  GRCh38.future
# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name
17\tassembled-molecule\t17\tChromosome\tCM000679.9\t=\tNC_000017.99\tPrimary Assembly\t83257441\tchr17
";

    #[test]
    fn layered_classifies_report_only_version() {
        use crate::hgvs::variant::Accession;
        let report = parse_assembly_report(GRCH38_FUTURE_SAMPLE);
        let table = ContigAliases::from_assembly_reports(&[(GenomeBuild::GRCh38, &report)]);

        let future = Accession::new("NC", "000017", Some(99));
        // Hardcoded table has never seen .99 → None.
        assert_eq!(infer_genome_build_from_accession(&future), None);
        // Data-driven layer classifies it from the report.
        assert_eq!(
            infer_genome_build_layered(Some(&table), &future),
            Some("GRCh38")
        );
    }

    #[test]
    fn layered_falls_back_to_hardcoded_when_absent_from_report() {
        use crate::hgvs::variant::Accession;
        let report = parse_assembly_report(GRCH38_FUTURE_SAMPLE);
        let table = ContigAliases::from_assembly_reports(&[(GenomeBuild::GRCh38, &report)]);

        // GRCh37 chr1 (.10) is not in this GRCh38-only report → hardcoded fallback.
        let grch37_chr1 = Accession::new("NC", "000001", Some(10));
        assert_eq!(
            infer_genome_build_layered(Some(&table), &grch37_chr1),
            Some("GRCh37")
        );
    }

    #[test]
    fn layered_with_no_table_matches_hardcoded() {
        use crate::hgvs::variant::Accession;
        let grch38_chr17 = Accession::new("NC", "000017", Some(11));
        assert_eq!(
            infer_genome_build_layered(None, &grch38_chr17),
            infer_genome_build_from_accession(&grch38_chr17)
        );
    }
}
