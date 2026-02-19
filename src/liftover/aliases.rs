//! Contig name aliases for different naming conventions.
//!
//! Handles mapping between different contig naming conventions:
//! - UCSC: chr1, chr2, ..., chrX, chrY, chrM
//! - RefSeq: NC_000001.10 (GRCh37), NC_000001.11 (GRCh38)
//! - Ensembl: 1, 2, ..., X, Y, MT

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

        // Mitochondrial
        aliases.add_alias("chrM", GenomeBuild::GRCh37, "NC_012920.1");
        aliases.add_alias("MT", GenomeBuild::GRCh37, "NC_012920.1");
        aliases.add_alias("NC_012920.1", GenomeBuild::GRCh37, "NC_012920.1");
        aliases.add_alias("chrM", GenomeBuild::GRCh38, "NC_012920.1");
        aliases.add_alias("MT", GenomeBuild::GRCh38, "NC_012920.1");
        aliases
            .refseq_to_ucsc
            .insert("NC_012920.1".to_string(), "chrM".to_string());
        aliases
            .refseq_to_ensembl
            .insert("NC_012920.1".to_string(), "MT".to_string());

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

#[cfg(test)]
mod tests {
    use super::*;

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
}
