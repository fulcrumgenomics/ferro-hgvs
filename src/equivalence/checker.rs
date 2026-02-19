//! Equivalence checker implementation.

use crate::error::FerroError;
use crate::hgvs::variant::HgvsVariant;
use crate::normalize::{NormalizeConfig, Normalizer};
use crate::reference::ReferenceProvider;

/// Level of equivalence between two variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EquivalenceLevel {
    /// Exactly identical (same string representation).
    Identical,
    /// Equivalent after normalization (e.g., shifted to same position).
    NormalizedMatch,
    /// Same variant but different accession versions (e.g., NM_000088.3 vs NM_000088.4).
    AccessionVersionDifference,
    /// Not equivalent - represent different changes.
    NotEquivalent,
}

impl EquivalenceLevel {
    /// Returns true if the variants are considered equivalent.
    pub fn is_equivalent(&self) -> bool {
        !matches!(self, EquivalenceLevel::NotEquivalent)
    }

    /// Returns a human-readable description of the equivalence level.
    pub fn description(&self) -> &'static str {
        match self {
            EquivalenceLevel::Identical => "Identical representation",
            EquivalenceLevel::NormalizedMatch => "Equivalent after normalization",
            EquivalenceLevel::AccessionVersionDifference => {
                "Same variant, different accession versions"
            }
            EquivalenceLevel::NotEquivalent => "Not equivalent",
        }
    }
}

impl std::fmt::Display for EquivalenceLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.description())
    }
}

/// Result of an equivalence check with additional details.
#[derive(Debug, Clone)]
pub struct EquivalenceResult {
    /// The determined equivalence level.
    pub level: EquivalenceLevel,
    /// The normalized form of the first variant (if normalization was performed).
    pub normalized_first: Option<String>,
    /// The normalized form of the second variant (if normalization was performed).
    pub normalized_second: Option<String>,
    /// Additional notes about the comparison.
    pub notes: Vec<String>,
}

impl EquivalenceResult {
    /// Create a new equivalence result.
    pub fn new(level: EquivalenceLevel) -> Self {
        Self {
            level,
            normalized_first: None,
            normalized_second: None,
            notes: Vec::new(),
        }
    }

    /// Add normalized forms to the result.
    pub fn with_normalized(mut self, first: String, second: String) -> Self {
        self.normalized_first = Some(first);
        self.normalized_second = Some(second);
        self
    }

    /// Add a note to the result.
    pub fn with_note(mut self, note: impl Into<String>) -> Self {
        self.notes.push(note.into());
        self
    }

    /// Returns true if the variants are considered equivalent.
    pub fn is_equivalent(&self) -> bool {
        self.level.is_equivalent()
    }
}

/// Checks equivalence between HGVS variants.
///
/// The checker determines if two variants represent the same genomic change,
/// even if they have different string representations.
pub struct EquivalenceChecker<P: ReferenceProvider> {
    normalizer: Normalizer<P>,
}

impl<P: ReferenceProvider> EquivalenceChecker<P> {
    /// Create a new equivalence checker with a reference provider.
    ///
    /// # Arguments
    ///
    /// * `provider` - Reference sequence provider for normalization
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::MockProvider;
    /// use ferro_hgvs::equivalence::EquivalenceChecker;
    ///
    /// let provider = MockProvider::with_test_data();
    /// let checker = EquivalenceChecker::new(provider);
    /// ```
    pub fn new(provider: P) -> Self {
        Self {
            normalizer: Normalizer::new(provider),
        }
    }

    /// Create a new equivalence checker with a custom normalizer configuration.
    ///
    /// # Arguments
    ///
    /// * `provider` - Reference sequence provider
    /// * `config` - Normalization configuration
    pub fn with_config(provider: P, config: NormalizeConfig) -> Self {
        Self {
            normalizer: Normalizer::with_config(provider, config),
        }
    }

    /// Check if two variants are equivalent.
    ///
    /// This method compares two variants and determines their level of equivalence.
    /// It first checks for string identity, then normalizes both variants and
    /// compares the normalized forms.
    ///
    /// # Arguments
    ///
    /// * `v1` - First variant
    /// * `v2` - Second variant
    ///
    /// # Returns
    ///
    /// * `Ok(EquivalenceResult)` - The equivalence result with details
    /// * `Err(FerroError)` - If normalization fails
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::{parse_hgvs, MockProvider};
    /// use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
    ///
    /// let provider = MockProvider::with_test_data();
    /// let checker = EquivalenceChecker::new(provider);
    ///
    /// let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
    /// let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
    ///
    /// let result = checker.check(&v1, &v2).unwrap();
    /// assert!(result.is_equivalent());
    /// ```
    pub fn check(
        &self,
        v1: &HgvsVariant,
        v2: &HgvsVariant,
    ) -> Result<EquivalenceResult, FerroError> {
        let str1 = v1.to_string();
        let str2 = v2.to_string();

        // Check for identical string representation
        if str1 == str2 {
            return Ok(EquivalenceResult::new(EquivalenceLevel::Identical)
                .with_normalized(str1.clone(), str2.clone()));
        }

        // Check for accession version differences
        if let Some(result) = self.check_accession_version_difference(v1, v2, &str1, &str2) {
            return Ok(result);
        }

        // Normalize both variants and compare
        let norm1 = self.normalizer.normalize(v1)?;
        let norm2 = self.normalizer.normalize(v2)?;

        let norm_str1 = norm1.to_string();
        let norm_str2 = norm2.to_string();

        if norm_str1 == norm_str2 {
            Ok(EquivalenceResult::new(EquivalenceLevel::NormalizedMatch)
                .with_normalized(norm_str1, norm_str2)
                .with_note("Variants normalize to the same form"))
        } else {
            // Check if they might still be equivalent with version differences
            if let Some(result) =
                self.check_accession_version_difference(&norm1, &norm2, &norm_str1, &norm_str2)
            {
                return Ok(result.with_note("Equivalent after normalization, different versions"));
            }

            Ok(EquivalenceResult::new(EquivalenceLevel::NotEquivalent)
                .with_normalized(norm_str1, norm_str2)
                .with_note("Variants do not normalize to the same form"))
        }
    }

    /// Check if two variants represent the same change but on different accession versions.
    fn check_accession_version_difference(
        &self,
        v1: &HgvsVariant,
        v2: &HgvsVariant,
        str1: &str,
        str2: &str,
    ) -> Option<EquivalenceResult> {
        // Get accessions (handle allele and null/unknown cases gracefully)
        let acc1 = match v1 {
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => return None,
            HgvsVariant::Allele(_) => return None, // Complex case, skip for now
            _ => v1.accession()?,
        };
        let acc2 = match v2 {
            HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => return None,
            HgvsVariant::Allele(_) => return None,
            _ => v2.accession()?,
        };

        // Check if base accessions (without version) match but versions differ
        if acc1.prefix == acc2.prefix && acc1.number == acc2.number && acc1.version != acc2.version
        {
            // Extract the part after the accession to compare the variant itself
            let variant_part1 = extract_variant_part(str1);
            let variant_part2 = extract_variant_part(str2);

            if variant_part1 == variant_part2 {
                return Some(
                    EquivalenceResult::new(EquivalenceLevel::AccessionVersionDifference)
                        .with_normalized(str1.to_string(), str2.to_string())
                        .with_note(format!(
                            "Same variant on different versions: {} vs {}",
                            acc1.version
                                .map(|v| v.to_string())
                                .unwrap_or_else(|| "no version".to_string()),
                            acc2.version
                                .map(|v| v.to_string())
                                .unwrap_or_else(|| "no version".to_string())
                        )),
                );
            }
        }

        None
    }

    /// Check if multiple variants are all equivalent to each other.
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of variants to compare
    ///
    /// # Returns
    ///
    /// * `true` if all variants are equivalent
    /// * `false` if any pair is not equivalent
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::{parse_hgvs, MockProvider};
    /// use ferro_hgvs::equivalence::EquivalenceChecker;
    ///
    /// let provider = MockProvider::with_test_data();
    /// let checker = EquivalenceChecker::new(provider);
    ///
    /// let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
    /// let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
    /// let v3 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
    ///
    /// assert!(checker.all_equivalent(&[v1, v2, v3]).unwrap());
    /// ```
    pub fn all_equivalent(&self, variants: &[HgvsVariant]) -> Result<bool, FerroError> {
        if variants.len() < 2 {
            return Ok(true);
        }

        let first = &variants[0];
        for variant in &variants[1..] {
            let result = self.check(first, variant)?;
            if !result.is_equivalent() {
                return Ok(false);
            }
        }

        Ok(true)
    }

    /// Group variants by equivalence.
    ///
    /// Returns groups of variants that are equivalent to each other.
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of variants to group
    ///
    /// # Returns
    ///
    /// Vector of variant groups, where variants in each group are equivalent.
    pub fn group_by_equivalence(
        &self,
        variants: &[HgvsVariant],
    ) -> Result<Vec<Vec<HgvsVariant>>, FerroError> {
        let mut groups: Vec<Vec<HgvsVariant>> = Vec::new();

        for variant in variants {
            let mut found_group = false;

            for group in &mut groups {
                if !group.is_empty() {
                    let result = self.check(&group[0], variant)?;
                    if result.is_equivalent() {
                        group.push(variant.clone());
                        found_group = true;
                        break;
                    }
                }
            }

            if !found_group {
                groups.push(vec![variant.clone()]);
            }
        }

        Ok(groups)
    }
}

/// Extract the variant part from an HGVS string (everything after the colon).
fn extract_variant_part(hgvs: &str) -> &str {
    if let Some(pos) = hgvs.find(':') {
        &hgvs[pos..]
    } else {
        hgvs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;
    use crate::reference::MockProvider;

    fn checker() -> EquivalenceChecker<MockProvider> {
        EquivalenceChecker::new(MockProvider::with_test_data())
    }

    #[test]
    fn test_identical_variants() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::Identical);
        assert!(result.is_equivalent());
    }

    #[test]
    fn test_different_variants() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.20A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
        assert!(!result.is_equivalent());
    }

    #[test]
    fn test_equivalence_level_is_equivalent() {
        assert!(EquivalenceLevel::Identical.is_equivalent());
        assert!(EquivalenceLevel::NormalizedMatch.is_equivalent());
        assert!(EquivalenceLevel::AccessionVersionDifference.is_equivalent());
        assert!(!EquivalenceLevel::NotEquivalent.is_equivalent());
    }

    #[test]
    fn test_equivalence_level_description() {
        assert!(!EquivalenceLevel::Identical.description().is_empty());
        assert!(!EquivalenceLevel::NormalizedMatch.description().is_empty());
        assert!(!EquivalenceLevel::AccessionVersionDifference
            .description()
            .is_empty());
        assert!(!EquivalenceLevel::NotEquivalent.description().is_empty());
    }

    #[test]
    fn test_equivalence_level_display() {
        let level = EquivalenceLevel::Identical;
        assert_eq!(level.to_string(), level.description());
    }

    #[test]
    fn test_all_equivalent_empty() {
        let checker = checker();
        assert!(checker.all_equivalent(&[]).unwrap());
    }

    #[test]
    fn test_all_equivalent_single() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        assert!(checker.all_equivalent(&[v1]).unwrap());
    }

    #[test]
    fn test_all_equivalent_same() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v3 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        assert!(checker.all_equivalent(&[v1, v2, v3]).unwrap());
    }

    #[test]
    fn test_all_equivalent_different() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.20A>G").unwrap();
        assert!(!checker.all_equivalent(&[v1, v2]).unwrap());
    }

    #[test]
    fn test_group_by_equivalence() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v3 = parse_hgvs("NM_000088.3:c.20A>G").unwrap();

        let groups = checker.group_by_equivalence(&[v1, v2, v3]).unwrap();
        assert_eq!(groups.len(), 2);
        // First group should have 2 variants (the identical ones)
        assert!(groups.iter().any(|g| g.len() == 2));
        // Second group should have 1 variant
        assert!(groups.iter().any(|g| g.len() == 1));
    }

    #[test]
    fn test_extract_variant_part() {
        assert_eq!(extract_variant_part("NM_000088.3:c.10A>G"), ":c.10A>G");
        assert_eq!(extract_variant_part("no_colon"), "no_colon");
    }

    #[test]
    fn test_equivalence_result_with_note() {
        let result = EquivalenceResult::new(EquivalenceLevel::NormalizedMatch)
            .with_note("Test note 1")
            .with_note("Test note 2");

        assert_eq!(result.notes.len(), 2);
        assert_eq!(result.notes[0], "Test note 1");
        assert_eq!(result.notes[1], "Test note 2");
    }

    #[test]
    fn test_substitution_at_different_positions() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.11A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
    }

    #[test]
    fn test_genomic_variants_identical() {
        let checker = checker();
        let v1 = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        let v2 = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::Identical);
    }

    #[test]
    fn test_genomic_variants_different() {
        let checker = checker();
        let v1 = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        let v2 = parse_hgvs("NC_000001.11:g.12346A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
    }

    #[test]
    fn test_protein_variants_identical() {
        let checker = checker();
        let v1 = parse_hgvs("NP_000079.2:p.Val600Glu").unwrap();
        let v2 = parse_hgvs("NP_000079.2:p.Val600Glu").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::Identical);
    }

    #[test]
    fn test_deletion_variants() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10del").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.10del").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::Identical);
    }

    // =========================================================================
    // P2: Cross-accession version tests
    // =========================================================================

    #[test]
    fn test_accession_version_difference_same_variant() {
        let checker = checker();
        // Same position and change, different accession versions
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.4:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
        assert!(result.is_equivalent());
        assert!(result.notes.iter().any(|n| n.contains("version")));
    }

    #[test]
    fn test_accession_version_difference_detected_before_normalize() {
        let checker = checker();
        // Versions differ, should detect before normalization
        let v1 = parse_hgvs("NC_000001.10:g.12345A>G").unwrap();
        let v2 = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
    }

    #[test]
    fn test_accession_version_difference_with_deletion() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10_12del").unwrap();
        let v2 = parse_hgvs("NM_000088.4:c.10_12del").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
    }

    #[test]
    fn test_accession_version_difference_with_insertion() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10_11insATG").unwrap();
        let v2 = parse_hgvs("NM_000088.4:c.10_11insATG").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
    }

    #[test]
    fn test_accession_version_different_variant_not_equivalent() {
        let checker = checker();
        // Different versions AND different variants = not equivalent
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.4:c.20A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
        assert!(!result.is_equivalent());
    }

    #[test]
    fn test_different_accessions_entirely_not_equivalent() {
        let checker = checker();
        // Different base accessions (not just version)
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_001234.1:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
    }

    #[test]
    fn test_protein_accession_version_difference() {
        let checker = checker();
        let v1 = parse_hgvs("NP_000079.1:p.Val600Glu").unwrap();
        let v2 = parse_hgvs("NP_000079.2:p.Val600Glu").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
    }

    #[test]
    fn test_no_version_vs_versioned() {
        let checker = checker();
        // One without version, one with
        let v1 = parse_hgvs("NM_000088:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        // Should still detect as version difference
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
    }

    #[test]
    fn test_version_difference_result_has_normalized_forms() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.4:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert!(result.normalized_first.is_some());
        assert!(result.normalized_second.is_some());
    }

    #[test]
    fn test_genomic_version_difference_grch37_vs_grch38() {
        let checker = checker();
        // GRCh37 (NC_000001.10) vs GRCh38 (NC_000001.11)
        let v1 = parse_hgvs("NC_000001.10:g.100000A>G").unwrap();
        let v2 = parse_hgvs("NC_000001.11:g.100000A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
        // Note: actual position may differ between builds, but this tests the detection
    }

    #[test]
    fn test_version_difference_note_contains_versions() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.5:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        // Notes should mention the specific versions
        let has_version_note = result
            .notes
            .iter()
            .any(|n| n.contains("3") && n.contains("5"));
        assert!(has_version_note, "Notes should contain version numbers");
    }

    #[test]
    fn test_multiple_variants_with_version_mix() {
        let checker = checker();
        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v3 = parse_hgvs("NM_000088.4:c.10A>G").unwrap();

        // all_equivalent should return true since they're all equivalent
        assert!(checker
            .all_equivalent(&[v1.clone(), v2.clone(), v3.clone()])
            .unwrap());

        // group_by_equivalence should put them all in one group
        let groups = checker.group_by_equivalence(&[v1, v2, v3]).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 3);
    }

    #[test]
    fn test_lrg_accession_version_difference() {
        let checker = checker();
        // LRG transcripts
        let v1 = parse_hgvs("LRG_1t1:c.10A>G").unwrap();
        let v2 = parse_hgvs("LRG_1t2:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        // These are different transcripts (t1 vs t2), not versions
        // Should be NotEquivalent
        assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
    }

    #[test]
    fn test_ensembl_accession_version_difference() {
        let checker = checker();
        let v1 = parse_hgvs("ENST00000123456.1:c.10A>G").unwrap();
        let v2 = parse_hgvs("ENST00000123456.2:c.10A>G").unwrap();

        let result = checker.check(&v1, &v2).unwrap();
        assert_eq!(result.level, EquivalenceLevel::AccessionVersionDifference);
    }

    // =========================================================================
    // P3: Equivalence grouping performance tests
    // =========================================================================

    #[test]
    fn test_grouping_100_identical_variants() {
        let checker = checker();

        // Generate 100 identical variants
        let variants: Vec<HgvsVariant> = (0..100)
            .map(|_| parse_hgvs("NM_000088.3:c.10A>G").unwrap())
            .collect();

        let groups = checker.group_by_equivalence(&variants).unwrap();
        // All should be in one group
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 100);
    }

    #[test]
    fn test_grouping_100_unique_variants() {
        let checker = checker();

        // Generate 100 unique variants (different positions)
        let variants: Vec<HgvsVariant> = (1..=100)
            .map(|i| parse_hgvs(&format!("NM_000088.3:c.{}A>G", i)).unwrap())
            .collect();

        let groups = checker.group_by_equivalence(&variants).unwrap();
        // Each should be in its own group
        assert_eq!(groups.len(), 100);
        for group in &groups {
            assert_eq!(group.len(), 1);
        }
    }

    #[test]
    fn test_grouping_mixed_identical_and_unique() {
        let checker = checker();

        // 50 copies of variant A, 50 copies of variant B
        let mut variants: Vec<HgvsVariant> = Vec::with_capacity(100);

        for _ in 0..50 {
            variants.push(parse_hgvs("NM_000088.3:c.10A>G").unwrap());
        }
        for _ in 0..50 {
            variants.push(parse_hgvs("NM_000088.3:c.20C>T").unwrap());
        }

        let groups = checker.group_by_equivalence(&variants).unwrap();
        assert_eq!(groups.len(), 2);

        // Order may vary, but sizes should be 50 each
        let sizes: Vec<usize> = groups.iter().map(|g| g.len()).collect();
        assert!(sizes.contains(&50));
        assert_eq!(sizes.iter().sum::<usize>(), 100);
    }

    #[test]
    fn test_grouping_interleaved_variants() {
        let checker = checker();

        // Interleave 3 different variants
        let mut variants: Vec<HgvsVariant> = Vec::with_capacity(60);

        for i in 0..20 {
            variants.push(parse_hgvs(&format!("NM_000088.3:c.{}A>G", 10)).unwrap());
            variants.push(parse_hgvs(&format!("NM_000088.3:c.{}A>G", 20)).unwrap());
            variants.push(parse_hgvs(&format!("NM_000088.3:c.{}A>G", 30)).unwrap());
            let _ = i; // silence unused warning
        }

        let groups = checker.group_by_equivalence(&variants).unwrap();
        assert_eq!(groups.len(), 3);

        // Each group should have 20 variants
        for group in &groups {
            assert_eq!(group.len(), 20);
        }
    }

    #[test]
    fn test_grouping_accession_version_differences() {
        let checker = checker();

        // Mix of version 3 and version 4 of same variant
        let mut variants: Vec<HgvsVariant> = Vec::with_capacity(40);

        for _ in 0..20 {
            variants.push(parse_hgvs("NM_000088.3:c.10A>G").unwrap());
        }
        for _ in 0..20 {
            variants.push(parse_hgvs("NM_000088.4:c.10A>G").unwrap());
        }

        let groups = checker.group_by_equivalence(&variants).unwrap();
        // They're considered equivalent (AccessionVersionDifference), so all in one group
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 40);
    }

    #[test]
    fn test_grouping_diverse_variant_types() {
        let checker = checker();

        // Different variant types at different positions
        let variant_strs = [
            "NM_000088.3:c.10A>G",       // substitution
            "NM_000088.3:c.20del",       // deletion
            "NM_000088.3:c.30_31insATG", // insertion
            "NM_000088.3:c.40dup",       // duplication
            "NC_000001.11:g.12345A>G",   // genomic
        ];

        let variants: Vec<HgvsVariant> = variant_strs
            .iter()
            .map(|s| parse_hgvs(s).unwrap())
            .collect();

        let groups = checker.group_by_equivalence(&variants).unwrap();
        // All different, so 5 groups
        assert_eq!(groups.len(), 5);
    }

    #[test]
    fn test_grouping_empty_input() {
        let checker = checker();
        let variants: Vec<HgvsVariant> = vec![];

        let groups = checker.group_by_equivalence(&variants).unwrap();
        assert!(groups.is_empty());
    }

    #[test]
    fn test_grouping_single_variant() {
        let checker = checker();
        let variants = vec![parse_hgvs("NM_000088.3:c.10A>G").unwrap()];

        let groups = checker.group_by_equivalence(&variants).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 1);
    }

    #[test]
    fn test_grouping_performance_200_variants() {
        use std::time::Instant;

        let checker = checker();

        // Generate 200 variants: 10 groups of 20 each
        let mut variants: Vec<HgvsVariant> = Vec::with_capacity(200);

        for pos in (1..=10).map(|x| x * 10) {
            for _ in 0..20 {
                variants.push(parse_hgvs(&format!("NM_000088.3:c.{}A>G", pos)).unwrap());
            }
        }

        let start = Instant::now();
        let groups = checker.group_by_equivalence(&variants).unwrap();
        let duration = start.elapsed();

        assert_eq!(groups.len(), 10);
        for group in &groups {
            assert_eq!(group.len(), 20);
        }

        // Should complete in reasonable time (< 5 seconds)
        assert!(
            duration.as_secs() < 5,
            "Grouping took too long: {:?}",
            duration
        );
    }

    #[test]
    fn test_all_equivalent_performance_100() {
        use std::time::Instant;

        let checker = checker();

        // 100 identical variants
        let variants: Vec<HgvsVariant> = (0..100)
            .map(|_| parse_hgvs("NM_000088.3:c.10A>G").unwrap())
            .collect();

        let start = Instant::now();
        let result = checker.all_equivalent(&variants).unwrap();
        let duration = start.elapsed();

        assert!(result);
        assert!(
            duration.as_secs() < 2,
            "all_equivalent took too long: {:?}",
            duration
        );
    }

    #[test]
    fn test_all_equivalent_early_exit_on_non_equivalent() {
        use std::time::Instant;

        let checker = checker();

        // First 2 are different, rest are same
        let mut variants: Vec<HgvsVariant> = Vec::with_capacity(100);
        variants.push(parse_hgvs("NM_000088.3:c.10A>G").unwrap());
        variants.push(parse_hgvs("NM_000088.3:c.20A>G").unwrap()); // different!

        for _ in 0..98 {
            variants.push(parse_hgvs("NM_000088.3:c.10A>G").unwrap());
        }

        let start = Instant::now();
        let result = checker.all_equivalent(&variants).unwrap();
        let duration = start.elapsed();

        assert!(!result);
        // Should exit early and be very fast
        assert!(
            duration.as_millis() < 1000,
            "all_equivalent should exit early"
        );
    }

    #[test]
    fn test_pairwise_check_performance() {
        use std::time::Instant;

        let checker = checker();

        let v1 = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let v2 = parse_hgvs("NM_000088.3:c.20A>G").unwrap();

        // Run 1000 pairwise checks
        let start = Instant::now();
        for _ in 0..1000 {
            let _ = checker.check(&v1, &v2).unwrap();
        }
        let duration = start.elapsed();

        // 1000 checks should complete in < 2 seconds
        assert!(
            duration.as_secs() < 2,
            "Pairwise checks too slow: {:?}",
            duration
        );
    }

    #[test]
    fn test_grouping_stability_same_results() {
        let checker = checker();

        // Same variants in same order should produce same groups
        let variants: Vec<HgvsVariant> = vec![
            parse_hgvs("NM_000088.3:c.10A>G").unwrap(),
            parse_hgvs("NM_000088.3:c.20A>G").unwrap(),
            parse_hgvs("NM_000088.3:c.10A>G").unwrap(),
            parse_hgvs("NM_000088.3:c.30A>G").unwrap(),
            parse_hgvs("NM_000088.3:c.20A>G").unwrap(),
        ];

        let groups1 = checker.group_by_equivalence(&variants).unwrap();
        let groups2 = checker.group_by_equivalence(&variants).unwrap();

        // Same number of groups
        assert_eq!(groups1.len(), groups2.len());

        // Groups have same sizes
        let mut sizes1: Vec<_> = groups1.iter().map(|g| g.len()).collect();
        let mut sizes2: Vec<_> = groups2.iter().map(|g| g.len()).collect();
        sizes1.sort();
        sizes2.sort();
        assert_eq!(sizes1, sizes2);
    }

    #[test]
    fn test_grouping_preserves_variant_order_in_groups() {
        let checker = checker();

        // Variants added in specific order
        let variants: Vec<HgvsVariant> = vec![
            parse_hgvs("NM_000088.3:c.10A>G").unwrap(), // Group A, first
            parse_hgvs("NM_000088.3:c.10A>G").unwrap(), // Group A, second
            parse_hgvs("NM_000088.3:c.10A>G").unwrap(), // Group A, third
        ];

        let groups = checker.group_by_equivalence(&variants).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 3);

        // All variants in the group should be identical
        for v in &groups[0] {
            assert_eq!(v.to_string(), "NM_000088.3:c.10A>G");
        }
    }
}
