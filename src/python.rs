//! Python bindings for ferro-hgvs using PyO3
//!
//! This module provides Python bindings for the HGVS parser and normalizer.
//! Enable with the `python` feature flag.

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::path::Path;

use crate::backtranslate::{Backtranslator, CodonChange, CodonTable};
use crate::batch::{BatchProcessor, BatchProgress, BatchResult};
use crate::convert::CoordinateMapper;
use crate::coords::{OneBasedPos, ZeroBasedPos};
use crate::effect::{Consequence, EffectPredictor, Impact, ProteinEffect};
use crate::equivalence::{EquivalenceChecker, EquivalenceLevel, EquivalenceResult};
use crate::error_handling::{CorrectionWarning, ErrorConfig, ErrorMode, ErrorOverride, ErrorType};
use crate::hgvs::location::{AminoAcid, CdsPos, TxPos};
use crate::hgvs::variant::HgvsVariant;
use crate::mave::{is_mave_short_form, parse_mave_hgvs, MaveContext};
use crate::prepare::{check_references, prepare_references, PrepareConfig, ReferenceManifest};
use crate::python_helpers::{
    get_variant_edit_type, get_variant_reference, parse_direction, variant_type_str,
};
use crate::reference::provider::ReferenceProvider;
use crate::reference::transcript::{GenomeBuild, Strand};
use crate::reference::MockProvider;
use crate::rsid::{format_rsid, parse_rsid as rust_parse_rsid, InMemoryRsIdLookup, RsIdResult};
use crate::spdi::{hgvs_to_spdi_simple, parse_spdi as rust_parse_spdi, spdi_to_hgvs, SpdiVariant};
use crate::vcf::{vcf_to_genomic_hgvs as rust_vcf_to_hgvs, VcfRecord};
use crate::{parse_hgvs, NormalizeConfig, Normalizer};

/// Parse an HGVS variant string
///
/// Args:
///     hgvs_string: The HGVS variant description to parse
///
/// Returns:
///     A PyHgvsVariant object representing the parsed variant
///
/// Raises:
///     ValueError: If the HGVS string cannot be parsed
#[pyfunction]
fn parse(hgvs_string: &str) -> PyResult<PyHgvsVariant> {
    match parse_hgvs(hgvs_string) {
        Ok(variant) => Ok(PyHgvsVariant { inner: variant }),
        Err(e) => Err(PyValueError::new_err(format!("Parse error: {}", e))),
    }
}

/// Normalize an HGVS variant string
///
/// **WARNING**: This function uses built-in test data for reference sequences.
/// For production use, create a Normalizer instance with your reference data.
///
/// Args:
///     hgvs_string: The HGVS variant description to normalize
///     direction: Shuffle direction - "3prime" (default) or "5prime"
///
/// Returns:
///     The normalized HGVS string
///
/// Raises:
///     ValueError: If the HGVS string cannot be parsed
///     RuntimeError: If normalization fails
#[pyfunction]
#[pyo3(signature = (hgvs_string, direction="3prime"))]
fn normalize(hgvs_string: &str, direction: &str) -> PyResult<String> {
    let variant = parse_hgvs(hgvs_string)
        .map_err(|e| PyValueError::new_err(format!("Parse error: {}", e)))?;

    let config = NormalizeConfig::default().with_direction(parse_direction(direction));
    // Note: Uses built-in test data. For production use, create a Normalizer with reference_json.
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::with_config(provider, config);

    let normalized = normalizer
        .normalize(&variant)
        .map_err(|e| PyRuntimeError::new_err(format!("Normalization error: {}", e)))?;

    Ok(normalized.to_string())
}

/// Python wrapper for HgvsVariant
#[pyclass(name = "HgvsVariant")]
#[derive(Clone)]
pub struct PyHgvsVariant {
    pub(crate) inner: HgvsVariant,
}

#[pymethods]
impl PyHgvsVariant {
    /// Get the string representation of the variant
    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    /// Get the repr of the variant
    fn __repr__(&self) -> String {
        format!("HgvsVariant('{}')", self.inner)
    }

    /// Equality comparison
    fn __eq__(&self, other: &Self) -> bool {
        self.inner.to_string() == other.inner.to_string()
    }

    /// Hash for use in sets and dict keys
    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.inner.to_string().hash(&mut hasher);
        hasher.finish()
    }

    /// Get the variant type as a string
    #[getter]
    fn variant_type(&self) -> &'static str {
        variant_type_str(&self.inner)
    }

    /// Get the reference (accession) of the variant
    #[getter]
    fn reference(&self) -> PyResult<String> {
        get_variant_reference(&self.inner).map_err(PyValueError::new_err)
    }

    /// Get the edit type as a string
    #[getter]
    fn edit_type(&self) -> &'static str {
        get_variant_edit_type(&self.inner)
    }

    /// Check if this is a genomic variant (g. prefix)
    fn is_genomic(&self) -> bool {
        matches!(self.inner, HgvsVariant::Genome(_))
    }

    /// Check if this is a coding variant (c. prefix)
    fn is_coding(&self) -> bool {
        matches!(self.inner, HgvsVariant::Cds(_))
    }

    /// Check if this is a non-coding variant (n. prefix)
    fn is_noncoding(&self) -> bool {
        matches!(self.inner, HgvsVariant::Tx(_))
    }

    /// Check if this is a protein variant (p. prefix)
    fn is_protein(&self) -> bool {
        matches!(self.inner, HgvsVariant::Protein(_))
    }

    /// Check if this is an RNA variant (r. prefix)
    fn is_rna(&self) -> bool {
        matches!(self.inner, HgvsVariant::Rna(_))
    }

    /// Check if this is a mitochondrial variant (m. prefix)
    fn is_mitochondrial(&self) -> bool {
        matches!(self.inner, HgvsVariant::Mt(_))
    }

    /// Check if this is a substitution
    fn is_substitution(&self) -> bool {
        get_variant_edit_type(&self.inner) == "substitution"
    }

    /// Check if this is a deletion
    fn is_deletion(&self) -> bool {
        get_variant_edit_type(&self.inner) == "deletion"
    }

    /// Check if this is an insertion
    fn is_insertion(&self) -> bool {
        get_variant_edit_type(&self.inner) == "insertion"
    }

    /// Check if this is a duplication
    fn is_duplication(&self) -> bool {
        get_variant_edit_type(&self.inner) == "duplication"
    }

    /// Check if this is a deletion-insertion (delins)
    fn is_delins(&self) -> bool {
        get_variant_edit_type(&self.inner) == "delins"
    }

    /// Normalize this variant
    ///
    /// **WARNING**: This method uses built-in test data for reference sequences.
    /// For production use, create a Normalizer instance with your reference data.
    ///
    /// Args:
    ///     direction: Shuffle direction - "3prime" (default) or "5prime"
    ///
    /// Returns:
    ///     A new PyHgvsVariant representing the normalized variant
    #[pyo3(signature = (direction="3prime"))]
    fn normalize(&self, direction: &str) -> PyResult<PyHgvsVariant> {
        let config = NormalizeConfig::default().with_direction(parse_direction(direction));
        // Note: Uses built-in test data. For production use, create a Normalizer with reference_json.
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::with_config(provider, config);

        let normalized = normalizer
            .normalize(&self.inner)
            .map_err(|e| PyRuntimeError::new_err(format!("Normalization error: {}", e)))?;

        Ok(PyHgvsVariant { inner: normalized })
    }

    /// Convert to a dictionary representation
    fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        let dict = PyDict::new(py);
        dict.set_item("string", self.inner.to_string())?;
        dict.set_item("variant_type", self.variant_type())?;
        if let Ok(ref_str) = self.reference() {
            dict.set_item("reference", ref_str)?;
        }
        dict.set_item("edit_type", self.edit_type())?;

        Ok(dict.into())
    }

    /// Convert to JSON string
    fn to_json(&self) -> PyResult<String> {
        serde_json::to_string(&self.inner)
            .map_err(|e| PyRuntimeError::new_err(format!("JSON serialization failed: {}", e)))
    }
}

/// HGVS variant normalizer using a reference provider
#[pyclass(name = "Normalizer")]
pub struct PyNormalizer {
    provider: MockProvider,
    config: NormalizeConfig,
}

#[pymethods]
impl PyNormalizer {
    /// Create a new normalizer
    ///
    /// Args:
    ///     reference_json: Path to a transcripts.json file for reference data.
    ///         If not provided, uses built-in test data (limited coverage,
    ///         **not suitable for production use**).
    ///     direction: Shuffle direction - "3prime" (default) or "5prime"
    #[new]
    #[pyo3(signature = (reference_json=None, direction="3prime"))]
    fn new(reference_json: Option<&str>, direction: &str) -> PyResult<Self> {
        let provider = match reference_json {
            Some(path) => MockProvider::from_json(Path::new(path))
                .map_err(|e| PyRuntimeError::new_err(format!("Failed to load reference: {}", e)))?,
            None => {
                // Using test data - this is explicitly documented in the API
                MockProvider::with_test_data()
            }
        };

        let config = NormalizeConfig::default().with_direction(parse_direction(direction));

        Ok(Self { provider, config })
    }

    /// Parse an HGVS string
    fn parse(&self, hgvs_string: &str) -> PyResult<PyHgvsVariant> {
        parse(hgvs_string)
    }

    /// Normalize an HGVS variant
    fn normalize_variant(&self, variant: &PyHgvsVariant) -> PyResult<PyHgvsVariant> {
        let normalizer = Normalizer::with_config(self.provider.clone(), self.config.clone());

        let normalized = normalizer
            .normalize(&variant.inner)
            .map_err(|e| PyRuntimeError::new_err(format!("Normalization error: {}", e)))?;

        Ok(PyHgvsVariant { inner: normalized })
    }

    /// Parse and normalize an HGVS string
    fn normalize(&self, hgvs_string: &str) -> PyResult<String> {
        let variant = parse_hgvs(hgvs_string)
            .map_err(|e| PyValueError::new_err(format!("Parse error: {}", e)))?;

        let normalizer = Normalizer::with_config(self.provider.clone(), self.config.clone());

        let normalized = normalizer
            .normalize(&variant)
            .map_err(|e| PyRuntimeError::new_err(format!("Normalization error: {}", e)))?;

        Ok(normalized.to_string())
    }
}

// ============================================================================
// SPDI Module
// ============================================================================

/// Parse an SPDI variant string
///
/// Args:
///     spdi_string: The SPDI variant description to parse (e.g., "NC_000001.11:12344:A:G")
///
/// Returns:
///     A PySpdiVariant object representing the parsed variant
///
/// Raises:
///     ValueError: If the SPDI string cannot be parsed
#[pyfunction]
fn parse_spdi(spdi_string: &str) -> PyResult<PySpdiVariant> {
    match rust_parse_spdi(spdi_string) {
        Ok(spdi) => Ok(PySpdiVariant { inner: spdi }),
        Err(e) => Err(PyValueError::new_err(format!("Parse error: {}", e))),
    }
}

/// Convert an HGVS variant to SPDI format (simple conversion without reference lookup)
///
/// Args:
///     variant: An HgvsVariant object
///
/// Returns:
///     A PySpdiVariant object
///
/// Raises:
///     ValueError: If the conversion fails
#[pyfunction]
fn hgvs_to_spdi(variant: &PyHgvsVariant) -> PyResult<PySpdiVariant> {
    match hgvs_to_spdi_simple(&variant.inner) {
        Ok(spdi) => Ok(PySpdiVariant { inner: spdi }),
        Err(e) => Err(PyValueError::new_err(format!("Conversion error: {}", e))),
    }
}

/// Convert an SPDI variant to HGVS format
///
/// Args:
///     spdi: A PySpdiVariant object
///
/// Returns:
///     An HgvsVariant object
///
/// Raises:
///     ValueError: If the conversion fails
#[pyfunction]
fn spdi_to_hgvs_variant(spdi: &PySpdiVariant) -> PyResult<PyHgvsVariant> {
    match spdi_to_hgvs(&spdi.inner) {
        Ok(variant) => Ok(PyHgvsVariant { inner: variant }),
        Err(e) => Err(PyValueError::new_err(format!("Conversion error: {}", e))),
    }
}

/// Python wrapper for SpdiVariant
#[pyclass(name = "SpdiVariant")]
#[derive(Clone)]
pub struct PySpdiVariant {
    inner: SpdiVariant,
}

#[pymethods]
impl PySpdiVariant {
    /// Create a new SPDI variant
    ///
    /// Args:
    ///     sequence: Reference sequence identifier (e.g., "NC_000001.11")
    ///     position: 0-based interbase position
    ///     deletion: Deleted sequence (can be empty for insertions)
    ///     insertion: Inserted sequence (can be empty for deletions)
    #[new]
    fn new(sequence: &str, position: u64, deletion: &str, insertion: &str) -> Self {
        Self {
            inner: SpdiVariant::new(sequence, position, deletion, insertion),
        }
    }

    /// Get the string representation of the variant
    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    /// Get the repr of the variant
    fn __repr__(&self) -> String {
        format!("SpdiVariant('{}')", self.inner)
    }

    /// Equality comparison
    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    /// Hash for use in sets and dict keys
    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }

    /// Reference sequence identifier
    #[getter]
    fn sequence(&self) -> &str {
        &self.inner.sequence
    }

    /// 0-based interbase position
    #[getter]
    fn position(&self) -> u64 {
        self.inner.position
    }

    /// Deleted sequence
    #[getter]
    fn deletion(&self) -> &str {
        &self.inner.deletion
    }

    /// Inserted sequence
    #[getter]
    fn insertion(&self) -> &str {
        &self.inner.insertion
    }

    /// Check if this is a substitution (SNV or MNV)
    fn is_substitution(&self) -> bool {
        self.inner.is_substitution()
    }

    /// Check if this is a single nucleotide variant
    fn is_snv(&self) -> bool {
        self.inner.is_substitution() && self.inner.deletion.len() == 1
    }

    /// Check if this is a pure deletion
    fn is_deletion(&self) -> bool {
        self.inner.is_deletion()
    }

    /// Check if this is a pure insertion
    fn is_insertion(&self) -> bool {
        self.inner.is_insertion()
    }

    /// Check if this is a deletion-insertion (delins)
    fn is_delins(&self) -> bool {
        self.inner.is_delins()
    }

    /// Check if this represents an identity (no change)
    fn is_identity(&self) -> bool {
        self.inner.is_identity()
    }

    /// Get the variant type as a string
    fn variant_type(&self) -> &'static str {
        self.inner.variant_type()
    }

    /// Convert 0-based SPDI position to 1-based HGVS position
    fn to_one_based_position(&self) -> u64 {
        self.inner.to_one_based_position()
    }

    /// Convert to a dictionary representation
    fn to_dict(&self, py: Python<'_>) -> PyResult<PyObject> {
        let dict = PyDict::new(py);
        dict.set_item("sequence", &self.inner.sequence)?;
        dict.set_item("position", self.inner.position)?;
        dict.set_item("deletion", &self.inner.deletion)?;
        dict.set_item("insertion", &self.inner.insertion)?;
        dict.set_item("variant_type", self.inner.variant_type())?;
        Ok(dict.into())
    }
}

// ============================================================================
// Coordinates Module
// ============================================================================

/// Python wrapper for ZeroBasedPos
#[pyclass(name = "ZeroBasedPos")]
#[derive(Clone)]
pub struct PyZeroBasedPos {
    inner: ZeroBasedPos,
}

#[pymethods]
impl PyZeroBasedPos {
    /// Create a new 0-based position
    #[new]
    fn new(pos: u64) -> Self {
        Self {
            inner: ZeroBasedPos::new(pos),
        }
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        format!("ZeroBasedPos({})", self.inner.value())
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }

    fn __lt__(&self, other: &Self) -> bool {
        self.inner < other.inner
    }

    fn __le__(&self, other: &Self) -> bool {
        self.inner <= other.inner
    }

    fn __gt__(&self, other: &Self) -> bool {
        self.inner > other.inner
    }

    fn __ge__(&self, other: &Self) -> bool {
        self.inner >= other.inner
    }

    /// Get the raw value
    #[getter]
    fn value(&self) -> u64 {
        self.inner.value()
    }

    /// Convert to 1-based position
    fn to_one_based(&self) -> PyOneBasedPos {
        PyOneBasedPos {
            inner: self.inner.to_one_based(),
        }
    }

    /// Use as array index (returns int)
    fn as_index(&self) -> usize {
        self.inner.as_index()
    }
}

/// Python wrapper for OneBasedPos
#[pyclass(name = "OneBasedPos")]
#[derive(Clone)]
pub struct PyOneBasedPos {
    inner: OneBasedPos,
}

#[pymethods]
impl PyOneBasedPos {
    /// Create a new 1-based position
    ///
    /// Raises:
    ///     ValueError: If pos is 0 (invalid for 1-based coordinates)
    #[new]
    fn new(pos: u64) -> PyResult<Self> {
        if pos == 0 {
            return Err(PyValueError::new_err("1-based position cannot be 0"));
        }
        Ok(Self {
            inner: OneBasedPos::new(pos),
        })
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        format!("OneBasedPos({})", self.inner.value())
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }

    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.inner.hash(&mut hasher);
        hasher.finish()
    }

    fn __lt__(&self, other: &Self) -> bool {
        self.inner < other.inner
    }

    fn __le__(&self, other: &Self) -> bool {
        self.inner <= other.inner
    }

    fn __gt__(&self, other: &Self) -> bool {
        self.inner > other.inner
    }

    fn __ge__(&self, other: &Self) -> bool {
        self.inner >= other.inner
    }

    /// Get the raw value
    #[getter]
    fn value(&self) -> u64 {
        self.inner.value()
    }

    /// Convert to 0-based position
    fn to_zero_based(&self) -> PyZeroBasedPos {
        PyZeroBasedPos {
            inner: self.inner.to_zero_based(),
        }
    }
}

/// Convert a 1-based HGVS position to a 0-based array index
#[pyfunction]
fn hgvs_pos_to_index(pos: u64) -> usize {
    crate::coords::hgvs_pos_to_index(pos)
}

/// Convert a 0-based array index to a 1-based HGVS position
#[pyfunction]
fn index_to_hgvs_pos(idx: usize) -> u64 {
    crate::coords::index_to_hgvs_pos(idx)
}

// ============================================================================
// Equivalence Module
// ============================================================================

/// Equivalence level between two variants
#[pyclass(name = "EquivalenceLevel", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum PyEquivalenceLevel {
    /// Exactly identical (same string representation)
    Identical = 0,
    /// Equivalent after normalization
    NormalizedMatch = 1,
    /// Same variant but different accession versions
    AccessionVersionDifference = 2,
    /// Not equivalent - represent different changes
    NotEquivalent = 3,
}

impl From<EquivalenceLevel> for PyEquivalenceLevel {
    fn from(level: EquivalenceLevel) -> Self {
        match level {
            EquivalenceLevel::Identical => PyEquivalenceLevel::Identical,
            EquivalenceLevel::NormalizedMatch => PyEquivalenceLevel::NormalizedMatch,
            EquivalenceLevel::AccessionVersionDifference => {
                PyEquivalenceLevel::AccessionVersionDifference
            }
            EquivalenceLevel::NotEquivalent => PyEquivalenceLevel::NotEquivalent,
        }
    }
}

#[pymethods]
impl PyEquivalenceLevel {
    fn __str__(&self) -> &'static str {
        match self {
            PyEquivalenceLevel::Identical => "Identical",
            PyEquivalenceLevel::NormalizedMatch => "NormalizedMatch",
            PyEquivalenceLevel::AccessionVersionDifference => "AccessionVersionDifference",
            PyEquivalenceLevel::NotEquivalent => "NotEquivalent",
        }
    }

    /// Returns true if the variants are considered equivalent
    fn is_equivalent(&self) -> bool {
        !matches!(self, PyEquivalenceLevel::NotEquivalent)
    }

    /// Get a human-readable description
    fn description(&self) -> &'static str {
        match self {
            PyEquivalenceLevel::Identical => "Identical representation",
            PyEquivalenceLevel::NormalizedMatch => "Equivalent after normalization",
            PyEquivalenceLevel::AccessionVersionDifference => {
                "Same variant, different accession versions"
            }
            PyEquivalenceLevel::NotEquivalent => "Not equivalent",
        }
    }
}

/// Result of an equivalence check
#[pyclass(name = "EquivalenceResult")]
#[derive(Clone)]
pub struct PyEquivalenceResult {
    /// The determined equivalence level
    #[pyo3(get)]
    pub level: PyEquivalenceLevel,
    /// The normalized form of the first variant (if normalization was performed)
    #[pyo3(get)]
    pub normalized_first: Option<String>,
    /// The normalized form of the second variant (if normalization was performed)
    #[pyo3(get)]
    pub normalized_second: Option<String>,
    /// Additional notes about the comparison
    #[pyo3(get)]
    pub notes: Vec<String>,
}

#[pymethods]
impl PyEquivalenceResult {
    fn __repr__(&self) -> String {
        format!(
            "EquivalenceResult(level={}, is_equivalent={})",
            self.level.__str__(),
            self.is_equivalent()
        )
    }

    /// Returns true if the variants are considered equivalent
    fn is_equivalent(&self) -> bool {
        self.level.is_equivalent()
    }
}

impl From<EquivalenceResult> for PyEquivalenceResult {
    fn from(result: EquivalenceResult) -> Self {
        Self {
            level: result.level.into(),
            normalized_first: result.normalized_first,
            normalized_second: result.normalized_second,
            notes: result.notes,
        }
    }
}

/// Equivalence checker for comparing HGVS variants
#[pyclass(name = "EquivalenceChecker")]
pub struct PyEquivalenceChecker {
    checker: EquivalenceChecker<MockProvider>,
}

#[pymethods]
impl PyEquivalenceChecker {
    /// Create a new equivalence checker
    ///
    /// Args:
    ///     reference_json: Optional path to a transcripts.json file for reference data.
    ///         If not provided, uses built-in test data.
    #[new]
    #[pyo3(signature = (reference_json=None))]
    fn new(reference_json: Option<&str>) -> PyResult<Self> {
        let provider = match reference_json {
            Some(path) => MockProvider::from_json(Path::new(path))
                .map_err(|e| PyRuntimeError::new_err(format!("Failed to load reference: {}", e)))?,
            None => MockProvider::with_test_data(),
        };

        Ok(Self {
            checker: EquivalenceChecker::new(provider),
        })
    }

    /// Check if two variants are equivalent
    ///
    /// Args:
    ///     v1: First variant
    ///     v2: Second variant
    ///
    /// Returns:
    ///     EquivalenceResult with details about the comparison
    fn check(&self, v1: &PyHgvsVariant, v2: &PyHgvsVariant) -> PyResult<PyEquivalenceResult> {
        self.checker
            .check(&v1.inner, &v2.inner)
            .map(|r| r.into())
            .map_err(|e| PyRuntimeError::new_err(format!("Equivalence check failed: {}", e)))
    }

    /// Check if multiple variants are all equivalent to each other
    ///
    /// Args:
    ///     variants: List of variants to compare
    ///
    /// Returns:
    ///     True if all variants are equivalent
    fn all_equivalent(&self, variants: Vec<PyRef<PyHgvsVariant>>) -> PyResult<bool> {
        let rust_variants: Vec<HgvsVariant> = variants.iter().map(|v| v.inner.clone()).collect();
        self.checker
            .all_equivalent(&rust_variants)
            .map_err(|e| PyRuntimeError::new_err(format!("Equivalence check failed: {}", e)))
    }
}

// ============================================================================
// Effect Prediction Module
// ============================================================================

/// Sequence Ontology consequence term
#[pyclass(name = "Consequence", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum PyConsequence {
    TranscriptAblation = 0,
    SpliceAcceptorVariant = 1,
    SpliceDonorVariant = 2,
    StopGained = 3,
    FrameshiftVariant = 4,
    StopLost = 5,
    StartLost = 6,
    MissenseVariant = 7,
    InframeInsertion = 8,
    InframeDeletion = 9,
    ProteinAlteringVariant = 10,
    SpliceRegionVariant = 11,
    SynonymousVariant = 12,
    StartRetainedVariant = 13,
    StopRetainedVariant = 14,
    FivePrimeUtrVariant = 15,
    ThreePrimeUtrVariant = 16,
    IntronVariant = 17,
    CodingSequenceVariant = 18,
}

impl From<Consequence> for PyConsequence {
    fn from(c: Consequence) -> Self {
        match c {
            Consequence::TranscriptAblation => PyConsequence::TranscriptAblation,
            Consequence::SpliceAcceptorVariant => PyConsequence::SpliceAcceptorVariant,
            Consequence::SpliceDonorVariant => PyConsequence::SpliceDonorVariant,
            Consequence::StopGained => PyConsequence::StopGained,
            Consequence::FrameshiftVariant => PyConsequence::FrameshiftVariant,
            Consequence::StopLost => PyConsequence::StopLost,
            Consequence::StartLost => PyConsequence::StartLost,
            Consequence::MissenseVariant => PyConsequence::MissenseVariant,
            Consequence::InframeInsertion => PyConsequence::InframeInsertion,
            Consequence::InframeDeletion => PyConsequence::InframeDeletion,
            Consequence::ProteinAlteringVariant => PyConsequence::ProteinAlteringVariant,
            Consequence::SpliceRegionVariant => PyConsequence::SpliceRegionVariant,
            Consequence::SynonymousVariant => PyConsequence::SynonymousVariant,
            Consequence::StartRetainedVariant => PyConsequence::StartRetainedVariant,
            Consequence::StopRetainedVariant => PyConsequence::StopRetainedVariant,
            Consequence::FivePrimeUtrVariant => PyConsequence::FivePrimeUtrVariant,
            Consequence::ThreePrimeUtrVariant => PyConsequence::ThreePrimeUtrVariant,
            Consequence::IntronVariant => PyConsequence::IntronVariant,
            Consequence::CodingSequenceVariant => PyConsequence::CodingSequenceVariant,
        }
    }
}

impl From<PyConsequence> for Consequence {
    fn from(c: PyConsequence) -> Self {
        match c {
            PyConsequence::TranscriptAblation => Consequence::TranscriptAblation,
            PyConsequence::SpliceAcceptorVariant => Consequence::SpliceAcceptorVariant,
            PyConsequence::SpliceDonorVariant => Consequence::SpliceDonorVariant,
            PyConsequence::StopGained => Consequence::StopGained,
            PyConsequence::FrameshiftVariant => Consequence::FrameshiftVariant,
            PyConsequence::StopLost => Consequence::StopLost,
            PyConsequence::StartLost => Consequence::StartLost,
            PyConsequence::MissenseVariant => Consequence::MissenseVariant,
            PyConsequence::InframeInsertion => Consequence::InframeInsertion,
            PyConsequence::InframeDeletion => Consequence::InframeDeletion,
            PyConsequence::ProteinAlteringVariant => Consequence::ProteinAlteringVariant,
            PyConsequence::SpliceRegionVariant => Consequence::SpliceRegionVariant,
            PyConsequence::SynonymousVariant => Consequence::SynonymousVariant,
            PyConsequence::StartRetainedVariant => Consequence::StartRetainedVariant,
            PyConsequence::StopRetainedVariant => Consequence::StopRetainedVariant,
            PyConsequence::FivePrimeUtrVariant => Consequence::FivePrimeUtrVariant,
            PyConsequence::ThreePrimeUtrVariant => Consequence::ThreePrimeUtrVariant,
            PyConsequence::IntronVariant => Consequence::IntronVariant,
            PyConsequence::CodingSequenceVariant => Consequence::CodingSequenceVariant,
        }
    }
}

#[pymethods]
impl PyConsequence {
    fn __str__(&self) -> &'static str {
        Consequence::from(*self).so_term()
    }

    /// Get the Sequence Ontology term
    fn so_term(&self) -> &'static str {
        Consequence::from(*self).so_term()
    }

    /// Get the Sequence Ontology ID
    fn so_id(&self) -> &'static str {
        Consequence::from(*self).so_id()
    }

    /// Get the impact level
    fn impact(&self) -> PyImpact {
        Consequence::from(*self).impact().into()
    }

    /// Get a human-readable description
    fn description(&self) -> &'static str {
        Consequence::from(*self).description()
    }
}

/// Variant impact level (VEP-style)
#[pyclass(name = "Impact", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum PyImpact {
    Modifier = 0,
    Low = 1,
    Moderate = 2,
    High = 3,
}

impl From<Impact> for PyImpact {
    fn from(i: Impact) -> Self {
        match i {
            Impact::Modifier => PyImpact::Modifier,
            Impact::Low => PyImpact::Low,
            Impact::Moderate => PyImpact::Moderate,
            Impact::High => PyImpact::High,
        }
    }
}

#[pymethods]
impl PyImpact {
    fn __str__(&self) -> &'static str {
        match self {
            PyImpact::High => "HIGH",
            PyImpact::Moderate => "MODERATE",
            PyImpact::Low => "LOW",
            PyImpact::Modifier => "MODIFIER",
        }
    }
}

/// Protein effect prediction result
#[pyclass(name = "ProteinEffect")]
#[derive(Clone)]
pub struct PyProteinEffect {
    inner: ProteinEffect,
}

#[pymethods]
impl PyProteinEffect {
    /// Get all applicable consequences
    #[getter]
    fn consequences(&self) -> Vec<PyConsequence> {
        self.inner
            .consequences
            .iter()
            .map(|c| (*c).into())
            .collect()
    }

    /// Get the highest impact level
    #[getter]
    fn impact(&self) -> PyImpact {
        self.inner.impact.into()
    }

    /// Check if this is a high-impact variant
    fn is_high_impact(&self) -> bool {
        self.inner.is_high_impact()
    }

    /// Check if this is a protein-altering variant
    fn is_protein_altering(&self) -> bool {
        self.inner.is_protein_altering()
    }

    fn __repr__(&self) -> String {
        format!(
            "ProteinEffect(impact={}, consequences={:?})",
            self.inner.impact,
            self.inner
                .consequences
                .iter()
                .map(|c| c.so_term())
                .collect::<Vec<_>>()
        )
    }
}

/// Protein effect predictor
#[pyclass(name = "EffectPredictor")]
pub struct PyEffectPredictor {
    inner: EffectPredictor,
}

#[pymethods]
impl PyEffectPredictor {
    /// Create a new effect predictor
    #[new]
    fn new() -> Self {
        Self {
            inner: EffectPredictor::new(),
        }
    }

    /// Classify an amino acid change
    ///
    /// Args:
    ///     ref_aa: Reference amino acid (3-letter code, e.g., "Val")
    ///     alt_aa: Alternate amino acid (3-letter code, e.g., "Glu")
    ///     position: Position in protein (1-based)
    ///
    /// Returns:
    ///     ProteinEffect with consequences and impact
    fn classify_amino_acid_change(
        &self,
        ref_aa: &str,
        alt_aa: &str,
        position: u64,
    ) -> PyResult<PyProteinEffect> {
        let ref_aa = parse_amino_acid(ref_aa)?;
        let alt_aa = parse_amino_acid(alt_aa)?;
        Ok(PyProteinEffect {
            inner: self
                .inner
                .classify_amino_acid_change(&ref_aa, &alt_aa, position),
        })
    }

    /// Classify an indel by frame effect
    ///
    /// Args:
    ///     ref_len: Length of reference sequence
    ///     alt_len: Length of alternate sequence
    ///
    /// Returns:
    ///     ProteinEffect with consequences and impact
    fn classify_indel(&self, ref_len: usize, alt_len: usize) -> PyProteinEffect {
        PyProteinEffect {
            inner: self.inner.classify_indel(ref_len, alt_len),
        }
    }

    /// Classify a splice site variant by distance from splice site
    ///
    /// Args:
    ///     offset: Distance from splice site (negative for acceptor, positive for donor)
    ///
    /// Returns:
    ///     ProteinEffect with consequences and impact
    fn classify_splice_variant(&self, offset: i64) -> PyProteinEffect {
        PyProteinEffect {
            inner: self.inner.classify_splice_variant(offset),
        }
    }

    /// Classify a UTR variant
    ///
    /// Args:
    ///     is_5_prime: True for 5' UTR, False for 3' UTR
    ///
    /// Returns:
    ///     ProteinEffect with consequences and impact
    fn classify_utr_variant(&self, is_5_prime: bool) -> PyProteinEffect {
        PyProteinEffect {
            inner: self.inner.classify_utr_variant(is_5_prime),
        }
    }
}

/// Helper function to parse amino acid from 3-letter code
fn parse_amino_acid(code: &str) -> PyResult<AminoAcid> {
    match code.to_uppercase().as_str() {
        "ALA" | "A" => Ok(AminoAcid::Ala),
        "ARG" | "R" => Ok(AminoAcid::Arg),
        "ASN" | "N" => Ok(AminoAcid::Asn),
        "ASP" | "D" => Ok(AminoAcid::Asp),
        "CYS" | "C" => Ok(AminoAcid::Cys),
        "GLN" | "Q" => Ok(AminoAcid::Gln),
        "GLU" | "E" => Ok(AminoAcid::Glu),
        "GLY" | "G" => Ok(AminoAcid::Gly),
        "HIS" | "H" => Ok(AminoAcid::His),
        "ILE" | "I" => Ok(AminoAcid::Ile),
        "LEU" | "L" => Ok(AminoAcid::Leu),
        "LYS" | "K" => Ok(AminoAcid::Lys),
        "MET" | "M" => Ok(AminoAcid::Met),
        "PHE" | "F" => Ok(AminoAcid::Phe),
        "PRO" | "P" => Ok(AminoAcid::Pro),
        "SER" | "S" => Ok(AminoAcid::Ser),
        "THR" | "T" => Ok(AminoAcid::Thr),
        "TRP" | "W" => Ok(AminoAcid::Trp),
        "TYR" | "Y" => Ok(AminoAcid::Tyr),
        "VAL" | "V" => Ok(AminoAcid::Val),
        "TER" | "*" | "X" => Ok(AminoAcid::Ter),
        "SEC" | "U" => Ok(AminoAcid::Sec),
        _ => Err(PyValueError::new_err(format!(
            "Unknown amino acid: {}",
            code
        ))),
    }
}

// ============================================================================
// MAVE Module
// ============================================================================

/// Context for parsing MAVE-HGVS short-form notation
#[pyclass(name = "MaveContext")]
#[derive(Clone)]
pub struct PyMaveContext {
    inner: MaveContext,
}

#[pymethods]
impl PyMaveContext {
    /// Create a new empty context
    #[new]
    fn new() -> Self {
        Self {
            inner: MaveContext::new(),
        }
    }

    /// Set the protein sequence accession for p. variants
    fn with_protein_accession(&self, accession: &str) -> Self {
        Self {
            inner: self.inner.clone().with_protein_accession(accession),
        }
    }

    /// Set the coding sequence accession for c. variants
    fn with_coding_accession(&self, accession: &str) -> Self {
        Self {
            inner: self.inner.clone().with_coding_accession(accession),
        }
    }

    /// Set the non-coding transcript accession for n. variants
    fn with_noncoding_accession(&self, accession: &str) -> Self {
        Self {
            inner: self.inner.clone().with_noncoding_accession(accession),
        }
    }

    /// Set the genomic sequence accession for g. variants
    fn with_genomic_accession(&self, accession: &str) -> Self {
        Self {
            inner: self.inner.clone().with_genomic_accession(accession),
        }
    }

    /// Set the gene symbol (informational)
    fn with_gene_symbol(&self, symbol: &str) -> Self {
        Self {
            inner: self.inner.clone().with_gene_symbol(symbol),
        }
    }

    /// Check if this context has any accessions defined
    fn has_accessions(&self) -> bool {
        self.inner.has_accessions()
    }

    /// Check if this context can handle a specific coordinate type
    ///
    /// Args:
    ///     coord_type: Single character coordinate type ('p', 'c', 'n', 'g', etc.)
    fn supports_coordinate_type(&self, coord_type: &str) -> PyResult<bool> {
        let c = coord_type
            .chars()
            .next()
            .ok_or_else(|| PyValueError::new_err("coord_type must be a non-empty string"))?;
        Ok(self.inner.supports_coordinate_type(c))
    }

    #[getter]
    fn noncoding_accession(&self) -> Option<String> {
        self.inner.noncoding_accession.clone()
    }

    #[getter]
    fn protein_accession(&self) -> Option<String> {
        self.inner.protein_accession.clone()
    }

    #[getter]
    fn coding_accession(&self) -> Option<String> {
        self.inner.coding_accession.clone()
    }

    #[getter]
    fn genomic_accession(&self) -> Option<String> {
        self.inner.genomic_accession.clone()
    }

    #[getter]
    fn gene_symbol(&self) -> Option<String> {
        self.inner.gene_symbol.clone()
    }
}

/// Parse a MAVE-HGVS variant string with context
///
/// Args:
///     hgvs_string: The MAVE-HGVS variant description (e.g., "p.Glu6Val")
///     context: MaveContext with accession information
///
/// Returns:
///     A HgvsVariant object with the accession filled in from context
///
/// Raises:
///     ValueError: If parsing fails or context doesn't support the coordinate type
#[pyfunction]
fn parse_mave_hgvs_variant(hgvs_string: &str, context: &PyMaveContext) -> PyResult<PyHgvsVariant> {
    match parse_mave_hgvs(hgvs_string, &context.inner) {
        Ok(variant) => Ok(PyHgvsVariant { inner: variant }),
        Err(e) => Err(PyValueError::new_err(format!("Parse error: {}", e))),
    }
}

/// Check if a string is in MAVE short-form notation (no accession)
///
/// Args:
///     hgvs_string: The HGVS string to check
///
/// Returns:
///     True if the string is in short form (e.g., "p.Val600Glu"), False otherwise
#[pyfunction]
fn is_mave_short_form_variant(hgvs_string: &str) -> bool {
    is_mave_short_form(hgvs_string)
}

// ============================================================================
// Batch Processing Module
// ============================================================================

/// Progress information for batch processing
#[pyclass(name = "BatchProgress")]
#[derive(Clone)]
pub struct PyBatchProgress {
    #[pyo3(get)]
    pub processed: usize,
    #[pyo3(get)]
    pub total: usize,
    #[pyo3(get)]
    pub successes: usize,
    #[pyo3(get)]
    pub failures: usize,
}

#[pymethods]
impl PyBatchProgress {
    /// Get the percentage complete
    fn percent(&self) -> f64 {
        if self.total == 0 {
            100.0
        } else {
            (self.processed as f64 / self.total as f64) * 100.0
        }
    }
}

impl From<&BatchProgress> for PyBatchProgress {
    fn from(p: &BatchProgress) -> Self {
        Self {
            processed: p.processed,
            total: p.total,
            successes: p.success,
            failures: p.errors,
        }
    }
}

/// Result of batch processing
#[pyclass(name = "BatchResult")]
pub struct PyBatchResult {
    inner: BatchResult<HgvsVariant>,
}

#[pymethods]
impl PyBatchResult {
    /// Total number of items processed
    fn total(&self) -> usize {
        self.inner.total()
    }

    /// Number of successful items
    fn success_count(&self) -> usize {
        self.inner.success_count()
    }

    /// Number of failed items
    fn error_count(&self) -> usize {
        self.inner.error_count()
    }

    /// Success rate as a percentage
    fn success_rate(&self) -> f64 {
        self.inner.success_rate()
    }

    /// Get successfully parsed/normalized variants
    fn successes(&self) -> Vec<PyHgvsVariant> {
        self.inner
            .results
            .iter()
            .filter_map(|r| match r {
                crate::batch::ItemResult::Ok(v) => Some(PyHgvsVariant { inner: v.clone() }),
                _ => None,
            })
            .collect()
    }

    /// Get errors as (index, error_message) tuples
    fn errors(&self) -> Vec<(usize, String)> {
        self.inner
            .results
            .iter()
            .enumerate()
            .filter_map(|(i, r)| match r {
                crate::batch::ItemResult::Err { error, .. } => Some((i, error.to_string())),
                _ => None,
            })
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "BatchResult(total={}, successes={}, errors={})",
            self.total(),
            self.success_count(),
            self.error_count()
        )
    }
}

/// Batch processor for parsing and normalizing multiple variants
#[pyclass(name = "BatchProcessor")]
pub struct PyBatchProcessor {
    processor: BatchProcessor<MockProvider>,
}

#[pymethods]
impl PyBatchProcessor {
    /// Create a new batch processor
    ///
    /// Args:
    ///     reference_json: Optional path to a transcripts.json file for reference data.
    ///         If not provided, uses built-in test data.
    #[new]
    #[pyo3(signature = (reference_json=None))]
    fn new(reference_json: Option<&str>) -> PyResult<Self> {
        let provider = match reference_json {
            Some(path) => MockProvider::from_json(Path::new(path))
                .map_err(|e| PyRuntimeError::new_err(format!("Failed to load reference: {}", e)))?,
            None => MockProvider::with_test_data(),
        };

        Ok(Self {
            processor: BatchProcessor::new(provider),
        })
    }

    /// Parse multiple HGVS strings
    ///
    /// Args:
    ///     variants: List of HGVS strings to parse
    ///
    /// Returns:
    ///     BatchResult with parsed variants and errors
    fn parse(&self, variants: Vec<String>) -> PyBatchResult {
        PyBatchResult {
            inner: self.processor.parse(&variants),
        }
    }

    /// Parse and normalize multiple HGVS strings
    ///
    /// Args:
    ///     variants: List of HGVS strings to parse and normalize
    ///
    /// Returns:
    ///     BatchResult with normalized variants and errors
    fn parse_and_normalize(&self, variants: Vec<String>) -> PyBatchResult {
        PyBatchResult {
            inner: self.processor.parse_and_normalize(&variants),
        }
    }

    /// Parse multiple HGVS strings with progress callback
    ///
    /// Args:
    ///     variants: List of HGVS strings to parse
    ///     callback: Python callable that receives BatchProgress
    ///
    /// Returns:
    ///     BatchResult with parsed variants and errors
    fn parse_with_progress(
        &self,
        variants: Vec<String>,
        callback: Bound<'_, pyo3::PyAny>,
    ) -> PyResult<PyBatchResult> {
        let result = self.processor.parse_with_progress(&variants, |progress| {
            let py_progress = PyBatchProgress::from(&progress);
            // Call the Python callback, ignoring errors
            let _ = callback.call1((py_progress,));
        });
        Ok(PyBatchResult { inner: result })
    }
}

// ============================================================================
// Error Handling Module
// ============================================================================

/// Error handling mode
#[pyclass(name = "ErrorMode", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum PyErrorMode {
    /// Strict mode - reject all non-standard input
    Strict = 0,
    /// Lenient mode - auto-correct with warnings
    Lenient = 1,
    /// Silent mode - auto-correct without warnings
    Silent = 2,
}

impl From<ErrorMode> for PyErrorMode {
    fn from(mode: ErrorMode) -> Self {
        match mode {
            ErrorMode::Strict => PyErrorMode::Strict,
            ErrorMode::Lenient => PyErrorMode::Lenient,
            ErrorMode::Silent => PyErrorMode::Silent,
        }
    }
}

impl From<PyErrorMode> for ErrorMode {
    fn from(mode: PyErrorMode) -> Self {
        match mode {
            PyErrorMode::Strict => ErrorMode::Strict,
            PyErrorMode::Lenient => ErrorMode::Lenient,
            PyErrorMode::Silent => ErrorMode::Silent,
        }
    }
}

/// Error type for configurable error handling
#[pyclass(name = "ErrorType", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum PyErrorType {
    LowercaseAminoAcid = 0,
    MissingVersion = 1,
    WrongDashCharacter = 2,
    ExtraWhitespace = 3,
    ProteinSubstitutionArrow = 4,
    PositionZero = 5,
    SingleLetterAminoAcid = 6,
    WrongQuoteCharacter = 7,
    LowercaseAccessionPrefix = 8,
    MixedCaseEditType = 9,
    OldSubstitutionSyntax = 10,
    InvalidUnicodeCharacter = 11,
    SwappedPositions = 12,
    TrailingAnnotation = 13,
    MissingCoordinatePrefix = 14,
    OldAlleleFormat = 15,
    RefSeqMismatch = 16,
}

impl From<ErrorType> for PyErrorType {
    fn from(t: ErrorType) -> Self {
        match t {
            ErrorType::LowercaseAminoAcid => PyErrorType::LowercaseAminoAcid,
            ErrorType::MissingVersion => PyErrorType::MissingVersion,
            ErrorType::WrongDashCharacter => PyErrorType::WrongDashCharacter,
            ErrorType::ExtraWhitespace => PyErrorType::ExtraWhitespace,
            ErrorType::ProteinSubstitutionArrow => PyErrorType::ProteinSubstitutionArrow,
            ErrorType::PositionZero => PyErrorType::PositionZero,
            ErrorType::SingleLetterAminoAcid => PyErrorType::SingleLetterAminoAcid,
            ErrorType::WrongQuoteCharacter => PyErrorType::WrongQuoteCharacter,
            ErrorType::LowercaseAccessionPrefix => PyErrorType::LowercaseAccessionPrefix,
            ErrorType::MixedCaseEditType => PyErrorType::MixedCaseEditType,
            ErrorType::OldSubstitutionSyntax => PyErrorType::OldSubstitutionSyntax,
            ErrorType::InvalidUnicodeCharacter => PyErrorType::InvalidUnicodeCharacter,
            ErrorType::SwappedPositions => PyErrorType::SwappedPositions,
            ErrorType::TrailingAnnotation => PyErrorType::TrailingAnnotation,
            ErrorType::MissingCoordinatePrefix => PyErrorType::MissingCoordinatePrefix,
            ErrorType::OldAlleleFormat => PyErrorType::OldAlleleFormat,
            ErrorType::RefSeqMismatch => PyErrorType::RefSeqMismatch,
        }
    }
}

impl From<PyErrorType> for ErrorType {
    fn from(t: PyErrorType) -> Self {
        match t {
            PyErrorType::LowercaseAminoAcid => ErrorType::LowercaseAminoAcid,
            PyErrorType::MissingVersion => ErrorType::MissingVersion,
            PyErrorType::WrongDashCharacter => ErrorType::WrongDashCharacter,
            PyErrorType::ExtraWhitespace => ErrorType::ExtraWhitespace,
            PyErrorType::ProteinSubstitutionArrow => ErrorType::ProteinSubstitutionArrow,
            PyErrorType::PositionZero => ErrorType::PositionZero,
            PyErrorType::SingleLetterAminoAcid => ErrorType::SingleLetterAminoAcid,
            PyErrorType::WrongQuoteCharacter => ErrorType::WrongQuoteCharacter,
            PyErrorType::LowercaseAccessionPrefix => ErrorType::LowercaseAccessionPrefix,
            PyErrorType::MixedCaseEditType => ErrorType::MixedCaseEditType,
            PyErrorType::OldSubstitutionSyntax => ErrorType::OldSubstitutionSyntax,
            PyErrorType::InvalidUnicodeCharacter => ErrorType::InvalidUnicodeCharacter,
            PyErrorType::SwappedPositions => ErrorType::SwappedPositions,
            PyErrorType::TrailingAnnotation => ErrorType::TrailingAnnotation,
            PyErrorType::MissingCoordinatePrefix => ErrorType::MissingCoordinatePrefix,
            PyErrorType::OldAlleleFormat => ErrorType::OldAlleleFormat,
            PyErrorType::RefSeqMismatch => ErrorType::RefSeqMismatch,
        }
    }
}

/// Override behavior for a specific error type
#[pyclass(name = "ErrorOverride", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum PyErrorOverride {
    /// Use the base mode's behavior
    Default = 0,
    /// Always return an error
    Reject = 1,
    /// Auto-correct and emit a warning
    WarnCorrect = 2,
    /// Auto-correct without warning
    SilentCorrect = 3,
    /// Accept the input as-is without correction
    Accept = 4,
}

impl From<ErrorOverride> for PyErrorOverride {
    fn from(o: ErrorOverride) -> Self {
        match o {
            ErrorOverride::Default => PyErrorOverride::Default,
            ErrorOverride::Reject => PyErrorOverride::Reject,
            ErrorOverride::WarnCorrect => PyErrorOverride::WarnCorrect,
            ErrorOverride::SilentCorrect => PyErrorOverride::SilentCorrect,
            ErrorOverride::Accept => PyErrorOverride::Accept,
        }
    }
}

impl From<PyErrorOverride> for ErrorOverride {
    fn from(o: PyErrorOverride) -> Self {
        match o {
            PyErrorOverride::Default => ErrorOverride::Default,
            PyErrorOverride::Reject => ErrorOverride::Reject,
            PyErrorOverride::WarnCorrect => ErrorOverride::WarnCorrect,
            PyErrorOverride::SilentCorrect => ErrorOverride::SilentCorrect,
            PyErrorOverride::Accept => ErrorOverride::Accept,
        }
    }
}

/// Warning generated during preprocessing
#[pyclass(name = "CorrectionWarning")]
#[derive(Clone)]
pub struct PyCorrectionWarning {
    #[pyo3(get)]
    pub error_type: PyErrorType,
    #[pyo3(get)]
    pub message: String,
    #[pyo3(get)]
    pub original: String,
    #[pyo3(get)]
    pub corrected: String,
}

impl From<&CorrectionWarning> for PyCorrectionWarning {
    fn from(w: &CorrectionWarning) -> Self {
        Self {
            error_type: w.error_type.into(),
            message: w.message.clone(),
            original: w.original.clone(),
            corrected: w.corrected.clone(),
        }
    }
}

#[pymethods]
impl PyCorrectionWarning {
    fn __repr__(&self) -> String {
        format!(
            "CorrectionWarning({:?}: '{}' -> '{}')",
            self.error_type, self.original, self.corrected
        )
    }
}

/// Error handling configuration
#[pyclass(name = "ErrorConfig")]
#[derive(Clone)]
pub struct PyErrorConfig {
    inner: ErrorConfig,
}

#[pymethods]
impl PyErrorConfig {
    /// Create a strict configuration (reject all non-standard input)
    #[staticmethod]
    fn strict() -> Self {
        Self {
            inner: ErrorConfig::strict(),
        }
    }

    /// Create a lenient configuration (auto-correct with warnings)
    #[staticmethod]
    fn lenient() -> Self {
        Self {
            inner: ErrorConfig::lenient(),
        }
    }

    /// Create a silent configuration (auto-correct without warnings)
    #[staticmethod]
    fn silent() -> Self {
        Self {
            inner: ErrorConfig::silent(),
        }
    }

    /// Add an override for a specific error type
    fn with_override(&self, error_type: PyErrorType, override_: PyErrorOverride) -> Self {
        Self {
            inner: self
                .inner
                .clone()
                .with_override(error_type.into(), override_.into()),
        }
    }

    /// Get the current error mode
    #[getter]
    fn mode(&self) -> PyErrorMode {
        self.inner.mode.into()
    }

    /// Check if the given error type should be rejected
    fn should_reject(&self, error_type: PyErrorType) -> bool {
        self.inner.should_reject(error_type.into())
    }

    /// Check if the given error type should be corrected
    fn should_correct(&self, error_type: PyErrorType) -> bool {
        self.inner.should_correct(error_type.into())
    }

    /// Check if the given error type should emit a warning
    fn should_warn(&self, error_type: PyErrorType) -> bool {
        self.inner.should_warn(error_type.into())
    }
}

/// Parse result with warnings
#[pyclass(name = "ParseResultWithWarnings")]
pub struct PyParseResultWithWarnings {
    #[pyo3(get)]
    pub variant: PyHgvsVariant,
    #[pyo3(get)]
    pub warnings: Vec<PyCorrectionWarning>,
    #[pyo3(get)]
    pub original_input: String,
    #[pyo3(get)]
    pub preprocessed_input: String,
}

#[pymethods]
impl PyParseResultWithWarnings {
    /// Returns true if there were any corrections made
    fn had_corrections(&self) -> bool {
        self.original_input != self.preprocessed_input
    }

    /// Returns true if there are any warnings
    fn has_warnings(&self) -> bool {
        !self.warnings.is_empty()
    }
}

/// Parse an HGVS string with lenient error handling
///
/// Args:
///     hgvs_string: The HGVS variant description to parse
///     config: Optional ErrorConfig (defaults to lenient mode)
///
/// Returns:
///     ParseResultWithWarnings containing the parsed variant and any warnings
///
/// Raises:
///     ValueError: If the HGVS string cannot be parsed even after corrections
#[pyfunction]
#[pyo3(signature = (hgvs_string, config=None))]
fn parse_lenient(
    hgvs_string: &str,
    config: Option<&PyErrorConfig>,
) -> PyResult<PyParseResultWithWarnings> {
    let config = config
        .map(|c| c.inner.clone())
        .unwrap_or_else(ErrorConfig::lenient);

    let preprocessor = config.preprocessor();
    let preprocess_result = preprocessor.preprocess(hgvs_string);

    if !preprocess_result.success {
        return Err(PyValueError::new_err(format!(
            "Preprocessing failed: {:?}",
            preprocess_result.warnings
        )));
    }

    let variant = parse_hgvs(&preprocess_result.preprocessed)
        .map_err(|e| PyValueError::new_err(format!("Parse error: {}", e)))?;

    Ok(PyParseResultWithWarnings {
        variant: PyHgvsVariant { inner: variant },
        warnings: preprocess_result
            .warnings
            .iter()
            .map(|w| w.into())
            .collect(),
        original_input: hgvs_string.to_string(),
        preprocessed_input: preprocess_result.preprocessed,
    })
}

// ============================================================================
// Backtranslation Module
// ============================================================================

/// A codon change representing a DNA variant
#[pyclass(name = "CodonChange")]
#[derive(Clone)]
pub struct PyCodonChange {
    inner: CodonChange,
}

#[pymethods]
impl PyCodonChange {
    /// Reference codon (e.g., "CTT")
    #[getter]
    fn ref_codon(&self) -> String {
        self.inner.ref_codon.to_string()
    }

    /// Alternate codon (e.g., "TTT")
    #[getter]
    fn alt_codon(&self) -> String {
        self.inner.alt_codon.to_string()
    }

    /// Position(s) in codon that changed (1-indexed)
    #[getter]
    fn changed_positions(&self) -> Vec<u8> {
        self.inner.changed_positions.clone()
    }

    /// Number of nucleotide changes
    fn num_changes(&self) -> usize {
        self.inner.nucleotide_changes.len()
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        format!(
            "CodonChange({} -> {}, positions={:?})",
            self.inner.ref_codon, self.inner.alt_codon, self.inner.changed_positions
        )
    }
}

/// Codon table (genetic code)
#[pyclass(name = "CodonTable")]
#[derive(Clone)]
pub struct PyCodonTable {
    inner: CodonTable,
}

#[pymethods]
impl PyCodonTable {
    /// Get the standard genetic code
    #[staticmethod]
    fn standard() -> Self {
        Self {
            inner: CodonTable::standard(),
        }
    }
}

/// Backtranslation engine for converting protein changes to DNA changes
#[pyclass(name = "Backtranslator")]
pub struct PyBacktranslator {
    inner: Backtranslator,
}

#[pymethods]
impl PyBacktranslator {
    /// Create a new backtranslator with the given codon table
    #[new]
    #[pyo3(signature = (codon_table=None))]
    fn new(codon_table: Option<&PyCodonTable>) -> Self {
        let table = codon_table
            .map(|t| t.inner.clone())
            .unwrap_or_else(CodonTable::standard);
        Self {
            inner: Backtranslator::new(table),
        }
    }

    /// Create a backtranslator with the standard genetic code
    #[staticmethod]
    fn standard() -> Self {
        Self {
            inner: Backtranslator::standard(),
        }
    }

    /// Backtranslate an amino acid substitution
    ///
    /// Args:
    ///     ref_aa: Reference amino acid (3-letter code, e.g., "Leu")
    ///     alt_aa: Alternate amino acid (3-letter code, e.g., "Phe")
    ///
    /// Returns:
    ///     List of possible CodonChange objects (single-nucleotide changes only)
    fn backtranslate_substitution(
        &self,
        ref_aa: &str,
        alt_aa: &str,
    ) -> PyResult<Vec<PyCodonChange>> {
        let ref_aa = parse_amino_acid(ref_aa)?;
        let alt_aa = parse_amino_acid(alt_aa)?;
        Ok(self
            .inner
            .backtranslate_substitution(&ref_aa, &alt_aa)
            .into_iter()
            .map(|c| PyCodonChange { inner: c })
            .collect())
    }

    /// Backtranslate a nonsense mutation (amino acid to stop codon)
    ///
    /// Args:
    ///     ref_aa: Reference amino acid (3-letter code)
    ///
    /// Returns:
    ///     List of possible CodonChange objects
    fn backtranslate_to_stop(&self, ref_aa: &str) -> PyResult<Vec<PyCodonChange>> {
        let ref_aa = parse_amino_acid(ref_aa)?;
        Ok(self
            .inner
            .backtranslate_to_stop(&ref_aa)
            .into_iter()
            .map(|c| PyCodonChange { inner: c })
            .collect())
    }

    /// Backtranslate a stop loss (stop codon to amino acid)
    ///
    /// Args:
    ///     alt_aa: Alternate amino acid (3-letter code)
    ///
    /// Returns:
    ///     List of possible CodonChange objects
    fn backtranslate_stop_loss(&self, alt_aa: &str) -> PyResult<Vec<PyCodonChange>> {
        let alt_aa = parse_amino_acid(alt_aa)?;
        Ok(self
            .inner
            .backtranslate_stop_loss(&alt_aa)
            .into_iter()
            .map(|c| PyCodonChange { inner: c })
            .collect())
    }
}

// ============================================================================
// rsID Module
// ============================================================================

/// Parse rsID string to numeric value
///
/// Args:
///     rsid: rsID string (e.g., "rs121913529" or "121913529")
///
/// Returns:
///     Numeric rsID value
///
/// Raises:
///     ValueError: If the rsID cannot be parsed
#[pyfunction]
fn parse_rsid_value(rsid: &str) -> PyResult<u64> {
    rust_parse_rsid(rsid).map_err(|e| PyValueError::new_err(format!("Invalid rsID: {}", e)))
}

/// Format numeric rsID to string with "rs" prefix
///
/// Args:
///     rsid_num: Numeric rsID value
///
/// Returns:
///     rsID string (e.g., "rs121913529")
#[pyfunction]
fn format_rsid_value(rsid_num: u64) -> String {
    format_rsid(rsid_num)
}

/// Result of rsID lookup
#[pyclass(name = "RsIdResult")]
#[derive(Clone)]
pub struct PyRsIdResult {
    inner: RsIdResult,
}

#[pymethods]
impl PyRsIdResult {
    #[getter]
    fn rsid(&self) -> &str {
        &self.inner.rsid
    }

    #[getter]
    fn contig(&self) -> &str {
        &self.inner.contig
    }

    #[getter]
    fn position(&self) -> u64 {
        self.inner.position
    }

    #[getter]
    fn reference(&self) -> &str {
        &self.inner.reference
    }

    #[getter]
    fn alternate(&self) -> &str {
        &self.inner.alternate
    }

    #[getter]
    fn hgvs(&self) -> Option<&str> {
        self.inner.hgvs.as_deref()
    }

    #[getter]
    fn allele_frequency(&self) -> Option<f64> {
        self.inner.allele_frequency
    }

    #[getter]
    fn clinical_significance(&self) -> Option<&str> {
        self.inner.clinical_significance.as_deref()
    }

    /// Check if this is a substitution (SNV)
    fn is_snv(&self) -> bool {
        self.inner.is_snv()
    }

    /// Check if this is a deletion
    fn is_deletion(&self) -> bool {
        self.inner.is_deletion()
    }

    /// Check if this is an insertion
    fn is_insertion(&self) -> bool {
        self.inner.is_insertion()
    }

    /// Generate simple genomic HGVS notation
    fn to_hgvs(&self) -> String {
        self.inner.to_hgvs()
    }

    fn __repr__(&self) -> String {
        format!(
            "RsIdResult({}, {}:{}, {}>{})",
            self.inner.rsid,
            self.inner.contig,
            self.inner.position,
            self.inner.reference,
            self.inner.alternate
        )
    }
}

/// Simple in-memory rsID lookup for testing
#[pyclass(name = "InMemoryRsIdLookup")]
pub struct PyInMemoryRsIdLookup {
    inner: InMemoryRsIdLookup,
}

#[pymethods]
impl PyInMemoryRsIdLookup {
    /// Create with common test variants (BRAF V600E, BRCA1 185delAG)
    #[staticmethod]
    fn with_test_data() -> Self {
        Self {
            inner: InMemoryRsIdLookup::with_test_data(),
        }
    }

    /// Look up rsID and return matching variants
    ///
    /// Args:
    ///     rsid: rsID string (e.g., "rs121913529")
    ///
    /// Returns:
    ///     List of RsIdResult objects
    ///
    /// Raises:
    ///     ValueError: If rsID not found
    fn lookup(&self, rsid: &str) -> PyResult<Vec<PyRsIdResult>> {
        use crate::rsid::RsIdLookup;
        self.inner
            .lookup(rsid)
            .map(|results| {
                results
                    .into_iter()
                    .map(|r| PyRsIdResult { inner: r })
                    .collect()
            })
            .map_err(|e| PyValueError::new_err(format!("Lookup failed: {}", e)))
    }

    /// Check if rsID exists
    fn contains(&self, rsid: &str) -> bool {
        use crate::rsid::RsIdLookup;
        self.inner.contains(rsid)
    }
}

// ============================================================================
// VCF Module
// ============================================================================

/// A VCF record
#[pyclass(name = "VcfRecord")]
#[derive(Clone)]
pub struct PyVcfRecord {
    inner: VcfRecord,
}

#[pymethods]
impl PyVcfRecord {
    /// Create a new VCF record
    ///
    /// Args:
    ///     chrom: Chromosome name (e.g., "chr7")
    ///     pos: 1-based position
    ///     reference: Reference allele
    ///     alternate: Alternate allele(s) - can be a single string or list of strings
    ///     id: Optional variant ID (e.g., rsID)
    #[new]
    #[pyo3(signature = (chrom, pos, reference, alternate, id=None))]
    fn new(chrom: &str, pos: u64, reference: &str, alternate: &str, id: Option<&str>) -> Self {
        Self {
            inner: VcfRecord {
                chrom: chrom.to_string(),
                pos,
                id: id.map(|s| s.to_string()),
                reference: reference.to_string(),
                alternate: vec![alternate.to_string()],
                quality: None,
                filter: None,
                info: std::collections::HashMap::new(),
                format: None,
                samples: Vec::new(),
                genome_build: GenomeBuild::default(),
            },
        }
    }

    /// Create a SNV record
    ///
    /// Args:
    ///     chrom: Chromosome name
    ///     pos: 1-based position
    ///     ref_base: Reference base (single character)
    ///     alt_base: Alternate base (single character)
    ///
    /// Raises:
    ///     ValueError: If ref_base or alt_base is empty
    #[staticmethod]
    fn snv(chrom: &str, pos: u64, ref_base: &str, alt_base: &str) -> PyResult<Self> {
        let ref_char = ref_base
            .chars()
            .next()
            .ok_or_else(|| PyValueError::new_err("ref_base must be a non-empty string"))?;
        let alt_char = alt_base
            .chars()
            .next()
            .ok_or_else(|| PyValueError::new_err("alt_base must be a non-empty string"))?;
        Ok(Self {
            inner: VcfRecord::snv(chrom, pos, ref_char, alt_char),
        })
    }

    #[getter]
    fn chrom(&self) -> &str {
        &self.inner.chrom
    }

    #[getter]
    fn pos(&self) -> u64 {
        self.inner.pos
    }

    #[getter]
    fn id(&self) -> Option<&str> {
        self.inner.id.as_deref()
    }

    #[getter]
    fn reference(&self) -> &str {
        &self.inner.reference
    }

    /// Get the first alternate allele (for simple variants)
    #[getter]
    fn alternate(&self) -> Option<&str> {
        self.inner.alternate.first().map(|s| s.as_str())
    }

    /// Get all alternate alleles
    #[getter]
    fn alternates(&self) -> Vec<String> {
        self.inner.alternate.clone()
    }

    fn __repr__(&self) -> String {
        let alt_str = self.inner.alternate.join(",");
        format!(
            "VcfRecord({}:{} {}>{})",
            self.inner.chrom, self.inner.pos, self.inner.reference, alt_str
        )
    }
}

/// Convert a VCF record to genomic HGVS notation
///
/// Args:
///     record: VcfRecord object
///     alt_index: Index of the alternate allele to convert (default 0)
///
/// Returns:
///     HgvsVariant object
///
/// Raises:
///     ValueError: If conversion fails
#[pyfunction]
#[pyo3(signature = (record, alt_index=0))]
fn vcf_to_genomic_hgvs(record: &PyVcfRecord, alt_index: usize) -> PyResult<PyHgvsVariant> {
    rust_vcf_to_hgvs(&record.inner, alt_index)
        .map(|v| PyHgvsVariant {
            inner: HgvsVariant::Genome(v),
        })
        .map_err(|e| PyValueError::new_err(format!("Conversion error: {}", e)))
}

// ============================================================================
// Prepare Module
// ============================================================================

/// Configuration for reference data preparation
#[pyclass(name = "PrepareConfig")]
#[derive(Clone)]
pub struct PyPrepareConfig {
    inner: PrepareConfig,
}

#[pymethods]
impl PyPrepareConfig {
    /// Create a new prepare configuration
    #[new]
    #[pyo3(signature = (
        output_dir="ferro-reference",
        download_transcripts=true,
        download_genome=false,
        download_genome_grch37=false,
        download_refseqgene=false,
        download_lrg=false,
        download_cdot=true,
        skip_existing=true,
        dry_run=false
    ))]
    fn new(
        output_dir: &str,
        download_transcripts: bool,
        download_genome: bool,
        download_genome_grch37: bool,
        download_refseqgene: bool,
        download_lrg: bool,
        download_cdot: bool,
        skip_existing: bool,
        dry_run: bool,
    ) -> Self {
        Self {
            inner: PrepareConfig {
                output_dir: std::path::PathBuf::from(output_dir),
                download_transcripts,
                download_genome,
                download_genome_grch37,
                download_refseqgene,
                download_lrg,
                download_cdot,
                skip_existing,
                clinvar_file: None,
                patterns_file: None,
                dry_run,
            },
        }
    }

    #[getter]
    fn output_dir(&self) -> String {
        self.inner.output_dir.display().to_string()
    }

    #[getter]
    fn download_transcripts(&self) -> bool {
        self.inner.download_transcripts
    }

    #[getter]
    fn download_genome(&self) -> bool {
        self.inner.download_genome
    }

    #[getter]
    fn download_cdot(&self) -> bool {
        self.inner.download_cdot
    }
}

/// Manifest of prepared reference data
#[pyclass(name = "ReferenceManifest")]
#[derive(Clone)]
pub struct PyReferenceManifest {
    inner: ReferenceManifest,
}

#[pymethods]
impl PyReferenceManifest {
    #[getter]
    fn prepared_at(&self) -> &str {
        &self.inner.prepared_at
    }

    #[getter]
    fn transcript_count(&self) -> usize {
        self.inner.transcript_count
    }

    #[getter]
    fn transcript_fastas(&self) -> Vec<String> {
        self.inner
            .transcript_fastas
            .iter()
            .map(|p| p.display().to_string())
            .collect()
    }

    #[getter]
    fn genome_fasta(&self) -> Option<String> {
        self.inner
            .genome_fasta
            .as_ref()
            .map(|p| p.display().to_string())
    }

    #[getter]
    fn cdot_json(&self) -> Option<String> {
        self.inner
            .cdot_json
            .as_ref()
            .map(|p| p.display().to_string())
    }

    #[getter]
    fn available_prefixes(&self) -> Vec<String> {
        self.inner.available_prefixes.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "ReferenceManifest(transcripts={}, prefixes={:?})",
            self.inner.transcript_count, self.inner.available_prefixes
        )
    }
}

/// Prepare reference data for normalization
///
/// Args:
///     config: PrepareConfig with download options
///
/// Returns:
///     ReferenceManifest describing the prepared data
///
/// Raises:
///     RuntimeError: If preparation fails
#[pyfunction]
fn prepare_reference_data(config: &PyPrepareConfig) -> PyResult<PyReferenceManifest> {
    prepare_references(&config.inner)
        .map(|m| PyReferenceManifest { inner: m })
        .map_err(|e| PyRuntimeError::new_err(format!("Prepare failed: {}", e)))
}

/// Check existing reference data
///
/// Args:
///     directory: Path to reference data directory
///
/// Returns:
///     ReferenceManifest if manifest.json exists
///
/// Raises:
///     RuntimeError: If check fails
#[pyfunction]
fn check_reference_data(directory: &str) -> PyResult<PyReferenceManifest> {
    check_references(Path::new(directory))
        .map(|m| PyReferenceManifest { inner: m })
        .map_err(|e| PyRuntimeError::new_err(format!("Check failed: {}", e)))
}

// ============================================================================
// Reference Module
// ============================================================================

/// Genome build (GRCh37 or GRCh38)
#[pyclass(name = "GenomeBuild", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum PyGenomeBuild {
    GRCh37 = 0,
    GRCh38 = 1,
    Unknown = 2,
}

impl From<GenomeBuild> for PyGenomeBuild {
    fn from(b: GenomeBuild) -> Self {
        match b {
            GenomeBuild::GRCh37 => PyGenomeBuild::GRCh37,
            GenomeBuild::GRCh38 => PyGenomeBuild::GRCh38,
            GenomeBuild::Unknown => PyGenomeBuild::Unknown,
        }
    }
}

#[pymethods]
impl PyGenomeBuild {
    fn __str__(&self) -> &'static str {
        match self {
            PyGenomeBuild::GRCh37 => "GRCh37",
            PyGenomeBuild::GRCh38 => "GRCh38",
            PyGenomeBuild::Unknown => "Unknown",
        }
    }
}

/// Strand direction
#[pyclass(name = "Strand", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum PyStrand {
    Plus = 0,
    Minus = 1,
}

impl From<Strand> for PyStrand {
    fn from(s: Strand) -> Self {
        match s {
            Strand::Plus => PyStrand::Plus,
            Strand::Minus => PyStrand::Minus,
        }
    }
}

#[pymethods]
impl PyStrand {
    fn __str__(&self) -> &'static str {
        match self {
            PyStrand::Plus => "+",
            PyStrand::Minus => "-",
        }
    }
}

// ============================================================================
// Coordinate Mapper Module
// ============================================================================

/// Coordinate mapper for converting between HGVS coordinate systems
///
/// This class provides methods for converting between genomic (g.), coding (c.),
/// transcript (n.), and protein (p.) coordinate systems.
#[pyclass(name = "CoordinateMapper")]
pub struct PyCoordinateMapper {
    provider: MockProvider,
}

#[pymethods]
impl PyCoordinateMapper {
    /// Create a new coordinate mapper
    ///
    /// Args:
    ///     reference_json: Optional path to a transcripts.json file for reference data.
    ///         If not provided, uses built-in test data.
    #[new]
    #[pyo3(signature = (reference_json=None))]
    fn new(reference_json: Option<&str>) -> PyResult<Self> {
        let provider = match reference_json {
            Some(path) => MockProvider::from_json(Path::new(path))
                .map_err(|e| PyRuntimeError::new_err(format!("Failed to load reference: {}", e)))?,
            None => MockProvider::with_test_data(),
        };

        Ok(Self { provider })
    }

    /// Convert a CDS position to genomic position
    ///
    /// Args:
    ///     transcript_id: Transcript accession (e.g., "NM_000088.3")
    ///     cds_position: CDS position (e.g., 100 for c.100)
    ///     offset: Optional intronic offset (e.g., 5 for c.100+5)
    ///
    /// Returns:
    ///     Tuple of (chromosome, genomic_position) or None if position is intronic
    ///     without offset support
    ///
    /// Raises:
    ///     ValueError: If transcript not found or has no genomic coordinates
    #[pyo3(signature = (transcript_id, cds_position, offset=None))]
    fn c_to_g(
        &self,
        transcript_id: &str,
        cds_position: i64,
        offset: Option<i64>,
    ) -> PyResult<Option<(String, u64)>> {
        let transcript = self.provider.get_transcript(transcript_id).map_err(|e| {
            PyValueError::new_err(format!("Transcript not found: {}: {}", transcript_id, e))
        })?;

        let mapper = CoordinateMapper::new(&transcript);

        // Treat offset=0 as no offset (semantically equivalent)
        let effective_offset = offset.filter(|&o| o != 0);

        let cds_pos = CdsPos {
            base: cds_position,
            offset: effective_offset,
            utr3: false,
        };

        if effective_offset.is_some() {
            // Use intronic-aware conversion
            match mapper.cds_to_genomic_with_intron(&cds_pos) {
                Ok(genomic_pos) => {
                    let chrom = mapper
                        .chromosome()
                        .ok_or_else(|| PyValueError::new_err("Transcript has no chromosome"))?
                        .to_string();
                    Ok(Some((chrom, genomic_pos)))
                }
                Err(e) => Err(PyValueError::new_err(format!("Conversion error: {}", e))),
            }
        } else {
            // Use standard conversion
            match mapper.cds_to_genomic(&cds_pos) {
                Ok(Some(genomic_pos)) => {
                    let chrom = mapper
                        .chromosome()
                        .ok_or_else(|| PyValueError::new_err("Transcript has no chromosome"))?
                        .to_string();
                    Ok(Some((chrom, genomic_pos)))
                }
                Ok(None) => Ok(None),
                Err(e) => Err(PyValueError::new_err(format!("Conversion error: {}", e))),
            }
        }
    }

    /// Convert a genomic position to CDS position
    ///
    /// Args:
    ///     transcript_id: Transcript accession (e.g., "NM_000088.3")
    ///     genomic_position: Genomic position (1-based)
    ///
    /// Returns:
    ///     Tuple of (cds_position, offset, is_utr3). For exonic positions,
    ///     offset will be None. For intronic positions, offset will contain
    ///     the distance from the nearest exon boundary.
    ///
    /// Raises:
    ///     ValueError: If transcript not found, position is outside transcript bounds,
    ///         or conversion fails
    fn g_to_c(
        &self,
        transcript_id: &str,
        genomic_position: u64,
    ) -> PyResult<(i64, Option<i64>, bool)> {
        let transcript = self.provider.get_transcript(transcript_id).map_err(|e| {
            PyValueError::new_err(format!("Transcript not found: {}: {}", transcript_id, e))
        })?;

        let mapper = CoordinateMapper::new(&transcript);

        match mapper.genomic_to_cds(genomic_position) {
            Ok(Some(cds_pos)) => Ok((cds_pos.base, cds_pos.offset, cds_pos.utr3)),
            Ok(None) => {
                // Position is intronic, try with intronic support
                mapper
                    .genomic_to_cds_intronic(genomic_position)
                    .map(|cds_pos| (cds_pos.base, cds_pos.offset, cds_pos.utr3))
                    .map_err(|e| {
                        PyValueError::new_err(format!(
                        "Position {} is outside transcript bounds or in an unsupported region: {}",
                        genomic_position, e
                    ))
                    })
            }
            Err(e) => Err(PyValueError::new_err(format!("Conversion error: {}", e))),
        }
    }

    /// Convert a CDS position to protein position
    ///
    /// Args:
    ///     transcript_id: Transcript accession (e.g., "NM_000088.3")
    ///     cds_position: CDS position (e.g., 100 for c.100)
    ///
    /// Returns:
    ///     Protein position (1-based codon number)
    ///
    /// Raises:
    ///     ValueError: If transcript not found or position is in UTR/intronic
    fn c_to_p(&self, transcript_id: &str, cds_position: i64) -> PyResult<u64> {
        let transcript = self.provider.get_transcript(transcript_id).map_err(|e| {
            PyValueError::new_err(format!("Transcript not found: {}: {}", transcript_id, e))
        })?;

        let mapper = CoordinateMapper::new(&transcript);

        let cds_pos = CdsPos::new(cds_position);

        mapper
            .cds_to_protein(&cds_pos)
            .map(|prot_pos| prot_pos.number)
            .map_err(|e| PyValueError::new_err(format!("Conversion error: {}", e)))
    }

    /// Convert a CDS position to transcript position
    ///
    /// Args:
    ///     transcript_id: Transcript accession (e.g., "NM_000088.3")
    ///     cds_position: CDS position (e.g., 100 for c.100)
    ///     offset: Optional intronic offset
    ///     utr3: Whether this is a 3' UTR position (c.*N notation)
    ///
    /// Returns:
    ///     Tuple of (transcript_position, offset)
    ///
    /// Raises:
    ///     ValueError: If transcript not found
    #[pyo3(signature = (transcript_id, cds_position, offset=None, utr3=false))]
    fn c_to_n(
        &self,
        transcript_id: &str,
        cds_position: i64,
        offset: Option<i64>,
        utr3: bool,
    ) -> PyResult<(i64, Option<i64>)> {
        let transcript = self.provider.get_transcript(transcript_id).map_err(|e| {
            PyValueError::new_err(format!("Transcript not found: {}: {}", transcript_id, e))
        })?;

        let mapper = CoordinateMapper::new(&transcript);

        let cds_pos = CdsPos {
            base: cds_position,
            offset,
            utr3,
        };

        mapper
            .cds_to_tx(&cds_pos)
            .map(|tx_pos| (tx_pos.base, tx_pos.offset))
            .map_err(|e| PyValueError::new_err(format!("Conversion error: {}", e)))
    }

    /// Convert a transcript position to CDS position
    ///
    /// Args:
    ///     transcript_id: Transcript accession (e.g., "NM_000088.3")
    ///     tx_position: Transcript position (1-based)
    ///     offset: Optional intronic offset
    ///     downstream: Whether this is a downstream position (n.*100 notation for
    ///         positions past the end of the transcript). Defaults to False.
    ///
    /// Returns:
    ///     Tuple of (cds_position, offset, is_utr3)
    ///
    /// Raises:
    ///     ValueError: If transcript not found or has no CDS
    #[pyo3(signature = (transcript_id, tx_position, offset=None, downstream=false))]
    fn n_to_c(
        &self,
        transcript_id: &str,
        tx_position: i64,
        offset: Option<i64>,
        downstream: bool,
    ) -> PyResult<(i64, Option<i64>, bool)> {
        let transcript = self.provider.get_transcript(transcript_id).map_err(|e| {
            PyValueError::new_err(format!("Transcript not found: {}: {}", transcript_id, e))
        })?;

        let mapper = CoordinateMapper::new(&transcript);

        let tx_pos = TxPos {
            base: tx_position,
            offset,
            downstream,
        };

        mapper
            .tx_to_cds(&tx_pos)
            .map(|cds_pos| (cds_pos.base, cds_pos.offset, cds_pos.utr3))
            .map_err(|e| PyValueError::new_err(format!("Conversion error: {}", e)))
    }

    /// Get the strand of a transcript
    ///
    /// Args:
    ///     transcript_id: Transcript accession
    ///
    /// Returns:
    ///     Strand (Plus or Minus)
    fn get_strand(&self, transcript_id: &str) -> PyResult<PyStrand> {
        let transcript = self.provider.get_transcript(transcript_id).map_err(|e| {
            PyValueError::new_err(format!("Transcript not found: {}: {}", transcript_id, e))
        })?;

        Ok(transcript.strand.into())
    }

    /// Get the chromosome of a transcript
    ///
    /// Args:
    ///     transcript_id: Transcript accession
    ///
    /// Returns:
    ///     Chromosome name or None if not set
    fn get_chromosome(&self, transcript_id: &str) -> PyResult<Option<String>> {
        let transcript = self.provider.get_transcript(transcript_id).map_err(|e| {
            PyValueError::new_err(format!("Transcript not found: {}: {}", transcript_id, e))
        })?;

        Ok(transcript.chromosome.clone())
    }

    /// Check if a transcript exists in the reference
    fn has_transcript(&self, transcript_id: &str) -> bool {
        self.provider.get_transcript(transcript_id).is_ok()
    }
}

// ============================================================================
// Enhanced PyHgvsVariant
// ============================================================================

/// ferro-hgvs Python module
#[pymodule]
fn ferro_hgvs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Core functions
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    m.add_function(wrap_pyfunction!(normalize, m)?)?;

    // SPDI functions
    m.add_function(wrap_pyfunction!(parse_spdi, m)?)?;
    m.add_function(wrap_pyfunction!(hgvs_to_spdi, m)?)?;
    m.add_function(wrap_pyfunction!(spdi_to_hgvs_variant, m)?)?;

    // Coordinate functions
    m.add_function(wrap_pyfunction!(hgvs_pos_to_index, m)?)?;
    m.add_function(wrap_pyfunction!(index_to_hgvs_pos, m)?)?;

    // MAVE functions
    m.add_function(wrap_pyfunction!(parse_mave_hgvs_variant, m)?)?;
    m.add_function(wrap_pyfunction!(is_mave_short_form_variant, m)?)?;

    // Error handling functions
    m.add_function(wrap_pyfunction!(parse_lenient, m)?)?;

    // rsID functions
    m.add_function(wrap_pyfunction!(parse_rsid_value, m)?)?;
    m.add_function(wrap_pyfunction!(format_rsid_value, m)?)?;

    // VCF functions
    m.add_function(wrap_pyfunction!(vcf_to_genomic_hgvs, m)?)?;

    // Prepare functions
    m.add_function(wrap_pyfunction!(prepare_reference_data, m)?)?;
    m.add_function(wrap_pyfunction!(check_reference_data, m)?)?;

    // Core classes
    m.add_class::<PyHgvsVariant>()?;
    m.add_class::<PyNormalizer>()?;

    // SPDI classes
    m.add_class::<PySpdiVariant>()?;

    // Coordinate classes
    m.add_class::<PyZeroBasedPos>()?;
    m.add_class::<PyOneBasedPos>()?;

    // Equivalence classes
    m.add_class::<PyEquivalenceLevel>()?;
    m.add_class::<PyEquivalenceResult>()?;
    m.add_class::<PyEquivalenceChecker>()?;

    // Effect prediction classes
    m.add_class::<PyConsequence>()?;
    m.add_class::<PyImpact>()?;
    m.add_class::<PyProteinEffect>()?;
    m.add_class::<PyEffectPredictor>()?;

    // MAVE classes
    m.add_class::<PyMaveContext>()?;

    // Batch processing classes
    m.add_class::<PyBatchProgress>()?;
    m.add_class::<PyBatchResult>()?;
    m.add_class::<PyBatchProcessor>()?;

    // Error handling classes
    m.add_class::<PyErrorMode>()?;
    m.add_class::<PyErrorType>()?;
    m.add_class::<PyErrorOverride>()?;
    m.add_class::<PyCorrectionWarning>()?;
    m.add_class::<PyErrorConfig>()?;
    m.add_class::<PyParseResultWithWarnings>()?;

    // Backtranslation classes
    m.add_class::<PyCodonChange>()?;
    m.add_class::<PyCodonTable>()?;
    m.add_class::<PyBacktranslator>()?;

    // rsID classes
    m.add_class::<PyRsIdResult>()?;
    m.add_class::<PyInMemoryRsIdLookup>()?;

    // VCF classes
    m.add_class::<PyVcfRecord>()?;

    // Prepare classes
    m.add_class::<PyPrepareConfig>()?;
    m.add_class::<PyReferenceManifest>()?;

    // Reference classes
    m.add_class::<PyGenomeBuild>()?;
    m.add_class::<PyStrand>()?;

    // Coordinate mapping classes
    m.add_class::<PyCoordinateMapper>()?;

    // Add version
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
