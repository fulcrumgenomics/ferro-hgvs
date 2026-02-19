"""ferro-hgvs: HGVS variant normalizer - Python bindings for the ferro bioinformatics toolkit.

This module provides Python bindings for HGVS variant parsing, normalization,
and various related utilities including SPDI conversion, effect prediction,
and batch processing.

Example:
    >>> import ferro_hgvs
    >>> variant = ferro_hgvs.parse("NM_000088.3:c.100A>G")
    >>> print(variant.variant_type)
    'coding'
"""

# Re-export everything from the native Rust extension module
from .ferro_hgvs import *

# Explicitly list the public API to help with IDE completion and type checking
__all__ = [
    # Version
    "__version__",
    # Core functions
    "parse",
    "normalize",
    # SPDI functions
    "parse_spdi",
    "hgvs_to_spdi",
    "spdi_to_hgvs_variant",
    # Coordinate functions
    "hgvs_pos_to_index",
    "index_to_hgvs_pos",
    # MAVE functions
    "parse_mave_hgvs_variant",
    "is_mave_short_form_variant",
    # Error handling functions
    "parse_lenient",
    # rsID functions
    "parse_rsid_value",
    "format_rsid_value",
    # VCF functions
    "vcf_to_genomic_hgvs",
    # Prepare functions
    "prepare_reference_data",
    "check_reference_data",
    # Core classes
    "HgvsVariant",
    "Normalizer",
    # SPDI classes
    "SpdiVariant",
    # Coordinate classes
    "ZeroBasedPos",
    "OneBasedPos",
    # Equivalence classes
    "EquivalenceLevel",
    "EquivalenceResult",
    "EquivalenceChecker",
    # Effect prediction classes
    "Consequence",
    "Impact",
    "ProteinEffect",
    "EffectPredictor",
    # MAVE classes
    "MaveContext",
    # Batch processing classes
    "BatchProgress",
    "BatchResult",
    "BatchProcessor",
    # Error handling classes
    "ErrorMode",
    "ErrorType",
    "ErrorOverride",
    "CorrectionWarning",
    "ErrorConfig",
    "ParseResultWithWarnings",
    # Backtranslation classes
    "CodonChange",
    "CodonTable",
    "Backtranslator",
    # rsID classes
    "RsIdResult",
    "InMemoryRsIdLookup",
    # VCF classes
    "VcfRecord",
    # Prepare classes
    "PrepareConfig",
    "ReferenceManifest",
    # Reference classes
    "GenomeBuild",
    "Strand",
    # Coordinate mapping classes
    "CoordinateMapper",
]
