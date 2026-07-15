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

from collections.abc import Iterable

# Re-export everything from the native Rust extension module
from .ferro_hgvs import *


class FerroError(Exception):
    """Base class for all ferro-hgvs variant-processing errors.

    Every ferro failure that stems from parsing, normalizing, projecting, or
    resolving reference data is an instance of ``FerroError`` (or one of its
    subclasses). Catch ``FerroError`` to handle any of them; catch a specific
    subclass to discriminate.

    Each subclass also inherits the built-in exception the corresponding
    operation historically raised (``ParseError`` is a ``ValueError``,
    ``NormalizationError`` is a ``RuntimeError``, and so on), so existing
    ``except ValueError`` / ``except RuntimeError`` callers keep working.

    Pure argument-validation failures (an empty string, an out-of-range index,
    an unrecognized enum spelling) are raised as a plain ``ValueError`` — they
    are programming errors, not variant-processing failures, and are outside
    this hierarchy.

    Attributes:
        code: The structured ``E####`` / ``W####`` error-code string, or
            ``None`` when the underlying failure carries no code.
        mutalyzer_codes: Equivalent mutalyzer diagnostic codes (e.g.
            ``"EINTRONIC"``); empty when there is no mutalyzer equivalent.
    """

    def __init__(
        self,
        message: str,
        code: str | None = None,
        mutalyzer_codes: Iterable[str] = (),
    ) -> None:
        super().__init__(message)
        self.code = code
        self.mutalyzer_codes: tuple[str, ...] = tuple(mutalyzer_codes)

    def __reduce__(
        self,
    ) -> tuple[type["FerroError"], tuple[str, str | None, tuple[str, ...]]]:
        # ``super().__init__(message)`` only seeds ``self.args``, so the default
        # ``BaseException.__reduce__`` reconstructs from ``args`` alone and drops
        # ``code`` / ``mutalyzer_codes`` on pickle / ``copy.deepcopy`` / cross-process
        # propagation (e.g. an exception raised in a ``multiprocessing`` worker) —
        # exactly the structured data this hierarchy exists to carry. Reconstruct
        # via the full 3-arg ``__init__`` instead. ``self.__class__`` preserves the
        # concrete subclass, and every subclass shares this ``__init__`` signature.
        message = self.args[0] if self.args else ""
        return (self.__class__, (message, self.code, self.mutalyzer_codes))


class ParseError(FerroError, ValueError):
    """Raised when an HGVS or SPDI string cannot be parsed."""


class NormalizationError(FerroError, RuntimeError):
    """Raised when normalizing a variant fails."""


# Named ``ReferenceDataError`` (not ``ReferenceError``) so it does not shadow the
# built-in ``ReferenceError`` (a weakref error) when callers ``from ferro_hgvs
# import *``.
class ReferenceDataError(FerroError, ValueError, RuntimeError):
    """Raised when reference data is unavailable or a reference lookup fails."""


class ProjectionError(FerroError, ValueError, RuntimeError):
    """Raised when projecting or converting a variant between coordinate systems fails."""


class EquivalenceError(FerroError, RuntimeError):
    """Raised when an equivalence check between variants fails."""


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
    # Convert-GFF functions
    "convert_gff",
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
    "NormalizationWarning",
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
    # Convert-GFF classes
    "ConvertGffConfig",
    "ConvertGffReport",
    # Reference classes
    "GenomeBuild",
    "Strand",
    # Coordinate mapping classes
    "CoordinateMapper",
    # Exception hierarchy
    "FerroError",
    "ParseError",
    "NormalizationError",
    "ReferenceDataError",
    "ProjectionError",
    "EquivalenceError",
]
