"""Type stubs for ferro-hgvs Python bindings."""

from enum import IntEnum
from typing import Any, Callable, Dict, List, Optional

__version__: str

# ============================================================================
# Core Functions
# ============================================================================

def parse(hgvs_string: str) -> HgvsVariant:
    """Parse an HGVS variant string.

    Args:
        hgvs_string: The HGVS variant description to parse

    Returns:
        A HgvsVariant object representing the parsed variant

    Raises:
        ValueError: If the HGVS string cannot be parsed

    Example:
        >>> variant = parse("NM_000088.3:c.100A>G")
        >>> variant.variant_type
        'coding'
    """
    ...

def normalize(hgvs_string: str, direction: str = "3prime") -> str:
    """Normalize an HGVS variant string.

    Args:
        hgvs_string: The HGVS variant description to normalize
        direction: Shuffle direction - "3prime" (default) or "5prime"

    Returns:
        The normalized HGVS string

    Raises:
        ValueError: If the HGVS string cannot be parsed
        RuntimeError: If normalization fails

    Example:
        >>> normalize("NM_000088.3:c.100delA")
        'NM_000088.3:c.100del'
    """
    ...

# ============================================================================
# SPDI Functions
# ============================================================================

def parse_spdi(spdi_string: str) -> SpdiVariant:
    """Parse an SPDI variant string.

    Args:
        spdi_string: The SPDI variant description (e.g., "NC_000001.11:12344:A:G")

    Returns:
        A SpdiVariant object representing the parsed variant

    Raises:
        ValueError: If the SPDI string cannot be parsed
    """
    ...

def hgvs_to_spdi(variant: HgvsVariant) -> SpdiVariant:
    """Convert an HGVS variant to SPDI format.

    Args:
        variant: An HgvsVariant object

    Returns:
        A SpdiVariant object

    Raises:
        ValueError: If the conversion fails
    """
    ...

def spdi_to_hgvs_variant(spdi: SpdiVariant) -> HgvsVariant:
    """Convert an SPDI variant to HGVS format.

    Args:
        spdi: A SpdiVariant object

    Returns:
        An HgvsVariant object

    Raises:
        ValueError: If the conversion fails
    """
    ...

# ============================================================================
# Coordinate Functions
# ============================================================================

def hgvs_pos_to_index(pos: int) -> int:
    """Convert a 1-based HGVS position to a 0-based array index."""
    ...

def index_to_hgvs_pos(idx: int) -> int:
    """Convert a 0-based array index to a 1-based HGVS position."""
    ...

# ============================================================================
# MAVE Functions
# ============================================================================

def parse_mave_hgvs_variant(hgvs_string: str, context: MaveContext) -> HgvsVariant:
    """Parse a MAVE-HGVS variant string with context.

    Args:
        hgvs_string: The MAVE-HGVS variant description (e.g., "p.Glu6Val")
        context: MaveContext with accession information

    Returns:
        A HgvsVariant object with the accession filled in from context

    Raises:
        ValueError: If parsing fails or context doesn't support the coordinate type
    """
    ...

def is_mave_short_form_variant(hgvs_string: str) -> bool:
    """Check if a string is in MAVE short-form notation (no accession).

    Args:
        hgvs_string: The HGVS string to check

    Returns:
        True if the string is in short form (e.g., "p.Val600Glu"), False otherwise
    """
    ...

# ============================================================================
# Error Handling Functions
# ============================================================================

def parse_lenient(
    hgvs_string: str, config: Optional[ErrorConfig] = None
) -> ParseResultWithWarnings:
    """Parse an HGVS string with lenient error handling.

    Args:
        hgvs_string: The HGVS variant description to parse
        config: Optional ErrorConfig (defaults to lenient mode)

    Returns:
        ParseResultWithWarnings containing the parsed variant and any warnings

    Raises:
        ValueError: If the HGVS string cannot be parsed even after corrections
    """
    ...

# ============================================================================
# rsID Functions
# ============================================================================

def parse_rsid_value(rsid: str) -> int:
    """Parse rsID string to numeric value.

    Args:
        rsid: rsID string (e.g., "rs121913529" or "121913529")

    Returns:
        Numeric rsID value

    Raises:
        ValueError: If the rsID cannot be parsed
    """
    ...

def format_rsid_value(rsid_num: int) -> str:
    """Format numeric rsID to string with "rs" prefix.

    Args:
        rsid_num: Numeric rsID value

    Returns:
        rsID string (e.g., "rs121913529")
    """
    ...

# ============================================================================
# VCF Functions
# ============================================================================

def vcf_to_genomic_hgvs(record: VcfRecord, alt_index: int = 0) -> HgvsVariant:
    """Convert a VCF record to genomic HGVS notation.

    Args:
        record: VcfRecord object
        alt_index: Index of the alternate allele to convert (default 0)

    Returns:
        HgvsVariant object

    Raises:
        ValueError: If conversion fails
    """
    ...

# ============================================================================
# Prepare Functions
# ============================================================================

def prepare_reference_data(config: PrepareConfig) -> ReferenceManifest:
    """Prepare reference data for normalization.

    Args:
        config: PrepareConfig with download options

    Returns:
        ReferenceManifest describing the prepared data

    Raises:
        RuntimeError: If preparation fails
    """
    ...

def check_reference_data(directory: str) -> ReferenceManifest:
    """Check existing reference data.

    Args:
        directory: Path to reference data directory

    Returns:
        ReferenceManifest if manifest.json exists

    Raises:
        RuntimeError: If check fails
    """
    ...

# ============================================================================
# Core Classes
# ============================================================================

class HgvsVariant:
    """A parsed HGVS variant.

    Attributes:
        variant_type: The type of variant (genomic, coding, non_coding, protein, rna, mitochondrial)
        reference: The reference accession (e.g., "NM_000088.3")
        edit_type: The type of edit (substitution, deletion, duplication, insertion, delins, etc.)
    """

    @property
    def variant_type(self) -> str:
        """Get the variant type as a string."""
        ...

    @property
    def reference(self) -> str:
        """Get the reference (accession) of the variant."""
        ...

    @property
    def edit_type(self) -> str:
        """Get the edit type as a string."""
        ...

    def normalize(self, direction: str = "3prime") -> HgvsVariant:
        """Normalize this variant."""
        ...

    def to_dict(self) -> Dict[str, Any]:
        """Convert to a dictionary representation."""
        ...

    def to_json(self) -> str:
        """Convert to JSON string."""
        ...

    def is_genomic(self) -> bool:
        """Check if this is a genomic variant (g. prefix)."""
        ...

    def is_coding(self) -> bool:
        """Check if this is a coding variant (c. prefix)."""
        ...

    def is_noncoding(self) -> bool:
        """Check if this is a non-coding variant (n. prefix)."""
        ...

    def is_protein(self) -> bool:
        """Check if this is a protein variant (p. prefix)."""
        ...

    def is_rna(self) -> bool:
        """Check if this is an RNA variant (r. prefix)."""
        ...

    def is_mitochondrial(self) -> bool:
        """Check if this is a mitochondrial variant (m. prefix)."""
        ...

    def is_substitution(self) -> bool:
        """Check if this is a substitution."""
        ...

    def is_deletion(self) -> bool:
        """Check if this is a deletion."""
        ...

    def is_insertion(self) -> bool:
        """Check if this is an insertion."""
        ...

    def is_duplication(self) -> bool:
        """Check if this is a duplication."""
        ...

    def is_delins(self) -> bool:
        """Check if this is a deletion-insertion (delins)."""
        ...

    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...

class Normalizer:
    """HGVS variant normalizer using a reference provider."""

    def __init__(self, reference_json: Optional[str] = None, direction: str = "3prime") -> None:
        """Create a new normalizer.

        Args:
            reference_json: Optional path to a transcripts.json file
            direction: Shuffle direction - "3prime" (default) or "5prime"
        """
        ...

    def parse(self, hgvs_string: str) -> HgvsVariant:
        """Parse an HGVS string."""
        ...

    def normalize_variant(self, variant: HgvsVariant) -> HgvsVariant:
        """Normalize an HGVS variant object."""
        ...

    def normalize(self, hgvs_string: str) -> str:
        """Parse and normalize an HGVS string."""
        ...

# ============================================================================
# SPDI Classes
# ============================================================================

class SpdiVariant:
    """SPDI (Sequence, Position, Deletion, Insertion) variant representation."""

    def __init__(self, sequence: str, position: int, deletion: str, insertion: str) -> None:
        """Create a new SPDI variant.

        Args:
            sequence: Reference sequence identifier (e.g., "NC_000001.11")
            position: 0-based interbase position
            deletion: Deleted sequence (can be empty for insertions)
            insertion: Inserted sequence (can be empty for deletions)
        """
        ...

    @property
    def sequence(self) -> str:
        """Reference sequence identifier."""
        ...

    @property
    def position(self) -> int:
        """0-based interbase position."""
        ...

    @property
    def deletion(self) -> str:
        """Deleted sequence."""
        ...

    @property
    def insertion(self) -> str:
        """Inserted sequence."""
        ...

    def is_substitution(self) -> bool:
        """Check if this is a substitution (SNV or MNV)."""
        ...

    def is_snv(self) -> bool:
        """Check if this is a single nucleotide variant."""
        ...

    def is_deletion(self) -> bool:
        """Check if this is a pure deletion."""
        ...

    def is_insertion(self) -> bool:
        """Check if this is a pure insertion."""
        ...

    def is_delins(self) -> bool:
        """Check if this is a deletion-insertion (delins)."""
        ...

    def is_identity(self) -> bool:
        """Check if this represents an identity (no change)."""
        ...

    def variant_type(self) -> str:
        """Get the variant type as a string."""
        ...

    def to_one_based_position(self) -> int:
        """Convert 0-based SPDI position to 1-based HGVS position."""
        ...

    def to_dict(self) -> Dict[str, Any]:
        """Convert to a dictionary representation."""
        ...

    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...

# ============================================================================
# Coordinate Classes
# ============================================================================

class ZeroBasedPos:
    """A 0-based position (array-style indexing, used by SPDI and BED)."""

    def __init__(self, pos: int) -> None:
        """Create a new 0-based position."""
        ...

    @property
    def value(self) -> int:
        """Get the raw value."""
        ...

    def to_one_based(self) -> OneBasedPos:
        """Convert to 1-based position."""
        ...

    def as_index(self) -> int:
        """Use as array index."""
        ...

    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __lt__(self, other: ZeroBasedPos) -> bool: ...
    def __le__(self, other: ZeroBasedPos) -> bool: ...
    def __gt__(self, other: ZeroBasedPos) -> bool: ...
    def __ge__(self, other: ZeroBasedPos) -> bool: ...

class OneBasedPos:
    """A 1-based position (human-readable, used by HGVS and VCF)."""

    def __init__(self, pos: int) -> None:
        """Create a new 1-based position.

        Raises:
            ValueError: If pos is 0 (invalid for 1-based coordinates)
        """
        ...

    @property
    def value(self) -> int:
        """Get the raw value."""
        ...

    def to_zero_based(self) -> ZeroBasedPos:
        """Convert to 0-based position."""
        ...

    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __lt__(self, other: OneBasedPos) -> bool: ...
    def __le__(self, other: OneBasedPos) -> bool: ...
    def __gt__(self, other: OneBasedPos) -> bool: ...
    def __ge__(self, other: OneBasedPos) -> bool: ...

# ============================================================================
# Equivalence Classes
# ============================================================================

class EquivalenceLevel(IntEnum):
    """Equivalence level between two variants."""

    Identical = 0
    NormalizedMatch = 1
    AccessionVersionDifference = 2
    NotEquivalent = 3

    def is_equivalent(self) -> bool:
        """Returns true if the variants are considered equivalent."""
        ...

    def description(self) -> str:
        """Get a human-readable description."""
        ...

class EquivalenceResult:
    """Result of an equivalence check."""

    @property
    def level(self) -> EquivalenceLevel:
        """The determined equivalence level."""
        ...

    @property
    def normalized_first(self) -> Optional[str]:
        """The normalized form of the first variant."""
        ...

    @property
    def normalized_second(self) -> Optional[str]:
        """The normalized form of the second variant."""
        ...

    @property
    def notes(self) -> List[str]:
        """Additional notes about the comparison."""
        ...

    def is_equivalent(self) -> bool:
        """Returns true if the variants are considered equivalent."""
        ...

class EquivalenceChecker:
    """Equivalence checker for comparing HGVS variants."""

    def __init__(self, reference_json: Optional[str] = None) -> None:
        """Create a new equivalence checker.

        Args:
            reference_json: Optional path to a transcripts.json file
        """
        ...

    def check(self, v1: HgvsVariant, v2: HgvsVariant) -> EquivalenceResult:
        """Check if two variants are equivalent."""
        ...

    def all_equivalent(self, variants: List[HgvsVariant]) -> bool:
        """Check if multiple variants are all equivalent to each other."""
        ...

# ============================================================================
# Effect Prediction Classes
# ============================================================================

class Consequence(IntEnum):
    """Sequence Ontology consequence term."""

    TranscriptAblation = 0
    SpliceAcceptorVariant = 1
    SpliceDonorVariant = 2
    StopGained = 3
    FrameshiftVariant = 4
    StopLost = 5
    StartLost = 6
    MissenseVariant = 7
    InframeInsertion = 8
    InframeDeletion = 9
    ProteinAlteringVariant = 10
    SpliceRegionVariant = 11
    SynonymousVariant = 12
    StartRetainedVariant = 13
    StopRetainedVariant = 14
    FivePrimeUtrVariant = 15
    ThreePrimeUtrVariant = 16
    IntronVariant = 17
    CodingSequenceVariant = 18

    def so_term(self) -> str:
        """Get the Sequence Ontology term."""
        ...

    def so_id(self) -> str:
        """Get the Sequence Ontology ID."""
        ...

    def impact(self) -> Impact:
        """Get the impact level."""
        ...

    def description(self) -> str:
        """Get a human-readable description."""
        ...

class Impact(IntEnum):
    """Variant impact level (VEP-style)."""

    Modifier = 0
    Low = 1
    Moderate = 2
    High = 3

class ProteinEffect:
    """Protein effect prediction result."""

    @property
    def consequences(self) -> List[Consequence]:
        """Get all applicable consequences."""
        ...

    @property
    def impact(self) -> Impact:
        """Get the highest impact level."""
        ...

    def is_high_impact(self) -> bool:
        """Check if this is a high-impact variant."""
        ...

    def is_protein_altering(self) -> bool:
        """Check if this is a protein-altering variant."""
        ...

class EffectPredictor:
    """Protein effect predictor."""

    def __init__(self) -> None:
        """Create a new effect predictor."""
        ...

    def classify_amino_acid_change(self, ref_aa: str, alt_aa: str, position: int) -> ProteinEffect:
        """Classify an amino acid change.

        Args:
            ref_aa: Reference amino acid (3-letter code, e.g., "Val")
            alt_aa: Alternate amino acid (3-letter code, e.g., "Glu")
            position: Position in protein (1-based)
        """
        ...

    def classify_indel(self, ref_len: int, alt_len: int) -> ProteinEffect:
        """Classify an indel by frame effect."""
        ...

    def classify_splice_variant(self, offset: int) -> ProteinEffect:
        """Classify a splice site variant by distance from splice site."""
        ...

    def classify_utr_variant(self, is_5_prime: bool) -> ProteinEffect:
        """Classify a UTR variant."""
        ...

# ============================================================================
# MAVE Classes
# ============================================================================

class MaveContext:
    """Context for parsing MAVE-HGVS short-form notation."""

    def __init__(self) -> None:
        """Create a new empty context."""
        ...

    @property
    def protein_accession(self) -> Optional[str]: ...
    @property
    def coding_accession(self) -> Optional[str]: ...
    @property
    def genomic_accession(self) -> Optional[str]: ...
    @property
    def gene_symbol(self) -> Optional[str]: ...
    def with_protein_accession(self, accession: str) -> MaveContext:
        """Set the protein sequence accession for p. variants."""
        ...

    def with_coding_accession(self, accession: str) -> MaveContext:
        """Set the coding sequence accession for c. variants."""
        ...

    def with_noncoding_accession(self, accession: str) -> MaveContext:
        """Set the non-coding transcript accession for n. variants."""
        ...

    def with_genomic_accession(self, accession: str) -> MaveContext:
        """Set the genomic sequence accession for g. variants."""
        ...

    def with_gene_symbol(self, symbol: str) -> MaveContext:
        """Set the gene symbol (informational)."""
        ...

    @property
    def noncoding_accession(self) -> Optional[str]: ...
    def has_accessions(self) -> bool:
        """Check if this context has any accessions defined."""
        ...

    def supports_coordinate_type(self, coord_type: str) -> bool:
        """Check if this context can handle a specific coordinate type.

        Args:
            coord_type: Single character coordinate type ('p', 'c', 'n', 'g', etc.)

        Raises:
            ValueError: If coord_type is empty
        """
        ...

# ============================================================================
# Batch Processing Classes
# ============================================================================

class BatchProgress:
    """Progress information for batch processing."""

    @property
    def processed(self) -> int: ...
    @property
    def total(self) -> int: ...
    @property
    def successes(self) -> int: ...
    @property
    def failures(self) -> int: ...
    def percent(self) -> float:
        """Get the percentage complete."""
        ...

class BatchResult:
    """Result of batch processing."""

    def total(self) -> int:
        """Total number of items processed."""
        ...

    def success_count(self) -> int:
        """Number of successful items."""
        ...

    def error_count(self) -> int:
        """Number of failed items."""
        ...

    def success_rate(self) -> float:
        """Success rate as a percentage."""
        ...

    def successes(self) -> List[HgvsVariant]:
        """Get successfully parsed/normalized variants."""
        ...

    def errors(self) -> List[tuple[int, str]]:
        """Get errors as (index, error_message) tuples."""
        ...

class BatchProcessor:
    """Batch processor for parsing and normalizing multiple variants."""

    def __init__(self, reference_json: Optional[str] = None) -> None:
        """Create a new batch processor.

        Args:
            reference_json: Optional path to a transcripts.json file
        """
        ...

    def parse(self, variants: List[str]) -> BatchResult:
        """Parse multiple HGVS strings."""
        ...

    def parse_and_normalize(self, variants: List[str]) -> BatchResult:
        """Parse and normalize multiple HGVS strings."""
        ...

    def parse_with_progress(
        self, variants: List[str], callback: Callable[[BatchProgress], None]
    ) -> BatchResult:
        """Parse multiple HGVS strings with progress callback."""
        ...

# ============================================================================
# Error Handling Classes
# ============================================================================

class ErrorMode(IntEnum):
    """Error handling mode."""

    Strict = 0
    Lenient = 1
    Silent = 2

class ErrorType(IntEnum):
    """Error type for configurable error handling."""

    LowercaseAminoAcid = 0
    MissingVersion = 1
    WrongDashCharacter = 2
    ExtraWhitespace = 3
    ProteinSubstitutionArrow = 4
    PositionZero = 5
    SingleLetterAminoAcid = 6
    WrongQuoteCharacter = 7
    LowercaseAccessionPrefix = 8
    MixedCaseEditType = 9
    OldSubstitutionSyntax = 10
    InvalidUnicodeCharacter = 11
    SwappedPositions = 12
    TrailingAnnotation = 13
    MissingCoordinatePrefix = 14
    OldAlleleFormat = 15
    RefSeqMismatch = 16

class ErrorOverride(IntEnum):
    """Override behavior for a specific error type."""

    Default = 0
    Reject = 1
    WarnCorrect = 2
    SilentCorrect = 3
    Accept = 4

class CorrectionWarning:
    """Warning generated during preprocessing."""

    @property
    def error_type(self) -> ErrorType: ...
    @property
    def message(self) -> str: ...
    @property
    def original(self) -> str: ...
    @property
    def corrected(self) -> str: ...

class ErrorConfig:
    """Error handling configuration."""

    @staticmethod
    def strict() -> ErrorConfig:
        """Create a strict configuration (reject all non-standard input)."""
        ...

    @staticmethod
    def lenient() -> ErrorConfig:
        """Create a lenient configuration (auto-correct with warnings)."""
        ...

    @staticmethod
    def silent() -> ErrorConfig:
        """Create a silent configuration (auto-correct without warnings)."""
        ...

    @property
    def mode(self) -> ErrorMode: ...
    def with_override(self, error_type: ErrorType, override: ErrorOverride) -> ErrorConfig:
        """Add an override for a specific error type."""
        ...

    def should_reject(self, error_type: ErrorType) -> bool:
        """Check if the given error type should be rejected."""
        ...

    def should_correct(self, error_type: ErrorType) -> bool:
        """Check if the given error type should be corrected."""
        ...

    def should_warn(self, error_type: ErrorType) -> bool:
        """Check if the given error type should emit a warning."""
        ...

class ParseResultWithWarnings:
    """Parse result with warnings."""

    @property
    def variant(self) -> HgvsVariant: ...
    @property
    def warnings(self) -> List[CorrectionWarning]: ...
    @property
    def original_input(self) -> str: ...
    @property
    def preprocessed_input(self) -> str: ...
    def had_corrections(self) -> bool:
        """Returns true if there were any corrections made."""
        ...

    def has_warnings(self) -> bool:
        """Returns true if there are any warnings."""
        ...

# ============================================================================
# Backtranslation Classes
# ============================================================================

class CodonChange:
    """A codon change representing a DNA variant."""

    @property
    def ref_codon(self) -> str:
        """Reference codon (e.g., "CTT")."""
        ...

    @property
    def alt_codon(self) -> str:
        """Alternate codon (e.g., "TTT")."""
        ...

    @property
    def changed_positions(self) -> List[int]:
        """Position(s) in codon that changed (1-indexed)."""
        ...

    def num_changes(self) -> int:
        """Number of nucleotide changes."""
        ...

class CodonTable:
    """Codon table (genetic code)."""

    @staticmethod
    def standard() -> CodonTable:
        """Get the standard genetic code."""
        ...

class Backtranslator:
    """Backtranslation engine for converting protein changes to DNA changes."""

    def __init__(self, codon_table: Optional[CodonTable] = None) -> None:
        """Create a new backtranslator with the given codon table."""
        ...

    @staticmethod
    def standard() -> Backtranslator:
        """Create a backtranslator with the standard genetic code."""
        ...

    def backtranslate_substitution(self, ref_aa: str, alt_aa: str) -> List[CodonChange]:
        """Backtranslate an amino acid substitution."""
        ...

    def backtranslate_to_stop(self, ref_aa: str) -> List[CodonChange]:
        """Backtranslate a nonsense mutation (amino acid to stop codon)."""
        ...

    def backtranslate_stop_loss(self, alt_aa: str) -> List[CodonChange]:
        """Backtranslate a stop loss (stop codon to amino acid)."""
        ...

# ============================================================================
# rsID Classes
# ============================================================================

class RsIdResult:
    """Result of rsID lookup."""

    @property
    def rsid(self) -> str: ...
    @property
    def contig(self) -> str: ...
    @property
    def position(self) -> int: ...
    @property
    def reference(self) -> str: ...
    @property
    def alternate(self) -> str: ...
    @property
    def hgvs(self) -> Optional[str]: ...
    @property
    def allele_frequency(self) -> Optional[float]: ...
    @property
    def clinical_significance(self) -> Optional[str]: ...
    def is_snv(self) -> bool:
        """Check if this is a substitution (SNV)."""
        ...

    def is_deletion(self) -> bool:
        """Check if this is a deletion."""
        ...

    def is_insertion(self) -> bool:
        """Check if this is an insertion."""
        ...

    def to_hgvs(self) -> str:
        """Generate simple genomic HGVS notation."""
        ...

class InMemoryRsIdLookup:
    """Simple in-memory rsID lookup for testing."""

    @staticmethod
    def with_test_data() -> InMemoryRsIdLookup:
        """Create with common test variants (BRAF V600E, BRCA1 185delAG)."""
        ...

    def lookup(self, rsid: str) -> List[RsIdResult]:
        """Look up rsID and return matching variants."""
        ...

    def contains(self, rsid: str) -> bool:
        """Check if rsID exists."""
        ...

# ============================================================================
# VCF Classes
# ============================================================================

class VcfRecord:
    """A VCF record."""

    def __init__(
        self, chrom: str, pos: int, reference: str, alternate: str, id: Optional[str] = None
    ) -> None:
        """Create a new VCF record.

        Args:
            chrom: Chromosome name (e.g., "chr7")
            pos: 1-based position
            reference: Reference allele
            alternate: Alternate allele
            id: Optional variant ID (e.g., rsID)
        """
        ...

    @staticmethod
    def snv(chrom: str, pos: int, ref_base: str, alt_base: str) -> VcfRecord:
        """Create a SNV record.

        Args:
            chrom: Chromosome name
            pos: 1-based position
            ref_base: Reference base (single character)
            alt_base: Alternate base (single character)

        Raises:
            ValueError: If ref_base or alt_base is empty
        """
        ...

    @property
    def chrom(self) -> str: ...
    @property
    def pos(self) -> int: ...
    @property
    def id(self) -> Optional[str]: ...
    @property
    def reference(self) -> str: ...
    @property
    def alternate(self) -> Optional[str]:
        """Get the first alternate allele (for simple variants)."""
        ...
    @property
    def alternates(self) -> List[str]:
        """Get all alternate alleles."""
        ...

# ============================================================================
# Prepare Classes
# ============================================================================

class PrepareConfig:
    """Configuration for reference data preparation."""

    def __init__(
        self,
        output_dir: str = "ferro-reference",
        download_transcripts: bool = True,
        download_genome: bool = False,
        download_genome_grch37: bool = False,
        download_refseqgene: bool = False,
        download_lrg: bool = False,
        download_cdot: bool = True,
        skip_existing: bool = True,
        dry_run: bool = False,
    ) -> None: ...
    @property
    def output_dir(self) -> str: ...
    @property
    def download_transcripts(self) -> bool: ...
    @property
    def download_genome(self) -> bool: ...
    @property
    def download_cdot(self) -> bool: ...

class ReferenceManifest:
    """Manifest of prepared reference data."""

    @property
    def prepared_at(self) -> str: ...
    @property
    def transcript_count(self) -> int: ...
    @property
    def transcript_fastas(self) -> List[str]: ...
    @property
    def genome_fasta(self) -> Optional[str]: ...
    @property
    def cdot_json(self) -> Optional[str]: ...
    @property
    def available_prefixes(self) -> List[str]: ...

# ============================================================================
# Reference Classes
# ============================================================================

class GenomeBuild(IntEnum):
    """Genome build (GRCh37 or GRCh38)."""

    GRCh37 = 0
    GRCh38 = 1
    Unknown = 2

class Strand(IntEnum):
    """Strand direction."""

    Plus = 0
    Minus = 1

# ============================================================================
# Coordinate Mapping Classes
# ============================================================================

class CoordinateMapper:
    """Coordinate mapper for converting between HGVS coordinate systems.

    This class provides methods for converting between genomic (g.), coding (c.),
    transcript (n.), and protein (p.) coordinate systems.
    """

    def __init__(self, reference_json: Optional[str] = None) -> None:
        """Create a new coordinate mapper.

        Args:
            reference_json: Optional path to a transcripts.json file for reference data.
                If not provided, uses built-in test data.
        """
        ...

    def c_to_g(
        self, transcript_id: str, cds_position: int, offset: Optional[int] = None
    ) -> Optional[tuple[str, int]]:
        """Convert a CDS position to genomic position.

        Args:
            transcript_id: Transcript accession (e.g., "NM_000088.3")
            cds_position: CDS position (e.g., 100 for c.100)
            offset: Optional intronic offset (e.g., 5 for c.100+5).
                Note: offset=0 is treated as no offset (semantically equivalent).

        Returns:
            Tuple of (chromosome, genomic_position) or None if position is intronic
            without offset support

        Raises:
            ValueError: If transcript not found or has no genomic coordinates
        """
        ...

    def g_to_c(self, transcript_id: str, genomic_position: int) -> tuple[int, Optional[int], bool]:
        """Convert a genomic position to CDS position.

        Args:
            transcript_id: Transcript accession (e.g., "NM_000088.3")
            genomic_position: Genomic position (1-based)

        Returns:
            Tuple of (cds_position, offset, is_utr3). For exonic positions,
            offset will be None. For intronic positions, offset will contain
            the distance from the nearest exon boundary.

        Raises:
            ValueError: If transcript not found, position is outside transcript bounds,
                or conversion fails
        """
        ...

    def c_to_p(self, transcript_id: str, cds_position: int) -> int:
        """Convert a CDS position to protein position.

        Args:
            transcript_id: Transcript accession (e.g., "NM_000088.3")
            cds_position: CDS position (e.g., 100 for c.100)

        Returns:
            Protein position (1-based codon number)

        Raises:
            ValueError: If transcript not found or position is in UTR/intronic
        """
        ...

    def c_to_n(
        self,
        transcript_id: str,
        cds_position: int,
        offset: Optional[int] = None,
        utr3: bool = False,
    ) -> tuple[int, Optional[int]]:
        """Convert a CDS position to transcript position.

        Args:
            transcript_id: Transcript accession (e.g., "NM_000088.3")
            cds_position: CDS position (e.g., 100 for c.100)
            offset: Optional intronic offset
            utr3: Whether this is a 3' UTR position (c.*N notation)

        Returns:
            Tuple of (transcript_position, offset)

        Raises:
            ValueError: If transcript not found
        """
        ...

    def n_to_c(
        self,
        transcript_id: str,
        tx_position: int,
        offset: Optional[int] = None,
        downstream: bool = False,
    ) -> tuple[int, Optional[int], bool]:
        """Convert a transcript position to CDS position.

        Args:
            transcript_id: Transcript accession (e.g., "NM_000088.3")
            tx_position: Transcript position (1-based)
            offset: Optional intronic offset
            downstream: Whether this is a downstream position (n.*100 notation for
                positions past the end of the transcript). Defaults to False.

        Returns:
            Tuple of (cds_position, offset, is_utr3)

        Raises:
            ValueError: If transcript not found or has no CDS
        """
        ...

    def get_strand(self, transcript_id: str) -> Strand:
        """Get the strand of a transcript.

        Args:
            transcript_id: Transcript accession

        Returns:
            Strand (Plus or Minus)
        """
        ...

    def get_chromosome(self, transcript_id: str) -> Optional[str]:
        """Get the chromosome of a transcript.

        Args:
            transcript_id: Transcript accession

        Returns:
            Chromosome name or None if not set
        """
        ...

    def has_transcript(self, transcript_id: str) -> bool:
        """Check if a transcript exists in the reference."""
        ...
