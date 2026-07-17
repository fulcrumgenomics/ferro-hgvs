"""Type stubs for ferro-hgvs Python bindings."""

from collections.abc import Iterable
from enum import IntEnum
from typing import Any, Callable, Literal, overload

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
        ParseError: If the HGVS string cannot be parsed (a subclass of
            ValueError, so ``except ValueError`` still catches it).

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
        ParseError: If the HGVS string cannot be parsed (a subclass of
            ValueError).
        ValueError: If ``direction`` is not one of
            "3prime"/"5prime"/"3"/"5"/"3'"/"5'" (case-insensitive).
        NormalizationError: If normalization fails (a subclass of RuntimeError).

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
        ParseError: If the SPDI string cannot be parsed (a subclass of
            ValueError).
    """
    ...

def hgvs_to_spdi(variant: HgvsVariant) -> SpdiVariant:
    """Convert an HGVS variant to SPDI format.

    Args:
        variant: An HgvsVariant object

    Returns:
        A SpdiVariant object

    Raises:
        ProjectionError: If the conversion fails (a subclass of ValueError and
            RuntimeError).
    """
    ...

def spdi_to_hgvs_variant(spdi: SpdiVariant) -> HgvsVariant:
    """Convert an SPDI variant to HGVS format.

    Args:
        spdi: A SpdiVariant object

    Returns:
        An HgvsVariant object

    Raises:
        ProjectionError: If the conversion fails (a subclass of ValueError and
            RuntimeError).
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
        ParseError: If parsing fails or the context doesn't support the
            coordinate type (a subclass of ValueError).
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

def parse_lenient(hgvs_string: str, config: ErrorConfig | None = None) -> ParseResultWithWarnings:
    """Parse an HGVS string with lenient error handling.

    Args:
        hgvs_string: The HGVS variant description to parse
        config: Optional ErrorConfig (defaults to lenient mode)

    Returns:
        ParseResultWithWarnings containing the parsed variant and any warnings

    Raises:
        ParseError: If the HGVS string cannot be parsed even after corrections
            (a subclass of ValueError).
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
        ParseError: If the rsID cannot be parsed (a subclass of ValueError).
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
        ProjectionError: If conversion fails (a subclass of ValueError and
            RuntimeError).
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
        ReferenceDataError: If preparation fails (a subclass of RuntimeError and
            ValueError).
    """
    ...

def check_reference_data(directory: str) -> ReferenceManifest:
    """Check existing reference data.

    Args:
        directory: Path to reference data directory

    Returns:
        ReferenceManifest if manifest.json exists

    Raises:
        ReferenceDataError: If check fails (a subclass of RuntimeError and
            ValueError).
    """
    ...

def convert_gff(config: ConvertGffConfig) -> ConvertGffReport:
    """Convert a GFF3/GTF annotation into the ``transcripts.json`` format.

    In-process equivalent of ``ferro convert-gff``; both call the same library
    serializer, so the output is byte-identical for the same inputs and flags.
    Use it to build a reference for ``Normalizer(reference_json=...)`` without
    shelling out to the CLI.

    Args:
        config: A ConvertGffConfig describing the inputs and options.

    Returns:
        ConvertGffReport with the loader summary, the written output path (or the
        JSON text when ``output`` is ``None``), and any warnings.

    Raises:
        ValueError: If ``error_mode`` is not recognized, or
            ``emit_genomic_sequences`` is set without a ``fasta``.
        ReferenceDataError: If the annotation or FASTA cannot be read/parsed, or
            (in strict mode) an error diagnostic is recorded.
    """
    ...

def build_transcript(config: BuildTranscriptConfig) -> BuildTranscriptReport:
    """Build a single-construct ``transcripts.json`` from a FASTA + CDS bounds.

    In-process equivalent of ``ferro build-transcript``; both call the same
    library builder, so the output is byte-identical for the same inputs and
    flags. Use it to wrap a synthetic construct's FASTA as a reference for
    ``Normalizer(reference_json=...)`` without shelling out to the CLI.

    Args:
        config: A BuildTranscriptConfig describing the inputs and options.

    Returns:
        BuildTranscriptReport with the transcript id, the written output path (or
        the JSON text when ``output`` is ``None``), and any warnings.

    Raises:
        ValueError: If ``strand`` is not "+"/"-", or the CDS bounds are invalid.
        ReferenceDataError: If the FASTA cannot be read or the contig cannot be
            resolved.
    """
    ...

# ============================================================================
# Core Classes
# ============================================================================

class HgvsVariant:
    """A parsed HGVS variant.

    Attributes:
        variant_type: The type of variant (genomic, coding, non_coding, protein, rna,
            mitochondrial, circular, rna_fusion, genome_ring, allele, null_allele,
            unknown_allele, supernumerary)
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

    @property
    def start(self) -> int | None:
        """Get the 1-based start position of the variant.

        Returns the base position (without intronic offset) for genomic, coding,
        non-coding, RNA, mitochondrial, and circular variants. For single-element
        alleles, delegates to the sub-variant. Returns None for protein variants,
        RNA fusions, genome ring variants, null/unknown alleles, and alleles with
        multiple sub-variants (whose start is ambiguous).

        Note: 5' UTR (``c.-5A>G``) and 3' UTR (``c.*5A>G``) positions are returned
        as raw base values and are indistinguishable from CDS positions at the same
        numeric value.
        """
        ...

    @property
    def end(self) -> int | None:
        """Get the 1-based end position (inclusive) of the variant.

        For point variants, end equals start. For single-element alleles,
        delegates to the sub-variant. Returns None for protein variants, RNA
        fusions, genome ring variants, null/unknown alleles, and alleles with
        multiple sub-variants.
        """
        ...

    @property
    def offset(self) -> int | None:
        """Get the intronic offset of the start position.

        Meaningful for coding (c.), non-coding (n.), and RNA (r.) variants with
        intronic positions. For ``c.93+1G>T``, returns 1. For exonic positions,
        returns None. Always returns None for variant types without intronic
        offsets (genomic, mitochondrial, circular, protein, fusion, allele).
        """
        ...

    @property
    def substitution_bases(self) -> tuple[str, str] | None:
        """Get the substitution reference and alternative bases.

        Returns a tuple (ref_base, alt_base) of single-character strings for
        substitution edits, e.g., ("A", "G") for A>G. Returns None for
        non-substitution edits.
        """
        ...

    @property
    def num_variants(self) -> int:
        """Get the number of sub-variants.

        Returns 1 for simple variants, N for alleles with N sub-variants.
        """
        ...

    @property
    def indel_length(self) -> int | None:
        """Get the net indel length (bases gained or lost).

        - Substitution/inversion/identity: 0
        - Deletion: negative (e.g., -3 for a 3bp deletion)
        - Insertion: positive (length of inserted sequence)
        - Delins: inserted_length - deleted_span
        - Duplication: positive (span of duplicated region)

        Returns None if the length cannot be determined.
        """
        ...

    def variants(self) -> list[HgvsVariant]:
        """Get sub-variants as a list.

        For alleles, returns the constituent variants. For simple variants,
        returns a single-element list containing self.
        """
        ...

    def is_identity(self) -> bool:
        """Check if this is an identity (no-change) variant."""
        ...

    def is_frameshift(self) -> bool:
        """Check if this variant causes a frameshift (indel_length % 3 != 0)."""
        ...

    def normalize(self, direction: str = "3prime") -> HgvsVariant:
        """Normalize this variant.

        Args:
            direction: Shuffle direction - "3prime" (default) or "5prime".

        Raises:
            ValueError: If ``direction`` is not one of
                "3prime"/"5prime"/"3"/"5"/"3'"/"5'" (case-insensitive).
        """
        ...

    def to_dict(self) -> dict[str, Any]:
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

    def __init__(self, reference_json: str | None = None, direction: str = "3prime") -> None:
        """Create a new normalizer.

        Args:
            reference_json: Optional path to a transcripts.json file
            direction: Shuffle direction - "3prime" (default) or "5prime"

        Raises:
            ValueError: If ``direction`` is not one of
                "3prime"/"5prime"/"3"/"5"/"3'"/"5'" (case-insensitive).
        """
        ...

    @staticmethod
    def from_manifest(manifest_path: str, direction: str = "3prime") -> Normalizer:
        """Create a normalizer from a reference manifest written by ``ferro prepare``.

        Args:
            manifest_path: Path to a manifest.json file (typically inside a
                directory produced by ``ferro prepare``).
            direction: Shuffle direction - "3prime" (default) or "5prime".

        Returns:
            A Normalizer backed by a MultiFastaProvider.

        Raises:
            ValueError: If ``direction`` is not one of
                "3prime"/"5prime"/"3"/"5"/"3'"/"5'" (case-insensitive).
            ReferenceDataError: If the manifest cannot be loaded (a subclass of
                RuntimeError and ValueError).
        """
        ...

    def has_genomic_data(self) -> bool:
        """Return True if the backing reference provides genomic sequence data.

        A build with no genomic data (built-in test data or a transcripts.json
        reference) cannot perform genome-dependent normalization; use
        ``Normalizer.from_manifest(...)`` for full capability.
        """
        ...

    def has_protein_data(self) -> bool:
        """Return True if the backing reference provides protein sequence data."""
        ...

    def reference_summary(self) -> dict[str, Any]:
        """Summarize the reference backend's capabilities.

        Returns a dict with keys ``provider_kind`` (one of ``"test_data"``,
        ``"json"``, or ``"manifest"``), ``has_genomic_data``, and
        ``has_protein_data``.
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

    def to_dict(self) -> dict[str, Any]:
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
    def normalized_first(self) -> str | None:
        """The normalized form of the first variant."""
        ...

    @property
    def normalized_second(self) -> str | None:
        """The normalized form of the second variant."""
        ...

    @property
    def notes(self) -> list[str]:
        """Additional notes about the comparison."""
        ...

    def is_equivalent(self) -> bool:
        """Returns true if the variants are considered equivalent."""
        ...

class EquivalenceChecker:
    """Equivalence checker for comparing HGVS variants."""

    def __init__(self, reference_json: str | None = None) -> None:
        """Create a new equivalence checker.

        Args:
            reference_json: Optional path to a transcripts.json file
        """
        ...

    @staticmethod
    def from_manifest(manifest_path: str) -> EquivalenceChecker:
        """Create an equivalence checker from a reference manifest.

        Args:
            manifest_path: Path to a manifest.json file produced by ``ferro prepare``.

        Returns:
            An EquivalenceChecker backed by a MultiFastaProvider.

        Raises:
            ReferenceDataError: If the manifest cannot be loaded (a subclass of
                RuntimeError and ValueError).
        """
        ...

    def has_genomic_data(self) -> bool:
        """Return True if the backing reference provides genomic sequence data.

        A checker with no genomic data (built-in test data or a transcript-only
        ``reference_json``) has limited genome-dependent capability; use
        ``EquivalenceChecker.from_manifest(...)`` for full capability.
        """
        ...

    def has_protein_data(self) -> bool:
        """Return True if the backing reference provides protein sequence data."""
        ...

    def reference_summary(self) -> dict[str, Any]:
        """Summarize the reference backend's capabilities.

        Returns a dict with keys ``provider_kind`` (one of ``"test_data"``,
        ``"json"``, or ``"manifest"``), ``has_genomic_data``, and
        ``has_protein_data``.
        """
        ...

    def check(self, v1: HgvsVariant, v2: HgvsVariant) -> EquivalenceResult:
        """Check if two variants are equivalent.

        Raises:
            EquivalenceError: If the check fails, e.g. an input cannot be
                normalized (a subclass of RuntimeError).
        """
        ...

    def all_equivalent(self, variants: list[HgvsVariant]) -> bool:
        """Check if multiple variants are all equivalent to each other.

        Raises:
            EquivalenceError: If the check fails, e.g. an input cannot be
                normalized (a subclass of RuntimeError).
        """
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
    def consequences(self) -> list[Consequence]:
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

    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...

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
    def protein_accession(self) -> str | None: ...
    @property
    def coding_accession(self) -> str | None: ...
    @property
    def genomic_accession(self) -> str | None: ...
    @property
    def gene_symbol(self) -> str | None: ...
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
    def noncoding_accession(self) -> str | None: ...
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

    def successes(self) -> list[HgvsVariant]:
        """Get successfully parsed/normalized variants."""
        ...

    def errors(self) -> list[tuple[int, str]]:
        """Get errors as (index, error_message) tuples."""
        ...

class BatchProcessor:
    """Batch processor for parsing and normalizing multiple variants."""

    def __init__(self, reference_json: str | None = None) -> None:
        """Create a new batch processor.

        Args:
            reference_json: Optional path to a transcripts.json file
        """
        ...

    @staticmethod
    def from_manifest(manifest_path: str) -> BatchProcessor:
        """Create a batch processor from a reference manifest.

        Args:
            manifest_path: Path to a manifest.json file produced by ``ferro prepare``.

        Returns:
            A BatchProcessor backed by a MultiFastaProvider.

        Raises:
            ReferenceDataError: If the manifest cannot be loaded (a subclass of
                RuntimeError and ValueError).
        """
        ...

    def has_genomic_data(self) -> bool:
        """Return True if the backing reference provides genomic sequence data.

        A processor with no genomic data (built-in test data or a transcript-only
        ``reference_json``) has limited genome-dependent capability; use
        ``BatchProcessor.from_manifest(...)`` for full capability.
        """
        ...

    def has_protein_data(self) -> bool:
        """Return True if the backing reference provides protein sequence data."""
        ...

    def reference_summary(self) -> dict[str, Any]:
        """Summarize the reference backend's capabilities.

        Returns a dict with keys ``provider_kind`` (one of ``"test_data"``,
        ``"json"``, or ``"manifest"``), ``has_genomic_data``, and
        ``has_protein_data``.
        """
        ...

    def parse(self, variants: list[str], workers: int = 0) -> BatchResult:
        """Parse multiple HGVS strings in parallel (GIL released).

        The returned BatchResult preserves the order of the input variants
        regardless of the worker count used.

        Args:
            variants: HGVS strings to parse.
            workers: Worker threads. 0 (default) uses all cores; 1 is serial; N uses N threads.
        """
        ...

    def parse_and_normalize(self, variants: list[str], workers: int = 0) -> BatchResult:
        """Parse and normalize multiple HGVS strings in parallel (GIL released).

        The returned BatchResult preserves the order of the input variants
        regardless of the worker count used.

        Args:
            variants: HGVS strings to parse and normalize.
            workers: Worker threads. 0 (default) uses all cores; 1 is serial; N uses N threads.
        """
        ...

    def parse_with_progress(
        self, variants: list[str], callback: Callable[[BatchProgress], None]
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
    # Discriminant 5 was `PositionZero` (W4002), retired in issue #269.
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
    DeprecatedStopCodonStar = 17
    DeprecatedStopCodonX = 18
    DeprecatedFrameshiftStar = 19
    DeprecatedFrameshiftX = 20
    DelSizeSuffix = 21
    EmptyDelinsInsert = 22
    RedundantRepeatLabel = 23
    SinglePositionRange = 24
    DeprecatedIvsNotation = 25
    DeprecatedConSyntax = 26
    LengthMismatch = 27
    AlleleFractionAnnotation = 28
    ClinVarProseMultiAllelic = 29
    RnaThymineCanonicalized = 30
    ProteinBracketedAaInsertion = 31
    PositionPastEnd = 32
    VariantExceedsReference = 33
    NonSpecMosaicForm = 34
    OverlapConflictingEdits = 35
    InitiatorMetCanonicalization = 36
    DupSizeSuffix = 37
    DupExplicitSeq = 38
    DelExplicitSeq = 39
    NonConformantBracketCardinality = 40
    UnresolvableCentromere = 41
    TranscriptFlankNotDescribable = 42
    IntronicOnBareTranscript = 43
    IncompleteCdsStartReference = 44

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
    def warnings(self) -> list[CorrectionWarning]: ...
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
    def changed_positions(self) -> list[int]:
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

    def __init__(self, codon_table: CodonTable | None = None) -> None:
        """Create a new backtranslator with the given codon table."""
        ...

    @staticmethod
    def standard() -> Backtranslator:
        """Create a backtranslator with the standard genetic code."""
        ...

    def backtranslate_substitution(self, ref_aa: str, alt_aa: str) -> list[CodonChange]:
        """Backtranslate an amino acid substitution."""
        ...

    def backtranslate_to_stop(self, ref_aa: str) -> list[CodonChange]:
        """Backtranslate a nonsense mutation (amino acid to stop codon)."""
        ...

    def backtranslate_stop_loss(self, alt_aa: str) -> list[CodonChange]:
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
    def hgvs(self) -> str | None: ...
    @property
    def allele_frequency(self) -> float | None: ...
    @property
    def clinical_significance(self) -> str | None: ...
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

    def lookup(self, rsid: str) -> list[RsIdResult]:
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
        self, chrom: str, pos: int, reference: str, alternate: str, id: str | None = None
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
    def id(self) -> str | None: ...
    @property
    def reference(self) -> str: ...
    @property
    def alternate(self) -> str | None:
        """Get the first alternate allele (for simple variants)."""
        ...
    @property
    def alternates(self) -> list[str]:
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
        download_cdot_grch37: bool = False,
        download_ensembl: bool = False,
    ) -> None: ...
    @property
    def output_dir(self) -> str: ...
    @property
    def download_transcripts(self) -> bool: ...
    @property
    def download_genome(self) -> bool: ...
    @property
    def download_cdot(self) -> bool: ...
    @property
    def download_cdot_grch37(self) -> bool: ...
    @property
    def download_ensembl(self) -> bool: ...

class ReferenceManifest:
    """Manifest of prepared reference data."""

    @property
    def prepared_at(self) -> str: ...
    @property
    def transcript_count(self) -> int: ...
    @property
    def transcript_fastas(self) -> list[str]: ...
    @property
    def genome_fasta(self) -> str | None: ...
    @property
    def cdot_json(self) -> str | None: ...
    @property
    def cdot_grch37_json(self) -> str | None: ...
    @property
    def available_prefixes(self) -> list[str]: ...

# ============================================================================
# Convert-GFF Classes
# ============================================================================

class ConvertGffConfig:
    """Configuration for :func:`convert_gff`, mirroring the ``ferro convert-gff``
    CLI flags."""

    def __init__(
        self,
        gff: str,
        fasta: str | None = None,
        output: str | None = None,
        build: str = "GRCh38",
        mane_only: bool = False,
        transcripts: list[str] | str | None = None,
        genes: list[str] | str | None = None,
        error_mode: str = "lenient",
        validate_fasta: bool = True,
        emit_genomic_sequences: bool = False,
        diagnostics_json: str | None = None,
    ) -> None:
        """Create a convert-gff configuration.

        Args:
            gff: Path to the input GFF3/GTF annotation file.
            fasta: Optional reference FASTA, used to extract exonic sequences and
                (with ``emit_genomic_sequences``) the embedded contig sequences.
            output: Optional output path for the ``transcripts.json``. If ``None``
                (default), the JSON is returned in
                ``ConvertGffReport.transcripts_json`` instead of written to disk.
            build: Genome build recorded in the output ("GRCh38" or "GRCh37").
            mane_only: Emit only MANE Select / MANE Plus Clinical transcripts.
            transcripts: Optional transcript-ID filter — a list of IDs or a single
                comma-separated string.
            genes: Optional gene-symbol filter — a list of symbols or a single
                comma-separated string.
            error_mode: "lenient" (default), "strict", or "silent".
            validate_fasta: Run CDS-length / start-codon FASTA validation when a
                FASTA is supplied (default True; inverse of ``--no-validate-fasta``).
            emit_genomic_sequences: Embed referenced contig sequences so the output
                is genome-capable. Requires ``fasta``.
            diagnostics_json: Optional path to write the sampled loader diagnostics.
        """
        ...

class ConvertGffReport:
    """Result of a :func:`convert_gff` run."""

    @property
    def summary(self) -> str:
        """One-line loader summary (records read, transcripts built, diagnostics)."""
        ...
    @property
    def transcripts_json(self) -> str | None:
        """The serialized ``transcripts.json`` text, present only when the config's
        ``output`` was ``None`` (otherwise the JSON was written to that path and
        this is ``None``)."""
        ...
    @property
    def output_path(self) -> str | None:
        """The path the ``transcripts.json`` was written to, or ``None`` if it was
        returned in ``transcripts_json`` instead."""
        ...
    @property
    def transcript_count(self) -> int:
        """Number of transcripts emitted."""
        ...
    @property
    def emitted_genomic_bytes(self) -> int:
        """Total bytes of genomic sequence embedded under ``genomic_sequences``
        (0 when ``emit_genomic_sequences`` was off or nothing was placed)."""
        ...
    @property
    def warnings(self) -> list[str]:
        """Non-fatal warnings raised during conversion (also emitted as
        ``UserWarning``)."""
        ...

# ============================================================================
# Build-Transcript Classes
# ============================================================================

class BuildTranscriptConfig:
    """Configuration for :func:`build_transcript`, mirroring the
    ``ferro build-transcript`` CLI flags."""

    def __init__(
        self,
        fasta: str,
        cds_start: int,
        cds_end: int,
        output: str | None = None,
        id: str | None = None,
        strand: str = "+",
        contig: str | None = None,
        gene: str | None = None,
        genome_build: str = "GRCh38",
        emit_genomic_sequences: bool = False,
    ) -> None:
        """Create a build-transcript configuration.

        Args:
            fasta: Path to the FASTA (indexed or plain; the index is built on the
                fly if absent).
            cds_start: CDS start position (1-based inclusive, transcript coordinates).
            cds_end: CDS end position (1-based inclusive, transcript coordinates).
            output: Optional output path for the ``transcripts.json``. If ``None``
                (default), the JSON is returned in
                ``BuildTranscriptReport.transcripts_json`` instead of written to disk.
            id: Transcript ID; defaults to the FASTA contig name.
            strand: "+" (default) or "-". Minus-strand sequences are reverse-
                complemented, so the CDS positions are relative to that sequence.
            contig: Contig to use when the FASTA has multiple contigs.
            gene: Optional gene symbol to embed in the transcript record.
            genome_build: Genome build name embedded verbatim (default "GRCh38").
            emit_genomic_sequences: Embed the contig's forward bytes so the output
                is genome-capable.
        """
        ...

class BuildTranscriptReport:
    """Result of a :func:`build_transcript` run."""

    @property
    def transcript_id(self) -> str:
        """The transcript ID that was emitted (the configured id or contig name)."""
        ...
    @property
    def transcripts_json(self) -> str | None:
        """The serialized ``transcripts.json`` text, present only when the config's
        ``output`` was ``None`` (otherwise the JSON was written to that path and
        this is ``None``)."""
        ...
    @property
    def output_path(self) -> str | None:
        """The path the ``transcripts.json`` was written to, or ``None`` if it was
        returned in ``transcripts_json`` instead."""
        ...
    @property
    def emitted_genomic_bytes(self) -> int:
        """Total bytes of genomic sequence embedded under ``genomic_sequences``
        (0 when ``emit_genomic_sequences`` was off)."""
        ...
    @property
    def warnings(self) -> list[str]:
        """Non-fatal warnings raised during the build (also emitted as
        ``UserWarning``)."""
        ...

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

    def __init__(self, reference_json: str | None = None) -> None:
        """Create a new coordinate mapper.

        Args:
            reference_json: Optional path to a transcripts.json file for reference data.
                If not provided, uses built-in test data.
        """
        ...

    @staticmethod
    def from_manifest(manifest_path: str) -> CoordinateMapper:
        """Create a coordinate mapper from a reference manifest.

        Args:
            manifest_path: Path to a manifest.json file produced by ``ferro prepare``.

        Returns:
            A CoordinateMapper backed by a MultiFastaProvider.

        Raises:
            ReferenceDataError: If the manifest cannot be loaded (a subclass of
                RuntimeError and ValueError).
        """
        ...

    def has_genomic_data(self) -> bool:
        """Return True if the backing reference provides genomic sequence data.

        A mapper with no genomic data (built-in test data or a transcript-only
        ``reference_json``) cannot resolve genomic coordinates; use
        ``CoordinateMapper.from_manifest(...)`` for full capability.
        """
        ...

    def has_protein_data(self) -> bool:
        """Return True if the backing reference provides protein sequence data."""
        ...

    def reference_summary(self) -> dict[str, Any]:
        """Summarize the reference backend's capabilities.

        Returns a dict with keys ``provider_kind`` (one of ``"test_data"``,
        ``"json"``, or ``"manifest"``), ``has_genomic_data``, and
        ``has_protein_data``.
        """
        ...

    def c_to_g(
        self, transcript_id: str, cds_position: int, offset: int | None = None
    ) -> tuple[str, int] | None:
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
            ReferenceDataError: If the transcript is not found.
            ProjectionError: If it has no genomic coordinates (both subclass
                ValueError).
        """
        ...

    def g_to_c(self, transcript_id: str, genomic_position: int) -> tuple[int, int | None, bool]:
        """Convert a genomic position to CDS position.

        Args:
            transcript_id: Transcript accession (e.g., "NM_000088.3")
            genomic_position: Genomic position (1-based)

        Returns:
            Tuple of (cds_position, offset, is_utr3). For exonic positions,
            offset will be None. For intronic positions, offset will contain
            the distance from the nearest exon boundary.

        Raises:
            ReferenceDataError: If the transcript is not found.
            ProjectionError: If the position is outside transcript bounds or
                conversion fails (both subclass ValueError).
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
            ReferenceDataError: If the transcript is not found.
            ProjectionError: If the position is in a UTR or intron (both
                subclass ValueError).
        """
        ...

    def c_to_n(
        self,
        transcript_id: str,
        cds_position: int,
        offset: int | None = None,
        utr3: bool = False,
    ) -> tuple[int, int | None]:
        """Convert a CDS position to transcript position.

        Args:
            transcript_id: Transcript accession (e.g., "NM_000088.3")
            cds_position: CDS position (e.g., 100 for c.100)
            offset: Optional intronic offset
            utr3: Whether this is a 3' UTR position (c.*N notation)

        Returns:
            Tuple of (transcript_position, offset)

        Raises:
            ReferenceDataError: If the transcript is not found (a subclass of
                ValueError).
        """
        ...

    def n_to_c(
        self,
        transcript_id: str,
        tx_position: int,
        offset: int | None = None,
        downstream: bool = False,
    ) -> tuple[int, int | None, bool]:
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
            ReferenceDataError: If the transcript is not found.
            ProjectionError: If it has no CDS (both subclass ValueError).
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

    def get_chromosome(self, transcript_id: str) -> str | None:
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

# ============================================================================
# Variant Projection
# ============================================================================

class NormalizationWarning:
    """A diagnostic emitted during normalization.

    Carried by ``VariantProjection.warnings`` and
    ``NormalizeResultWithWarnings.warnings``.
    """

    @property
    def code(self) -> str:
        """Stable SVA code, e.g. ``"REFSEQ_MISMATCH"`` / ``"POSITION_PAST_END"``."""
        ...

    @property
    def message(self) -> str:
        """Human-readable description of the warning."""
        ...

class VariantProjection:
    """The result of projecting a g. variant onto a transcript."""

    @property
    def g_name(self) -> str | None:
        """The normalized g. variant as an HGVS string.

        None when no genomic representation is available (e.g. a bare-NM_
        coding input with no genome alignment; see #498).
        """
        ...

    @property
    def c_name(self) -> str | None:
        """The c./n. variant as an HGVS string, or None when projection skipped."""
        ...

    @property
    def n_name(self) -> str | None:
        """The n. (transcript-relative) variant as an HGVS string.

        Populated for both coding transcripts (derived genome-free from the c.
        form) and non-coding transcripts; None when no transcript coordinate is
        available (e.g. an empty allele).
        """
        ...

    @property
    def p_name(self) -> str | None:
        """The p. variant as an HGVS string, or None for intronic/UTR/non-coding variants.

        Rendered in the projector's configured style (#1050; default `Ter` /
        three-letter). Note that one-letter mode always spells the stop codon as
        `*` (there is no one-letter `Ter`), overriding `protein_stop="ter"`.
        """
        ...

    def p_name_styled(
        self,
        protein_stop: str | None = None,
        amino_acid_code: str | None = None,
    ) -> str | None:
        """The p. variant rendered with an explicit style (#1050).

        Overrides the projector's configured default per call. ``None`` for
        either argument keeps the projector's stored value for that axis, so a
        single axis can be overridden. Returns ``None`` in exactly the cases
        ``p_name`` does.

        Args:
            protein_stop: "ter" or "star", or None to keep the projector default.
            amino_acid_code: "three" or "one", or None to keep the default.

        Raises:
            ValueError: If a supplied value is not recognized.
        """
        ...

    @property
    def r_name(self) -> str | None:
        """The predicted r. (RNA) consequence as an HGVS string, or None.

        CDS-relative numbering (matches c.); None when not representable as a
        concrete exonic RNA edit (no transcript sequence, non-c./n. input, or an
        unresolved payload).
        """
        ...

    @property
    def warnings(self) -> list[NormalizationWarning]:
        """Warnings emitted while normalizing the input for this projection.

        For example an auto-corrected reference-sequence mismatch. Empty when the
        input normalized cleanly. Mirrors ``Normalizer.normalize_with_warnings``,
        making the same diagnostics reachable off a projection.
        """
        ...

    def has_warnings(self) -> bool:
        """True if any normalization warnings were emitted for this projection."""
        ...

    @property
    def transcript_id(self) -> str:
        """The transcript accession (e.g., "NM_000088.3")."""
        ...

    @property
    def gene_symbol(self) -> str | None:
        """Gene symbol, if available."""
        ...

    @property
    def is_frameshift(self) -> bool:
        """True if the variant causes a frameshift."""
        ...

    @property
    def is_intronic(self) -> bool:
        """True if the variant falls entirely within an intron."""
        ...

    @property
    def is_utr(self) -> bool:
        """True if the variant falls in a UTR."""
        ...

    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...

class VariantProjector:
    """Projects g. HGVS variants onto transcripts to produce c./p. equivalents."""

    def __init__(
        self,
        reference_json: str | None = None,
        direction: str = "3prime",
        assembly: str | None = None,
        protein_stop: str = "ter",
        amino_acid_code: str = "three",
    ) -> None:
        """Create a variant projector.

        Args:
            reference_json: Path to a transcripts.json file. If None, uses
                built-in test data (limited; not for production).
            direction: Shuffle direction passed to the internal normalizer
                ("3prime" or "5prime").
            assembly: Optional genome-build override ("GRCh37"/"GRCh38", or the
                aliases "hg19"/"hg38") for build-agnostic inputs. A bare NG_/LRG_
                input carries no build; this fills one in. An input whose
                accession already encodes a build (NC_*.10/.11) keeps it.
            protein_stop: Stop-codon spelling for rendered p. names — "ter"
                (default) or "star" (`*`).
            amino_acid_code: Amino-acid code width for rendered p. names —
                "three" (default) or "one".

        Raises:
            ValueError: If ``direction``, ``assembly``, ``protein_stop``, or
                ``amino_acid_code`` is not a recognized value.
        """
        ...

    @staticmethod
    def from_manifest(
        manifest_path: str,
        direction: str = "3prime",
        assembly: str | None = None,
        protein_stop: str = "ter",
        amino_acid_code: str = "three",
    ) -> "VariantProjector":
        """Create a projector from a ferro-prepare manifest.

        Args:
            manifest_path: Path to manifest.json produced by `ferro prepare`.
            direction: Shuffle direction ("3prime" or "5prime").
            assembly: Optional genome-build override ("GRCh37"/"GRCh38", or the
                aliases "hg19"/"hg38") for build-agnostic inputs.
            protein_stop: Stop-codon spelling for rendered p. names — "ter"
                (default) or "star" (`*`).
            amino_acid_code: Amino-acid code width for rendered p. names —
                "three" (default) or "one".

        Returns:
            A VariantProjector backed by MultiFastaProvider with cdot data.

        Raises:
            ValueError: If ``direction``, ``assembly``, ``protein_stop``, or
                ``amino_acid_code`` is not a recognized value.
            RuntimeError: If the manifest cannot be loaded.
        """
        ...

    def project(self, hgvs_string: str, transcript: str) -> VariantProjection:
        """Project a g. HGVS string onto a target transcript.

        Args:
            hgvs_string: A g. HGVS variant string (e.g., "NC_000017.11:g.48275363C>T").
            transcript: Transcript accession (e.g., "NM_000088.3").

        Returns:
            VariantProjection with g./c./p. representations and flags.

        Raises:
            ProjectionError: If parsing or projection fails (a subclass of
                ValueError and RuntimeError).
        """
        ...

    def project_all(self, hgvs_string: str) -> list[VariantProjection]:
        """Parse, normalize, and project a g. HGVS string onto the curated set of overlapping transcripts.

        Returns the curated enumerated set (#656) rather than every overlapping record: superseded
        transcript versions are collapsed (only the highest version per base accession is kept) and
        predicted XM_/XR_ models are dropped when a curated NM_/NR_ transcript covers the same locus;
        predicted models are kept only when they are the sole coverage.

        Returns projections in clinical priority order (MANE Select first, then Plus Clinical,
        then canonical, then longest CDS). Individual transcript failures are silently skipped.

        Args:
            hgvs_string: A g. HGVS variant string.

        Returns:
            List of VariantProjection objects, one per curated overlapping transcript.
            Empty list when no transcripts overlap the variant.

        Raises:
            ProjectionError: If parsing or normalization fails (a subclass of
                ValueError and RuntimeError).
        """
        ...

    def project_variant(self, variant: HgvsVariant, transcript: str) -> VariantProjection:
        """Normalize and project an already-parsed g. variant onto a transcript.

        Equivalent to ``project(str(variant), transcript)`` but skips the re-parse.
        Useful when the caller already holds an HgvsVariant (e.g. from
        ``ferro_hgvs.parse(...)``) and wants to project without going back through a
        string.

        Args:
            variant: An HgvsVariant (g.). Will be normalized before projection.
            transcript: Transcript accession (e.g., "NM_000088.3").

        Returns:
            VariantProjection with g./c./p. representations and flags.

        Raises:
            ProjectionError: If normalization or projection fails (a subclass of
                ValueError and RuntimeError).
        """
        ...

    def project_variant_all(self, variant: HgvsVariant) -> list[VariantProjection]:
        """Normalize and project an already-parsed g. variant onto the curated set of overlapping transcripts.

        Equivalent to ``project_all(str(variant))`` but skips the re-parse; see ``project_all`` for the
        curated enumeration policy (#656).

        Args:
            variant: An HgvsVariant (g.). Will be normalized before projection.

        Returns:
            List of VariantProjection objects in clinical priority order.
            Empty list when no transcripts overlap the variant.

        Raises:
            ProjectionError: If normalization or projection fails (a subclass of
                ValueError and RuntimeError).
        """
        ...

    def project_normalized(self, variant: HgvsVariant, transcript: str) -> VariantProjection:
        """Project an already-normalized g. variant onto a single transcript, skipping normalization.

        More efficient than `project` when the caller has pre-normalized the variant and wants to
        project it against many transcripts.

        Warning: passing a non-normalized variant will produce coordinates that may not match
        other tools' canonical form.

        Args:
            variant: An already-normalized HgvsVariant (g.).
            transcript: Transcript accession (e.g., "NM_000088.3").

        Returns:
            VariantProjection with g./c./p. representations and flags.

        Raises:
            ProjectionError: If projection fails (a subclass of ValueError and
                RuntimeError).
        """
        ...

    def project_normalized_all(self, variant: HgvsVariant) -> list[VariantProjection]:
        """Project an already-normalized g. variant onto the curated set of overlapping transcripts, skipping re-normalization.

        Returns the curated enumerated set (#656) rather than every overlapping record: superseded
        transcript versions are collapsed (only the highest version per base accession is kept) and
        predicted XM_/XR_ models are dropped when a curated NM_/NR_ transcript covers the same locus;
        predicted models are kept only when they are the sole coverage.

        Callers that pre-normalize once and then fan-out across transcripts should use this method
        to avoid repeated normalization overhead.

        Args:
            variant: An already-normalized HgvsVariant (g.).

        Returns:
            List of VariantProjection objects in clinical priority order.
            Empty list when no transcripts overlap the variant.

        Raises:
            ProjectionError: If projection fails (a subclass of ValueError and
                RuntimeError).
        """
        ...

    @overload
    def project_many(
        self, hgvs_strings: list[str], return_exceptions: Literal[False] = False
    ) -> list[list[VariantProjection]]:
        """Batched parse + normalize + project_all over a list of g. HGVS strings.

        Equivalent to ``[self.project_all(s) for s in hgvs_strings]``, but takes a
        single Python→Rust call, releases the GIL for the entire batch, and reuses
        the projector's internal transcript / ref-protein caches across all inputs.
        Use this for fan-out workloads (thousands+ of variants).

        Args:
            hgvs_strings: List of g. HGVS variant strings.
            return_exceptions: When False (default) the first failing input aborts
                the batch by raising ``ProjectionError``. When True nothing is
                raised: each output element is either that input's list of
                projections or the ``ProjectionError`` for its failure, so one bad
                input does not discard the rest.

        Returns:
            List of result lists — one inner list per input, in the same order.
            Each inner list holds the per-transcript projections for that variant
            in clinical priority order. With ``return_exceptions=True`` an element
            may instead be the ``ProjectionError`` for that input.

        Raises:
            ProjectionError: On the first parse / projection error, unless
                ``return_exceptions=True`` (then failures are returned in place).
        """
        ...

    @overload
    def project_many(
        self, hgvs_strings: list[str], return_exceptions: Literal[True]
    ) -> list[list[VariantProjection] | ProjectionError]: ...
    def project_to_genomic(self, variant: HgvsVariant, normalize: bool = True) -> HgvsVariant:
        """Project a transcript-coordinate variant (c./n./r.) onto its parent
        genomic reference and return a Genome-kind HgvsVariant.

        The output g. variant carries the parent NG/NC accession from the input's
        ``Accession.genomic_context``.

        By default (``normalize=True``) the result is the projector-normalized
        genomic form — what most callers want; a Genome input is canonicalized
        too. The normalizer follows this projector's configured shuffle
        direction (``VariantProjector(direction=...)``): with the default
        ``direction="3prime"`` that is the spec-canonical, 3'-shifted form (e.g.
        a non-3'-most ``g.1003del`` in a poly-A run → ``g.1007del``); a
        ``direction="5prime"`` projector instead returns the 5'-anchored form.
        Pass ``normalize=False`` for the raw pivot, which intentionally does not
        normalize its input (#785): a non-canonical input then yields a
        non-canonical genomic output, and a Genome input passes through
        unchanged (idempotent).

        Limitations:
            - ``Allele`` inputs are rejected pending #328.
            - Plus-strand ``Base::U`` in r. inputs is forwarded verbatim into the
              g. output; callers should pre-translate U→T.
            - Intronic offsets on non-coding transcripts are rejected pending
              #332.

        Args:
            variant: A c./n./r./g. HgvsVariant. Transcript-coordinate variants
                must carry a genomic_context (NG/NC parent) on their Accession.
            normalize: When True (default), return the normalized genomic form;
                when False, return the raw un-normalized pivot.

        Returns:
            The Genome-kind HgvsVariant for the requested projection.

        Raises:
            ProjectionError: If the input lacks a parent reference, carries an
                unknown (`?`) position, or is otherwise unsupported (a subclass
                of ValueError and RuntimeError).
        """
        ...

    @overload
    def project_normalized_many(
        self, variants: list[HgvsVariant], return_exceptions: Literal[False] = False
    ) -> list[list[VariantProjection]]:
        """Batched ``project_normalized_all`` over a list of already-normalized
        g. variants.

        Same batching semantics as ``project_many`` but skips renormalization on
        each input.

        Args:
            variants: List of HgvsVariant objects (must already be normalized).
            return_exceptions: When False (default) the first failing input aborts
                the batch by raising ``ProjectionError``. When True nothing is
                raised: each output element is either that input's list of
                projections or the ``ProjectionError`` for its failure.

        Returns:
            List of result lists, one per input. With ``return_exceptions=True`` an
            element may instead be the ``ProjectionError`` for that input.

        Raises:
            ProjectionError: On the first projection error, unless
                ``return_exceptions=True`` (then failures are returned in place).
        """
        ...

    @overload
    def project_normalized_many(
        self, variants: list[HgvsVariant], return_exceptions: Literal[True]
    ) -> list[list[VariantProjection] | ProjectionError]: ...
    def has_genomic_data(self) -> bool:
        """Return True if the backing reference provides genomic sequence data.

        A projector with no genomic data (built-in test data or a transcript-only
        ``reference_json``) cannot ``project_to_genomic``; use
        ``VariantProjector.from_manifest(...)`` for full capability.
        """
        ...

    def has_protein_data(self) -> bool:
        """Return True if the backing reference provides protein sequence data."""
        ...

    def reference_summary(self) -> dict[str, Any]:
        """Summarize the reference backend's capabilities.

        Returns a dict with keys ``provider_kind`` (one of ``"test_data"``,
        ``"json"``, or ``"manifest"``), ``has_genomic_data``, and
        ``has_protein_data``.
        """
        ...

# ============================================================================
# Exception hierarchy
# ============================================================================

class FerroError(Exception):
    """Base class for all ferro-hgvs variant-processing errors.

    Every ferro failure that stems from parsing, normalizing, projecting, or
    resolving reference data is an instance of ``FerroError`` (or one of its
    subclasses). Catch ``FerroError`` to handle any of them; catch a specific
    subclass to discriminate.

    Pure argument-validation failures (an empty string, an out-of-range index,
    an unrecognized enum spelling) are raised as plain ``ValueError`` — they are
    programming errors, not variant-processing failures, and are not part of
    this hierarchy.

    Attributes:
        code: The structured ``E####`` / ``W####`` error-code string, or
            ``None`` when the underlying failure carries no code.
        mutalyzer_codes: Equivalent mutalyzer diagnostic codes (e.g.
            ``"EINTRONIC"``); empty when there is no mutalyzer equivalent.
    """

    code: str | None
    mutalyzer_codes: tuple[str, ...]
    def __init__(
        self,
        message: str,
        code: str | None = ...,
        mutalyzer_codes: Iterable[str] = ...,
    ) -> None: ...

class ParseError(FerroError, ValueError):
    """Raised when an HGVS or SPDI string cannot be parsed."""

class NormalizationError(FerroError, RuntimeError):
    """Raised when normalizing a variant fails."""

class ReferenceDataError(FerroError, ValueError, RuntimeError):
    """Raised when reference data is unavailable or a reference lookup fails."""

class ProjectionError(FerroError, ValueError, RuntimeError):
    """Raised when projecting or converting a variant between coordinate systems fails."""

class EquivalenceError(FerroError, RuntimeError):
    """Raised when an equivalence check between variants fails."""
