"""Type stubs for ferro-hgvs Python bindings."""

from collections.abc import Iterable, Iterator
from typing import Any, Callable, Literal, overload

from typing_extensions import Self

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

def normalize(hgvs_string: str) -> str:
    """Normalize an HGVS variant string.

    Args:
        hgvs_string: The HGVS variant description to normalize

    Returns:
        The normalized HGVS string

    Raises:
        ParseError: If the HGVS string cannot be parsed (a subclass of
            ValueError).
        NormalizationError: If normalization fails (a subclass of RuntimeError).

    Example:
        >>> normalize("NM_000088.3:c.100delA")
        'NM_000088.3:c.100del'
    """
    ...

def normalize_with_warnings(hgvs_string: str) -> NormalizeResultWithWarnings:
    """Normalize an HGVS variant string, returning the warnings with it.

    The warning-bearing sibling of :func:`normalize`, which returns only the
    string and so cannot tell a caller that normalization *repaired* something.
    Several repairs are lossy in a way the returned string does not record —
    ``MEMBERS_COALESCED_FROM_REPORTED_FORM`` (separately reported cis members
    were merged), ``INSERTED_SEQUENCE_EXPANDED`` (a bracketed/range ``ins``
    payload became a literal), ``REFSEQ_MISMATCH`` (a stated reference base was
    corrected) — so on :func:`normalize` they look exactly like a clean pass.

    The normalized string is identical to what :func:`normalize` returns; only
    the diagnostics are added.

    Args:
        hgvs_string: The HGVS variant description to normalize

    Returns:
        A NormalizeResultWithWarnings.

    Raises:
        ParseError: If the HGVS string cannot be parsed (a subclass of
            ValueError).
        NormalizationError: If normalization fails (a subclass of RuntimeError).

    Example:
        >>> result = normalize_with_warnings("NM_000088.3:c.100delA")
        >>> str(result.result)
        'NM_000088.3:c.100del'
        >>> [w.code for w in result.warnings]
        []
    """
    ...

# ============================================================================
# Sequence-pair Functions
# ============================================================================

def from_sequences(
    accession: str,
    position: int,
    reference: str,
    alternate: str,
    *,
    max_grid_cells: int | None = None,
) -> HgvsVariant:
    """Derive an HGVS description from a reference/alternate sequence pair.

    The caller supplies the bases; this supplies the description. It reads **no
    reference sequence**, so its output is a pure function of its arguments —
    the same bases give the same description on any machine, against any
    reference build, with no hidden input.

    **Case is not a refusal.** Both sequences are upper-cased before anything
    reads them, so a soft-masked (lower-case) window derives exactly what its
    upper-case twin derives.

    This delivers the two normalization rules that are always achievable —
    conformant and deterministic — and deliberately not the two that need the
    reference, recommended form and confluent. So an output may be 3'-shiftable
    further than a window-local function could shift it, which is not a defect.
    Run ``Normalizer.normalize`` afterwards if you want it, or
    ``Normalizer.from_sequences(..., recommended_form=True)`` to do both in one
    call.

    Args:
        accession: The sequence the window is on. A transcript or protein
            accession is refused — a ``g.`` description on an ``NM_`` denotes
            nothing.
        position: 1-based position of the window's first base.
        reference: The reference bases over the window. Taken on trust and not
            verified — verifying it would need the reference and would make the
            provider a hidden input, costing exactly the determinism this
            function exists to provide.
        alternate: The observed bases over the same window.
        max_grid_cells: Largest alignment grid, in cells, the partitioner will
            build. Defaults to ``(4096 + 1) ** 2``, about 310 MB at roughly 18
            bytes a cell. Exceeding it raises rather than answering with a
            weaker rule.

    Returns:
        The derived HgvsVariant.

    Raises:
        ValueError: for a ``max_grid_cells`` of 0, checked at the binding
            before any derivation.
        NormalizationError: for everything the derivation itself refuses — a
            zero position, an empty reference, a non-nucleotide symbol, a
            transcript or protein accession, a grid over budget, or an inserted
            payload resting against the window's 5' edge with no anchor inside
            it. A subclass of ``FerroError`` and ``RuntimeError``.

    The two classes are **not** one hierarchy, and no single class catches both:
    the binding's ``ValueError``s are outside ``FerroError`` by design, which is
    the repository-wide split between an argument-shape mistake and a
    variant-processing failure. Catch ``(ValueError, ferro_hgvs.FerroError)`` —
    or simply ``Exception`` — to handle everything this function can raise.

    Example:
        >>> str(from_sequences("NC_000001.11", 1000, "AGCG", "AG"))
        'NC_000001.11:g.1002_1003del'
    """
    ...

def from_sequences_detailed(
    accession: str,
    position: int,
    reference: str,
    alternate: str,
    *,
    max_grid_cells: int | None = None,
) -> DerivedDescription:
    """:func:`from_sequences`, reporting also whether the derivation reached a window edge.

    Same arguments and same refusals; see :func:`from_sequences`.
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
    """Convert a 1-based HGVS position to a 0-based array index.

    Raises:
        ValueError: If ``pos`` is 0, which is not a valid 1-based position.
    """
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

class _NativeEnum:
    """Shared behaviour of ferro's native enums.

    These are pyo3 classes, not :class:`enum.Enum` subclasses, and the
    difference is visible: ``isinstance(member, enum.Enum)`` is ``False`` and
    the class itself is not iterable (no ``list(Cls)``, no ``__members__``).
    pyo3 cannot provide either for a ``#[pyclass]`` enum.

    Everything else an ``IntEnum`` offers is available: members hash, order,
    compare equal to their integer, expose :attr:`name` and :attr:`value`, and
    can be looked up by value with ``Cls(value)`` (raising ``ValueError`` for an
    unknown one). ``Cls(value)`` returns the *interned* member — the same object
    bound to the class attribute — so ``is`` behaves as it does on a standard
    library enum, and members hash like the integers they equal.

    The stubs previously declared these as ``IntEnum``, which type-checkers
    honoured while the runtime did not — see issue #1245.
    """

    def __init__(self, value: int) -> None: ...
    @property
    def name(self) -> str:
        """The member's name, matching the attribute it is bound to."""
        ...

    @property
    def value(self) -> int:
        """The member's integer value."""
        ...

    def __int__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __eq__(self, other: object) -> bool: ...
    # Ordering accepts a member of the *same* enum or a plain int, mirroring
    # what pyo3's `ord` + `eq_int` actually implement. `Self` narrows to the
    # concrete subclass, so `Axis.Genomic < Strand.Plus` is rejected the way the
    # runtime rejects it (`TypeError`) rather than silently type-checking. `int`
    # stays in the union because `Axis.Genomic < 3` genuinely works.
    #
    # `__eq__` keeps `object`: equality is total, answering `False` across enums
    # instead of raising.
    def __lt__(self, other: Self | int) -> bool: ...
    def __le__(self, other: Self | int) -> bool: ...
    def __gt__(self, other: Self | int) -> bool: ...
    def __ge__(self, other: Self | int) -> bool: ...

class Axis(_NativeEnum):
    """The coordinate axis (reference molecule / coordinate system) of a variant.

    Returned by :attr:`HgvsVariant.axis`. Each member carries its single-letter
    HGVS coordinate code (:meth:`code`) and a molecule-type grouping
    (:meth:`is_dna` / :meth:`is_rna` / :meth:`is_protein`).
    """

    Genomic: Axis  # 0
    Coding: Axis  # 1
    NonCoding: Axis  # 2
    Rna: Axis  # 3
    Protein: Axis  # 4
    Mitochondrial: Axis  # 5
    Circular: Axis  # 6

    def code(self) -> str:
        """The single-letter HGVS coordinate code (``g``/``c``/``n``/``r``/``p``/``m``/``o``)."""
        ...

    def is_dna(self) -> bool:
        """Whether this axis addresses a DNA molecule (``g``/``c``/``n``/``m``/``o``).

        Treats mitochondrial and circular DNA as DNA, unlike the deprecated
        ``is_genomic()`` predicate (which is ``g.``-only).
        """
        ...

    def is_rna(self) -> bool:
        """Whether this axis addresses an RNA molecule (``r.``)."""
        ...

    def is_protein(self) -> bool:
        """Whether this axis addresses a protein (``p.``)."""
        ...

    def __str__(self) -> str: ...

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
    def axis(self) -> Axis | None:
        """The coordinate axis (reference molecule / coordinate system) of this variant.

        Returns an :class:`Axis`, or ``None`` if there is no single well-defined
        axis. For a leaf variant this is its coordinate kind; for an allele
        (``ACC:c.[...]``) it is the axis shared by every member, so - unlike the
        deprecated ``is_*`` predicates - it works consistently whether or not the
        edit is wrapped in an allele.

        Returns ``None`` for an empty or mixed-axis allele, a bare ``[0]``/``[?]``
        marker, or an RNA-fusion ``::`` construct joining two different
        transcripts. A genome ring (``ACC:g.[seg1::seg2]``) is a single genomic
        accession, so it resolves to :attr:`Axis.Genomic`.
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

    def normalize(self) -> HgvsVariant:
        """Normalize this variant.

        Args:

        Raises:
        """
        ...

    def to_dict(self) -> dict[str, Any]:
        """Convert to a dictionary representation."""
        ...

    def to_json(self) -> str:
        """Convert to JSON string."""
        ...

    def is_genomic(self) -> bool:
        """Check if this is a genomic variant (g. prefix).

        .. deprecated::
            Use the :attr:`axis` property instead (``variant.axis == Axis.Genomic``).
            Unlike this predicate, ``axis`` resolves an allele's shared axis rather
            than returning ``False`` for every allele parent. Emits a
            ``DeprecationWarning``.
        """
        ...

    def is_coding(self) -> bool:
        """Check if this is a coding variant (c. prefix).

        .. deprecated::
            Use the :attr:`axis` property instead (``variant.axis == Axis.Coding``).
            Emits a ``DeprecationWarning``.
        """
        ...

    def is_noncoding(self) -> bool:
        """Check if this is a non-coding variant (n. prefix).

        .. deprecated::
            Use the :attr:`axis` property instead (``variant.axis == Axis.NonCoding``).
            Emits a ``DeprecationWarning``.
        """
        ...

    def is_protein(self) -> bool:
        """Check if this is a protein variant (p. prefix).

        .. deprecated::
            Use the :attr:`axis` property instead (``variant.axis == Axis.Protein``).
            Emits a ``DeprecationWarning``.
        """
        ...

    def is_rna(self) -> bool:
        """Check if this is an RNA variant (r. prefix).

        .. deprecated::
            Use the :attr:`axis` property instead (``variant.axis == Axis.Rna``).
            Emits a ``DeprecationWarning``.
        """
        ...

    def is_mitochondrial(self) -> bool:
        """Check if this is a mitochondrial variant (m. prefix).

        .. deprecated::
            Use the :attr:`axis` property instead (``variant.axis == Axis.Mitochondrial``).
            Emits a ``DeprecationWarning``.
        """
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
    """HGVS variant normalizer using a reference provider.

    Every method that reads the reference — :meth:`normalize`,
    :meth:`normalize_variant`, :meth:`normalize_with_warnings`,
    :meth:`to_spdi`, :meth:`canonical_spdi` and :meth:`apply_to_reference` —
    releases the GIL for the duration of that work, so calls from several
    ``threading.Thread`` workers overlap instead of serialising (#1455). A
    ``Normalizer`` is safe to share between threads.
    """

    def __init__(
        self,
        reference_json: str | None = None,
        error_config: ErrorConfig | None = None,
    ) -> None:
        """Create a new normalizer.

        Args:
            reference_json: Optional path to a transcripts.json file
                error_config: Optional ErrorConfig controlling reference-mismatch
                handling (e.g. ``ErrorConfig.strict()`` to reject
                wrong-reference variants). Defaults to lenient.

        Raises:
        """
        ...

    @staticmethod
    def from_manifest(
        manifest_path: str,
        error_config: ErrorConfig | None = None,
    ) -> Normalizer:
        """Create a normalizer from a reference manifest written by ``ferro prepare``.

        Args:
            manifest_path: Path to a manifest.json file (typically inside a
                directory produced by ``ferro prepare``).
            error_config: Optional ErrorConfig controlling reference-mismatch
                handling (e.g. ``ErrorConfig.strict()`` to reject
                wrong-reference variants). Defaults to lenient.

        Returns:
            A Normalizer backed by a MultiFastaProvider.

        Raises:
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

    def to_spdi(self, variant: HgvsVariant) -> SpdiVariant:
        """Convert to SPDI, resolving edits that need the reference (#1159).

        The module-level `hgvs_to_spdi` does no reference lookup and so fails on
        `del`, `delins`, `inv` and `dup`; this uses the reference this Normalizer
        holds and resolves them.

        Preserves how the variant was written, so it is a transliteration and NOT
        an encoding-invariant key — use `canonical_spdi` for that. A multi-member
        allele is not one triple and is refused.
        """
        ...

    def canonical_spdi(self, variant: HgvsVariant) -> SpdiVariant:
        """An encoding-invariant SPDI key, derived from the resulting bases (#1159).

        Two descriptions on one accession denoting the same edit give the same
        key, whatever their spelling, member count or member order — so a
        spanning `delins` and its decomposed allele match. Normalized HGVS
        strings cannot serve for this, because normalization is not
        encoding-invariant for complex indels (#1157).

        Blunt-trimmed to the block that actually changed: deterministic and
        independent of how wide a span the input named, but not claimed to be
        "maximally shifted" in the SPDI specification's sense.

        Raises:
            ProjectionError: if there is no single well-defined resulting
                sequence — a trans/mosaic/chimeric or null allele, members on
                different accessions, members that overlap, an edit SPDI cannot
                represent, or a span wider than 100 000 bases.
        """
        ...

    def apply_to_reference(self, variant: HgvsVariant) -> AppliedVariant:
        """Apply the variant and return the window before and after (#1159).

        The ground truth for equivalence, and the more general of the two
        primitives: two descriptions denote the same edit exactly when their
        `resulting` bases match. `canonical_spdi` is this reduced to a key.

        Raises:
            ProjectionError: as `canonical_spdi`.
        """
        ...

    def normalize(self, hgvs_string: str) -> str:
        """Parse and normalize an HGVS string.

        Applies 3'/5' shifting and local re-spelling rules to the *given*
        spelling, on any axis (``g/c/n/r/p/m``). This is the entry point for a
        caller who already has a description and wants it normalized.

        Contrast :meth:`rederive`, which discards the spelling and re-derives a
        description from the denoted bases — the confluent, genomic-only path,
        where two spellings of one variant reach one result. See
        :meth:`rederive` for when to reach for each.
        """
        ...

    def to_sequences(self, variant: HgvsVariant, pad: int = 128) -> SequencePair:
        """Apply the variant and return the padded window as a sequence pair.

        The inverse of :func:`from_sequences`, and what turns any HGVS
        description into derivation input — so a caller that already has
        descriptions needs no new plumbing to reach it. Differs from
        `apply_to_reference` in three ways that all matter to that use: the
        position is 1-based, both sides are padded, and it reports whether the
        window's 3' edge is settled (``window_is_final``).

        That flag is **not** "the full pad was served" — see its own doc. A
        window clipped by the sequence end served less than the full pad and is
        still settled, because there is nothing further to read.

        The pad is not decoration. `dup` typing reads the reference bases
        immediately 5' of an insertion point, so a member flush with the
        window's 5' edge comes back as an `ins` rather than a `dup`.

        Args:
            variant: The variant to express as sequences.
            pad: Flank in bases, on each side. Defaults to 128.

        Raises:
            ProjectionError: as `apply_to_reference`.
        """
        ...

    def reanchor(
        self,
        pair: SequencePair,
        start: int | None = None,
        end: int | None = None,
    ) -> SequencePair:
        """Move a window to [start, end], padding from the reference or trimming.

        The reference-holding half of re-anchoring; ``SequencePair.trim_to`` is
        the pure half and can only narrow. ``None`` leaves that edge where it
        is; both bounds are 1-based inclusive.

        Use it to hold a derivation inside a region it must not leave — a target
        region, an amplicon, a tiling window — provided every raw window
        overlaps that region.

        **It moves a window's edges; it does not relocate the window.** Each
        edge may go outwards (padded from the reference) or inwards (trimmed),
        in any combination — but the requested window must overlap the pair's
        own, and the overlap must still hold the bases the two sequences
        disagree on. The bases come back upper-cased, as from
        ``to_sequences``; ``SequencePair.trim_to`` fetches nothing and leaves
        case alone.

        **Not the tool for making heterogeneous inputs agree in general.** That
        is already available and better:
        ``from_sequences(..., recommended_form=True)`` or a round trip through
        ``to_sequences``, both of which reach the
        reference-anchored placement.

        Raises:
            NormalizationError: for every refusal — a bound outside the sequence
                (refused rather than clamped back to the contig), a requested
                window disjoint from the pair's own, and anything
                ``SequencePair.trim_to`` refuses. This method raises one class;
                it is not ``ReferenceDataError``, which this stub named until
                2026-08-12 and which no input could produce.
        """
        ...

    def from_sequences(
        self,
        accession: str,
        position: int,
        reference: str,
        alternate: str,
        *,
        max_grid_cells: int | None = None,
        recommended_form: bool = False,
    ) -> HgvsVariant:
        """:func:`from_sequences`, against this normalizer's reference.

        The free function is a pure function of its arguments and so cannot
        range-check ``position`` — doing that needs the reference, which would
        make the provider a hidden input. This one holds a provider, so it does
        check, and refuses an interval running past the end of the sequence. It
        also offers ``normalize``, which the free function cannot.

        Args:
            accession: The sequence the window is on.
            position: 1-based position of the window's first base.
            reference: The reference bases over the window.
            alternate: The observed bases over the same window.
            max_grid_cells: As the free function.
            recommended_form: Route the derived description through ``normalize``
                to reach ferro's recommended, reference-anchored form before
                returning it. Defaults to False, but prefer True unless you have
                a reason not to: in an internal sweep of many synthetic shapes
                ``normalize`` moved a meaningful share of derived descriptions
                (repeat notation, reference-anchored
                member re-derivation, and inversions spread across several
                members). All three are the recommended-form and confluence rules
                this design assigns to ``normalize``, so False still yields a
                conformant, deterministic description.

        Raises:
            ValueError: for a ``max_grid_cells`` of 0. Checked at the binding
                and outside the ``FerroError`` hierarchy.
            NormalizationError: everything else, including the two refusals this
                method adds over the free function — an unknown accession, and
                an interval running past the end of the sequence — plus any
                refusal from ``normalize`` when ``recommended_form=True``. Those
                two were documented as ``ReferenceDataError`` until 2026-08-12
                and never raised it.
        """
        ...

    def rederive(
        self,
        description: str,
        *,
        max_grid_cells: int | None = None,
        recommended_form: bool = False,
    ) -> str:
        """Re-derive a description from the bases it denotes.

        The "one canonical description per variant" round trip, as a single
        call: express the variant as a padded reference/alternate window
        (``to_sequences``), derive a description from those bases alone
        (``from_sequences``), and — while a member still rests on a window edge
        that can still move — double the pad and retry. Two spellings of one
        variant therefore reach one description, decided by the observed bases
        rather than by how either was written.

        Prefer this over :meth:`normalize` when the input's spelling should
        carry no weight. The two are different operations, not two names for
        one:

        - :meth:`normalize` shifts and re-spells the *given* description in
          place (3'/5' shifting and local re-spelling rules). It works on every
          axis (``g/c/n/r/p/m``).
        - :meth:`rederive` throws the spelling away and re-derives from the
          denoted bases, so two spellings of one variant reach one result
          (confluence). Because it re-derives from a genomic sequence window it
          is genomic-only — ``g.`` (and ``m.`` on the two rCRS accessions);
          any other axis is refused.

        The loop reads ``DerivedDescription``'s two per-side flags apart, so a
        placement pinned to the sequence's own start or end is recognised as
        settled rather than chased.

        Args:
            description: The HGVS description to re-derive. Its axis must be one
                ``from_sequences`` emits — genomic ``g.`` (and ``m.`` on the two
                rCRS mitochondrial accessions); any other is refused.
            max_grid_cells: As ``from_sequences``.
            recommended_form: Route the re-derived description through
                ``normalize`` before returning it. Defaults to False, which
                yields the alignment-derived form (conformant + deterministic)
                rather than ferro's recommended form. Set True for the
                recommended, reference-anchored form.

        Returns:
            The canonical HGVS description, as a string.

        Raises:
            ParseError: ``description`` does not parse.
            ValueError: for a ``max_grid_cells`` of 0, as ``from_sequences``.
            NormalizationError: everything else — a variant with no single
                resulting sequence, a non-genomic axis, a grid over
                ``max_grid_cells``, or a placement still resting on the edge of
                the widest window the loop reads (a repeat tract whose shift
                depends on how much reference is read). Plus any refusal from
                ``normalize`` when ``recommended_form=True``.
        """
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

class EquivalenceLevel(_NativeEnum):
    """Equivalence level between two variants."""

    Identical: EquivalenceLevel  # 0
    NormalizedMatch: EquivalenceLevel  # 1
    AccessionVersionDifference: EquivalenceLevel  # 2
    NotEquivalent: EquivalenceLevel  # 3
    SequenceMatch: EquivalenceLevel  # 4
    CrossAxisSequenceMatch: EquivalenceLevel  # 5
    Indeterminate: EquivalenceLevel  # 6

    def is_equivalent(self) -> bool:
        """Returns true if the variants are considered equivalent.

        False for ``Indeterminate`` too: undecidable is not a positive verdict.
        Use ``is_decided()`` to tell "no" from "cannot tell".
        """
        ...

    def is_decided(self) -> bool:
        """Returns true unless the checker could not decide.

        False only for ``Indeterminate``.
        """
        ...

    def is_at_least(self, floor: EquivalenceLevel) -> bool:
        """Whether this verdict is at least as strong as ``floor``.

        The denotational order is ``Identical`` > ``CrossAxisSequenceMatch`` >
        ``SequenceMatch``. ``NormalizedMatch`` and
        ``AccessionVersionDifference`` are off it and answer false in both
        directions — gating on ``NormalizedMatch`` would define the equivalence
        relation in terms of the normalizer the gate is about.
        """
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
    """Equivalence checker for comparing HGVS variants.

    :meth:`check` and :meth:`all_equivalent` release the GIL while they read the
    reference, so calls from several ``threading.Thread`` workers overlap
    instead of serialising (#1455). An ``EquivalenceChecker`` is safe to share
    between threads.
    """

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

    def check_confluence(
        self,
        variants: list[HgvsVariant],
        relation: Literal["cross_axis", "spdi"] = "cross_axis",
    ) -> dict[str, Any]:
        """Run the opt-in confluence self-check over a corpus of variants.

        Groups the inputs into equivalence classes under ``relation`` and reports
        every class whose members normalize to more than one distinct output — a
        non-confluence witness. A diagnostic: it reports, and never emits a
        pass/fail release verdict.

        Args:
            variants: The variants to check.
            relation: ``"cross_axis"`` (default) is apply-equality on every
                determined axis — the relation that establishes variant identity.
                ``"spdi"`` groups by same denoted bases on the description's own
                axis; it is weaker and insufficient for identity. Typed as a
                ``Literal`` so a misspelling is a type error rather than a
                runtime ``ValueError``; the binding still validates at runtime,
                for callers that are not type-checked.

                ``"cross_axis"`` has no key and so compares every pair, which is
                quadratic in the corpus size; ``"spdi"`` keys each input once and
                is linear. Batch by accession or sample for a large call set —
                classes never span accessions under either relation.

        Returns:
            A dict with keys ``relation`` (str), ``is_confluent`` (bool),
            ``is_complete`` (bool), ``classes_checked`` (int),
            ``undecided_pairs`` (int), ``violations`` (list of
            ``{"inputs": list[str], "outputs": list[str]}``), and ``skipped``
            (list of ``{"input": str, "kind": str, "reason": str,
            "decline_class": str | None, "decline_site": str | None}``).

            ``is_confluent`` is an observation about this corpus under this
            relation, not a release gate, and it is a statement about the *whole*
            corpus only when ``is_complete`` is true — nothing skipped, and no
            comparison undecidable. ``undecided_pairs`` counts comparisons that
            came back undecidable rather than deciding a pair apart; such a pair
            fails to merge exactly as a decided-apart pair does, so a non-zero
            count means the classes are coarser than reported and violations
            inside them are invisible. It is always 0 for ``"spdi"``, which
            compares keys rather than pairs.

            A skip's ``kind`` is ``"unplaceable"`` (the relation never placed it,
            so it is in no class) or ``"normalization_declined"`` (it was placed,
            and its class is counted, but it contributed no output). Coverage is
            therefore ``classes_checked`` plus the ``"unplaceable"`` skips —
            never ``classes_checked + len(skipped)``, which double-counts.

            ``decline_class`` and ``decline_site`` are the skip's refusal as
            machine-readable names rather than prose: the class is
            ``"reference_unavailable"`` (the reference could not serve the bases,
            so the description is not at fault) or ``"unrepresentable"`` (it
            denotes no single sequence SPDI can address), and the site names
            which check refused — ``"multi_molecule_allele"``,
            ``"null_or_unknown_allele"``, ``"empty_allele"``,
            ``"member_conversion"``, ``"cross_accession"`` or
            ``"unresolved_accession"``. Compare these rather than searching
            ``reason`` for a keyword; ``reason`` is a rendering of them and its
            wording is not API. Both are ``None`` where there is no such refusal
            to report: a ``"normalization_declined"`` skip, and the ``"spdi"``
            first-pass mismatch (keyless here, yet its triples project).

        Raises:
            ValueError: If ``relation`` is not ``"cross_axis"`` or ``"spdi"``.
        """
        ...

# ============================================================================
# Effect Prediction Classes
# ============================================================================

class Consequence(_NativeEnum):
    """Sequence Ontology consequence term."""

    TranscriptAblation: Consequence  # 0
    SpliceAcceptorVariant: Consequence  # 1
    SpliceDonorVariant: Consequence  # 2
    StopGained: Consequence  # 3
    FrameshiftVariant: Consequence  # 4
    StopLost: Consequence  # 5
    StartLost: Consequence  # 6
    MissenseVariant: Consequence  # 7
    InframeInsertion: Consequence  # 8
    InframeDeletion: Consequence  # 9
    ProteinAlteringVariant: Consequence  # 10
    SpliceRegionVariant: Consequence  # 11
    SynonymousVariant: Consequence  # 12
    StartRetainedVariant: Consequence  # 13
    StopRetainedVariant: Consequence  # 14
    FivePrimeUtrVariant: Consequence  # 15
    ThreePrimeUtrVariant: Consequence  # 16
    IntronVariant: Consequence  # 17
    CodingSequenceVariant: Consequence  # 18

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

class Impact(_NativeEnum):
    """Variant impact level (VEP-style)."""

    Modifier: Impact  # 0
    Low: Impact  # 1
    Moderate: Impact  # 2
    High: Impact  # 3

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

    def classify_splice_variant(self, offset: int) -> ProteinEffect | None:
        """Classify a splice site variant by distance from splice site.

        Returns ``None`` when the offset is unknown -- ``+?`` / ``-?``, carried
        as the sentinels ``2**63 - 1`` and ``-2**63``. Such a description states
        no distance from the splice site, so there is nothing to classify.

        Args:
            offset: Distance from splice site (negative for acceptor, positive
                for donor)
        """
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

class AppliedVariant:
    """A variant applied to its reference: the same window before and after (#1159)."""

    @property
    def accession(self) -> str:
        """The accession every member acts on."""
        ...

    @property
    def start(self) -> int:
        """0-based start of the window both sequences cover."""
        ...

    @property
    def reference(self) -> str:
        """The reference bases over the window."""
        ...

    @property
    def resulting(self) -> str:
        """The same window after every member has been applied."""
        ...

    def __repr__(self) -> str: ...

class SequencePair:
    """A reference/alternate pair over one window.

    The input `from_sequences` takes and the output `to_sequences` produces.
    """

    def __init__(self, accession: str, position: int, reference: str, alternate: str) -> None:
        """Build a pair from bases you already hold.

        A window out of a BAM, a VCF row, an aligner's output. ``position`` is
        1-based and names the first base of ``reference``.

        ``window_is_final`` is False on a pair built this way: caller-supplied
        bases carry no evidence about whether their 3' edge is where the
        sequence stops or where the read did.

        Raises:
            NormalizationError: for a zero position, an empty reference, or a
                symbol outside the IUPAC-IUBMB nucleotide set (``general.md:48``),
                which ``U`` is excluded from here because this surface's axis is
                DNA — exactly what :func:`from_sequences` refuses, so a pair that
                constructs is a pair that derives. The ambiguity codes (``Y``,
                ``R``, ``S``, …) are admitted; real ClinVar rows carry them.
                ``X`` and ``-`` are not, being alignment-only.
        """
        ...

    def trim_to(self, start: int | None = None, end: int | None = None) -> SequencePair:
        """Narrow this window to [start, end], trimming matching bases only.

        The reference-free half of re-anchoring: use it to hold a derivation
        inside a region it must not leave. ``None`` leaves that edge where it
        is; both bounds are 1-based inclusive.

        It can only narrow — widening needs bases this object does not hold, so
        that is ``Normalizer.reanchor``. Narrowing is not free: a bound that
        cuts an ambiguous run pulls the placement back to the bound, which is
        the point when the bound is a requirement and a footgun when it is
        arbitrary. The derivation reports it as ``placement_bounded_by_window``.

        A narrowing that removes nothing is always accepted, which
        ``Normalizer.reanchor`` depends on: it narrows to the intersection
        before it widens, so a request that contains this window narrows by a
        no-op.

        Raises:
            NormalizationError: for a bound that would widen the window, cut
                into a base the two sequences disagree on, empty the reference,
                or put start past end. Refuses rather than clamping.
        """
        ...

    def derive(self, *, max_grid_cells: int | None = None) -> DerivedDescription:
        """Derive an HGVS description from this window.

        :func:`from_sequences` over the four values this pair already carries,
        so a caller holding a pair — especially one just returned by
        ``trim_to`` or ``reanchor`` — does not re-spread them and cannot
        accidentally pair a pre-trim position with post-trim bases.

        Raises:
            ValueError: for a ``max_grid_cells`` of 0.
            NormalizationError: as the free :func:`from_sequences`.
        """
        ...

    @property
    def end(self) -> int:
        """1-based position of the last base of ``reference``."""
        ...

    @property
    def accession(self) -> str:
        """The accession the window is on."""
        ...

    @property
    def position(self) -> int:
        """1-based position of the window's first base.

        Note this differs from `AppliedVariant.start`, which is 0-based.
        """
        ...

    @property
    def reference(self) -> str:
        """The reference bases over the window."""
        ...

    @property
    def alternate(self) -> str:
        """The same window after the variant has been applied."""
        ...

    @property
    def window_is_final(self) -> bool:
        """Whether the window is as wide 3' as it can usefully be.

        True when it ends where the pad asked *or* where the sequence itself
        ends — so True is the ordinary answer, including near a sequence end
        where the pad was clipped. False means the provider stopped short of
        both and the window may have cut an ambiguous run in half.

        Not "the full pad was served", and 3'-only: the 5' flank is prepended
        separately and is not reflected here.
        """
        ...

    def __repr__(self) -> str: ...
    def __eq__(self, other: object) -> bool:
        """Equality over all five fields, not by identity."""
        ...

    def __hash__(self) -> int:
        """Hash consistent with :meth:`__eq__`."""
        ...

class DerivedDescription:
    """A derived description, plus the one caveat a window-local derivation owes its caller."""

    @property
    def variant(self) -> HgvsVariant:
        """The description."""
        ...

    @property
    def placement_bounded_by_window(self) -> bool:
        """Whether any member rests on an edge of the supplied window.

        A "could move" flag, not an "is wrong" flag, and conservative in that
        direction deliberately: the placement may lie against the edge of the
        bases supplied, so a wider window — or a `Normalizer.normalize` pass,
        which holds the reference — could move it.

        It is *not* "exactly and only when the answer is read-dependent". That
        wording was here and is measurably false: over a 4-base run a window
        flush with the tract is flagged and still gives the same answer a
        whole-sequence derivation gives. A flagged answer is never wrong — it
        denotes the same bases and carries the same canonical SPDI.
        """
        ...

    @property
    def bounded_at_start(self) -> bool:
        """Whether a member rests on the window's 5' edge.

        ``placement_bounded_by_window`` is ``bounded_at_start or
        bounded_at_end``; this says whether the 5' side specifically is the one
        on an edge. A caller widening a window one side at a time reads it to
        decide whether widening 5' could move the answer — and, together with
        the window already starting at base 1, to recognise a placement pinned
        to the sequence start that no widening can move.
        """
        ...

    @property
    def bounded_at_end(self) -> bool:
        """Whether a member rests on the window's 3' edge.

        The 3' counterpart of ``bounded_at_start``; see it. A placement on the
        3' edge of a window that already ends at the sequence's last base is
        settled by the sequence, not the window.
        """
        ...

    def __repr__(self) -> str: ...

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

    def parse_and_rederive(
        self,
        variants: list[str],
        *,
        max_grid_cells: int | None = None,
        recommended_form: bool = False,
        workers: int = 0,
    ) -> BatchResult:
        """Parse and re-derive multiple HGVS strings in parallel (GIL released).

        The batch, parallel counterpart of `Normalizer.rederive`: each variant is
        expressed as the bases it denotes and one canonical description is derived
        from those bases, so two spellings of one variant converge. This is a
        different operation from `parse_and_normalize` (which 3'-shifts the input's
        spelling), mirroring how `Normalizer.rederive` differs from
        `Normalizer.normalize`. The returned BatchResult preserves input order
        regardless of the worker count.

        Args:
            variants: HGVS strings to parse and re-derive.
            max_grid_cells: Largest alignment grid, in cells (as
                `Normalizer.rederive`). None uses the default.
            recommended_form: When True, also apply `normalize`'s recommended-form
                rules to each re-derived description.
            workers: Worker threads. 0 (default) uses all cores; 1 is serial; N uses N threads.

        Raises:
            ValueError: for a `max_grid_cells` of 0.
        """
        ...

    def parse_streaming(self, variants: Iterable[str], workers: int = 0) -> BatchStream:
        """Parse an iterable of HGVS strings lazily, yielding results in input order.

        The streaming counterpart to `parse`. `parse` takes a fully materialized
        list and returns a fully materialized BatchResult, so both are resident at
        once; this consumes the iterable a chunk at a time and yields, so peak
        memory does not grow with the input's length. Measured on 1M inputs:
        ~1782 MiB materializing through `parse` against ~38 MiB streaming, at
        about 1.8x the wall clock (the cost of one result object per item).

        Accepts anything iterable of `str` — a list, a generator, or an open file.
        The GIL is released for the parsing itself, as with `parse`.

        Args:
            variants: Iterable of HGVS strings.
            workers: Worker threads. 0 (default) uses all cores; 1 is serial; N uses N threads.
        """
        ...

    def parse_and_normalize_streaming(
        self, variants: Iterable[str], workers: int = 0
    ) -> BatchStream:
        """Parse and normalize an iterable of HGVS strings lazily, in input order.

        See `parse_streaming` for the rationale and the measured numbers.

        Args:
            variants: Iterable of HGVS strings.
            workers: Worker threads. 0 (default) uses all cores; 1 is serial; N uses N threads.
        """
        ...

    def parse_and_rederive_streaming(
        self,
        variants: Iterable[str],
        *,
        max_grid_cells: int | None = None,
        recommended_form: bool = False,
        workers: int = 0,
    ) -> BatchStream:
        """Parse and re-derive an iterable of HGVS strings lazily, in input order.

        The streaming counterpart of `parse_and_rederive`: use it for a large
        corpus, where materializing one result object per input is the dominant
        cost. See `parse_streaming` for the rationale and memory numbers, and
        `parse_and_rederive` for how re-derivation differs from normalization.

        Args:
            variants: Iterable of HGVS strings.
            max_grid_cells: Largest alignment grid, in cells (as
                `Normalizer.rederive`). None uses the default.
            recommended_form: When True, also apply `normalize`'s recommended-form
                rules to each re-derived description.
            workers: Worker threads. 0 (default) uses all cores; 1 is serial; N uses N threads.

        Raises:
            ValueError: for a `max_grid_cells` of 0.
        """
        ...

    def parse_with_progress(
        self, variants: list[str], callback: Callable[[BatchProgress], None]
    ) -> BatchResult:
        """Parse multiple HGVS strings with progress callback."""
        ...

class BatchItem:
    """One result from a BatchStream.

    A streaming API cannot hand back a BatchResult, whose purpose is to hold
    every result at once. This is the per-item equivalent: exactly one of
    `variant` / `error` is set.
    """

    @property
    def variant(self) -> HgvsVariant | None:
        """The parsed (and, for the normalize stream, normalized) variant, or None."""
        ...

    @property
    def input(self) -> str | None:
        """The input string, present only on failure."""
        ...

    @property
    def error(self) -> str | None:
        """The rendered error, or None on success."""
        ...

    @property
    def ok(self) -> bool:
        """Whether this item parsed (and normalized) successfully."""
        ...

    def __repr__(self) -> str: ...

class BatchStream(Iterator[BatchItem]):
    """Lazy, order-preserving stream of BatchItem, one per input.

    Its own iterator: `iter(stream) is stream`. Once exhausted it stays
    exhausted.
    """

    def __iter__(self) -> BatchStream: ...
    def __next__(self) -> BatchItem: ...

# ============================================================================
# Error Handling Classes
# ============================================================================

class ErrorMode(_NativeEnum):
    """Error handling mode."""

    Strict: ErrorMode  # 0
    Lenient: ErrorMode  # 1
    Silent: ErrorMode  # 2

class ErrorType(_NativeEnum):
    """Error type for configurable error handling."""

    LowercaseAminoAcid: ErrorType  # 0
    MissingVersion: ErrorType  # 1
    WrongDashCharacter: ErrorType  # 2
    ExtraWhitespace: ErrorType  # 3
    ProteinSubstitutionArrow: ErrorType  # 4
    # Discriminant 5 was `PositionZero` (W4002), retired in issue #269.
    SingleLetterAminoAcid: ErrorType  # 6
    WrongQuoteCharacter: ErrorType  # 7
    LowercaseAccessionPrefix: ErrorType  # 8
    MixedCaseEditType: ErrorType  # 9
    OldSubstitutionSyntax: ErrorType  # 10
    InvalidUnicodeCharacter: ErrorType  # 11
    SwappedPositions: ErrorType  # 12
    TrailingAnnotation: ErrorType  # 13
    MissingCoordinatePrefix: ErrorType  # 14
    OldAlleleFormat: ErrorType  # 15
    RefSeqMismatch: ErrorType  # 16
    DeprecatedStopCodonStar: ErrorType  # 17
    DeprecatedStopCodonX: ErrorType  # 18
    DeprecatedFrameshiftStar: ErrorType  # 19
    DeprecatedFrameshiftX: ErrorType  # 20
    DelSizeSuffix: ErrorType  # 21
    EmptyDelinsInsert: ErrorType  # 22
    RedundantRepeatLabel: ErrorType  # 23
    SinglePositionRange: ErrorType  # 24
    DeprecatedIvsNotation: ErrorType  # 25
    DeprecatedConSyntax: ErrorType  # 26
    LengthMismatch: ErrorType  # 27
    AlleleFractionAnnotation: ErrorType  # 28
    ClinVarProseMultiAllelic: ErrorType  # 29
    RnaThymineCanonicalized: ErrorType  # 30
    ProteinBracketedAaInsertion: ErrorType  # 31
    PositionPastEnd: ErrorType  # 32
    VariantExceedsReference: ErrorType  # 33
    NonSpecMosaicForm: ErrorType  # 34
    OverlapConflictingEdits: ErrorType  # 35
    InitiatorMetCanonicalization: ErrorType  # 36
    DupSizeSuffix: ErrorType  # 37
    DupExplicitSeq: ErrorType  # 38
    DelExplicitSeq: ErrorType  # 39
    NonConformantBracketCardinality: ErrorType  # 40
    UnresolvableCentromere: ErrorType  # 41
    TranscriptFlankNotDescribable: ErrorType  # 42
    IntronicOnBareTranscript: ErrorType  # 43
    IncompleteCdsStartReference: ErrorType  # 44
    InsertionWithoutInsertedSequence: ErrorType  # 45
    MembersCoalescedFromReportedForm: ErrorType  # 46
    AlignmentOnlySymbolInDescription: ErrorType  # 47
    NonCodingPositionOutsideTranscript: ErrorType  # 48
    GenomicPositionOffset: ErrorType  # 49
    InsertionSizeCountWithoutSequence: ErrorType  # 50

class ErrorOverride(_NativeEnum):
    """Override behavior for a specific error type."""

    Default: ErrorOverride  # 0
    Reject: ErrorOverride  # 1
    WarnCorrect: ErrorOverride  # 2
    SilentCorrect: ErrorOverride  # 3
    Accept: ErrorOverride  # 4

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

class GenomeBuild(_NativeEnum):
    """Genome build (GRCh37 or GRCh38)."""

    GRCh37: GenomeBuild  # 0
    GRCh38: GenomeBuild  # 1
    Unknown: GenomeBuild  # 2

class Strand(_NativeEnum):
    """Strand direction."""

    Plus: Strand  # 0
    Minus: Strand  # 1
    Unknown: Strand  # 2

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
                A downstream position is REFUSED: it names a nucleotide beyond
                the transcript's last base, which the n. axis cannot number and
                which has no CDS position (HGVS background/numbering.md:52, :54).
                It previously returned a position derived from `tx_position`
                alone, i.e. a different nucleotide, with no error.

        Returns:
            Tuple of (cds_position, offset, is_utr3)

        Raises:
            ReferenceDataError: If the transcript is not found.
            ProjectionError: If it has no CDS, or if `downstream` is True
                (both subclass ValueError).
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

class NormalizeResultWithWarnings:
    """A normalized variant together with the warnings normalization raised.

    Returned by :func:`normalize_with_warnings` and by
    :meth:`Normalizer.normalize_with_warnings`.
    """

    @property
    def result(self) -> HgvsVariant:
        """The normalized variant.

        Identical to what the warning-free entry points return — collecting the
        diagnostics does not move the normalized form.
        """
        ...

    @property
    def warnings(self) -> list[NormalizationWarning]:
        """Warnings raised during normalization; empty when it was clean."""
        ...

    def has_warnings(self) -> bool:
        """True when at least one warning was raised."""
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
        assembly: str | None = None,
        protein_stop: str = "ter",
        amino_acid_code: str = "three",
        error_config: ErrorConfig | None = None,
    ) -> None:
        """Create a variant projector.

        Args:
            reference_json: Path to a transcripts.json file. If None, uses
                built-in test data (limited; not for production).
                ("3prime" or "5prime").
            assembly: Optional genome-build override ("GRCh37"/"GRCh38", or the
                aliases "hg19"/"hg38") for build-agnostic inputs. A bare NG_/LRG_
                input carries no build; this fills one in. An input whose
                accession already encodes a build (NC_*.10/.11) keeps it.
            protein_stop: Stop-codon spelling for rendered p. names — "ter"
                (default) or "star" (`*`).
            amino_acid_code: Amino-acid code width for rendered p. names —
                "three" (default) or "one".
            error_config: Optional ErrorConfig controlling reference-mismatch
                handling (e.g. ``ErrorConfig.strict()`` to reject
                wrong-reference variants). Defaults to lenient.

        Raises:
            ValueError: If ``assembly``, ``protein_stop``, or
                ``amino_acid_code`` is not a recognized value.
        """
        ...

    @staticmethod
    def from_manifest(
        manifest_path: str,
        assembly: str | None = None,
        protein_stop: str = "ter",
        amino_acid_code: str = "three",
        error_config: ErrorConfig | None = None,
    ) -> "VariantProjector":
        """Create a projector from a ferro-prepare manifest.

        Args:
            manifest_path: Path to manifest.json produced by `ferro prepare`.
            assembly: Optional genome-build override ("GRCh37"/"GRCh38", or the
                aliases "hg19"/"hg38") for build-agnostic inputs.
            protein_stop: Stop-codon spelling for rendered p. names — "ter"
                (default) or "star" (`*`).
            amino_acid_code: Amino-acid code width for rendered p. names —
                "three" (default) or "one".
            error_config: Optional ErrorConfig controlling reference-mismatch
                handling (e.g. ``ErrorConfig.strict()`` to reject
                wrong-reference variants). Defaults to lenient.

        Returns:
            A VariantProjector backed by MultiFastaProvider with cdot data.

        Raises:
            ValueError: If ``assembly``, ``protein_stop``, or
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
        too. That is the spec-canonical, 3'-shifted form (e.g. a non-3'-most
        ``g.1003del`` in a poly-A run → ``g.1007del``); 3' is the only direction
        ferro shifts. Pass ``normalize=False`` for the raw pivot, which intentionally does not
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
