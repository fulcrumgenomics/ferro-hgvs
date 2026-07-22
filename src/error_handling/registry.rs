//! Registry of all error and warning codes with documentation.
//!
//! This module provides a comprehensive registry of all error (E-codes) and
//! warning (W-codes) in ferro-hgvs, enabling the `ferro explain` command.

use super::codes::{CodeCategory, CodeInfo, ModeBehavior};
use std::collections::HashMap;
use std::sync::OnceLock;

/// Static registry of all error and warning codes.
static CODE_REGISTRY: OnceLock<HashMap<&'static str, CodeInfo>> = OnceLock::new();

/// Get the code registry, initializing it if needed.
fn get_registry() -> &'static HashMap<&'static str, CodeInfo> {
    CODE_REGISTRY.get_or_init(build_registry)
}

/// Build the code registry.
fn build_registry() -> HashMap<&'static str, CodeInfo> {
    let mut map = HashMap::new();

    // =========================================================================
    // ERROR CODES (E-prefix)
    // =========================================================================

    // --- Parse Errors (E1xxx) ---

    map.insert(
        "E1001",
        CodeInfo {
            code: "E1001",
            name: "InvalidAccession",
            summary: "The accession identifier does not match any recognized format.",
            explanation: "HGVS expressions require a valid reference sequence accession number. \
                Valid prefixes include NM_ (mRNA), NP_ (protein), NC_ (chromosome), \
                NG_ (genomic region), NR_ (non-coding RNA), and LRG_ (Locus Reference Genomic).",
            category: CodeCategory::Parse,
            bad_examples: &["XY_000088.3:c.459A>G", "000088.3:c.459A>G"],
            good_examples: &["NM_000088.3:c.459A>G", "NC_000001.11:g.12345A>G"],
            mode_behavior: None,
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/general/"),
            related_codes: &["E1002", "W1003"],
        },
    );

    map.insert(
        "E1002",
        CodeInfo {
            code: "E1002",
            name: "UnknownVariantType",
            summary: "Unknown variant type prefix.",
            explanation: "HGVS expressions must specify a coordinate type prefix: \
                g. (genomic), c. (coding DNA), n. (non-coding DNA), r. (RNA), \
                p. (protein), m. (mitochondrial), or o. (circular).",
            category: CodeCategory::Parse,
            bad_examples: &["NM_000088.3:x.459A>G", "NM_000088.3:459A>G"],
            good_examples: &["NM_000088.3:c.459A>G", "NC_000001.11:g.12345A>G"],
            mode_behavior: None,
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/general/"),
            related_codes: &["E1001", "W3006"],
        },
    );

    map.insert(
        "E1003",
        CodeInfo {
            code: "E1003",
            name: "InvalidPosition",
            summary: "The position format is invalid.",
            explanation:
                "Positions must be valid integers or offset positions (for intronic variants). \
                Position 0 is never valid on numeric coordinate axes (c., g., n., r., m., o.); \
                numbering starts at 1. Negative positions indicate upstream of the start codon, \
                and *N positions indicate downstream of the stop codon. The protein axis is the \
                deliberate exception — `p.0` (\"no protein product\") and `p.0?` (\"possibly no \
                protein product\") are spec-defined and accepted.",
            category: CodeCategory::Parse,
            bad_examples: &["NM_000088.3:c.0A>G", "NM_000088.3:c.A>G"],
            good_examples: &["NM_000088.3:c.1A>G", "NM_000088.3:c.-10A>G"],
            mode_behavior: None,
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/numbering/",
            ),
            related_codes: &["E3001", "W4001"],
        },
    );

    map.insert(
        "E1004",
        CodeInfo {
            code: "E1004",
            name: "InvalidEdit",
            summary: "The edit type or format is invalid.",
            explanation: "Valid edit types include substitutions (A>G), deletions (del), \
                insertions (ins), duplications (dup), inversions (inv), and \
                deletion-insertions (delins). The edit format must match HGVS syntax.",
            category: CodeCategory::Parse,
            bad_examples: &["NM_000088.3:c.459A->G", "NM_000088.3:c.459mutation"],
            good_examples: &["NM_000088.3:c.459A>G", "NM_000088.3:c.459del"],
            mode_behavior: None,
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/",
            ),
            related_codes: &["E1005", "E1007"],
        },
    );

    map.insert(
        "E1005",
        CodeInfo {
            code: "E1005",
            name: "UnexpectedEnd",
            summary: "Unexpected end of input.",
            explanation: "The HGVS expression ended before it was complete. \
                This usually means a required component is missing.",
            category: CodeCategory::Parse,
            bad_examples: &["NM_000088.3:c.", "NM_000088.3:c.459"],
            good_examples: &["NM_000088.3:c.459A>G", "NM_000088.3:c.459del"],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E1004", "E1006"],
        },
    );

    map.insert(
        "E1006",
        CodeInfo {
            code: "E1006",
            name: "UnexpectedChar",
            summary: "Unexpected character in input.",
            explanation: "The parser encountered a character that is not valid at this position. \
                Check for typos, extra characters, or incorrect syntax.",
            category: CodeCategory::Parse,
            bad_examples: &["NM_000088.3:c.459A>>G", "NM_000088.3:c.459A>G;extra"],
            good_examples: &["NM_000088.3:c.459A>G"],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E1005", "W2004"],
        },
    );

    map.insert(
        "E1007",
        CodeInfo {
            code: "E1007",
            name: "InvalidBase",
            summary: "Invalid nucleotide base.",
            explanation: "Nucleotide bases must be one of A, C, G, T (or U for RNA). \
                For protein variants, use three-letter amino acid codes.",
            category: CodeCategory::Parse,
            bad_examples: &["NM_000088.3:c.459X>G", "NM_000088.3:c.459A>X"],
            good_examples: &["NM_000088.3:c.459A>G", "NM_000088.3:c.459A>T"],
            mode_behavior: None,
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/variant/substitution/",
            ),
            related_codes: &["E1008"],
        },
    );

    map.insert(
        "E1008",
        CodeInfo {
            code: "E1008",
            name: "InvalidAminoAcid",
            summary: "Invalid amino acid code.",
            explanation:
                "Amino acids must be specified using three-letter codes (e.g., Val, Glu, Ala). \
                Single-letter codes may be accepted in lenient mode.",
            category: CodeCategory::Parse,
            bad_examples: &["NP_000079.2:p.Xxx600Yyy", "NP_000079.2:p.Val600Xxx"],
            good_examples: &["NP_000079.2:p.Val600Glu", "NP_000079.2:p.Arg342Ter"],
            mode_behavior: None,
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/protein/variant/",
            ),
            related_codes: &["W1001", "W1002"],
        },
    );

    // --- Reference Errors (E2xxx) ---

    map.insert(
        "E2001",
        CodeInfo {
            code: "E2001",
            name: "ReferenceNotFound",
            summary: "Reference sequence or transcript not found.",
            explanation: "The specified accession could not be found in the reference data. \
                Ensure the accession is correct and that reference data has been loaded.",
            category: CodeCategory::Reference,
            bad_examples: &["NM_999999.1:c.100A>G"],
            good_examples: &["NM_000088.3:c.100A>G"],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E2002", "E2003", "E2004", "E2005"],
        },
    );

    map.insert(
        "E2002",
        CodeInfo {
            code: "E2002",
            name: "SequenceNotFound",
            summary: "Sequence data not available.",
            explanation: "The sequence for this accession is not available. \
                This may occur when transcript metadata exists but sequence data is missing.",
            category: CodeCategory::Reference,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E2001"],
        },
    );

    map.insert(
        "E2003",
        CodeInfo {
            code: "E2003",
            name: "ChromosomeNotFound",
            summary: "Chromosome or contig not found.",
            explanation:
                "The specified chromosome or contig is not available in the reference data. \
                Check that the genome build matches your reference data.",
            category: CodeCategory::Reference,
            bad_examples: &["NC_999999.1:g.100A>G"],
            good_examples: &["NC_000001.11:g.100A>G"],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E2001"],
        },
    );

    map.insert(
        "E2004",
        CodeInfo {
            code: "E2004",
            name: "TranscriptVersionNotExact",
            summary: "Exact transcript version not available under strict resolution.",
            explanation:
                "The requested transcript was found only at a different version (or via cdot \
                reconstruction); strict resolution refuses to substitute. Request an available \
                version or relax strict version pinning.",
            category: CodeCategory::Reference,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E2001"],
        },
    );

    map.insert(
        "E2005",
        CodeInfo {
            code: "E2005",
            name: "TranscriptSequenceUnreconstructable",
            summary: "Transcript sequence cannot be reconstructed from cdot.",
            explanation:
                "Base synthesis fell back to cdot's exon alignment against the genome (the \
                transcript FASTA lacks this accession), but the cdot CIGAR encodes an insertion \
                — transcript bases with no genome counterpart. cdot records only the insertion \
                length, not the inserted bases, so the sequence cannot be reconstructed; ferro \
                declines rather than serve a knowingly-divergent sequence.",
            category: CodeCategory::Reference,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E2001"],
        },
    );

    // --- Validation Errors (E3xxx) ---

    map.insert(
        "E3001",
        CodeInfo {
            code: "E3001",
            name: "PositionOutOfBounds",
            summary: "Position is outside the valid range.",
            explanation: "The specified position is beyond the length of the reference sequence. \
                Check that the position is correct for this transcript or chromosome.",
            category: CodeCategory::Validation,
            bad_examples: &["NM_000088.3:c.99999999A>G"],
            good_examples: &["NM_000088.3:c.459A>G"],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E1003", "E3003"],
        },
    );

    map.insert(
        "E3002",
        CodeInfo {
            code: "E3002",
            name: "ReferenceMismatch",
            summary: "Reference sequence mismatch.",
            explanation: "The reference base(s) stated in the HGVS expression do not match \
                the actual reference sequence at that position. This may indicate a wrong \
                transcript version or a typo in the variant description.",
            category: CodeCategory::Validation,
            bad_examples: &["NM_000088.3:c.459G>A (when ref is actually A)"],
            good_examples: &["NM_000088.3:c.459A>G"],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["W5001", "W3001"],
        },
    );

    map.insert(
        "E3003",
        CodeInfo {
            code: "E3003",
            name: "InvalidRange",
            summary: "Invalid coordinate range.",
            explanation: "The start position is greater than the end position, \
                or the range is otherwise invalid.",
            category: CodeCategory::Validation,
            bad_examples: &["NM_000088.3:c.200_100del"],
            good_examples: &["NM_000088.3:c.100_200del"],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["W4001"],
        },
    );

    map.insert(
        "E3004",
        CodeInfo {
            code: "E3004",
            name: "ExonIntronBoundary",
            summary: "Variant spans exon-intron boundary.",
            explanation:
                "The variant spans an exon-intron junction, which cannot be properly normalized. \
                Split into separate exonic and intronic components if needed.",
            category: CodeCategory::Validation,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E3005", "E4001"],
        },
    );

    map.insert(
        "E3005",
        CodeInfo {
            code: "E3005",
            name: "UtrCdsBoundary",
            summary: "Variant spans UTR-CDS boundary.",
            explanation: "The variant spans the boundary between UTR and coding sequence. \
                This requires special handling for normalization.",
            category: CodeCategory::Validation,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E3004"],
        },
    );

    map.insert(
        "E3006",
        CodeInfo {
            code: "E3006",
            name: "SelfCancellingAllele",
            summary: "Self-cancelling allele (overlapping del + dup or equivalent).",
            explanation: "HGVS does not allow descriptions that remove part of a reference \
                sequence and re-insert (part of) the same sequence in the same allele. \
                For example, `c.[762_768del;767_774dup]` deletes [762..=768] and re-introduces \
                bases [767..=774] from the reference, partially undoing the deletion. \
                Drop the redundant edit (or rewrite as a single `delins`).",
            category: CodeCategory::Validation,
            bad_examples: &["NM_004006.2:c.[762_768del;767_774dup]"],
            good_examples: &["NM_004006.2:c.[100_110del;200_210dup]"],
            mode_behavior: None,
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/general/"),
            related_codes: &["E3003", "W3003"],
        },
    );

    map.insert(
        "E3007",
        CodeInfo {
            code: "E3007",
            name: "AlignmentGap",
            summary:
                "Variant position falls within a transcript-genome alignment gap (CIGAR indel).",
            explanation: "The variant position lands strictly inside a cdot transcript-genome \
                CIGAR indel (an `Insertion` is a transcript base with no genome counterpart; a \
                `Deletion` is a genome base with no transcript counterpart). The position has no \
                well-defined counterpart on the other axis, so projecting through it is refused \
                instead of emitting a silently-wrong coordinate.",
            category: CodeCategory::Validation,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E3004", "E5001"],
        },
    );

    // --- Normalization Errors (E4xxx) ---

    map.insert(
        "E4001",
        CodeInfo {
            code: "E4001",
            name: "IntronicVariant",
            summary: "Intronic variant cannot be normalized.",
            explanation: "Intronic variants (those with +/- offset positions) cannot be normalized \
                without genomic context. Use genomic coordinates for intronic variant normalization.",
            category: CodeCategory::Normalization,
            bad_examples: &["NM_000088.3:c.459+10del"],
            good_examples: &["NC_000017.11:g.12345del"],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E3004", "W4007"],
        },
    );

    map.insert(
        "E4002",
        CodeInfo {
            code: "E4002",
            name: "UnsupportedVariant",
            summary: "Unsupported variant type for this operation.",
            explanation: "This variant type is not supported for the requested operation. \
                Not all HGVS variant types can be normalized or converted.",
            category: CodeCategory::Normalization,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &[],
        },
    );

    // --- Conversion Errors (E5xxx) ---

    map.insert(
        "E5001",
        CodeInfo {
            code: "E5001",
            name: "ConversionFailed",
            summary: "Coordinate conversion failed.",
            explanation: "The coordinate conversion between reference sequences failed. \
                This may occur when converting between coding and genomic coordinates.",
            category: CodeCategory::Conversion,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E5002"],
        },
    );

    map.insert(
        "E5002",
        CodeInfo {
            code: "E5002",
            name: "NoOverlappingTranscript",
            summary: "No overlapping transcript found.",
            explanation:
                "No transcript was found that overlaps with the specified genomic position. \
                The position may be intergenic or the transcript data may be incomplete.",
            category: CodeCategory::Conversion,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E5001", "E2001"],
        },
    );

    // --- IO Errors (E9xxx) ---

    map.insert(
        "E9001",
        CodeInfo {
            code: "E9001",
            name: "IoError",
            summary: "File I/O error.",
            explanation: "An error occurred while reading or writing a file. \
                Check that the file exists and you have the necessary permissions.",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E9002"],
        },
    );

    map.insert(
        "E9002",
        CodeInfo {
            code: "E9002",
            name: "JsonError",
            summary: "JSON parsing error.",
            explanation: "An error occurred while parsing JSON data. \
                Check that the JSON is valid and matches the expected format.",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E9001"],
        },
    );

    // =========================================================================
    // WARNING CODES (W-prefix)
    // =========================================================================

    // --- Case/Capitalization Warnings (W1xxx) ---

    map.insert(
        "W1001",
        CodeInfo {
            code: "W1001",
            name: "LowercaseAminoAcid",
            summary: "Lowercase amino acid code.",
            explanation: "Three-letter amino acid codes should use standard capitalization \
                with the first letter uppercase (e.g., Val, Glu, Ala). Lowercase codes like \
                'val' or 'glu' are auto-corrected in lenient/silent modes.",
            category: CodeCategory::Case,
            bad_examples: &["p.val600glu", "p.VAL600GLU"],
            good_examples: &["p.Val600Glu"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/protein/variant/",
            ),
            related_codes: &["W1002", "E1008"],
        },
    );

    map.insert(
        "W1002",
        CodeInfo {
            code: "W1002",
            name: "SingleLetterAminoAcid",
            summary: "Single-letter amino acid code.",
            explanation: "HGVS recommends three-letter amino acid codes. Single-letter codes \
                like 'V' for Valine may be auto-expanded in lenient/silent modes.",
            category: CodeCategory::Case,
            bad_examples: &["p.V600E"],
            good_examples: &["p.Val600Glu"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/protein/variant/",
            ),
            related_codes: &["W1001"],
        },
    );

    map.insert(
        "W1003",
        CodeInfo {
            code: "W1003",
            name: "LowercaseAccessionPrefix",
            summary: "Lowercase accession prefix.",
            explanation: "Accession prefixes should be uppercase (NM_, NC_, NP_, etc.). \
                Lowercase prefixes like 'nm_' are auto-corrected in lenient/silent modes.",
            category: CodeCategory::Case,
            bad_examples: &["nm_000088.3:c.459A>G", "nc_000001.11:g.12345A>G"],
            good_examples: &["NM_000088.3:c.459A>G", "NC_000001.11:g.12345A>G"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: None,
            related_codes: &["E1001"],
        },
    );

    map.insert(
        "W1004",
        CodeInfo {
            code: "W1004",
            name: "MixedCaseEditType",
            summary: "Mixed case in edit type keyword.",
            explanation: "Edit type keywords should be lowercase (del, ins, dup, etc.). \
                Mixed case like 'Del' or 'INS' is auto-corrected in lenient/silent modes.",
            category: CodeCategory::Case,
            bad_examples: &["c.100Del", "c.100_101INS"],
            good_examples: &["c.100del", "c.100_101ins"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: None,
            related_codes: &[],
        },
    );

    // --- Character Warnings (W2xxx) ---

    map.insert(
        "W2001",
        CodeInfo {
            code: "W2001",
            name: "WrongDashCharacter",
            summary: "Wrong dash character (en-dash or em-dash).",
            explanation: "HGVS requires ASCII hyphen-minus (-). Word processors often \
                substitute en-dash (–) or em-dash (—). These are auto-corrected in \
                lenient/silent modes.",
            category: CodeCategory::Character,
            bad_examples: &["c.100–200del", "c.100—200del"],
            good_examples: &["c.100-200del"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: None,
            related_codes: &["W2002", "W2004"],
        },
    );

    map.insert(
        "W2002",
        CodeInfo {
            code: "W2002",
            name: "WrongQuoteCharacter",
            summary: "Smart quotes instead of ASCII quotes.",
            explanation: "HGVS requires ASCII quotes. Word processors often substitute \
                curly/smart quotes. These are auto-corrected in lenient/silent modes.",
            category: CodeCategory::Character,
            bad_examples: &["c.100_101ins\u{201C}ATG\u{201D}"],
            good_examples: &["c.100_101ins\"ATG\""],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: None,
            related_codes: &["W2001", "W2004"],
        },
    );

    map.insert(
        "W2003",
        CodeInfo {
            code: "W2003",
            name: "ExtraWhitespace",
            summary: "Extra whitespace in description.",
            explanation: "HGVS expressions should not contain spaces (the spec is explicit \
                that the format `reference:description` has spaces \"added for clarity only\"). \
                Leading, trailing, and embedded whitespace — and zero-width invisible \
                characters (U+200B ZERO WIDTH SPACE, U+200C ZERO WIDTH NON-JOINER, \
                U+200D ZERO WIDTH JOINER, U+FEFF BOM/ZWNBSP) — are stripped in \
                lenient/silent modes.",
            category: CodeCategory::Character,
            bad_examples: &[
                "c.100 A>G",
                "NM_000088.3 : c.459A>G",
                "c.[100A>G ; 200T>C]",
                "NM_000088.3\u{200B}:c.100A>G",
            ],
            good_examples: &["c.100A>G", "NM_000088.3:c.459A>G"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/general/"),
            related_codes: &["W2001", "W2004"],
        },
    );

    map.insert(
        "W2004",
        CodeInfo {
            code: "W2004",
            name: "InvalidUnicodeCharacter",
            summary: "Non-ASCII Unicode character.",
            explanation: "HGVS expressions should use ASCII characters only. \
                Common Unicode look-alikes are auto-corrected in lenient/silent modes.",
            category: CodeCategory::Character,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: None,
            related_codes: &["W2001", "W2002"],
        },
    );

    // --- Format/Syntax Warnings (W3xxx) ---

    map.insert(
        "W3001",
        CodeInfo {
            code: "W3001",
            name: "MissingVersion",
            summary: "Missing transcript version.",
            explanation: "Reference sequence accessions should include a version number \
                (e.g., NM_000088.3 not NM_000088). Without a version, the correct reference \
                sequence cannot be determined. In lenient/silent modes, parsing succeeds but \
                the variant is flagged for downstream handling.",
            category: CodeCategory::Format,
            bad_examples: &["NM_000088:c.459A>G"],
            good_examples: &["NM_000088.3:c.459A>G"],
            mode_behavior: Some(ModeBehavior::warn_accept()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/general/"),
            related_codes: &["E3002"],
        },
    );

    map.insert(
        "W3002",
        CodeInfo {
            code: "W3002",
            name: "ProteinSubstitutionArrow",
            summary: "Arrow in protein substitution.",
            explanation: "Protein substitutions should not use '>' between amino acids. \
                Use the format 'p.Val600Glu' not 'p.Val600>Glu'.",
            category: CodeCategory::Format,
            bad_examples: &["p.Val600>Glu"],
            good_examples: &["p.Val600Glu"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/protein/variant/substitution/"),
            related_codes: &[],
        },
    );

    map.insert(
        "W3003",
        CodeInfo {
            code: "W3003",
            name: "OldSubstitutionSyntax",
            summary: "Old-style multi-base substitution syntax.",
            explanation: "Multi-base substitutions should use 'delins' not '>'. \
                For example, 'c.100_102delinsATG' not 'c.100_102>ATG'.",
            category: CodeCategory::Format,
            bad_examples: &["c.100_102>ATG"],
            good_examples: &["c.100_102delinsATG"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: None,
            related_codes: &[],
        },
    );

    map.insert(
        "W3004",
        CodeInfo {
            code: "W3004",
            name: "OldAlleleFormat",
            summary: "Old/deprecated allele bracket format.",
            explanation: "In compound allele notation, the coordinate type should appear before \
                the brackets. Use 'NM_000088.3:c.[100A>G;200C>T]' not 'NM_000088.3:[c.100A>G;c.200C>T]'.",
            category: CodeCategory::Format,
            bad_examples: &["NM_000088.3:[c.100A>G;c.200C>T]"],
            good_examples: &["NM_000088.3:c.[100A>G;200C>T]"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: None,
            related_codes: &[],
        },
    );

    map.insert(
        "W3005",
        CodeInfo {
            code: "W3005",
            name: "TrailingAnnotation",
            summary: "Trailing protein annotation.",
            explanation: "ClinVar and other databases often append protein consequence annotations \
                to HGVS expressions (e.g., '(p.Lys236=)'). These are stripped in lenient/silent modes.",
            category: CodeCategory::Format,
            bad_examples: &["NM_000088.3:c.459A>G (p.Lys153=)", "NM_000088.3:c.100del (p.Val34fs)"],
            good_examples: &["NM_000088.3:c.459A>G", "NM_000088.3:c.100del"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: None,
            related_codes: &[],
        },
    );

    map.insert(
        "W3006",
        CodeInfo {
            code: "W3006",
            name: "MissingCoordinatePrefix",
            summary: "Missing coordinate type prefix.",
            explanation: "HGVS expressions require a coordinate type prefix (g., c., p., etc.). \
                For chromosomal accessions (NC_), 'g.' can be inferred in lenient/silent modes.",
            category: CodeCategory::Format,
            bad_examples: &["NC_000001.11:12345A>G"],
            good_examples: &["NC_000001.11:g.12345A>G"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: None,
            related_codes: &["E1002"],
        },
    );

    map.insert(
        "W3007",
        CodeInfo {
            code: "W3007",
            name: "DeprecatedStopCodonStar",
            summary: "'*' stop-codon glyph (spec-valid; retained but no longer emitted).",
            explanation: "Per the HGVS checklist ('Ter' or '*' should be used for a \
                translation stop codon; 'X' should not), '*' is a spec-valid glyph \
                co-equal with 'Ter' — not a deprecated form. The parser accepts '*' \
                natively in every mode (e.g. 'p.Arg97*') and canonicalizes it to the \
                preferred three-letter 'Ter' on display, so no correction or warning is \
                produced. This code is retained (python-exposed, stable discriminant) but \
                is never emitted (#1114).",
            category: CodeCategory::Format,
            bad_examples: &[],
            good_examples: &["NP_000079.2:p.Arg97Ter", "NP_000079.2:p.Trp10Ter"],
            mode_behavior: Some(ModeBehavior::accepted()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/checklist/"),
            related_codes: &["W3008", "W3009", "W3010"],
        },
    );

    map.insert(
        "W3008",
        CodeInfo {
            code: "W3008",
            name: "DeprecatedStopCodonX",
            summary: "Deprecated 'X' for stop codon in protein substitution.",
            explanation: "Per the HGVS checklist, 'the X should not be used' to indicate a \
                translation stop codon. The single-letter 'X' is also reserved for the \
                'any amino acid' symbol Xaa, making 'p.Arg97X' ambiguous. In lenient/silent \
                modes the input is rewritten to 'p.Arg97Ter'.",
            category: CodeCategory::Format,
            bad_examples: &["NP_000079.2:p.Arg97X", "NP_000079.2:p.Trp10X"],
            good_examples: &["NP_000079.2:p.Arg97Ter", "NP_000079.2:p.Trp10Ter"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/protein/variant/substitution/",
            ),
            related_codes: &["W3007", "W3009", "W3010"],
        },
    );

    map.insert(
        "W3009",
        CodeInfo {
            code: "W3009",
            name: "DeprecatedFrameshiftStar",
            summary: "'fs*N' frameshift-termination glyph (spec-valid; retained but no longer \
                emitted).",
            explanation: "Per HGVS (the checklist's 'Ter' or '*' guidance, and \
                'protein/frameshift.md', which lists 'p.Arg97Profs*23' as a valid example), \
                'fs*N' is a spec-valid frameshift-termination glyph — not a deprecated form. \
                The parser accepts it natively in every mode (e.g. 'p.Arg97fs*23') and \
                canonicalizes it to the preferred 'fsTerN' on display, so no correction or \
                warning is produced. This code is retained (python-exposed, stable \
                discriminant) but is never emitted (#1114).",
            category: CodeCategory::Format,
            bad_examples: &[],
            good_examples: &[
                "NP_000079.2:p.Arg97fsTer23",
                "NP_000079.2:p.Ser539AlafsTer110",
            ],
            mode_behavior: Some(ModeBehavior::accepted()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/protein/variant/frameshift/",
            ),
            related_codes: &["W3007", "W3008", "W3010"],
        },
    );

    map.insert(
        "W3010",
        CodeInfo {
            code: "W3010",
            name: "DeprecatedFrameshiftX",
            summary: "Deprecated 'fsXN' frameshift termination notation.",
            explanation: "Neither 'X' nor 'Xaa' is used in HGVS frameshift termination; the \
                canonical form is 'fsTerN'. The strict parser rejects 'fsXN' outright; in \
                lenient/silent modes the preprocessor rewrites 'p.Arg97fsX23' to \
                'p.Arg97fsTer23' and emits this warning.",
            category: CodeCategory::Format,
            bad_examples: &["NP_000079.2:p.Arg97fsX23"],
            good_examples: &["NP_000079.2:p.Arg97fsTer23"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/protein/variant/frameshift/",
            ),
            related_codes: &["W3007", "W3008", "W3009"],
        },
    );

    map.insert(
        "W3011",
        CodeInfo {
            code: "W3011",
            name: "DelSizeSuffix",
            summary: "Deletion described with a size-count suffix instead of a position range.",
            explanation: "HGVS requires both endpoints of a multi-residue deletion to be named \
                (`checklist.md:49`: `g.123del3` is \"not allowed, correct is `g.123_125del`\"). \
                When the anchor is a plain point position ferro treats this as MUST-level: the \
                default parse rejects, and lenient/silent modes rewrite `g.123del6` to \
                `g.123_128del`. When the anchor is intronic/offset the end position cannot be \
                synthesized safely, so there is no repair to offer and the form is rejected \
                outright. A range anchor that merely repeats the size (`c.100_102del3`) already \
                names both endpoints and only warns.",
            category: CodeCategory::Format,
            bad_examples: &["NG_012232.1:g.123del6"],
            good_examples: &["NG_012232.1:g.123_128del"],
            mode_behavior: Some(ModeBehavior::warn_accept()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/",
            ),
            related_codes: &["W4003"],
        },
    );

    map.insert(
        "W3012",
        CodeInfo {
            code: "W3012",
            name: "EmptyDelinsInsert",
            summary: "Deletion-insertion with empty inserted sequence.",
            explanation:
                "A `delins` with no inserted sequence is semantically equivalent to a plain \
                deletion. Per HGVS spec the canonical form is `del`. In lenient/silent modes, \
                ferro rewrites `delins` to `del` and warns once.",
            category: CodeCategory::Format,
            bad_examples: &["NC_000001.11:g.100_102delins"],
            good_examples: &["NC_000001.11:g.100_102del"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/"),
            related_codes: &[],
        },
    );

    map.insert(
        "W3013",
        CodeInfo {
            code: "W3013",
            name: "RedundantRepeatLabel",
            summary: "Repeat description with redundant base label.",
            explanation:
                "HGVS spec discourages including the repeat-unit base label when the positions \
                already define the unit, e.g. `r.-125_-123cug[4]` should be written as \
                `r.-125_-123[4]`. In lenient/silent modes, ferro strips the redundant base \
                segment and warns once.",
            category: CodeCategory::Format,
            bad_examples: &["NM_000088.3:r.100_102cug[4]"],
            good_examples: &["NM_000088.3:r.100_102[4]"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/RNA/repeated/",
            ),
            related_codes: &[],
        },
    );

    map.insert(
        "W3014",
        CodeInfo {
            code: "W3014",
            name: "DeprecatedIvsNotation",
            summary: "Retracted c.IVS intronic notation.",
            explanation: "The c.IVSn+offset / c.IVSn-offset notation has been retracted by \
                HGVS because the IVS number is not unique across transcripts and cannot be \
                resolved to a unique genomic position without metadata. Use the canonical \
                intronic-offset form, e.g. c.88+2T>G instead of c.IVS2+2T>G. ferro cannot \
                auto-rewrite this form without genomic context, so all modes reject the input \
                and emit an actionable hint.",
            category: CodeCategory::Format,
            bad_examples: &["c.IVS2+2T>G", "c.IVS5-1G>T"],
            good_examples: &["c.88+2T>G", "c.123-1G>T"],
            mode_behavior: Some(ModeBehavior::always_reject()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/background/numbering/"),
            related_codes: &["E1006", "E4001"],
        },
    );

    map.insert(
        "W3015",
        CodeInfo {
            code: "W3015",
            name: "DeprecatedConSyntax",
            summary: "Deprecated 'con' (sequence conversion) edit syntax.",
            explanation:
                "The HGVS spec retired the `con` edit type in favour of `delins`. Conversions \
                (a range of nucleotides replaced by a sequence from elsewhere) should be \
                described as deletion-insertions: `c.100_200conNM_001.1:c.5_105` becomes \
                `c.100_200delinsNM_001.1:c.5_105`. In lenient/silent modes ferro auto-rewrites; \
                in strict mode the input is rejected.",
            category: CodeCategory::Format,
            bad_examples: &["c.100_200conNM_001.1:c.5_105"],
            good_examples: &["c.100_200delinsNM_001.1:c.5_105"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/"),
            related_codes: &["W3003"],
        },
    );

    map.insert(
        "W3016",
        CodeInfo {
            code: "W3016",
            name: "LengthMismatch",
            summary: "Explicit reference sequence length does not match position range.",
            explanation:
                "For `del` / `dup` / `inv` / `delins` with an explicit reference sequence, the \
                sequence length must equal the position range length (`end - start + 1`). \
                `g.100_110delAAAATTTGCC` has range 11 but sequence length 10. Lenient mode emits \
                W3016 without rewriting (the principled correction depends on which endpoint is \
                wrong); strict mode rejects.",
            category: CodeCategory::Format,
            bad_examples: &["g.100_110delAAAATTTGCC", "g.100_105invTAGCA"],
            good_examples: &["g.100_109delAAAATTTGCC", "g.100_104invTAGCA"],
            mode_behavior: Some(ModeBehavior::warn_accept()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/general/"),
            related_codes: &["W5001"],
        },
    );

    map.insert(
        "W3017",
        CodeInfo {
            code: "W3017",
            name: "AlleleFractionAnnotation",
            summary: "Allele-fraction / heteroplasmy annotation appended to an HGVS expression.",
            explanation:
                "Annotations such as `[level=70%]`, `[heteroplasmy=70%]`, `[mosaic=80%]`, or \
                `(80%)` are sometimes appended to mitochondrial or somatic HGVS expressions to \
                record heteroplasmy / variant allele fraction. The HGVS spec does not encode \
                allele fraction inside the variant string; fraction belongs in accompanying \
                metadata (VCF `FORMAT/AF`, ClinVar's heteroplasmy field, etc.). ferro rejects \
                this annotation in all modes and emits this code so downstream tooling can \
                recognize the intent rather than seeing a generic trailing-character parse \
                error. The remediation is to strip the annotation from the HGVS string and \
                record the fraction in the appropriate metadata channel (e.g. VCF \
                `FORMAT/AF`, ClinVar's heteroplasmy field, or a sibling column in the \
                caller's data model).",
            category: CodeCategory::Format,
            bad_examples: &[
                "NC_012920.1:m.3243A>G[level=70%]",
                "NC_012920.1:m.3243A>G(80%)",
                "NC_012920.1:m.3243A>G[heteroplasmy=70%]",
            ],
            good_examples: &["NC_012920.1:m.3243A>G"],
            mode_behavior: Some(ModeBehavior::always_reject()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/",
            ),
            related_codes: &[],
        },
    );

    map.insert(
        "W3018",
        CodeInfo {
            code: "W3018",
            name: "ClinVarProseMultiAllelic",
            summary: "ClinVar prose multi-allelic shorthand m.<pos><ref>><alt>/<alt2>.",
            explanation: "Some ClinVar submissions report two alternative alleles at the same \
                mitochondrial position with the prose shorthand `m.<pos><ref>><alt>/<alt2>` \
                (e.g. `m.3243A>G/T`), where the RHS of the slash is a bare nucleotide letter — \
                no edit operator, no accession, no `=` reference marker. The HGVS spec does \
                not define this form: the mosaic `/` separator requires either two fully \
                qualified variants, or the spec compact mosaic `<pos>=/<alt-edit>` form (`=` \
                stands in for the reference, edit names the alternative). All modes reject \
                this shape and emit this code so downstream tooling can recognize the intent. \
                Three spec-supported alternatives: (a) compound brackets \
                `<acc>:m.[3243A>G;3243A>T]`; (b) dual fully qualified slash \
                `<acc>:m.3243A>G/<acc>:m.3243A>T`; (c) spec compact mosaic \
                `<acc>:m.3243=/A>T`.",
            category: CodeCategory::Format,
            bad_examples: &["NC_012920.1:m.3243A>G/T", "NC_012920.1:m.3243A>G/C"],
            good_examples: &[
                "NC_012920.1:m.[3243A>G;3243A>T]",
                "NC_012920.1:m.3243A>G/NC_012920.1:m.3243A>T",
                "NC_012920.1:m.3243=/A>T",
            ],
            mode_behavior: Some(ModeBehavior::always_reject()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/",
            ),
            related_codes: &["W3017"],
        },
    );

    map.insert(
        "W3020",
        CodeInfo {
            code: "W3020",
            name: "RnaThymineCanonicalized",
            summary: "Thymine ('t'/'T') used as a base in an r. (RNA) description.",
            explanation:
                "Per HGVS v21.0 RNA nomenclature, the RNA alphabet is `a/c/g/u`. Thymine is \
                non-canonical input inside an `r.` description. In lenient/silent modes ferro \
                rewrites each `t`/`T` to `u` (lowercased) and emits one warning per occurrence. \
                In strict mode the input is rejected. This applies to every RNA edit shape \
                where a base byte appears: substitution (ref/alt), deletion (stated ref), \
                insertion (alt), duplication (stated ref), inversion (stated ref), delins \
                (stated ref and alt), and repeat unit.",
            category: CodeCategory::Format,
            bad_examples: &[
                "NM_000088.3:r.123a>t",
                "NM_000088.3:r.123_125delAUT",
                "NM_000088.3:r.123aut[5]",
            ],
            good_examples: &[
                "NM_000088.3:r.123a>u",
                "NM_000088.3:r.123_125delauu",
                "NM_000088.3:r.123auu[5]",
            ],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/RNA/substitution/",
            ),
            related_codes: &[],
        },
    );

    map.insert(
        "W3019",
        CodeInfo {
            code: "W3019",
            name: "NonSpecMosaicForm",
            summary: "Non-spec mosaic / chimeric form (nested `/`+`//` or bracketed `[a/b]`).",
            explanation:
                "HGVS v21 does not define two related mosaic/chimeric shapes that ferro rejects \
                under one diagnostic family: (1) **nesting** `/` and `//` at the same level \
                (e.g. `m.[3243A>G/T]//[3243A>C/G]` — chimeric-of-mosaic, or the mirror image \
                mosaic-of-chimeric), and (2) `[a/b]` — a bracketed mosaic group (e.g. \
                `c.[100A>G/200T>C]`). Both rejections are upgraded to this targeted code so \
                downstream tools can key off a single family. Cannot be auto-corrected (no \
                canonical alternative exists); use compound brackets `[a;b]`, dual \
                fully-qualified slash `acc:c.X/acc:c.Y` (mosaic) or `acc:c.X//acc:c.Y` \
                (chimeric), or the compact short-hands `acc:c.<pos>=/<edit>` (mosaic) and \
                `acc:c.<pos>=//<edit>` (chimeric) instead.",
            category: CodeCategory::Format,
            bad_examples: &[
                "m.[3243A>G/T]//[3243A>C/G]",
                "m.[3243A>G//T]/[3243A>C//G]",
                "NM_000088.3:c.[100A>G/200T>C]",
                "NM_000088.3:c.[100A>G//200T>C]",
            ],
            good_examples: &[
                "NM_000088.3:c.[100A>G;200T>C]",
                "NM_000088.3:c.100A>G/NM_000088.3:c.200T>C",
                "NM_000088.3:c.100A>G//NM_000088.3:c.200T>C",
                "NM_000088.3:c.123=/A>G",
                "NM_000088.3:c.123=//A>G",
            ],
            mode_behavior: Some(ModeBehavior::always_reject()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/",
            ),
            related_codes: &["W3004"],
        },
    );

    map.insert(
        "W3021",
        CodeInfo {
            code: "W3021",
            name: "ProteinBracketedAaInsertion",
            summary: "Bracketed amino-acid list inside a protein insertion edit.",
            explanation:
                "HGVS v21's protein insertion notation concatenates 3-letter amino-acid codes \
                with no separators (`p.Arg97_Trp98insAlaPro`). The bracketed form \
                `p.Arg97_Trp98ins[Ala;Pro]` is not in the spec: brackets are reserved for \
                alleles at the variant level, not for amino-acid lists inside an edit. \
                ferro cannot auto-rewrite the bracketed form because mixing 1-letter and \
                3-letter codes inside `[...]` is ambiguous; all modes reject with a hint \
                pointing at the canonical `insAlaPro` shape.",
            category: CodeCategory::Format,
            bad_examples: &[
                "NP_000088.3:p.Arg97_Trp98ins[Ala;Pro]",
                "NP_000088.3:p.Arg97_Trp98ins[A;P]",
            ],
            good_examples: &[
                "NP_000088.3:p.Arg97_Trp98insAlaPro",
                "NP_000088.3:p.Arg97_Trp98insAla",
            ],
            mode_behavior: Some(ModeBehavior::always_reject()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/protein/variant/insertion/",
            ),
            related_codes: &["E1004"],
        },
    );

    map.insert(
        "W3022",
        CodeInfo {
            code: "W3022",
            name: "InitiatorMetCanonicalization",
            summary: "Canonical p.dup form includes the initiator methionine (p.Met1).",
            explanation: "HGVS Prioritization (general.md §3) requires duplications to be \
                preferred over insertions when both can describe the same change. \
                When ferro canonicalizes a protein insertion (or affix-trimmed \
                delins) into a duplication whose interval includes position 1, \
                the result touches the initiator methionine. The spec uses \
                Met1-inclusive ranges elsewhere (deletion.md §63-65: \
                p.(Met1_Leu46del)), so the duplication form is permitted; however \
                the predicted protein consequence may also be described per the \
                substitution rule for start-codon variants (substitution.md §45-65: \
                p.0, p.0?, or p.(Met1?)). Consumers should consider both \
                interpretations.",
            category: CodeCategory::Format,
            bad_examples: &["NP_000088.3:p.Met1_Lys2insMet"],
            good_examples: &["NP_000088.3:p.Met1dup"],
            mode_behavior: Some(ModeBehavior::warn_accept()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/general/"),
            related_codes: &["W3021"],
        },
    );

    map.insert(
        "W3023",
        CodeInfo {
            code: "W3023",
            name: "DupSizeSuffix",
            summary: "Duplication described with a size-count suffix instead of a position range.",
            explanation: "HGVS requires both endpoints of a multi-residue duplication to be named \
                (`DNA/duplication.md:140`). The legacy form `g.123dup6` names only one, and the \
                spec explicitly declines to disambiguate it (`g.123_128dup` vs `g.124_129dup`), \
                so there is no safe repair: a point-anchored `dup<N>` is rejected in every mode. \
                A range anchor that merely repeats the size (`c.20_21dup2`) already names both \
                endpoints and only warns.",
            category: CodeCategory::Format,
            bad_examples: &["NM_004006.2:c.20_21dup2", "NC_000023.11:g.123dup6"],
            good_examples: &["NM_004006.2:c.20_21dup", "NC_000023.11:g.123_124dup"],
            mode_behavior: Some(ModeBehavior::warn_accept()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/duplication/",
            ),
            related_codes: &["W3011", "W3016", "W3024"],
        },
    );

    map.insert(
        "W3024",
        CodeInfo {
            code: "W3024",
            name: "DupExplicitSeq",
            summary: "Duplication described with the explicit duplicated nucleotide sequence.",
            explanation:
                "Per HGVS recommendations (`DNA/duplication.md`), the duplicated sequence should \
                NOT be appended after `dup` because (1) the description is longer, (2) it contains \
                redundant information already implied by the position range, and (3) the chances \
                to make an error increase (e.g. writing `dupG` when the reference is `T`). The \
                canonical form is `c.20_23dup`, not `c.20_23dupTAGA`. ferro rewrites the explicit \
                form to the canonical form in lenient and silent modes; strict mode rejects.",
            category: CodeCategory::Format,
            bad_examples: &["NM_004006.2:c.20_23dupTAGA", "NM_004006.2:c.20dupT"],
            good_examples: &["NM_004006.2:c.20_23dup", "NM_004006.2:c.20dup"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/duplication/",
            ),
            related_codes: &["W3011", "W3016", "W3023", "W3025"],
        },
    );

    map.insert(
        "W3025",
        CodeInfo {
            code: "W3025",
            name: "DelExplicitSeq",
            summary: "Deletion described with the explicit deleted nucleotide sequence.",
            explanation:
                "Per HGVS recommendations (`DNA/deletion.md`), the deleted sequence should NOT \
                be appended after `del` because (1) the description is longer, (2) it contains \
                redundant information already implied by the position range, and (3) the chances \
                to make an error increase (e.g. writing `delG` when the reference is `A`). The \
                canonical form is `g.33344590_33344592del`, not `g.33344590_33344592delGAT`. \
                ferro rewrites the explicit form to the canonical form in lenient and silent \
                modes; strict mode rejects. This warning applies only to pure deletions; \
                deletion-insertions (`delins`) are excluded from this check.",
            category: CodeCategory::Format,
            bad_examples: &[
                "NC_000023.11:g.33344590_33344592delGAT",
                "NC_000023.11:g.33344591delA",
            ],
            good_examples: &[
                "NC_000023.11:g.33344590_33344592del",
                "NC_000023.11:g.33344591del",
            ],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/",
            ),
            related_codes: &["W3011", "W3016", "W3023", "W3024"],
        },
    );

    map.insert(
        "W3026",
        CodeInfo {
            code: "W3026",
            name: "NonConformantBracketCardinality",
            summary: "Standalone single-member allele bracket (e.g. c.[76A>C], p.[=]).",
            explanation: "`[ ]` is allele syntax. The HGVS spec admits only two conformant \
                shapes, identically across every coordinate system (c/g/n/m/o/r/p): one \
                bracket group with two or more cis members (`c.[76A>C;88G>T]`), or two or \
                more trans groups (`c.[76A>C];[88G>T]`). A single variant wrapped in \
                brackets standalone is explicitly invalid (DNA/alleles.md: the notation \
                `c.[76A>C]` without describing the second allele is misleading; the \
                recommended form is `c.[76A>C];[76=]`), and the pure markers `=`/`?`/`0` \
                are valid only as a whole group inside a multi-group trans construct \
                (`c.[2376G>C];[=]`), never standalone. The canonical repair drops the \
                redundant brackets (`c.[76A>C]` -> `c.76A>C`). Strict mode rejects; lenient \
                unwraps and warns; silent unwraps without a warning.",
            category: CodeCategory::Format,
            bad_examples: &[
                "NM_000088.3:c.[76A>C]",
                "NC_000001.11:g.[1000G>A]",
                "NP_000079.2:p.[=]",
            ],
            good_examples: &[
                "NM_000088.3:c.76A>C",
                "NM_000088.3:c.[76A>C;88G>T]",
                "NM_000088.3:c.[76A>C];[76=]",
            ],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/",
            ),
            related_codes: &["W3004", "W3021"],
        },
    );

    // --- Position/Range Warnings (W4xxx) ---

    map.insert(
        "W4001",
        CodeInfo {
            code: "W4001",
            name: "SwappedPositions",
            summary: "Swapped/inverted range positions.",
            explanation: "In a range, the start position should be less than the end position. \
                Ranges like 'c.200_100del' are corrected to 'c.100_200del' in lenient/silent modes.",
            category: CodeCategory::Position,
            bad_examples: &["c.200_100del", "g.50000_10000dup"],
            good_examples: &["c.100_200del", "g.10000_50000dup"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: None,
            related_codes: &["E3003"],
        },
    );

    // W4002 PositionZero was retired in issue #269: `c.0…` is rejected
    // by `detect_position_zero` as `E1003 InvalidPosition` (a hard error),
    // never as a soft-validation warning. Keeping a duplicate identity in
    // the registry was a documentation artefact only.

    map.insert(
        "W4003",
        CodeInfo {
            code: "W4003",
            name: "SinglePositionRange",
            summary: "Single-position range used where a single position is canonical.",
            explanation: "HGVS spec requires `position(s)_deleted` (and the analogous \
                `positions_duplicated` / `positions_inverted`) to contain two different positions. \
                Forms like `c.123_123del`, `c.123_123dup`, and `c.100_100inv` should be written \
                with the single-position form `c.123del` / `c.123dup` / `c.100inv` respectively. \
                In lenient/silent modes, ferro collapses the redundant range and warns once.",
            category: CodeCategory::Position,
            bad_examples: &[
                "NM_000088.3:c.123_123del",
                "NM_000088.3:c.123_123dup",
                "NM_000088.3:c.100_100inv",
            ],
            good_examples: &[
                "NM_000088.3:c.123del",
                "NM_000088.3:c.123dup",
                "NM_000088.3:c.100inv",
            ],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/",
            ),
            related_codes: &["W3011"],
        },
    );

    map.insert(
        "W4004",
        CodeInfo {
            code: "W4004",
            name: "PositionPastEnd",
            summary: "Position lies past its coordinate sub-axis bound.",
            explanation: "A `c.`, `c.*N`, `c.-N`, or `n.` position must reference a position \
                that exists in the resolved transcript. Plain integer `c.<N>` positions cannot \
                exceed the CDS length; 3'UTR `c.*N` positions cannot exceed the post-CDS \
                transcript length; 5'UTR `c.-N` positions cannot exceed the pre-CDS transcript \
                length; `n.<N>` positions cannot exceed the transcript length. Strict mode \
                rejects past-end inputs; lenient mode emits W4004 and short-circuits \
                normalize() to the canonical variant (no further shifting). Has no safe \
                auto-correction — the user could have meant a different position or a \
                different reference. Intronic-offset magnitude checks remain out of scope \
                (their bounds depend on intron size, which requires alignment data this \
                helper does not consult); bare-transcript intronic positions (where the \
                reference form itself is invalid) are rejected by W4007 (#486).",
            category: CodeCategory::Position,
            bad_examples: &[
                "NM_001001656.1:c.946G>C (CDS-end = 945)",
                "NM_001001656.1:c.946dup",
                "NM_001001656.1:c.935_946del",
            ],
            good_examples: &["NM_001001656.1:c.945G>C", "NM_001001656.1:c.*1G>C"],
            mode_behavior: Some(ModeBehavior::warn_accept()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/general/"),
            related_codes: &["E3001", "W3016", "W4007"],
        },
    );

    map.insert(
        "W4005",
        CodeInfo {
            code: "W4005",
            name: "UnresolvableCentromere",
            summary: "Centromere position cannot be resolved to a coordinate.",
            explanation: "The `cen` marker denotes a chromosome's centromere — an \
                assembly-annotated region, not a nucleotide derivable from the reference \
                sequence — so ferro cannot place it on a concrete coordinate to normalize. \
                (`pter`/`qter` resolve to the first/last nucleotide and do normalize.) \
                Strict mode rejects with FerroError::InvalidCoordinates; lenient/silent \
                modes preserve the input, with the UnresolvableSpecialPosition warning \
                surfaced by normalize_with_diagnostics. See #488.",
            category: CodeCategory::Position,
            bad_examples: &["NC_000001.11:g.cendel"],
            good_examples: &["NC_000001.11:g.pterdel (pter resolves to base 1)"],
            mode_behavior: Some(ModeBehavior::warn_accept()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/general/"),
            related_codes: &["W4004"],
        },
    );

    map.insert(
        "W4006",
        CodeInfo {
            code: "W4006",
            name: "TranscriptFlankNotDescribable",
            summary:
                "Telomere marker on a genomic-reference c. denotes a transcript-flank position.",
            explanation: "On a genomic reference (NG_/NC_/LRG_) c. description, `pter`/`qter` \
                denote the genomic parent's terminus, which lies in the 5'/3' transcript flank — \
                beyond the transcript's terminal exons. HGVS does not permit numbering flanking \
                nucleotides in c. coordinates (background/numbering.md transcript-flanking; the \
                flank-numbering proposal was rejected, see consultation/open-issues.md). Strict \
                mode rejects; lenient/silent preserve the input. Use the genomic g. form instead.",
            category: CodeCategory::Position,
            bad_examples: &["NG_012337.1(NM_003002.2):c.pterdel"],
            good_examples: &["NG_012337.1:g.pterdel (genomic form)"],
            mode_behavior: Some(ModeBehavior::warn_accept()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/background/numbering/"),
            related_codes: &["W4005"],
        },
    );

    map.insert(
        "W4007",
        CodeInfo {
            code: "W4007",
            name: "IntronicOnBareTranscript",
            summary: "Intronic offset on a bare transcript reference (no genomic context).",
            explanation: "A coding (c.) or non-coding (n.) DNA reference sequence does not \
                contain introns, so it cannot describe an intronic position (an offset `+M`/`-M` \
                from a splice site). The bare forms `NM_004006.2:c.357+1G>A` and \
                `NR_038420.1:n.100+10del` are spec-invalid (background/refseq.md): an intronic \
                description must supply a genomic reference, e.g. `NG_012337.1(NM_004006.2):c.357+1G>A`, \
                `NC_…(NM_…):c.357+1`, or `LRG_199t1:c.357+1`. Strict mode rejects (mapped to \
                EINTRONIC); lenient mode emits W4007 and returns the existing value unchanged; \
                silent accepts. Has no safe auto-correction — ferro cannot synthesize the \
                genomic reference.",
            category: CodeCategory::Position,
            bad_examples: &["NM_004006.2:c.357+1G>A", "NR_038420.1:n.100+10del"],
            good_examples: &[
                "NG_012337.1(NM_004006.2):c.357+1G>A",
                "LRG_199t1:c.357+1G>A",
            ],
            mode_behavior: Some(ModeBehavior::warn_accept()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/background/refseq/"),
            related_codes: &["W4004", "E4001"],
        },
    );

    // --- Semantic Warnings (W5xxx) ---

    map.insert(
        "W5001",
        CodeInfo {
            code: "W5001",
            name: "RefSeqMismatch",
            summary: "Reference sequence mismatch.",
            explanation: "The reference base(s) stated in the HGVS expression do not match \
                the actual reference sequence. In lenient/silent modes, normalization proceeds \
                using the actual reference sequence, but a warning is always emitted.",
            category: CodeCategory::Semantic,
            bad_examples: &["c.100G>A (when reference is T at position 100)"],
            good_examples: &["c.100T>A (matches actual reference)"],
            mode_behavior: Some(ModeBehavior::always_warn_if_not_rejected()),
            hgvs_spec_url: None,
            related_codes: &["E3002", "W3001"],
        },
    );

    map.insert(
        "W5002",
        CodeInfo {
            code: "W5002",
            name: "OverlapConflictingEdits",
            summary: "Two or more cis-allele edits share identical reference bounds.",
            explanation: "Two or more edits in a cis allele point at the same single \
                base or range with identical (start, end) bounds. The HGVS spec does \
                not define a canonical form for this case; ferro preserves the input \
                verbatim and emits this warning so callers can decide whether to flag \
                or reject. If the HGVS Variant Working Group rules these inputs invalid, \
                change the mode_behavior here to promote this code to an error in Strict.",
            category: CodeCategory::Semantic,
            bad_examples: &[
                "g.[100G>A;100A>C]",
                "g.[100del;100A>C]",
                "c.[100_103del;100_103inv]",
            ],
            good_examples: &["g.[100A>C;101A>C]", "c.[762_768del;767_774dup]"],
            mode_behavior: Some(ModeBehavior::always_warn_if_not_rejected()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/general/"),
            related_codes: &["W5001"],
        },
    );

    map.insert(
        "W5003",
        CodeInfo {
            code: "W5003",
            name: "VariantExceedsReference",
            summary: "Variant position range exceeds the available reference sequence.",
            explanation: "The provider returned fewer bytes than the HGVS interval span. \
                Per HGVS spec (recommendations/background/refseq.md \u{00A7}43), \"the entirety \
                of the variant sequence must be encompassed by the selected reference \
                sequence.\" Biocommons hgvs raises HGVSInvalidVariantError for this shape. \
                Strict mode rejects with FerroError::VariantExceedsReference; lenient mode \
                emits this warning (also published as `CanonicalSplitSkipped` in the \
                NormalizationWarning enum) and preserves the input verbatim; silent mode \
                preserves the input without warning. Common cause: a version-fallback \
                substituted a shorter reference (e.g. NG_032871.1 \u{2192} NG_032871.2), \
                or the variant positions are genuinely past the reference end.",
            category: CodeCategory::Semantic,
            bad_examples: &[
                "NG_032871.1:g.32476_53457delinsAATTAAGGTATA (when reference is shorter than 53457)",
            ],
            good_examples: &[
                "(no auto-correct; use a complete reference sequence per refseq.md \u{00A7}43)",
            ],
            mode_behavior: Some(ModeBehavior::warn_accept()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/background/refseq/#reference-sequences",
            ),
            related_codes: &["E3001", "W3016"],
        },
    );

    map.insert(
        "W5004",
        CodeInfo {
            code: "W5004",
            name: "IncompleteCdsStartReference",
            summary: "Transcript has a 5\u{2032}-incomplete CDS (cds_start_NF); not an \
                HGVS-recommended reference for c./p. description.",
            explanation: "The transcript's 5\u{2032} CDS is annotated incomplete (Ensembl \
                `cds_start_NF`): the ATG start codon is not present in the transcript record, \
                so `c.1`/`p.1` are undefined relative to it. Per HGVS recommendations, a \
                reference sequence used for coding (`c.`) or protein (`p.`) description must \
                have a complete, defined CDS; a 5\u{2032}-incomplete transcript is not an \
                HGVS-recommended reference for that purpose. Strict mode rejects `c.`/`p.` \
                variants against such a transcript; lenient mode emits this warning and \
                accepts the input unchanged; silent mode accepts without warning. Use the \
                genomic (`g.`) or non-coding (`n.`) representation instead, which do not \
                depend on a defined CDS start.",
            category: CodeCategory::Semantic,
            bad_examples: &["ENST00000381176.3:c.10A>G (transcript is cds_start_NF)"],
            good_examples: &[
                "NC_000001.11:g.94129A>G (genomic representation)",
                "ENST00000381176.3:n.100A>G (non-coding representation)",
            ],
            mode_behavior: Some(ModeBehavior::warn_accept()),
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/background/refseq/"),
            related_codes: &[],
        },
    );

    // =========================================================================
    // LOADER WARNINGS / ERRORS (E-LOAD-* / W-LOAD-*)
    // =========================================================================

    // --- Loader Errors (E-LOAD-*) ---

    map.insert(
        "E-LOAD-001",
        CodeInfo {
            code: "E-LOAD-001",
            name: "MalformedRecord",
            summary: "A required GFF/GTF column failed to parse and the offending row was dropped.",
            explanation: "During GFF3 or GTF loading, a record could not be interpreted because a \
                mandatory column was malformed: for example, a non-integer coordinate, a missing \
                tab delimiter, an invalid strand character, or an out-of-range phase value (must \
                be 0, 1, or 2). The row is silently discarded and the loader continues. \
                See the GFF/GTF loader design (§6 Stage 2 — record parsing).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E-LOAD-002", "W-LOAD-010"],
        },
    );

    map.insert(
        "E-LOAD-002",
        CodeInfo {
            code: "E-LOAD-002",
            name: "UnsupportedFormat",
            summary: "The input file could not be classified as GFF3 or GTF after extension and content sniffing.",
            explanation: "ferro attempts to detect the annotation format by examining the file \
                extension and scanning for format-specific header directives (e.g. `##gff-version 3` \
                for GFF3, `gene_id` attributes for GTF). When neither heuristic succeeds the file \
                is treated as unsupported. Currently reserved — Phase 1 falls back to GFF3 when \
                detection is ambiguous (see W-LOAD-001). \
                See the GFF/GTF loader design (§6 Stage 1 — format detection).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &[],
        },
    );

    map.insert(
        "E-LOAD-103",
        CodeInfo {
            code: "E-LOAD-103",
            name: "StrandRequired",
            summary: "A transcript has an unspecified or unknown strand and was dropped.",
            explanation: "HGVS coordinate arithmetic requires a defined strand (+ or −) for every \
                transcript. GFF3 allows `.` (unspecified) and `?` (unknown) as strand values; \
                transcripts carrying either of these values cannot be represented in ferro's \
                data model and are discarded. If the strand is determinable from context, \
                fix the source annotation before loading. \
                See the GFF/GTF loader design (§6 Stage 4 — transcript assembly).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["W-LOAD-100"],
        },
    );

    // --- Loader Warnings (W-LOAD-*) ---

    map.insert(
        "W-LOAD-001",
        CodeInfo {
            code: "W-LOAD-001",
            name: "UnknownFormat",
            summary: "Format could not be determined with high confidence; defaulted to GFF3.",
            explanation: "ferro examines the file extension and the presence of `##gff-version` \
                or GTF-style attribute syntax to decide which parser to invoke. When neither \
                indicator is present or they conflict, ferro falls back to GFF3 and emits this \
                warning. Inspect the file extension and ensure a `##gff-version 3` (or `2`) \
                directive is present on the first line. \
                See the GFF/GTF loader design (§6 Stage 1 — format detection).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &[],
        },
    );

    map.insert(
        "W-LOAD-002",
        CodeInfo {
            code: "W-LOAD-002",
            name: "InlineFastaIgnored",
            summary: "A `##FASTA` directive was encountered; inline sequences are not parsed.",
            explanation: "GFF3 files may embed reference sequences after a `##FASTA` line. \
                ferro's loader does not parse inline FASTA sections — if reference sequences are \
                needed for CDS validation, supply them via `--fasta` separately. \
                Currently reserved — Phase 1 does not yet detect the `##FASTA` directive. \
                See the GFF/GTF loader design (§6 Stage 2 — record parsing).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &[],
        },
    );

    map.insert(
        "W-LOAD-010",
        CodeInfo {
            code: "W-LOAD-010",
            name: "OrphanFeature",
            summary: "A feature references a parent ID that is not present in the file; the orphan is dropped.",
            explanation: "In GFF3, child features reference their parent via the `Parent=` \
                attribute; in GTF, children reference `transcript_id`. When the referenced \
                parent record is absent — common with truncated or region-filtered input files — \
                the child feature cannot be assembled into a transcript and is discarded. \
                Check that the source file is complete and that all parent features are present. \
                See the GFF/GTF loader design (§6 Stage 3 — parent-child linking).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["E-LOAD-001"],
        },
    );

    map.insert(
        "W-LOAD-100",
        CodeInfo {
            code: "W-LOAD-100",
            name: "TranscriptWithoutExons",
            summary: "A transcript-like feature has no exon, UTR, or CDS children and was dropped.",
            explanation: "After parent-child linking, a transcript record (mRNA, transcript, \
                ncRNA, lincRNA, etc.) was found with no child features from which an exon \
                structure could be derived. Without at least one exon, UTR, or CDS interval \
                there is nothing to build the transcript model from, so the transcript is \
                discarded. Verify that the source file contains matching child records with the \
                correct `Parent=` or `transcript_id` attributes. \
                See the GFF/GTF loader design (§6 Stage 4 — transcript assembly).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["W-LOAD-101", "E-LOAD-103"],
        },
    );

    map.insert(
        "W-LOAD-101",
        CodeInfo {
            code: "W-LOAD-101",
            name: "GeneAsTranscript",
            summary: "A `gene` feature with direct CDS children (no mRNA/transcript) is treated as a transcript.",
            explanation: "Some prokaryotic GFF3 files attach CDS records directly to a `gene` \
                feature without an intervening mRNA or transcript record. ferro detects this \
                pattern and promotes the gene to act as its own transcript, preserving the CDS \
                coordinates. If this was unintentional, restructure the annotation to include an \
                explicit transcript-level feature. \
                See the GFF/GTF loader design (§6 Stage 4 — transcript assembly).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["W-LOAD-100"],
        },
    );

    map.insert(
        "W-LOAD-110",
        CodeInfo {
            code: "W-LOAD-110",
            name: "PhaseAppliedToCdsStart",
            summary: "The leading CDS feature has a non-zero phase; `cds_start_genomic` was shifted by the phase value.",
            explanation: "GFF3 phase values indicate how many bases at the start of a CDS feature \
                should be skipped to reach the first complete in-frame codon. When the 5'-most CDS \
                record on a transcript carries a phase of 1 or 2, the loader advances \
                `cds_start_genomic` by that number of bases so that ferro's CDS coordinates \
                begin at a codon boundary. This is expected for fragmented or partial CDS \
                annotations but may indicate a truncated transcript. \
                See the GFF/GTF loader design (§6 Stage 4 — CDS start adjustment).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["W-LOAD-111"],
        },
    );

    map.insert(
        "W-LOAD-111",
        CodeInfo {
            code: "W-LOAD-111",
            name: "PhaseUnavailable",
            summary:
                "The leading CDS feature has no phase recorded; the unshifted CDS start is used.",
            explanation: "When the 5'-most CDS record on a transcript has a missing or `.` phase \
                field, ferro cannot determine whether the genomic start is offset from a codon \
                boundary. The loader falls back to using the unshifted coordinate, which may be \
                incorrect if the annotation was intended to encode a non-zero phase. Inspect the \
                source annotation and add an explicit phase value where possible. \
                See the GFF/GTF loader design (§6 Stage 4 — CDS start adjustment).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["W-LOAD-110"],
        },
    );

    map.insert(
        "W-LOAD-112",
        CodeInfo {
            code: "W-LOAD-112",
            name: "StopCodonAssumed",
            summary: "No explicit `stop_codon` record was found; `cds_end_genomic` was extended by 3 bp on the 3' side.",
            explanation: "ferro's CDS model is stop-codon-inclusive: the 3' boundary of the CDS \
                encompasses the stop codon. GFF3 files often encode the stop codon as a separate \
                `stop_codon` feature rather than including it in the CDS interval. When no \
                `stop_codon` record is present, the loader extends `cds_end_genomic` by 3 bp \
                toward the 3' end, clipped at the transcript boundary. If the extension would \
                overshoot the last exon, the CDS end may be slightly incorrect. Verify that the \
                transcript end is consistent with the annotated stop codon. \
                See the GFF/GTF loader design (§6 Stage 4 — stop codon handling).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["W-LOAD-110"],
        },
    );

    map.insert(
        "W-LOAD-120",
        CodeInfo {
            code: "W-LOAD-120",
            name: "ConflictingProteinId",
            summary: "CDS fragments of a single transcript reported disagreeing `protein_id` values; `protein_id` was left unset.",
            explanation: "GFF3 and GTF convention is that every CDS fragment of a single \
                transcript carries the same `protein_id` (they all encode parts of the same \
                polypeptide). When the loader observes two or more distinct `protein_id` \
                values on the CDS records of one transcript, it cannot pick one safely \
                without making the output depend on input record order. The transcript is \
                still loaded, but its `protein_id` is left unset so the projector's \
                transcript-id-based protein-accession fallback (see issue #310) takes over. \
                Inspect the source annotation and reconcile the disagreement upstream. \
                See the GFF/GTF loader design (§6 Stage 4 — transcript assembly).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &[],
        },
    );

    map.insert(
        "W-LOAD-200",
        CodeInfo {
            code: "W-LOAD-200",
            name: "CdsLengthNotMod3",
            summary: "FASTA validation found a CDS whose length is not divisible by 3.",
            explanation: "When FASTA-based validation is enabled (Phase 4), ferro checks that the \
                assembled CDS length (3' end inclusive of the stop codon) is a multiple of 3. A \
                non-multiple-of-3 length suggests a frameshift, an incomplete CDS record, or an \
                annotation error in the source file. Currently reserved — FASTA validation is not \
                yet shipped (Phase 4). \
                See the GFF/GTF loader design (§6 Stage 5 — FASTA validation).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["W-LOAD-201"],
        },
    );

    map.insert(
        "W-LOAD-201",
        CodeInfo {
            code: "W-LOAD-201",
            name: "NonCanonicalStartCodon",
            summary:
                "FASTA validation found that the first CDS codon is not in {ATG, CTG, GTG, TTG}.",
            explanation: "When FASTA-based validation is enabled (Phase 4), ferro reads the first \
                three bases of the CDS from the reference FASTA and checks for a canonical start \
                codon (ATG) or one of the known alternative initiators (CTG, GTG, TTG). Any other \
                triplet suggests a mis-annotated CDS start or a sequencing artefact. Currently \
                reserved — FASTA validation is not yet shipped (Phase 4). \
                See the GFF/GTF loader design (§6 Stage 5 — FASTA validation).",
            category: CodeCategory::Io,
            bad_examples: &[],
            good_examples: &[],
            mode_behavior: None,
            hgvs_spec_url: None,
            related_codes: &["W-LOAD-200"],
        },
    );

    map
}

/// Get information about a specific code.
/// Accepts both uppercase and lowercase codes (e.g., "W1001" or "w1001").
pub fn get_code_info(code: &str) -> Option<&'static CodeInfo> {
    let uppercase = code.to_uppercase();
    get_registry().get(uppercase.as_str())
}

/// List all codes.
pub fn list_all_codes() -> Vec<&'static CodeInfo> {
    let mut codes: Vec<_> = get_registry().values().collect();
    codes.sort_by_key(|c| c.code);
    codes
}

/// List all error codes (E-prefix).
pub fn list_error_codes() -> Vec<&'static CodeInfo> {
    let mut codes: Vec<_> = get_registry().values().filter(|c| c.is_error()).collect();
    codes.sort_by_key(|c| c.code);
    codes
}

/// List all warning codes (W-prefix).
pub fn list_warning_codes() -> Vec<&'static CodeInfo> {
    let mut codes: Vec<_> = get_registry().values().filter(|c| c.is_warning()).collect();
    codes.sort_by_key(|c| c.code);
    codes
}

/// List codes by category.
pub fn list_codes_by_category(category: CodeCategory) -> Vec<&'static CodeInfo> {
    let mut codes: Vec<_> = get_registry()
        .values()
        .filter(|c| c.category == category)
        .collect();
    codes.sort_by_key(|c| c.code);
    codes
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error_handling::codes::ModeAction;

    #[test]
    fn test_registry_has_codes() {
        assert!(!get_registry().is_empty());
    }

    #[test]
    fn test_get_code_info() {
        let e1001 = get_code_info("E1001");
        assert!(e1001.is_some());
        assert_eq!(e1001.unwrap().name, "InvalidAccession");

        let w1001 = get_code_info("W1001");
        assert!(w1001.is_some());
        assert_eq!(w1001.unwrap().name, "LowercaseAminoAcid");

        assert!(get_code_info("NOTEXIST").is_none());
    }

    #[test]
    fn test_list_error_codes() {
        let errors = list_error_codes();
        assert!(!errors.is_empty());
        for code in errors {
            assert!(code.code.starts_with('E'));
        }
    }

    #[test]
    fn test_list_warning_codes() {
        let warnings = list_warning_codes();
        assert!(!warnings.is_empty());
        for code in warnings {
            assert!(code.code.starts_with('W'));
        }
    }

    #[test]
    fn test_list_all_codes_sorted() {
        let codes = list_all_codes();
        let sorted: Vec<_> = codes.iter().map(|c| c.code).collect();
        let mut check = sorted.clone();
        check.sort();
        assert_eq!(sorted, check);
    }

    #[test]
    fn test_warning_codes_have_mode_behavior() {
        let warnings = list_warning_codes();
        for code in warnings {
            // Loader codes (W-LOAD-*) do not participate in the parser's
            // strict/lenient/silent configuration, so they may omit mode_behavior.
            if code.code.starts_with("W-LOAD-") {
                continue;
            }
            assert!(
                code.mode_behavior.is_some(),
                "Warning {} should have mode_behavior",
                code.code
            );
        }
    }

    #[test]
    fn test_error_codes_no_mode_behavior() {
        let errors = list_error_codes();
        for code in errors {
            assert!(
                code.mode_behavior.is_none(),
                "Error {} should not have mode_behavior",
                code.code
            );
        }
    }

    #[test]
    fn test_format_terminal() {
        let code = get_code_info("W1001").unwrap();
        let output = code.format_terminal(false);
        assert!(output.contains("W1001"));
        assert!(output.contains("LowercaseAminoAcid"));
        assert!(output.contains("ferro parse --ignore"));
    }

    #[test]
    fn test_format_json() {
        let code = get_code_info("E1001").unwrap();
        let json = code.format_json();
        assert!(json.contains("\"code\":\"E1001\""));
        assert!(json.contains("\"name\":\"InvalidAccession\""));
    }

    #[test]
    fn test_format_markdown() {
        let code = get_code_info("W1001").unwrap();
        let md = code.format_markdown();
        assert!(md.contains("## W1001:"));
        assert!(md.contains("| Mode | Action |"));
    }

    #[test]
    fn test_list_codes_by_category() {
        let parse_codes = list_codes_by_category(CodeCategory::Parse);
        assert!(!parse_codes.is_empty());
        for code in parse_codes {
            assert_eq!(code.category, CodeCategory::Parse);
        }
    }

    #[test]
    fn loader_error_codes_are_registered() {
        assert!(get_code_info("E-LOAD-001").is_some());
        assert!(get_code_info("E-LOAD-002").is_some());
        assert!(get_code_info("E-LOAD-103").is_some());
    }

    #[test]
    fn transcript_version_not_exact_is_registered() {
        // #809: E2004 must resolve so `ferro explain E2004` works (it shares
        // the engine's `code()` mapping but is a separate reference-family code).
        let info = get_code_info("E2004").expect("E2004 registered");
        assert_eq!(info.name, "TranscriptVersionNotExact");
        assert_eq!(info.category, CodeCategory::Reference);
    }

    #[test]
    fn transcript_sequence_unreconstructable_is_registered() {
        // #807: E2005 must resolve so `ferro explain E2005` works (the enum→registry
        // edge is not test-enforced, so a missing registry entry is otherwise silent).
        let info = get_code_info("E2005").expect("E2005 registered");
        assert_eq!(info.name, "TranscriptSequenceUnreconstructable");
        assert_eq!(info.category, CodeCategory::Reference);
    }

    #[test]
    fn incomplete_cds_start_reference_is_registered() {
        // #972: W5004 must resolve so `ferro explain W5004` works and Task 5's
        // c./p.-over-cds_start_NF gating can look up the mode behavior.
        let info = get_code_info("W5004").expect("W5004 registered");
        assert_eq!(info.name, "IncompleteCdsStartReference");
        assert_eq!(info.category, CodeCategory::Semantic);

        let behavior = info.mode_behavior.expect("W5004 has mode_behavior");
        assert_eq!(behavior, ModeBehavior::warn_accept());
        assert_eq!(behavior.strict, ModeAction::Reject);
        assert_eq!(behavior.lenient, ModeAction::WarnAndAccept);
        assert_eq!(behavior.silent, ModeAction::Accept);
    }

    #[test]
    fn loader_warning_codes_are_registered() {
        for code in &[
            "W-LOAD-001",
            "W-LOAD-002",
            "W-LOAD-010",
            "W-LOAD-100",
            "W-LOAD-101",
            "W-LOAD-110",
            "W-LOAD-111",
            "W-LOAD-112",
            "W-LOAD-120",
            "W-LOAD-200",
            "W-LOAD-201",
        ] {
            assert!(
                get_code_info(code).is_some(),
                "loader warning code {} not registered",
                code
            );
        }
    }

    #[test]
    fn loader_codes_have_consistent_metadata() {
        let info = get_code_info("W-LOAD-100").unwrap();
        assert_eq!(info.name, "TranscriptWithoutExons");
        assert!(!info.summary.is_empty());
        assert!(!info.explanation.is_empty());
    }
}
