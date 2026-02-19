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
            explanation: "Positions must be valid integers or offset positions (for intronic variants). \
                Position 0 is never valid in HGVS. Negative positions indicate upstream of the start codon, \
                and *N positions indicate downstream of the stop codon.",
            category: CodeCategory::Parse,
            bad_examples: &["NM_000088.3:c.0A>G", "NM_000088.3:c.A>G"],
            good_examples: &["NM_000088.3:c.1A>G", "NM_000088.3:c.-10A>G"],
            mode_behavior: None,
            hgvs_spec_url: Some("https://hgvs-nomenclature.org/stable/recommendations/DNA/numbering/"),
            related_codes: &["E3001", "W4001", "W4002"],
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
            related_codes: &["E2002", "E2003"],
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
            related_codes: &["E3004"],
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
            explanation: "HGVS expressions should not contain spaces (except in gene symbols). \
                Extra spaces are removed in lenient/silent modes.",
            category: CodeCategory::Character,
            bad_examples: &["c.100 A>G", "NM_000088.3 : c.459A>G"],
            good_examples: &["c.100A>G", "NM_000088.3:c.459A>G"],
            mode_behavior: Some(ModeBehavior::standard_correctable()),
            hgvs_spec_url: None,
            related_codes: &[],
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

    map.insert(
        "W4002",
        CodeInfo {
            code: "W4002",
            name: "PositionZero",
            summary: "Invalid position zero.",
            explanation: "Position 0 does not exist in HGVS notation. Coding sequences start at \
                position 1 (ATG start codon), with upstream positions numbered as -1, -2, etc. \
                This error cannot be auto-corrected.",
            category: CodeCategory::Position,
            bad_examples: &["c.0A>G", "g.0del"],
            good_examples: &["c.1A>G", "c.-1A>G"],
            mode_behavior: Some(ModeBehavior::always_reject()),
            hgvs_spec_url: Some(
                "https://hgvs-nomenclature.org/stable/recommendations/DNA/numbering/",
            ),
            related_codes: &["E1003"],
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
}
