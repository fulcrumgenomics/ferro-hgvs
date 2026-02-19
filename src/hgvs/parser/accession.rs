//! Accession parsing
//!
//! Parses reference sequence accessions like NC_000001.11, NM_000088.3,
//! ENST00000012345.1, etc.

use crate::hgvs::variant::Accession;
use memchr::memchr;
use nom::{
    branch::alt,
    bytes::complete::{tag, take_while, take_while1},
    character::complete::digit1,
    combinator::{map, opt},
    sequence::preceded,
    IResult, Parser,
};

/// Parse an accession (e.g., "NM_000088.3" or "ENST00000012345.1" or "GRCh37(chr23)" or "P54802")
/// Optimized with byte-based dispatch to avoid trying all 5 alternatives
#[inline]
pub fn parse_accession(input: &str) -> IResult<&str, Accession> {
    let bytes = input.as_bytes();
    if bytes.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Eof,
        )));
    }

    // Dispatch based on first character(s)
    match bytes[0] {
        // Standard RefSeq accessions: NC_, NM_, NP_, NG_, NR_, NW_, NT_, XM_, XR_, XP_
        b'N' | b'X' => {
            // Check for underscore pattern (standard accession)
            if bytes.len() > 2 && bytes[2] == b'_' {
                if let Ok(result) = parse_standard_accession(input) {
                    return Ok(result);
                }
            }
            // Fall through to simple accession
            parse_simple_accession(input)
        }
        // Ensembl accessions: ENST, ENSG, ENSP, ENSE, ENSR
        b'E' => {
            if bytes.len() >= 4 && bytes[1] == b'N' && bytes[2] == b'S' {
                if let Ok(result) = parse_ensembl_accession(input) {
                    return Ok(result);
                }
            }
            parse_simple_accession(input)
        }
        // Assembly notation: GRCh37, GRCh38
        b'G' => {
            if bytes.len() >= 4 && bytes[1] == b'R' && bytes[2] == b'C' && bytes[3] == b'h' {
                if let Ok(result) = parse_assembly_accession(input) {
                    return Ok(result);
                }
            }
            parse_simple_accession(input)
        }
        // Assembly notation: hg18, hg19, hg38
        b'h' => {
            if bytes.len() >= 2 && bytes[1] == b'g' {
                if let Ok(result) = parse_assembly_accession(input) {
                    return Ok(result);
                }
            }
            parse_simple_accession(input)
        }
        // LRG accessions: LRG_XXX
        b'L' => {
            if bytes.len() >= 4 && bytes[1] == b'R' && bytes[2] == b'G' && bytes[3] == b'_' {
                if let Ok(result) = parse_standard_accession(input) {
                    return Ok(result);
                }
            }
            parse_simple_accession(input)
        }
        // Potential UniProt: uppercase letter followed by digits (P54802, Q8TAM1, A0A...)
        b'A'..=b'Z' => {
            // Try UniProt first if second char is digit or looks like UniProt pattern
            if bytes.len() >= 6 && bytes[1].is_ascii_digit() {
                if let Ok(result) = parse_uniprot_accession(input) {
                    return Ok(result);
                }
            }
            // Try standard accession for other letter prefixes (like AC_)
            if bytes.len() > 2 && bytes[2] == b'_' {
                if let Ok(result) = parse_standard_accession(input) {
                    return Ok(result);
                }
            }
            parse_simple_accession(input)
        }
        // Anything else goes to simple accession
        _ => parse_simple_accession(input),
    }
}

// Note: HGVS type prefixes (c., g., n., p., r., m., o.) are checked inline
// using byte matching for better performance - see parse_simple_accession()

/// Check if a character is valid in a SAM reference sequence name.
///
/// Per SAM spec v1, reference names may contain printable ASCII `[!-~]` except:
/// `\ , " ' ( ) [ ] { } < >`
/// Additionally, they cannot start with `*` or `=` (handled separately).
fn is_sam_refname_char(c: char) -> bool {
    // Printable ASCII range: ! (33) to ~ (126)
    let code = c as u32;
    if !(33..=126).contains(&code) {
        return false;
    }
    // Disallowed characters: \ , " ' ( ) [ ] { } < >
    !matches!(
        c,
        '\\' | ',' | '"' | '\'' | '(' | ')' | '[' | ']' | '{' | '}' | '<' | '>'
    )
}

/// Parse a simple chromosome/contig accession using SAM-compatible character restrictions.
///
/// This handles bare reference sequence names like "chr1", "contig", "1", etc.
/// Uses look-ahead to correctly identify the HGVS separator (`:` followed by type prefix).
///
/// SAM spec: Reference names may contain printable ASCII `[!-~]` except `\ , " ' ( ) [ ] { } < >`
/// and cannot start with `*` or `=`.
fn parse_simple_accession(input: &str) -> IResult<&str, Accession> {
    // Check first character: cannot be empty, and cannot start with * or =
    let first_char = input.chars().next().ok_or_else(|| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Alpha))
    })?;

    if first_char == '*' || first_char == '=' || !is_sam_refname_char(first_char) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Alpha,
        )));
    }

    // Find the HGVS separator: a colon followed by a valid type prefix
    // Use memchr for faster byte scanning
    let mut separator_pos = None;
    let mut search_start = 0;
    let input_bytes = input.as_bytes();

    while let Some(colon_offset) = memchr(b':', &input_bytes[search_start..]) {
        let colon_pos = search_start + colon_offset;
        let after_colon = &input_bytes[colon_pos + 1..];

        // Check if this colon is followed by a valid HGVS type prefix (X. format)
        // All prefixes are single letter + dot: c., g., n., p., r., m., o.
        let is_hgvs_prefix = after_colon.len() >= 2
            && after_colon[1] == b'.'
            && matches!(
                after_colon[0],
                b'c' | b'g' | b'n' | b'p' | b'r' | b'm' | b'o'
            );

        if is_hgvs_prefix {
            separator_pos = Some(colon_pos);
            break;
        }

        // Also accept colon followed by digit, '(', or '?' for type-less variants
        // e.g., AL513220.9:40902_43653del, NG_011403.2:(80027_96047)_(99154_121150)del
        // But only if there's no valid type prefix later in the string
        // (to avoid breaking HLA-style accessions like HLA-A*01:01:p.Arg100)
        let is_typeless_variant =
            !after_colon.is_empty() && matches!(after_colon[0], b'0'..=b'9' | b'(' | b'?');

        if is_typeless_variant {
            // Check if there's a type prefix later - if so, skip this colon
            let has_later_type_prefix = after_colon.windows(2).any(|w| {
                w[1] == b'.' && matches!(w[0], b'c' | b'g' | b'n' | b'p' | b'r' | b'm' | b'o')
            });
            if !has_later_type_prefix {
                separator_pos = Some(colon_pos);
                break;
            }
        }

        // This colon is part of the accession name, keep searching
        search_start = colon_pos + 1;
    }

    let separator_pos = separator_pos.ok_or_else(|| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Tag))
    })?;

    // Validate all characters in the accession name are SAM-compatible
    let accession_str = &input[..separator_pos];
    if accession_str.is_empty() || !accession_str.chars().all(is_sam_refname_char) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Alpha,
        )));
    }

    // Check for optional version suffix: ".N" at the end where N is digits
    let (name, version) = if let Some(dot_pos) = accession_str.rfind('.') {
        let potential_version = &accession_str[dot_pos + 1..];
        if !potential_version.is_empty() && potential_version.chars().all(|c| c.is_ascii_digit()) {
            let version_num = potential_version.parse::<u32>().ok();
            (&accession_str[..dot_pos], version_num)
        } else {
            (accession_str, None)
        }
    } else {
        (accession_str, None)
    };

    Ok((
        &input[separator_pos..],
        Accession::with_style(name.to_string(), String::new(), version, true),
    ))
}

/// Parse a UniProt-style accession (e.g., "P54802", "Q8TAM1", "A0A024R1R8")
///
/// Formats:
/// - 6-char: single letter + 5 alphanumeric (e.g., P54802)
/// - 10-char: letter + digit + letter + 2 alphanumeric + digit + letter + 2 alphanumeric + digit (e.g., A0A024R1R8)
///
/// UniProt accessions don't have versions or underscores
fn parse_uniprot_accession(input: &str) -> IResult<&str, Accession> {
    // Try 10-char format first (longer match)
    if let Ok(result) = parse_uniprot_10char(input) {
        return Ok(result);
    }
    // Fall back to 6-char format
    parse_uniprot_6char(input)
}

/// Parse 6-character UniProt accession (e.g., "P54802", "Q8TAM1")
fn parse_uniprot_6char(input: &str) -> IResult<&str, Accession> {
    let bytes = input.as_bytes();

    // Need at least 6 chars for accession + 3 for ":p." or ":p("
    if bytes.len() < 9 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Eof,
        )));
    }

    // First char must be uppercase letter
    if !bytes[0].is_ascii_uppercase() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Alpha,
        )));
    }

    // Next 5 chars must be alphanumeric
    for &b in &bytes[1..6] {
        if !b.is_ascii_alphanumeric() {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::AlphaNumeric,
            )));
        }
    }

    // Check what follows (must be ':p.' or ':p(')
    if bytes[6] != b':' || bytes[7] != b'p' || (bytes[8] != b'.' && bytes[8] != b'(') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Safe: all bytes verified as ASCII
    let first = &input[0..1];
    let number_part = &input[1..6];
    Ok((
        &input[6..],
        Accession::with_style(first.to_string(), number_part.to_string(), None, true),
    ))
}

/// Parse 10-character UniProt accession (e.g., "A0A024R1R8")
/// Format: [A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]
fn parse_uniprot_10char(input: &str) -> IResult<&str, Accession> {
    let bytes = input.as_bytes();

    // Need at least 10 chars for accession + 3 for ":p." or ":p("
    if bytes.len() < 13 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Eof,
        )));
    }

    // First char must be uppercase letter
    if !bytes[0].is_ascii_uppercase() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Alpha,
        )));
    }

    // Next 9 chars must be alphanumeric
    for &b in &bytes[1..10] {
        if !b.is_ascii_alphanumeric() {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::AlphaNumeric,
            )));
        }
    }

    // Check what follows (must be ':p.' or ':p(')
    if bytes[10] != b':' || bytes[11] != b'p' || (bytes[12] != b'.' && bytes[12] != b'(') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Safe: all bytes verified as ASCII
    let first = &input[0..1];
    let number_part = &input[1..10];
    Ok((
        &input[10..],
        Accession::with_style(first.to_string(), number_part.to_string(), None, true),
    ))
}

/// Parse assembly/chromosome notation (e.g., "GRCh37(chr23)")
fn parse_assembly_accession(input: &str) -> IResult<&str, Accession> {
    // Parse assembly name (alphanumeric, e.g., GRCh37, GRCh38, hg19, hg38)
    let (input, assembly) = take_while1(|c: char| c.is_ascii_alphanumeric()).parse(input)?;

    // Must be a recognized assembly name
    if !matches!(
        assembly,
        "GRCh37" | "GRCh38" | "hg19" | "hg38" | "hg18" | "GRCh36"
    ) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    let (input, _) = tag("(").parse(input)?;
    let (input, chromosome) = take_while1(|c: char| c.is_ascii_alphanumeric()).parse(input)?;
    let (input, _) = tag(")").parse(input)?;

    Ok((
        input,
        Accession::from_assembly(assembly.to_string(), chromosome.to_string()),
    ))
}

/// Parse a standard RefSeq-style accession with underscore (e.g., "NM_000088.3")
fn parse_standard_accession(input: &str) -> IResult<&str, Accession> {
    let (input, prefix) = parse_prefix(input)?;
    let (input, _) = tag("_").parse(input)?;
    let (input, number) = parse_number(input)?;
    let (input, version) = opt(preceded(tag("."), parse_version)).parse(input)?;

    Ok((
        input,
        Accession::with_style(prefix.to_string(), number.to_string(), version, false),
    ))
}

/// Parse an Ensembl-style accession without underscore (e.g., "ENST00000012345.1")
fn parse_ensembl_accession(input: &str) -> IResult<&str, Accession> {
    let (input, prefix) = parse_ensembl_prefix(input)?;
    let (input, number) = digit1.parse(input)?;
    let (input, version) = opt(preceded(tag("."), parse_version)).parse(input)?;

    Ok((
        input,
        Accession::with_style(prefix.to_string(), number.to_string(), version, true),
    ))
}

/// Parse Ensembl prefix (ENST, ENSG, ENSP, ENSE, ENSR)
#[inline]
fn parse_ensembl_prefix(input: &str) -> IResult<&str, &str> {
    alt((
        tag("ENST"),
        tag("ENSG"),
        tag("ENSP"),
        tag("ENSE"),
        tag("ENSR"),
    ))
    .parse(input)
}

/// Parse accession prefix (letters like NC, NM, NP, LRG)
#[inline]
fn parse_prefix(input: &str) -> IResult<&str, &str> {
    take_while1(|c: char| c.is_ascii_alphabetic()).parse(input)
}

/// Parse accession number (digits, possibly with embedded letters)
#[inline]
fn parse_number(input: &str) -> IResult<&str, &str> {
    take_while1(|c: char| c.is_ascii_alphanumeric()).parse(input)
}

/// Parse version number
#[inline]
fn parse_version(input: &str) -> IResult<&str, u32> {
    let (remaining, s) = digit1.parse(input)?;
    // Use checked parsing to return error on overflow instead of silent 0
    let version: u32 = s.parse().map_err(|_| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
    })?;
    Ok((remaining, version))
}

/// Parse optional gene symbol in parentheses, e.g., "(COL1A1)"
///
/// Returns None for missing gene symbols or empty parentheses "()"
pub fn parse_gene_symbol(input: &str) -> IResult<&str, Option<String>> {
    opt(map(
        (tag("("), take_while(|c: char| c != ')'), tag(")")),
        |(_, symbol, _): (&str, &str, &str)| {
            // Treat empty strings as None (no gene symbol)
            if symbol.is_empty() {
                None
            } else {
                Some(symbol.to_string())
            }
        },
    ))
    .parse(input)
    .map(|(remaining, opt_opt)| (remaining, opt_opt.flatten()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_nm_accession() {
        let (remaining, acc) = parse_accession("NM_000088.3:c").unwrap();
        assert_eq!(remaining, ":c");
        assert_eq!(&*acc.prefix, "NM");
        assert_eq!(&*acc.number, "000088");
        assert_eq!(acc.version, Some(3));
        assert!(!acc.ensembl_style);
    }

    #[test]
    fn test_parse_nc_accession() {
        let (remaining, acc) = parse_accession("NC_000001.11:g").unwrap();
        assert_eq!(remaining, ":g");
        assert_eq!(&*acc.prefix, "NC");
        assert_eq!(&*acc.number, "000001");
        assert_eq!(acc.version, Some(11));
        assert!(!acc.ensembl_style);
    }

    #[test]
    fn test_parse_accession_no_version() {
        let (remaining, acc) = parse_accession("NM_000088:c").unwrap();
        assert_eq!(remaining, ":c");
        assert_eq!(acc.version, None);
    }

    #[test]
    fn test_parse_ensembl_transcript() {
        // Standard Ensembl format without underscore
        let (remaining, acc) = parse_accession("ENST00000012345.1:c").unwrap();
        assert_eq!(remaining, ":c");
        assert_eq!(&*acc.prefix, "ENST");
        assert_eq!(&*acc.number, "00000012345");
        assert_eq!(acc.version, Some(1));
        assert!(acc.ensembl_style);
        assert!(acc.is_ensembl());
        // Check display format - should be without underscore
        assert_eq!(acc.full(), "ENST00000012345.1");
    }

    #[test]
    fn test_parse_ensembl_gene() {
        let (remaining, acc) = parse_accession("ENSG00000141510.5:g").unwrap();
        assert_eq!(remaining, ":g");
        assert_eq!(&*acc.prefix, "ENSG");
        assert!(acc.ensembl_style);
        assert_eq!(acc.full(), "ENSG00000141510.5");
    }

    #[test]
    fn test_parse_ensembl_protein() {
        let (remaining, acc) = parse_accession("ENSP00000012345.2:p").unwrap();
        assert_eq!(remaining, ":p");
        assert_eq!(&*acc.prefix, "ENSP");
        assert!(acc.ensembl_style);
    }

    #[test]
    fn test_parse_ensembl_no_version() {
        let (remaining, acc) = parse_accession("ENST00000012345:c").unwrap();
        assert_eq!(remaining, ":c");
        assert_eq!(acc.version, None);
        assert!(acc.ensembl_style);
        assert_eq!(acc.full(), "ENST00000012345");
    }

    #[test]
    fn test_ensembl_validation() {
        // Valid: 11 digits
        let acc =
            Accession::with_style("ENST".to_string(), "00000012345".to_string(), Some(1), true);
        assert!(acc.validate_ensembl());

        // Valid: 15 digits
        let acc = Accession::with_style(
            "ENST".to_string(),
            "000000123456789".to_string(),
            Some(1),
            true,
        );
        assert!(acc.validate_ensembl());

        // Invalid: too few digits (10)
        let acc =
            Accession::with_style("ENST".to_string(), "0000001234".to_string(), Some(1), true);
        assert!(!acc.validate_ensembl());

        // Invalid: too many digits (16)
        let acc = Accession::with_style(
            "ENST".to_string(),
            "0000001234567890".to_string(),
            Some(1),
            true,
        );
        assert!(!acc.validate_ensembl());
    }

    #[test]
    fn test_inferred_variant_type() {
        let acc = Accession::new("NC".to_string(), "000001".to_string(), Some(11));
        assert_eq!(acc.inferred_variant_type(), Some("g"));

        let acc = Accession::new("NM".to_string(), "000088".to_string(), Some(3));
        assert_eq!(acc.inferred_variant_type(), Some("c"));

        let acc = Accession::new("NR".to_string(), "000001".to_string(), Some(1));
        assert_eq!(acc.inferred_variant_type(), Some("n"));

        let acc = Accession::new("NP".to_string(), "000079".to_string(), Some(2));
        assert_eq!(acc.inferred_variant_type(), Some("p"));

        let acc =
            Accession::with_style("ENST".to_string(), "00000012345".to_string(), Some(1), true);
        assert_eq!(acc.inferred_variant_type(), Some("c"));

        let acc =
            Accession::with_style("ENSG".to_string(), "00000141510".to_string(), Some(5), true);
        assert_eq!(acc.inferred_variant_type(), Some("g"));

        let acc =
            Accession::with_style("ENSP".to_string(), "00000012345".to_string(), Some(2), true);
        assert_eq!(acc.inferred_variant_type(), Some("p"));
    }

    #[test]
    fn test_parse_gene_symbol() {
        let (remaining, symbol) = parse_gene_symbol("(COL1A1)c.100").unwrap();
        assert_eq!(remaining, "c.100");
        assert_eq!(symbol, Some("COL1A1".to_string()));
    }

    #[test]
    fn test_parse_no_gene_symbol() {
        let (remaining, symbol) = parse_gene_symbol("c.100").unwrap();
        assert_eq!(remaining, "c.100");
        assert_eq!(symbol, None);
    }

    #[test]
    fn test_parse_simple_accession_chr() {
        // Basic chromosome accession
        let (remaining, acc) = parse_accession("chr1:g.123").unwrap();
        assert_eq!(remaining, ":g.123");
        assert_eq!(&*acc.prefix, "chr1");
        assert_eq!(&*acc.number, "");
        assert_eq!(acc.version, None);
        assert!(acc.ensembl_style); // Uses ensembl_style for no-underscore formatting
        assert_eq!(acc.full(), "chr1");
    }

    #[test]
    fn test_parse_simple_accession_contig() {
        // Bare contig name
        let (remaining, acc) = parse_accession("contig:g.123").unwrap();
        assert_eq!(remaining, ":g.123");
        assert_eq!(&*acc.prefix, "contig");
        assert_eq!(acc.full(), "contig");
    }

    #[test]
    fn test_parse_simple_accession_numeric() {
        // Purely numeric chromosome
        let (remaining, acc) = parse_accession("1:g.123").unwrap();
        assert_eq!(remaining, ":g.123");
        assert_eq!(&*acc.prefix, "1");
        assert_eq!(acc.full(), "1");
    }

    #[test]
    fn test_parse_simple_accession_versioned() {
        // Simple accession with version
        let (remaining, acc) = parse_accession("chr1.2:g.123").unwrap();
        assert_eq!(remaining, ":g.123");
        assert_eq!(&*acc.prefix, "chr1");
        assert_eq!(acc.version, Some(2));
        assert_eq!(acc.full(), "chr1.2");
    }

    #[test]
    fn test_parse_simple_accession_hla_style() {
        // HLA-style accession with internal colons
        let (remaining, acc) = parse_accession("HLA-A*01:01:p.Arg100").unwrap();
        assert_eq!(remaining, ":p.Arg100");
        assert_eq!(&*acc.prefix, "HLA-A*01:01");
        assert_eq!(acc.full(), "HLA-A*01:01");
    }

    #[test]
    fn test_parse_simple_accession_all_types() {
        // Test all HGVS type prefixes
        for type_prefix in ["c.", "g.", "n.", "p.", "r.", "m.", "o."] {
            let input = format!("chr1:{type_prefix}123");
            let (remaining, acc) = parse_accession(&input).unwrap();
            assert_eq!(remaining, format!(":{type_prefix}123"));
            assert_eq!(&*acc.prefix, "chr1");
        }
    }

    #[test]
    fn test_parse_simple_accession_special_chars() {
        // SAM-compatible special characters (no underscore to avoid matching RefSeq pattern)
        let (remaining, acc) = parse_accession("scaffold#1:g.123").unwrap();
        assert_eq!(remaining, ":g.123");
        assert_eq!(&*acc.prefix, "scaffold#1");

        // Pipe character (common in NCBI FASTA headers)
        let (remaining, acc) = parse_accession("gi|12345|ref:g.123").unwrap();
        assert_eq!(remaining, ":g.123");
        assert_eq!(&*acc.prefix, "gi|12345|ref");
    }

    #[test]
    fn test_parse_simple_accession_invalid_start() {
        // Cannot start with * or =
        assert!(parse_accession("*chr1:g.123").is_err());
        assert!(parse_accession("=chr1:g.123").is_err());
    }

    #[test]
    fn test_parse_simple_accession_refseq_still_works() {
        // Ensure RefSeq-style accessions are still parsed by parse_standard_accession
        let (remaining, acc) = parse_accession("NC_000001.11:g.123").unwrap();
        assert_eq!(remaining, ":g.123");
        assert_eq!(&*acc.prefix, "NC");
        assert_eq!(&*acc.number, "000001");
        assert_eq!(acc.version, Some(11));
        assert!(!acc.ensembl_style); // RefSeq uses underscore style
    }
}
