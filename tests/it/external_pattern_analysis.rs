//! Detailed analysis of external HGVS patterns
//!
//! This module separates patterns into:
//! - Valid HGVS that we support
//! - Valid HGVS that we don't support (gaps)
//! - Invalid HGVS that we correctly reject
//! - Invalid HGVS that we incorrectly accept (bugs)

use ferro_hgvs::parse_hgvs;

/// Patterns that are semantically INVALID per HGVS specification
/// These SHOULD fail to parse (or be rejected by validation)
const INVALID_HGVS_PATTERNS: &[(&str, &str)] = &[
    // Position 0 is invalid (genomic starts at 1)
    ("NG_007485.1:g.0del", "position 0 invalid"),
    ("NG_012337.1:g.0_1insAAA", "position 0 invalid"),
    // Negative genomic positions invalid
    ("NG_007485.1:g.-1del", "negative genomic position"),
    // Inverted ranges (end < start)
    ("NG_007485.1:g.4_3del", "inverted range"),
    ("NG_007485.1(NM_000077.4):c.135_130insT", "inverted range"),
    ("NG_012337.1(NM_012459.2):c.135_130insT", "inverted range"),
    (
        "NG_012337.1(NM_012459.2):c.-11025_11020inv",
        "inverted range",
    ),
    // Intronic offsets on genomic (genomic has no introns)
    ("NG_012337.1:g.7125+1G>T", "intronic offset on genomic"),
    ("NG_012337.1:g.5+3del", "intronic offset on genomic"),
    ("NG_012337.1:g.4_5+3del", "intronic offset on genomic"),
    ("NG_012337.1:g.5+3_7del", "intronic offset on genomic"),
    (
        "NG_012337.1:g.10_11insLRG_24:g.2+5",
        "intronic offset on genomic in reference",
    ),
    // Invalid syntax
    (
        "NG_012337.1:g.20>60",
        "substitution with number instead of base",
    ),
    (
        "NG_012337.1:g.20_21ins[30_?]",
        "incomplete range in insertion",
    ),
    (
        "NG_012337.1:g.20_21ins[?_40]",
        "incomplete range in insertion",
    ),
    (
        "NG_012337.1:g.20_21ins[30_40;50_?]",
        "incomplete range in insertion",
    ),
    (
        "NG_012337.1:g.20_21ins[30_40;?_60]",
        "incomplete range in insertion",
    ),
    ("NG_012337.1:g.20_21AA[30_40]", "invalid AA syntax"),
    ("NG_012337.1:g.20_21AA[?]", "invalid AA syntax"),
    // Protein positions require amino acid
    (
        "NM_004152.3:p.205_685del",
        "protein position without amino acid",
    ),
];

/// Patterns that are VALID HGVS but use non-standard/extended notation
/// These are valid per HGVS but may not be widely supported
const EXTENDED_VALID_PATTERNS: &[(&str, &str)] = &[
    // Assembly/chromosome notation (valid but non-RefSeq)
    ("GRCh37(chr23):g.3250del", "assembly notation"),
    ("GRCh37(chrX):g.3250del", "assembly notation"),
    ("GRCh37(chrY):g.3250del", "assembly notation"),
    ("GRCh38(chr23):g.3250del", "assembly notation"),
    ("GRCh38(chrX):g.3250del", "assembly notation"),
    ("GRCh38(chrY):g.3250del", "assembly notation"),
    // Complex insertions with references
    (
        "NG_012337.1(NM_012459.2):c.1_2ins[LRG_1:g.100000;AAA]",
        "complex insertion with reference",
    ),
    ("NG_012337.1:g.274dup[1;A]", "complex dup notation"),
    (
        "NG_012337.1:g.20_21ins[(123);30_40]",
        "complex insertion with ranges",
    ),
    (
        "NG_012337.1:g.20_21ins[123;30_40]",
        "complex insertion with ranges",
    ),
    (
        "NG_012337.1:g.20_21ins[30_40;123]",
        "complex insertion with ranges",
    ),
    ("NG_012337.1:g.20_21ins[?]", "uncertain insertion"),
    ("NG_012337.1:g.20_21ins[30_40;?]", "uncertain insertion"),
    // Uncertain dup ranges
    ("NG_012337.1:g.722_742dup(731_741)", "uncertain dup range"),
    ("NG_012337.1:g.722_742dup(731_?)", "uncertain dup range"),
    ("NG_012337.1:g.722_742dup(?_731)", "uncertain dup range"),
    ("NG_012337.1:g.722_742dup(?)", "uncertain dup range"),
    ("NG_012337.1:g.722_742dup?", "uncertain dup"),
    // Delins with reference
    (
        "NG_007485.1(NR_003529.3):n.10delinsNG_007485.1:g.0",
        "delins with reference",
    ),
    (
        "NG_007485.1(NR_003529.3):n.10delinsNG_007485.1:g.100000",
        "delins with reference",
    ),
];

/// Patterns that are legacy/lax notation (should parse, normalize on output)
const LEGACY_NOTATION_PATTERNS: &[(&str, &str)] = &[
    // Multi-base substitution (should be delins)
    (
        "NG_012337.1(NM_003002.2):c.274ATT>T",
        "multi-base substitution",
    ),
    // Inversion with length/sequence specifier
    (
        "NG_012337.1(NM_003002.2):c.274inv3",
        "inversion with length",
    ),
    (
        "NG_012337.1(NM_003002.2):c.274_279inv3",
        "inversion with length",
    ),
    (
        "NG_012337.1(NM_003002.2):c.274_279invAA",
        "inversion with sequence",
    ),
    (
        "NG_012337.1(NM_003002.2):c.274invAA",
        "inversion with sequence",
    ),
    ("NG_012337.1:g.274inv1", "inversion with length"),
    ("NG_012337.1:g.274invT", "inversion with sequence"),
    ("NG_012337.1:g.274_276inv3", "inversion with length"),
    ("NG_012337.1:g.274_276invTAC", "inversion with sequence"),
    // Dup with length/sequence (redundant but valid)
    ("NG_012337.1:g.274dup1", "dup with length"),
    ("NG_012337.1:g.274dupT", "dup with sequence"),
    // Del with explicit length
    ("NG_012337.1(NM_003002.2):c.274del1", "del with length"),
    ("NG_012337.1(NM_003002.2):c.274del3", "del with length"),
    ("NG_012337.1(NM_003002.2):c.274_275del2", "del with length"),
    ("NG_012337.1(NM_003002.2):c.274_279del3", "del with length"),
    // Delins with explicit deleted sequence
    ("NG_012337.1:g.7125delGACinsT", "delins with deleted seq"),
    ("NG_012337.1:g.7125delGinsT", "delins with deleted seq"),
];

/// Copy number notation (Invitae-specific, not standard HGVS)
const NONSTANDARD_PATTERNS: &[(&str, &str)] = &[
    ("NM_000352.3:c.2557_2694copy2", "copy number notation"),
    (
        "NC_000014.8:g.88401076_88459508copy4",
        "copy number notation",
    ),
    ("NM_000352.3:n.2557_2694copy2", "copy number notation"),
];

#[test]
fn test_invalid_patterns_correctly_rejected() {
    println!("\n=== Invalid HGVS Patterns (Should Reject) ===\n");

    let mut correctly_rejected = 0;
    let mut incorrectly_accepted = Vec::new();

    for (pattern, reason) in INVALID_HGVS_PATTERNS {
        let result = parse_hgvs(pattern);
        if result.is_err() {
            correctly_rejected += 1;
            println!("✓ Rejected: {} ({})", pattern, reason);
        } else {
            incorrectly_accepted.push((pattern, reason));
            println!("✗ ACCEPTED (should reject): {} ({})", pattern, reason);
        }
    }

    println!(
        "\nCorrectly rejected: {}/{}",
        correctly_rejected,
        INVALID_HGVS_PATTERNS.len()
    );

    if !incorrectly_accepted.is_empty() {
        println!("\nBUGS - Incorrectly accepted invalid patterns:");
        for (p, r) in &incorrectly_accepted {
            println!("  - {} ({})", p, r);
        }
    }
}

#[test]
fn test_extended_valid_patterns() {
    println!("\n=== Extended/Complex Valid HGVS Patterns ===\n");

    let mut supported = 0;
    let mut unsupported = Vec::new();

    for (pattern, description) in EXTENDED_VALID_PATTERNS {
        let result = parse_hgvs(pattern);
        if result.is_ok() {
            supported += 1;
            println!("✓ Supported: {} ({})", pattern, description);
        } else {
            unsupported.push((pattern, description));
            println!("✗ Unsupported: {} ({})", pattern, description);
        }
    }

    println!(
        "\nSupported: {}/{}",
        supported,
        EXTENDED_VALID_PATTERNS.len()
    );
}

#[test]
fn test_legacy_notation_patterns() {
    println!("\n=== Legacy/Lax Notation Patterns ===\n");

    let mut supported = 0;
    let mut unsupported = Vec::new();

    for (pattern, description) in LEGACY_NOTATION_PATTERNS {
        let result = parse_hgvs(pattern);
        if result.is_ok() {
            supported += 1;
            println!("✓ Supported: {} ({})", pattern, description);
        } else {
            unsupported.push((pattern, description));
            println!("✗ Unsupported: {} ({})", pattern, description);
        }
    }

    println!(
        "\nSupported: {}/{}",
        supported,
        LEGACY_NOTATION_PATTERNS.len()
    );
}

#[test]
fn test_nonstandard_patterns() {
    println!("\n=== Non-standard Patterns (Invitae-specific) ===\n");

    let mut supported = 0;

    for (pattern, description) in NONSTANDARD_PATTERNS {
        let result = parse_hgvs(pattern);
        if result.is_ok() {
            supported += 1;
            println!("✓ Supported: {} ({})", pattern, description);
        } else {
            println!("✗ Unsupported: {} ({})", pattern, description);
        }
    }

    println!("\nSupported: {}/{}", supported, NONSTANDARD_PATTERNS.len());
}

#[test]
fn test_recalculated_support_summary() {
    println!("\n============================================================");
    println!("       RECALCULATED HGVS PATTERN SUPPORT SUMMARY");
    println!("============================================================\n");

    // Count invalid patterns we correctly reject
    let invalid_rejected = INVALID_HGVS_PATTERNS
        .iter()
        .filter(|(p, _)| parse_hgvs(p).is_err())
        .count();
    let invalid_accepted = INVALID_HGVS_PATTERNS.len() - invalid_rejected;

    // Count extended valid patterns we support
    let extended_supported = EXTENDED_VALID_PATTERNS
        .iter()
        .filter(|(p, _)| parse_hgvs(p).is_ok())
        .count();
    let extended_unsupported = EXTENDED_VALID_PATTERNS.len() - extended_supported;

    // Count legacy patterns we support
    let legacy_supported = LEGACY_NOTATION_PATTERNS
        .iter()
        .filter(|(p, _)| parse_hgvs(p).is_ok())
        .count();
    let legacy_unsupported = LEGACY_NOTATION_PATTERNS.len() - legacy_supported;

    // Count non-standard patterns
    let nonstandard_supported = NONSTANDARD_PATTERNS
        .iter()
        .filter(|(p, _)| parse_hgvs(p).is_ok())
        .count();

    println!("Invalid HGVS (error cases):");
    println!(
        "  Correctly rejected: {}/{}",
        invalid_rejected,
        INVALID_HGVS_PATTERNS.len()
    );
    println!("  Incorrectly accepted (bugs): {}", invalid_accepted);

    println!("\nExtended Valid HGVS (complex notation):");
    println!(
        "  Supported: {}/{}",
        extended_supported,
        EXTENDED_VALID_PATTERNS.len()
    );
    println!("  Unsupported (gaps): {}", extended_unsupported);

    println!("\nLegacy Notation (should normalize):");
    println!(
        "  Supported: {}/{}",
        legacy_supported,
        LEGACY_NOTATION_PATTERNS.len()
    );
    println!("  Unsupported (gaps): {}", legacy_unsupported);

    println!("\nNon-standard (Invitae-specific):");
    println!(
        "  Supported: {}/{}",
        nonstandard_supported,
        NONSTANDARD_PATTERNS.len()
    );

    // Calculate adjusted Mutalyzer stats
    let total_valid = EXTENDED_VALID_PATTERNS.len() + LEGACY_NOTATION_PATTERNS.len();
    let total_supported = extended_supported + legacy_supported;
    let mutalyzer_valid_rate = total_supported as f64 / total_valid as f64 * 100.0;

    println!("\n------------------------------------------------------------");
    println!("MUTALYZER ADJUSTED STATS (excluding invalid test cases):");
    println!("  Valid patterns tested: {}", total_valid);
    println!("  Supported: {}", total_supported);
    println!("  Unsupported: {}", total_valid - total_supported);
    println!("  Support rate: {:.1}%", mutalyzer_valid_rate);
    println!("------------------------------------------------------------");
}
