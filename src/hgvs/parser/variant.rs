//! Full variant parsing
//!
//! Parses complete HGVS variant strings into the HgvsVariant type.

use crate::error::FerroError;
use crate::hgvs::edit::ProteinEdit;
use crate::hgvs::interval::{
    CdsInterval, GenomeInterval, Interval, ProtInterval, RnaInterval, TxInterval, UncertainBoundary,
};
use crate::hgvs::location::{AminoAcid, CdsPos, ProtPos};
use crate::hgvs::parser::accession::{parse_accession, parse_gene_symbol};
use crate::hgvs::parser::edit::{parse_na_edit, parse_protein_edit};
use crate::hgvs::parser::position::{
    parse_amino_acid, parse_cds_pos, parse_genome_pos, parse_prot_pos, parse_rna_pos, parse_tx_pos,
};
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{
    Accession, AllelePhase, AlleleVariant, CdsVariant, CircularVariant, GenomeVariant, HgvsVariant,
    LocEdit, MtVariant, ProteinVariant, RnaFusionBreakpoint, RnaFusionVariant, RnaVariant,
    TxVariant,
};
use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::{char, digit1},
    combinator::map,
    sequence::{delimited, preceded},
    IResult, Parser,
};

fn parse_uncertain_tx_pos(input: &str) -> IResult<&str, Mu<crate::hgvs::location::TxPos>> {
    alt((
        // Unknown: ?
        map(tag("?"), |_| Mu::Unknown),
        map(delimited(char('('), parse_tx_pos, char(')')), Mu::uncertain),
        map(parse_tx_pos, Mu::certain),
    ))
    .parse(input)
}

fn parse_uncertain_rna_pos(input: &str) -> IResult<&str, Mu<crate::hgvs::location::RnaPos>> {
    alt((
        // Unknown: ?
        map(tag("?"), |_| Mu::Unknown),
        map(
            delimited(char('('), parse_rna_pos, char(')')),
            Mu::uncertain,
        ),
        map(parse_rna_pos, Mu::certain),
    ))
    .parse(input)
}

/// Parse (pos_pos) with intronic offsets as a range boundary inner
fn parse_cds_intronic_range_inner(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::CdsPos>> {
    let (remaining, (start, _, end)) = (parse_cds_pos, tag("_"), parse_cds_pos).parse(input)?;
    // Only treat as range boundary if at least one has intronic offset
    if start.offset.is_some() || end.offset.is_some() {
        Ok((
            remaining,
            UncertainBoundary::range(Mu::certain(start), Mu::certain(end)),
        ))
    } else {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )))
    }
}

/// Parse unparenthesized ?_pos pattern (e.g., ?_-540 in c.?_-540_3117+?del)
/// Only succeeds if the pattern is followed by another _ (indicating a range continuation)
fn parse_cds_unparenthesized_range_start(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::CdsPos>> {
    let (after_question, _) = tag::<_, _, nom::error::Error<&str>>("?_").parse(input)?;
    let (remaining, pos) = parse_cds_pos(after_question)?;
    // Only use this parser if followed by _ (indicating another position follows)
    if !remaining.starts_with('_') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }
    Ok((
        remaining,
        UncertainBoundary::range(Mu::Unknown, Mu::certain(pos)),
    ))
}

/// Parse a CDS range boundary: (?_pos), (pos_?), (pos_pos), or (pos+off_pos-off)
fn parse_cds_range_boundary(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::CdsPos>> {
    delimited(
        char('('),
        alt((
            // (?_pos)
            map((tag("?"), tag("_"), parse_cds_pos), |(_, _, end)| {
                UncertainBoundary::range(Mu::Unknown, Mu::certain(end))
            }),
            // (pos_?)
            map((parse_cds_pos, tag("_"), tag("?")), |(start, _, _)| {
                UncertainBoundary::range(Mu::certain(start), Mu::Unknown)
            }),
            // (pos_pos) with intronic offsets
            parse_cds_intronic_range_inner,
            // Simple range (pos_pos) without offsets - for uncertain boundaries like (4_6)
            map(
                (parse_cds_pos, tag("_"), parse_cds_pos),
                |(start, _, end)| UncertainBoundary::range(Mu::certain(start), Mu::certain(end)),
            ),
        )),
        char(')'),
    )
    .parse(input)
}

/// Parse a complex CDS boundary that can be a single position or a range
/// Handles: pos, (pos), ?, (?_pos), (pos_?), (pos+off_pos-off)
/// Note: (pos_pos) without offsets/unknowns is NOT parsed here - it's handled at the interval level.
/// Parse a CDS position boundary with optional uncertainty
/// Optimized: fast path for positions starting with digit, *, or -
#[inline]
fn parse_cds_boundary(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::CdsPos>> {
    let bytes = input.as_bytes();
    if bytes.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Eof,
        )));
    }

    // Fast path: if starts with digit, *, or - it's likely a certain position
    if bytes[0].is_ascii_digit() || bytes[0] == b'*' || bytes[0] == b'-' {
        return map(parse_cds_pos, UncertainBoundary::certain).parse(input);
    }

    // Otherwise, try complex patterns
    alt((
        // Range boundary: (?_pos), (pos_?), or (pos+off_pos-off) with intronic offsets
        parse_cds_range_boundary,
        // Single uncertain position: (pos)
        map(
            delimited(char('('), parse_cds_pos, char(')')),
            UncertainBoundary::uncertain,
        ),
        // Unparenthesized ?_pos pattern (e.g., ?_-540 in ?_-540_3117+?del)
        // Only treat ?_pos as a range boundary if it's followed by another _
        // Otherwise ?_pos should be parsed as ? (start) and pos (end) of an interval
        parse_cds_unparenthesized_range_start,
        // Unknown: ?
        map(tag("?"), |_| UncertainBoundary::unknown()),
    ))
    .parse(input)
}

/// Parse a genome range boundary: (?_pos), (pos_?), or (pos_pos)
fn parse_genome_range_boundary(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::GenomePos>> {
    delimited(
        char('('),
        alt((
            // (?_pos)
            map((tag("?"), tag("_"), parse_genome_pos), |(_, _, end)| {
                UncertainBoundary::range(Mu::Unknown, Mu::certain(end))
            }),
            // (pos_?)
            map((parse_genome_pos, tag("_"), tag("?")), |(start, _, _)| {
                UncertainBoundary::range(Mu::certain(start), Mu::Unknown)
            }),
            // (pos_pos) - range with both positions known
            map(
                (parse_genome_pos, tag("_"), parse_genome_pos),
                |(start, _, end)| UncertainBoundary::range(Mu::certain(start), Mu::certain(end)),
            ),
        )),
        char(')'),
    )
    .parse(input)
}

/// Parse a complex genome boundary that can be a single position or a range
/// Only parses (?_pos) and (pos_?) as range boundaries - plain (pos_pos) is handled at interval level
/// Parse a genome position boundary with optional uncertainty
/// Optimized: fast path for simple digit positions (most common case)
#[inline]
fn parse_genome_boundary(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::GenomePos>> {
    let bytes = input.as_bytes();
    if bytes.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Eof,
        )));
    }

    // Fast path: if starts with digit, it's a simple certain position
    if bytes[0].is_ascii_digit() {
        return map(parse_genome_pos, UncertainBoundary::certain).parse(input);
    }

    // Otherwise, try complex patterns
    alt((
        // Range boundary with unknown: (?_pos) or (pos_?)
        parse_genome_range_boundary,
        // Single uncertain position: (pos)
        map(
            delimited(char('('), parse_genome_pos, char(')')),
            UncertainBoundary::uncertain,
        ),
        // Unknown: ?
        map(tag("?"), |_| UncertainBoundary::unknown()),
    ))
    .parse(input)
}

/// Check if input starts with a special position marker (pter, qter, cen)
#[inline]
fn starts_with_special_pos(bytes: &[u8]) -> bool {
    (bytes.len() >= 4 && (bytes[0..4] == *b"pter" || bytes[0..4] == *b"qter"))
        || (bytes.len() >= 3 && bytes[0..3] == *b"cen")
}

/// Parse an interval (single position or range) with optional uncertainty
/// Supports complex boundaries like (pos_pos) ranges
/// Note: Inverted ranges (start > end) are rejected as invalid HGVS
/// Optimized: fast path for simple positions and ranges (starting with digit)
#[inline]
fn parse_genome_interval(input: &str) -> IResult<&str, GenomeInterval> {
    let bytes = input.as_bytes();
    if bytes.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Eof,
        )));
    }

    // Fast path: if starts with digit or special position, it's a simple position or range
    if bytes[0].is_ascii_digit() || starts_with_special_pos(bytes) {
        let (remaining, start) = parse_genome_pos(input)?;
        // Check if this is a range (followed by _)
        if let Some(after_underscore) = remaining.strip_prefix('_') {
            // Use fast path if followed by another digit or special position
            let after_bytes = after_underscore.as_bytes();
            if !after_bytes.is_empty()
                && (after_bytes[0].is_ascii_digit() || starts_with_special_pos(after_bytes))
            {
                let (remaining, end) = parse_genome_pos(after_underscore)?;
                // Skip inverted range check if:
                // - either position is special (pter/qter have no comparable base values)
                // - the edit is a structural variant type that may use inverted coordinates
                //   (insertions, duplications, inversions, and complex delins)
                let allows_inverted = remaining.starts_with("dupins")
                    || remaining.starts_with("inv")
                    || remaining.starts_with("ins")
                    || remaining.starts_with("dup");
                if !start.is_special()
                    && !end.is_special()
                    && !allows_inverted
                    && start.base > end.base
                {
                    return Err(nom::Err::Error(nom::error::Error::new(
                        input,
                        nom::error::ErrorKind::Verify,
                    )));
                }
                return Ok((
                    remaining,
                    Interval::with_complex_boundaries(
                        UncertainBoundary::certain(start),
                        UncertainBoundary::certain(end),
                    ),
                ));
            }
            // Fall through to complex parser
        } else {
            // Single position (no underscore follows)
            return Ok((remaining, GenomeInterval::point(start)));
        }
    }

    // Complex cases: starts with ( or ?
    alt((
        // Simple uncertain interval: (pos_pos) NOT followed by _
        // This produces (pos)_(pos) in output
        |input| {
            let (remaining, (start, _, end)) = delimited(
                char('('),
                (parse_genome_pos, tag("_"), parse_genome_pos),
                char(')'),
            )
            .parse(input)?;
            // If followed by _, this is a complex boundary, not a simple uncertain interval
            if remaining.starts_with('_') {
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Verify,
                )));
            }
            // Reject inverted ranges
            if start.base > end.base {
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Verify,
                )));
            }
            Ok((
                remaining,
                Interval::with_uncertainty(Mu::uncertain(start), Mu::uncertain(end)),
            ))
        },
        // Complex interval with range boundaries: (start_boundary)_(end_boundary)
        // This handles patterns like: (?_100)_(200_?)del, (100_200)_(300_400)del
        |input| {
            let (remaining, (start, _, end)) =
                (parse_genome_boundary, tag("_"), parse_genome_boundary).parse(input)?;
            // For simple ranges with certain positions, reject inverted ranges
            if let (
                UncertainBoundary::Single(Mu::Certain(s)),
                UncertainBoundary::Single(Mu::Certain(e)),
            ) = (&start, &end)
            {
                if s.base > e.base {
                    return Err(nom::Err::Error(nom::error::Error::new(
                        input,
                        nom::error::ErrorKind::Verify,
                    )));
                }
            }
            Ok((remaining, Interval::with_complex_boundaries(start, end)))
        },
        // Single position (possibly uncertain, unknown, or a range boundary)
        map(parse_genome_boundary, |boundary| match boundary {
            UncertainBoundary::Single(Mu::Certain(inner)) => GenomeInterval::point(inner),
            UncertainBoundary::Single(Mu::Uncertain(inner)) => {
                Interval::with_uncertainty(Mu::uncertain(inner), Mu::uncertain(inner))
            }
            UncertainBoundary::Single(Mu::Unknown) => {
                Interval::with_uncertainty(Mu::Unknown, Mu::Unknown)
            }
            UncertainBoundary::Range { start, end } => {
                // A single range boundary used as a point - unusual but valid
                Interval::with_complex_boundaries(
                    UncertainBoundary::Range {
                        start: start.clone(),
                        end: end.clone(),
                    },
                    UncertainBoundary::Range { start, end },
                )
            }
        }),
    ))
    .parse(input)
}

/// Parse a simple uncertain CDS interval (pos_pos) NOT followed by _
/// This prevents eating the first boundary of a complex interval
fn parse_simple_uncertain_cds_interval(input: &str) -> IResult<&str, CdsInterval> {
    let (remaining, _) = char('(').parse(input)?;
    let (remaining, start) = parse_cds_pos(remaining)?;
    let (remaining, _) = tag("_").parse(remaining)?;
    let (remaining, end) = parse_cds_pos(remaining)?;
    let (remaining, _) = char(')').parse(remaining)?;

    // Check that this is NOT followed by _ (which would indicate a complex interval)
    if remaining.starts_with('_') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    // Only if NEITHER has an offset, treat as simple uncertain interval
    if start.offset.is_none() && end.offset.is_none() {
        Ok((
            remaining,
            Interval::with_uncertainty(Mu::uncertain(start), Mu::uncertain(end)),
        ))
    } else {
        // Check if followed by another parenthesized range without underscore
        // This handles patterns like (2727+1_2728-1)(2858+1_2859-1)
        if remaining.starts_with('(') {
            let (remaining2, _) = char('(').parse(remaining)?;
            let (remaining2, end_start) = parse_cds_pos(remaining2)?;
            let (remaining2, _) = tag("_").parse(remaining2)?;
            let (remaining2, end_end) = parse_cds_pos(remaining2)?;
            let (remaining2, _) = char(')').parse(remaining2)?;

            return Ok((
                remaining2,
                Interval::with_complex_boundaries(
                    UncertainBoundary::range(Mu::certain(start), Mu::certain(end)),
                    UncertainBoundary::range(Mu::certain(end_start), Mu::certain(end_end)),
                ),
            ));
        }

        // Has offsets - this is a single range boundary used as an interval
        // (unusual but valid - e.g., just "(100+1_101-1)del" as a point)
        Ok((
            remaining,
            Interval::with_complex_boundaries(
                UncertainBoundary::range(Mu::certain(start), Mu::certain(end)),
                UncertainBoundary::range(Mu::certain(start), Mu::certain(end)),
            ),
        ))
    }
}

/// Compare two CDS positions for ordering
/// Returns true if start > end (inverted range)
///
/// Note: When either position has an intronic offset, we cannot determine
/// the true genomic order without knowing intron lengths. In such cases,
/// we allow the range (return false) since patterns like c.85-47_84+48 are
/// valid intronic ranges that appear in databases like ExAC.
fn cds_pos_is_inverted(start: &CdsPos, end: &CdsPos) -> bool {
    // If either position is unknown (?), we can't determine order. Allow these ranges.
    if start.is_unknown() || end.is_unknown() {
        return false;
    }

    // If either position has an intronic offset, we can't reliably determine
    // the genomic order without knowing the intron length. Allow these ranges.
    if start.offset.is_some() || end.offset.is_some() {
        return false;
    }

    // UTR ordering: 5' UTR < CDS < 3' UTR
    // 5' UTR has utr3=false and negative base
    // CDS has utr3=false and positive base
    // 3' UTR has utr3=true

    // Different regions
    if start.utr3 != end.utr3 {
        // 3' UTR (utr3=true) is always after 5' UTR or CDS
        // So start=3'UTR, end=5'UTR/CDS is inverted
        return start.utr3;
    }

    // Same UTR region - compare bases
    if start.utr3 {
        // Both in 3' UTR: *1 < *100
        start.base > end.base
    } else {
        // 5' UTR or CDS: -100 < -1 < 1 < 100
        start.base > end.base
    }
}

/// Parse a CDS interval with support for complex boundaries
/// Note: Inverted ranges (start > end) are rejected as invalid HGVS
/// Parse a CDS interval
/// Optimized: fast path for simple positions (starting with digit, *, or -)
#[inline]
fn parse_cds_interval(input: &str) -> IResult<&str, CdsInterval> {
    let bytes = input.as_bytes();
    if bytes.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Eof,
        )));
    }

    // Fast path: if starts with digit, *, or - it's a simple position or range
    // Also handle ? when followed by +/- offset (e.g., ?-232) but NOT when followed by _
    // (e.g., ?_-540 needs the complex parser for patterns like ?_-540_3117+?del)
    let use_fast_path = bytes[0].is_ascii_digit()
        || bytes[0] == b'*'
        || bytes[0] == b'-'
        || (bytes[0] == b'?' && bytes.len() > 1 && (bytes[1] == b'+' || bytes[1] == b'-'));
    if use_fast_path {
        let (remaining, start) = parse_cds_pos(input)?;
        // Check if this is a range (followed by _)
        if let Some(after_underscore) = remaining.strip_prefix('_') {
            // Only use fast path if followed by a position start character (simple range)
            // Fall through to complex parser for patterns like 100_(200)del
            let after_bytes = after_underscore.as_bytes();
            // Note: ? is NOT included - it needs the complex parser
            let is_simple_range = !after_bytes.is_empty()
                && (after_bytes[0].is_ascii_digit()
                    || after_bytes[0] == b'*'
                    || after_bytes[0] == b'-');
            if is_simple_range {
                let (remaining, end) = parse_cds_pos(after_underscore)?;
                // Reject inverted ranges
                if cds_pos_is_inverted(&start, &end) {
                    return Err(nom::Err::Error(nom::error::Error::new(
                        input,
                        nom::error::ErrorKind::Verify,
                    )));
                }
                return Ok((
                    remaining,
                    Interval::with_complex_boundaries(
                        UncertainBoundary::certain(start),
                        UncertainBoundary::certain(end),
                    ),
                ));
            }
            // Fall through to complex parser
        } else {
            // Single position (no underscore follows)
            return Ok((remaining, CdsInterval::point(start)));
        }
    }

    // Complex cases: starts with ( or ?
    alt((
        // Simple uncertain interval: (pos_pos) NOT followed by _
        // This produces (pos)_(pos) in output - both positions are uncertain
        parse_simple_uncertain_cds_interval,
        // Complex interval with range boundaries: (start_boundary)_(end_boundary)
        // This handles patterns like: (4185+1_4186-1)_(4357+1_4358-1)del or simple ranges
        |input| {
            let (remaining, (start, _, end)) =
                (parse_cds_boundary, tag("_"), parse_cds_boundary).parse(input)?;
            // For simple ranges with certain positions, reject inverted ranges
            if let (
                UncertainBoundary::Single(Mu::Certain(s)),
                UncertainBoundary::Single(Mu::Certain(e)),
            ) = (&start, &end)
            {
                if cds_pos_is_inverted(s, e) {
                    return Err(nom::Err::Error(nom::error::Error::new(
                        input,
                        nom::error::ErrorKind::Verify,
                    )));
                }
            }
            Ok((remaining, Interval::with_complex_boundaries(start, end)))
        },
        // Single position (possibly uncertain, unknown, or a range boundary)
        map(parse_cds_boundary, |boundary| match boundary {
            UncertainBoundary::Single(Mu::Certain(inner)) => CdsInterval::point(inner),
            UncertainBoundary::Single(Mu::Uncertain(inner)) => {
                Interval::with_uncertainty(Mu::uncertain(inner), Mu::uncertain(inner))
            }
            UncertainBoundary::Single(Mu::Unknown) => {
                Interval::with_uncertainty(Mu::Unknown, Mu::Unknown)
            }
            UncertainBoundary::Range { start, end } => {
                // A single range boundary used as a point - unusual but valid
                Interval::with_complex_boundaries(
                    UncertainBoundary::Range {
                        start: start.clone(),
                        end: end.clone(),
                    },
                    UncertainBoundary::Range { start, end },
                )
            }
        }),
    ))
    .parse(input)
}

/// Parse a transcript range boundary: (pos_pos)
fn parse_tx_range_boundary(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::TxPos>> {
    delimited(
        char('('),
        map((parse_tx_pos, tag("_"), parse_tx_pos), |(start, _, end)| {
            UncertainBoundary::range(Mu::certain(start), Mu::certain(end))
        }),
        char(')'),
    )
    .parse(input)
}

/// Parse a transcript boundary: (pos_pos), (pos), or pos
fn parse_tx_boundary(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::TxPos>> {
    alt((
        // Range boundary: (pos_pos)
        parse_tx_range_boundary,
        // Single uncertain position: (pos)
        map(
            delimited(char('('), parse_tx_pos, char(')')),
            UncertainBoundary::uncertain,
        ),
        // Certain position: pos
        map(parse_tx_pos, UncertainBoundary::certain),
    ))
    .parse(input)
}

fn parse_tx_interval(input: &str) -> IResult<&str, TxInterval> {
    alt((
        // Whole interval uncertain: (100_200) NOT followed by _
        |input| {
            let (remaining, (start, _, end)) =
                delimited(char('('), (parse_tx_pos, tag("_"), parse_tx_pos), char(')'))
                    .parse(input)?;
            // If followed by _, this is a complex boundary, not a simple uncertain interval
            if remaining.starts_with('_') {
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Verify,
                )));
            }
            Ok((
                remaining,
                Interval::with_uncertainty(Mu::uncertain(start), Mu::uncertain(end)),
            ))
        },
        // Complex interval with range boundaries: (pos_pos)_(pos_pos)
        |input| {
            let (remaining, (start, _, end)) =
                (parse_tx_boundary, tag("_"), parse_tx_boundary).parse(input)?;
            Ok((remaining, Interval::with_complex_boundaries(start, end)))
        },
        // Range with individual uncertainty
        map(
            (parse_uncertain_tx_pos, tag("_"), parse_uncertain_tx_pos),
            |(start, _, end)| Interval::with_uncertainty(start, end),
        ),
        // Single position (possibly uncertain or unknown)
        map(parse_uncertain_tx_pos, |pos| match pos {
            Mu::Unknown => Interval::with_uncertainty(Mu::Unknown, Mu::Unknown),
            Mu::Uncertain(inner) => {
                Interval::with_uncertainty(Mu::uncertain(inner), Mu::uncertain(inner))
            }
            Mu::Certain(inner) => TxInterval::point(inner),
        }),
    ))
    .parse(input)
}

/// Parse a protein range boundary: (?_pos), (pos_?), or (pos_pos)
fn parse_prot_range_boundary(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::ProtPos>> {
    delimited(
        char('('),
        alt((
            // (?_pos) - unknown to known position
            map((tag("?"), tag("_"), parse_prot_pos), |(_, _, end)| {
                UncertainBoundary::range(Mu::Unknown, Mu::certain(end))
            }),
            // (pos_?) - known position to unknown
            map((parse_prot_pos, tag("_"), tag("?")), |(start, _, _)| {
                UncertainBoundary::range(Mu::certain(start), Mu::Unknown)
            }),
            // (pos_pos) - range with both positions known
            map(
                (parse_prot_pos, tag("_"), parse_prot_pos),
                |(start, _, end)| UncertainBoundary::range(Mu::certain(start), Mu::certain(end)),
            ),
        )),
        char(')'),
    )
    .parse(input)
}

/// Parse amino acid followed by position range: Met(?_1), Glu(1_?), or Glu2(?)
/// This handles patterns like p.Met(?_1)_Glu2(?) where the amino acid has an uncertain position range
fn parse_prot_aa_with_range(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::ProtPos>> {
    use crate::hgvs::location::ProtPos;
    use nom::character::complete::digit1;

    let (remaining, aa) = parse_amino_acid(input)?;

    // Check for AminoAcid + Number + (?) pattern: e.g., Glu2(?)
    if let Ok((number_remaining, num_str)) =
        digit1::<&str, nom::error::Error<&str>>.parse(remaining)
    {
        if number_remaining.starts_with("(?)") {
            let pos = num_str.parse::<u64>().unwrap_or(0);
            let (final_remaining, _) =
                tag::<&str, &str, nom::error::Error<&str>>("(?)").parse(number_remaining)?;
            return Ok((
                final_remaining,
                UncertainBoundary::uncertain(ProtPos::new(aa, pos)),
            ));
        }
    }

    // Must be followed by (
    if !remaining.starts_with('(') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    let (remaining, _) = char('(').parse(remaining)?;

    // Parse the range inside parens
    let (remaining, boundary) = alt((
        // (?_NUMBER) - unknown to known number (e.g., (?_1))
        map(
            (tag("?"), tag("_"), digit1),
            |(_, _, end): (&str, &str, &str)| {
                let pos = end.parse::<u64>().unwrap_or(0);
                UncertainBoundary::range(Mu::Unknown, Mu::certain(ProtPos::new(aa, pos)))
            },
        ),
        // (NUMBER_?) - known number to unknown (e.g., (1_?))
        map(
            (digit1, tag("_"), tag("?")),
            |(start, _, _): (&str, &str, &str)| {
                let pos = start.parse::<u64>().unwrap_or(0);
                UncertainBoundary::range(Mu::certain(ProtPos::new(aa, pos)), Mu::Unknown)
            },
        ),
        // (?) - fully uncertain position (no number before)
        map(tag("?"), |_| {
            UncertainBoundary::uncertain(ProtPos::new(aa, 0))
        }),
    ))
    .parse(remaining)?;

    let (remaining, _) = char(')').parse(remaining)?;

    Ok((remaining, boundary))
}

/// Parse a protein boundary: (?_pos), (pos_?), (pos_pos), (pos), ?, pos, or AminoAcid(?_N)
fn parse_prot_boundary(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::ProtPos>> {
    alt((
        // Amino acid with position range: Met(?_1), Glu(1_?), Ter(?)
        // Must be checked BEFORE certain position to avoid consuming just "Met"
        parse_prot_aa_with_range,
        // Range boundary: (?_pos), (pos_?), or (pos_pos)
        parse_prot_range_boundary,
        // Single uncertain position: (pos)
        map(
            delimited(char('('), parse_prot_pos, char(')')),
            UncertainBoundary::uncertain,
        ),
        // Unknown amino acid with position: ?NUMBER (e.g., ?4 means unknown amino acid at position 4)
        map(preceded(char('?'), digit1), |n: &str| {
            UncertainBoundary::certain(ProtPos::new(AminoAcid::Xaa, n.parse().unwrap_or(0)))
        }),
        // Unknown: ?
        map(tag("?"), |_| UncertainBoundary::unknown()),
        // Certain position: pos
        map(parse_prot_pos, UncertainBoundary::certain),
    ))
    .parse(input)
}

fn parse_prot_interval(input: &str) -> IResult<&str, ProtInterval> {
    alt((
        // Simple uncertain interval: (pos_pos) NOT followed by _
        // This is for patterns like (Lys23_Leu24)del where the whole interval is uncertain
        |input| {
            let (remaining, (start, _, end)) = delimited(
                char('('),
                (parse_prot_pos, tag("_"), parse_prot_pos),
                char(')'),
            )
            .parse(input)?;
            // If followed by _, this is a complex boundary, not a simple uncertain interval
            if remaining.starts_with('_') {
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Verify,
                )));
            }
            Ok((
                remaining,
                Interval::with_uncertainty(Mu::uncertain(start), Mu::uncertain(end)),
            ))
        },
        // Complex interval with range boundaries: (?_pos)_pos, pos_(pos_?), etc.
        |input| {
            let (remaining, (start, _, end)) =
                (parse_prot_boundary, tag("_"), parse_prot_boundary).parse(input)?;
            Ok((remaining, Interval::with_complex_boundaries(start, end)))
        },
        // Single position (possibly uncertain, unknown, or a range boundary)
        map(parse_prot_boundary, |boundary| match boundary {
            UncertainBoundary::Single(Mu::Certain(inner)) => ProtInterval::point(inner),
            UncertainBoundary::Single(Mu::Uncertain(inner)) => {
                Interval::with_uncertainty(Mu::uncertain(inner), Mu::uncertain(inner))
            }
            UncertainBoundary::Single(Mu::Unknown) => {
                Interval::with_uncertainty(Mu::Unknown, Mu::Unknown)
            }
            UncertainBoundary::Range { start, end } => {
                // A single range boundary used as a point - unusual but valid
                Interval::with_complex_boundaries(
                    UncertainBoundary::Range {
                        start: start.clone(),
                        end: end.clone(),
                    },
                    UncertainBoundary::Range { start, end },
                )
            }
        }),
    ))
    .parse(input)
}

fn parse_rna_interval(input: &str) -> IResult<&str, RnaInterval> {
    alt((
        // Whole interval uncertain: (100_200)
        map(
            delimited(
                char('('),
                (parse_rna_pos, tag("_"), parse_rna_pos),
                char(')'),
            ),
            |(start, _, end)| Interval::with_uncertainty(Mu::uncertain(start), Mu::uncertain(end)),
        ),
        // Range with individual uncertainty
        map(
            (parse_uncertain_rna_pos, tag("_"), parse_uncertain_rna_pos),
            |(start, _, end)| Interval::with_uncertainty(start, end),
        ),
        // Single position (possibly uncertain or unknown)
        map(parse_uncertain_rna_pos, |pos| match pos {
            Mu::Unknown => Interval::with_uncertainty(Mu::Unknown, Mu::Unknown),
            Mu::Uncertain(inner) => {
                Interval::with_uncertainty(Mu::uncertain(inner), Mu::uncertain(inner))
            }
            Mu::Certain(inner) => RnaInterval::point(inner),
        }),
    ))
    .parse(input)
}

/// Parse genomic variant (g.)
fn parse_genome_variant(
    accession: Accession,
    gene_symbol: Option<String>,
) -> impl FnMut(&str) -> IResult<&str, HgvsVariant> {
    move |input: &str| {
        let (input, _) = tag("g.").parse(input)?;

        // First try whole-genome identity (g.=) - no position
        if let Ok((remaining, edit)) = parse_whole_genome_identity(input) {
            let dummy_pos = crate::hgvs::location::GenomePos::new(1);
            let dummy_interval = GenomeInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Genome(GenomeVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(dummy_interval, edit),
                }),
            ));
        }

        // Try whole-genome unknown (g.?) - no position
        if let Ok((remaining, edit)) = parse_whole_genome_unknown(input) {
            let dummy_pos = crate::hgvs::location::GenomePos::new(1);
            let dummy_interval = GenomeInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Genome(GenomeVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(dummy_interval, edit),
                }),
            ));
        }

        // Try compound allele notation: g.[pos1edit1;pos2edit2;...]
        if input.starts_with('[') {
            if let Ok((remaining, allele)) =
                parse_genome_compound_allele(input, accession.clone(), gene_symbol.clone())
            {
                return Ok((remaining, HgvsVariant::Allele(allele)));
            }
        }

        let (input, interval) = parse_genome_interval(input)?;

        // Try to parse an edit; if no edit found, use position-only identity
        let (input, edit) = if let Ok((remaining, edit)) = parse_na_edit(input) {
            (remaining, edit)
        } else if input.is_empty() || input.starts_with(']') || input.starts_with(';') {
            // Position-only pattern (e.g., g.12345 or g.100_200)
            // Use identity with no sequence to represent "position reference"
            (
                input,
                crate::hgvs::edit::NaEdit::Identity {
                    sequence: None,
                    whole_entity: false,
                },
            )
        } else {
            // Unknown pattern - fail
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )));
        };

        Ok((
            input,
            HgvsVariant::Genome(GenomeVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, edit),
            }),
        ))
    }
}

/// Parse a compound allele with genomic variants: [pos1edit1;pos2edit2;...]
fn parse_genome_compound_allele(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, AlleleVariant> {
    let (remaining, _) = char('[').parse(input)?;

    let mut variants = Vec::new();
    let mut current = remaining;

    loop {
        // Parse interval and edit
        let (rest, interval) = parse_genome_interval(current)?;
        // Edit is optional - if not present, treat as identity (no change)
        let (rest, edit) = if let Ok((r, e)) = parse_na_edit(rest) {
            (r, e)
        } else {
            // No edit found - use identity (reference)
            (
                rest,
                crate::hgvs::edit::NaEdit::Identity {
                    sequence: None,
                    whole_entity: false,
                },
            )
        };

        variants.push(HgvsVariant::Genome(GenomeVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit: LocEdit::new(interval, edit),
        }));

        // Check for more variants (semicolon separator) or end of allele
        if let Some(after_semi) = rest.strip_prefix(';') {
            current = after_semi;
        } else if let Some(after_bracket) = rest.strip_prefix(']') {
            current = after_bracket;
            break;
        } else {
            return Err(nom::Err::Error(nom::error::Error::new(
                rest,
                nom::error::ErrorKind::Tag,
            )));
        }
    }

    if variants.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    Ok((current, AlleleVariant::cis(variants)))
}

/// Parse whole-genome identity: g.=
fn parse_whole_genome_identity(input: &str) -> IResult<&str, crate::hgvs::edit::NaEdit> {
    map(tag("="), |_| {
        crate::hgvs::edit::NaEdit::whole_entity_identity()
    })
    .parse(input)
}

/// Parse whole-genome unknown: g.?
fn parse_whole_genome_unknown(input: &str) -> IResult<&str, crate::hgvs::edit::NaEdit> {
    // Only match terminal ? (not ? followed by _ which indicates uncertain position)
    let (remaining, _) = tag("?").parse(input)?;
    if remaining.starts_with('_') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }
    Ok((remaining, crate::hgvs::edit::NaEdit::whole_entity_unknown()))
}

/// Parse whole-CDS identity: c.=
fn parse_whole_cds_identity(input: &str) -> IResult<&str, crate::hgvs::edit::NaEdit> {
    map(tag("="), |_| {
        crate::hgvs::edit::NaEdit::whole_entity_identity()
    })
    .parse(input)
}

/// Parse whole-CDS unknown: c.?
/// Only matches if ? is at end (not followed by _, +, -, or edit keywords like dup/del/ins)
fn parse_whole_cds_unknown(input: &str) -> IResult<&str, crate::hgvs::edit::NaEdit> {
    let (remaining, _) = tag("?").parse(input)?;
    // If followed by '_', this is an interval with unknown start, not whole entity unknown
    // If followed by '+' or '-', this has an offset (e.g., ?-232)
    // If followed by edit keywords, it's an unknown position with edit (e.g., ?dup)
    if remaining.starts_with('_')
        || remaining.starts_with('+')
        || remaining.starts_with('-')
        || remaining.starts_with("dup")
        || remaining.starts_with("del")
        || remaining.starts_with("ins")
        || remaining.starts_with("inv")
    {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }
    Ok((remaining, crate::hgvs::edit::NaEdit::whole_entity_unknown()))
}

/// Parse CDS variant (c.)
fn parse_cds_variant(
    accession: Accession,
    gene_symbol: Option<String>,
) -> impl FnMut(&str) -> IResult<&str, HgvsVariant> {
    move |input: &str| {
        let (input, _) = tag("c.").parse(input)?;

        // Check for cis allele shorthand: c.[145C>T;147C>G]
        if input.starts_with('[') {
            return parse_cds_allele_shorthand(input, accession.clone(), gene_symbol.clone());
        }

        // Check for position-based unknown phase shorthand: c.2376G>C(;)3103del
        // The (;) separator indicates unknown phase between variants
        if input.contains("(;)") && !input.starts_with('[') {
            return parse_cds_position_unknown_phase(input, accession.clone(), gene_symbol.clone());
        }

        // First try whole-CDS identity (c.=) - no position
        if let Ok((remaining, edit)) = parse_whole_cds_identity(input) {
            // Use a dummy interval for whole-entity identity
            let dummy_pos = crate::hgvs::location::CdsPos {
                base: 1,
                offset: None,
                utr3: false,
            };
            let dummy_interval = CdsInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Cds(CdsVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(dummy_interval, edit),
                }),
            ));
        }

        // Try whole-CDS unknown (c.?) - no position
        if let Ok((remaining, edit)) = parse_whole_cds_unknown(input) {
            // Use a dummy interval for whole-entity unknown
            let dummy_pos = crate::hgvs::location::CdsPos {
                base: 1,
                offset: None,
                utr3: false,
            };
            let dummy_interval = CdsInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Cds(CdsVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(dummy_interval, edit),
                }),
            ));
        }

        // Try predicted variant: c.(positionEdit) where the whole position+edit is in parens
        // e.g., c.(9740C>A) or c.(742G>A)
        if input.starts_with('(') && !input.starts_with("(?") {
            // Check if this is a predicted variant pattern (simple position+edit in parens)
            // vs an uncertain interval like (100_200)
            if let Some(close_paren) = input.find(')') {
                let inner = &input[1..close_paren];
                // If inner contains a simple substitution pattern (no underscore for range)
                // and doesn't start with a digit followed by underscore, it's likely a predicted variant
                if inner.contains('>') && !inner.contains('_') {
                    // Parse the inner content as position + edit
                    if let Ok((after_edit, interval)) = parse_cds_interval(inner) {
                        if let Ok((remaining_inner, edit)) = parse_na_edit(after_edit) {
                            if remaining_inner.is_empty() {
                                // Mark as predicted by wrapping in uncertainty
                                let remaining = &input[close_paren + 1..];
                                return Ok((
                                    remaining,
                                    HgvsVariant::Cds(CdsVariant {
                                        accession: accession.clone(),
                                        gene_symbol: gene_symbol.clone(),
                                        loc_edit: LocEdit::new_predicted(interval, edit),
                                    }),
                                ));
                            }
                        }
                    }
                }
            }
        }

        let (input, interval) = parse_cds_interval(input)?;

        // Try to parse an edit - if none found, check if position-only is allowed
        let (input, edit) = if let Ok((remaining, edit)) = parse_na_edit(input) {
            (remaining, edit)
        } else if input.is_empty() || !input.chars().next().unwrap_or(' ').is_alphanumeric() {
            // No edit found - position-only is only allowed for:
            // 1. Ranges (different start and end positions)
            // 2. UTR positions (3' or 5' UTR)
            let is_range = match (&interval.start, &interval.end) {
                (
                    crate::hgvs::interval::UncertainBoundary::Single(
                        crate::hgvs::uncertainty::Mu::Certain(start),
                    ),
                    crate::hgvs::interval::UncertainBoundary::Single(
                        crate::hgvs::uncertainty::Mu::Certain(end),
                    ),
                ) => start != end,
                _ => true, // Complex boundaries count as ranges
            };

            let is_utr = match &interval.start {
                crate::hgvs::interval::UncertainBoundary::Single(
                    crate::hgvs::uncertainty::Mu::Certain(pos),
                ) => pos.utr3 || pos.base < 0, // 3' UTR or 5' UTR
                _ => false,
            };

            if is_range || is_utr {
                (input, crate::hgvs::edit::NaEdit::PositionOnly)
            } else {
                // Simple single position without edit - reject it
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Tag,
                )));
            }
        } else {
            // There's more input but it doesn't parse as an edit - this is an error
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )));
        };

        Ok((
            input,
            HgvsVariant::Cds(CdsVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, edit),
            }),
        ))
    }
}

/// Parse position-based unknown phase shorthand: 2376G>C(;)3103del
fn parse_cds_position_unknown_phase(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    // Pre-allocate with estimated capacity (most alleles have 2-4 variants)
    let mut variants = Vec::with_capacity(4);

    // Split by (;) and parse each part
    for part in input.split("(;)") {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }

        // Parse interval + edit
        let (edit_remaining, interval) = parse_cds_interval(part)?;
        let (final_remaining, edit) = parse_na_edit(edit_remaining)?;

        if !final_remaining.trim().is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                final_remaining,
                nom::error::ErrorKind::Tag,
            )));
        }

        variants.push(HgvsVariant::Cds(CdsVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit: LocEdit::new(interval, edit),
        }));
    }

    if variants.len() < 2 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    Ok((
        "",
        HgvsVariant::Allele(AlleleVariant::new(variants, AllelePhase::Unknown)),
    ))
}

/// Parse CDS allele shorthand: [145C>T;147C>G], [145C>T(;)147C>G], or [76A>C];[0]
/// Also handles mixed phase: [123A>G;456C>T(;)789G>A] where some variants are cis
/// and others have unknown phase relationship.
fn parse_cds_allele_shorthand(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    // Must start with [
    if !input.starts_with('[') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Check for trans allele pattern: [var];[var] or [var];[0] or [var];[?]
    if input.contains("];[") {
        return parse_cds_trans_allele_shorthand(input, accession, gene_symbol);
    }

    let close_bracket = input.rfind(']').ok_or_else(|| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Tag))
    })?;

    let content = &input[1..close_bracket];
    let remaining = &input[close_bracket + 1..];

    // Check for mixed phase: contains both ; and (;)
    // We need to check for ; that's NOT part of (;)
    let has_unknown_phase = content.contains("(;)");
    let has_cis_separator = {
        // Replace (;) temporarily to check for standalone ;
        let temp = content.replace("(;)", "\x00");
        temp.contains(';')
    };

    // Determine phase and parsing strategy
    let phase = if has_unknown_phase {
        // If there's any unknown phase, the overall phase is unknown
        AllelePhase::Unknown
    } else {
        AllelePhase::Cis
    };

    // Parse variants - handle mixed separators
    // Pre-allocate with estimated capacity (most alleles have 2-4 variants)
    let mut variants = Vec::with_capacity(4);

    if has_unknown_phase && has_cis_separator {
        // Mixed phase: split on (;) first, then on ;
        for unknown_group in content.split("(;)") {
            for part in unknown_group.split(';') {
                let part = part.trim();
                if part.is_empty() {
                    continue;
                }

                let (edit_remaining, interval) = parse_cds_interval(part)?;
                // Edit is optional - if not present, treat as identity (no change)
                let (final_remaining, edit) = if let Ok((r, e)) = parse_na_edit(edit_remaining) {
                    (r, e)
                } else {
                    (
                        edit_remaining,
                        crate::hgvs::edit::NaEdit::Identity {
                            sequence: None,
                            whole_entity: false,
                        },
                    )
                };

                if !final_remaining.trim().is_empty() {
                    return Err(nom::Err::Error(nom::error::Error::new(
                        final_remaining,
                        nom::error::ErrorKind::Tag,
                    )));
                }

                variants.push(HgvsVariant::Cds(CdsVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(interval, edit),
                }));
            }
        }
    } else {
        // Single separator type
        let separator = if has_unknown_phase { "(;)" } else { ";" };

        for part in content.split(separator) {
            let part = part.trim();
            if part.is_empty() {
                continue;
            }

            // Parse interval + edit (e.g., "145C>T" or "100_200del")
            let (edit_remaining, interval) = parse_cds_interval(part)?;
            // Edit is optional - if not present, treat as identity (no change)
            let (final_remaining, edit) = if let Ok((r, e)) = parse_na_edit(edit_remaining) {
                (r, e)
            } else {
                (
                    edit_remaining,
                    crate::hgvs::edit::NaEdit::Identity {
                        sequence: None,
                        whole_entity: false,
                    },
                )
            };

            if !final_remaining.trim().is_empty() {
                return Err(nom::Err::Error(nom::error::Error::new(
                    final_remaining,
                    nom::error::ErrorKind::Tag,
                )));
            }

            variants.push(HgvsVariant::Cds(CdsVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, edit),
            }));
        }
    }

    if variants.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    Ok((
        remaining,
        HgvsVariant::Allele(AlleleVariant::new(variants, phase)),
    ))
}

/// Parse CDS trans allele shorthand: [76A>C];[0], [76A>C];[?], [var1];[var2]
fn parse_cds_trans_allele_shorthand(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    // Pre-allocate with capacity for trans alleles (typically 2)
    let mut variants = Vec::with_capacity(2);
    let mut remaining = input;

    while !remaining.is_empty() {
        // Must start with [
        if !remaining.starts_with('[') {
            break;
        }

        // Find closing bracket
        let close_bracket = remaining.find(']').ok_or_else(|| {
            nom::Err::Error(nom::error::Error::new(
                remaining,
                nom::error::ErrorKind::Tag,
            ))
        })?;

        let content = &remaining[1..close_bracket].trim();

        // Check for special markers
        let variant = if *content == "0" {
            HgvsVariant::NullAllele
        } else if *content == "?" {
            HgvsVariant::UnknownAllele
        } else {
            // Parse as CDS variant (interval + edit)
            let (edit_remaining, interval) = parse_cds_interval(content)?;
            let (final_remaining, edit) = parse_na_edit(edit_remaining)?;

            if !final_remaining.trim().is_empty() {
                return Err(nom::Err::Error(nom::error::Error::new(
                    final_remaining,
                    nom::error::ErrorKind::Tag,
                )));
            }

            HgvsVariant::Cds(CdsVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, edit),
            })
        };

        variants.push(variant);
        remaining = &remaining[close_bracket + 1..];

        // Check for more variants (semicolon separator)
        if remaining.starts_with(';') {
            remaining = &remaining[1..];
        }
    }

    if variants.len() < 2 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    Ok((
        remaining,
        HgvsVariant::Allele(AlleleVariant::new(variants, AllelePhase::Trans)),
    ))
}

/// Parse transcript variant (n.)
fn parse_tx_variant(
    accession: Accession,
    gene_symbol: Option<String>,
) -> impl FnMut(&str) -> IResult<&str, HgvsVariant> {
    move |input: &str| {
        let (input, _) = tag("n.").parse(input)?;
        let (input, interval) = parse_tx_interval(input)?;
        let (input, edit) = parse_na_edit(input)?;

        Ok((
            input,
            HgvsVariant::Tx(TxVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, edit),
            }),
        ))
    }
}

/// Parse mitochondrial variant (m.)
fn parse_mt_variant(
    accession: Accession,
    gene_symbol: Option<String>,
) -> impl FnMut(&str) -> IResult<&str, HgvsVariant> {
    move |input: &str| {
        let (input, _) = tag("m.").parse(input)?;
        let (input, interval) = parse_genome_interval(input)?;
        let (input, edit) = parse_na_edit(input)?;

        Ok((
            input,
            HgvsVariant::Mt(MtVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, edit),
            }),
        ))
    }
}

/// Parse whole-protein identity: p.= or p.(=)
fn parse_whole_protein_identity(input: &str) -> IResult<&str, ProteinEdit> {
    alt((
        // Predicted identity: (=)
        map(tag("(=)"), |_| {
            ProteinEdit::whole_protein_identity_predicted()
        }),
        // Certain identity: =
        map(tag("="), |_| ProteinEdit::whole_protein_identity()),
    ))
    .parse(input)
}

/// Parse no protein production: p.0 or p.0?
fn parse_whole_protein_no_protein(input: &str) -> IResult<&str, ProteinEdit> {
    alt((
        // Predicted no-protein: 0?
        map(tag("0?"), |_| ProteinEdit::no_protein_predicted()),
        // Certain no-protein: 0
        map(tag("0"), |_| ProteinEdit::no_protein()),
    ))
    .parse(input)
}

/// Parse whole-protein unknown: p.? or p.(?)
fn parse_whole_protein_unknown(input: &str) -> IResult<&str, ProteinEdit> {
    alt((
        // Predicted unknown: (?)
        map(tag("(?)"), |_| {
            ProteinEdit::whole_protein_unknown_predicted()
        }),
        // Certain unknown: ?
        map(tag("?"), |_| ProteinEdit::whole_protein_unknown()),
    ))
    .parse(input)
}

/// Parse a protein allele shorthand: p.[edit1;edit2;...]
fn parse_protein_allele_shorthand(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, AlleleVariant> {
    let (remaining, _) = char('[').parse(input)?;

    let mut variants = Vec::new();
    let mut current = remaining;

    loop {
        // Parse interval and edit
        let (rest, interval) = parse_prot_interval(current)?;
        let (rest, edit) = parse_protein_edit(rest)?;

        variants.push(HgvsVariant::Protein(ProteinVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit: LocEdit::new(interval, edit),
        }));

        // Check for more variants (semicolon separator) or end of allele
        if let Some(after_semi) = rest.strip_prefix(';') {
            current = after_semi;
        } else if let Some(after_bracket) = rest.strip_prefix(']') {
            current = after_bracket;
            break;
        } else {
            return Err(nom::Err::Error(nom::error::Error::new(
                rest,
                nom::error::ErrorKind::Tag,
            )));
        }
    }

    if variants.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    Ok((current, AlleleVariant::cis(variants)))
}

/// Parse a numeric-only protein repeat (e.g., p.223PA[3], p.179_180AP[9])
/// These patterns have a numeric position (no amino acid) followed by single-letter AA repeat
fn parse_protein_numeric_repeat(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    // Parse the start position (just a number)
    let (remaining, start_num) = digit1.parse(input)?;
    let start: u64 = start_num.parse().map_err(|_| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
    })?;

    // Check for range (pos_pos)
    let (remaining, end_opt) = if remaining.starts_with('_') {
        let (rest, _) = tag("_").parse(remaining)?;
        let (rest, end_num) = digit1.parse(rest)?;
        let end: u64 = end_num.parse().map_err(|_| {
            nom::Err::Error(nom::error::Error::new(
                remaining,
                nom::error::ErrorKind::Digit,
            ))
        })?;
        (rest, Some(end))
    } else {
        (remaining, None)
    };

    // Parse the protein repeat edit (single-letter AAs followed by [count])
    let (remaining, edit) = parse_protein_edit(remaining)?;

    // Verify this is actually a repeat edit (single or multi)
    if !matches!(
        edit,
        ProteinEdit::Repeat { .. } | ProteinEdit::MultiRepeat { .. }
    ) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    // Create protein positions using a placeholder amino acid (the actual AA is in the edit)
    let start_pos = ProtPos::new(AminoAcid::Xaa, start);
    let interval = match end_opt {
        Some(end) => {
            let end_pos = ProtPos::new(AminoAcid::Xaa, end);
            ProtInterval::new(start_pos, end_pos)
        }
        None => ProtInterval::point(start_pos),
    };

    Ok((
        remaining,
        HgvsVariant::Protein(ProteinVariant {
            accession,
            gene_symbol,
            loc_edit: LocEdit::new(interval, edit),
        }),
    ))
}

/// Parse numeric-only protein position (e.g., p.78, p.4894Q, p.42_43insAspAla)
/// These are non-standard notations that appear in ClinVar
/// - p.78 = position 78, unknown effect
/// - p.4894Q = position 4894, substitution to Gln (single-letter code)
/// - p.42_43insAspAla = insertion at positions 42-43
fn parse_protein_numeric_position(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    // Parse the start position number
    let (remaining, start_num) = digit1.parse(input)?;
    let start: u64 = start_num.parse().map_err(|_| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
    })?;

    if start == 0 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    // Check for range (e.g., 42_43)
    let (remaining, end_opt) = if remaining.starts_with('_') {
        let (rest, _) = tag("_").parse(remaining)?;
        let (rest, end_num) = digit1.parse(rest)?;
        let end: u64 = end_num.parse().map_err(|_| {
            nom::Err::Error(nom::error::Error::new(
                remaining,
                nom::error::ErrorKind::Digit,
            ))
        })?;
        (rest, Some(end))
    } else {
        (remaining, None)
    };

    // Create the interval
    let start_pos = ProtPos::new(AminoAcid::Xaa, start);
    let interval = match end_opt {
        Some(end) => {
            let end_pos = ProtPos::new(AminoAcid::Xaa, end);
            ProtInterval::new(start_pos, end_pos)
        }
        None => ProtInterval::point(start_pos),
    };

    // Check for edit (ins, del, dup) or single letter amino acid
    let (remaining, edit) = if remaining.starts_with("ins")
        || remaining.starts_with("del")
        || remaining.starts_with("dup")
    {
        // Parse the protein edit
        parse_protein_edit(remaining)?
    } else if remaining.starts_with(|c: char| c.is_ascii_uppercase()) {
        // Try to parse a single-letter amino acid
        let aa_char = remaining.chars().next().unwrap();
        if let Some(aa) = AminoAcid::from_one_letter(aa_char) {
            // This is a substitution to the given amino acid
            // We use Xaa as the reference since we don't know the original
            (
                &remaining[1..],
                ProteinEdit::Substitution {
                    reference: AminoAcid::Xaa,
                    alternative: aa,
                },
            )
        } else {
            // Not a valid amino acid letter - reject
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }
    } else if remaining.is_empty() || remaining.starts_with(' ') || remaining.starts_with(',') {
        // Just a position, no amino acid - treat as position-only unknown
        (remaining, ProteinEdit::position_unknown())
    } else {
        // Something else follows - reject and let other parsers try
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    };

    Ok((
        remaining,
        HgvsVariant::Protein(ProteinVariant {
            accession,
            gene_symbol,
            loc_edit: LocEdit::new(interval, edit),
        }),
    ))
}

/// Parse protein variant (p.)
fn parse_protein_variant(
    accession: Accession,
    gene_symbol: Option<String>,
) -> impl FnMut(&str) -> IResult<&str, HgvsVariant> {
    move |input: &str| {
        let (input, _) = tag("p.").parse(input)?;

        // First try whole-protein identity (p.= or p.(=)) - no position
        if let Ok((remaining, edit)) = parse_whole_protein_identity(input) {
            let dummy_pos = ProtPos::new(AminoAcid::Met, 1);
            let dummy_interval = ProtInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Protein(ProteinVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(dummy_interval, edit),
                }),
            ));
        }

        // Try no protein production (p.0 or p.0?) - no position
        if let Ok((remaining, edit)) = parse_whole_protein_no_protein(input) {
            let dummy_pos = ProtPos::new(AminoAcid::Met, 1);
            let dummy_interval = ProtInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Protein(ProteinVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(dummy_interval, edit),
                }),
            ));
        }

        // Try whole-protein unknown (p.?) - no position
        if let Ok((remaining, edit)) = parse_whole_protein_unknown(input) {
            let dummy_pos = ProtPos::new(AminoAcid::Met, 1);
            let dummy_interval = ProtInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Protein(ProteinVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(dummy_interval, edit),
                }),
            ));
        }

        // Try protein allele: p.[edit1;edit2;...]
        if input.starts_with('[') {
            if let Ok((remaining, allele)) =
                parse_protein_allele_shorthand(input, accession.clone(), gene_symbol.clone())
            {
                return Ok((remaining, HgvsVariant::Allele(allele)));
            }
        }

        // Try predicted protein change: p.(position+edit) - wrapped in parentheses
        if input.starts_with('(') && !input.starts_with("(=)") {
            // Check for predicted change pattern
            if let Ok((remaining, (interval, edit))) = delimited(
                char('('),
                (parse_prot_interval, parse_protein_edit),
                char(')'),
            )
            .parse(input)
            {
                return Ok((
                    remaining,
                    HgvsVariant::Protein(ProteinVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::with_uncertainty(interval, Mu::Uncertain(edit)),
                    }),
                ));
            }
        }

        // Try numeric-only protein repeat (e.g., p.223PA[3], p.179_180AP[9])
        // These have just a number position (no amino acid) followed by single-letter AA repeat
        if let Ok((remaining, variant)) =
            parse_protein_numeric_repeat(input, accession.clone(), gene_symbol.clone())
        {
            return Ok((remaining, variant));
        }

        // Try numeric-only protein position (e.g., p.78, p.4894Q)
        // These are non-standard but appear in ClinVar
        if let Ok((remaining, variant)) =
            parse_protein_numeric_position(input, accession.clone(), gene_symbol.clone())
        {
            return Ok((remaining, variant));
        }

        // Normal case: position + edit
        let (input, interval) = parse_prot_interval(input)?;

        // Try to parse an edit; if no edit is present and input is exhausted,
        // treat as position-specific unknown (common in VEP output like p.Met1_?4)
        let (input, edit) = match parse_protein_edit(input) {
            Ok((remaining, edit)) => (remaining, edit),
            Err(_) if input.is_empty() || input.starts_with(' ') || input.starts_with(',') => {
                // No edit specified - treat as position-specific unknown
                (input, ProteinEdit::position_unknown())
            }
            Err(e) => return Err(e),
        };

        Ok((
            input,
            HgvsVariant::Protein(ProteinVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, edit),
            }),
        ))
    }
}

/// Parse RNA variant (r.)
fn parse_rna_variant(
    accession: Accession,
    gene_symbol: Option<String>,
) -> impl FnMut(&str) -> IResult<&str, HgvsVariant> {
    move |input: &str| {
        let (input, _) = tag("r.").parse(input)?;

        // Check for cis allele shorthand: r.[118_261del;118_373del]
        if input.starts_with('[') {
            return parse_rna_allele_shorthand(input, accession.clone(), gene_symbol.clone());
        }

        // Check for position-based unknown phase shorthand: r.100A>U(;)200del
        if input.contains("(;)") && !input.starts_with('[') {
            return parse_rna_position_unknown_phase(input, accession.clone(), gene_symbol.clone());
        }

        // First try whole-RNA identity (r.=) - no position
        if let Ok((remaining, edit)) = parse_whole_rna_identity(input) {
            let dummy_pos = crate::hgvs::location::RnaPos {
                base: 1,
                offset: None,
                utr3: false,
            };
            let dummy_interval = RnaInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Rna(RnaVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(dummy_interval, edit),
                }),
            ));
        }

        // Try whole-RNA unknown (r.?) - no position
        if let Ok((remaining, edit)) = parse_whole_rna_unknown(input) {
            let dummy_pos = crate::hgvs::location::RnaPos {
                base: 1,
                offset: None,
                utr3: false,
            };
            let dummy_interval = RnaInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Rna(RnaVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(dummy_interval, edit),
                }),
            ));
        }

        // Try RNA-specific splice pattern (r.spl)
        if let Ok((remaining, edit)) = parse_rna_splice(input) {
            let dummy_pos = crate::hgvs::location::RnaPos {
                base: 1,
                offset: None,
                utr3: false,
            };
            let dummy_interval = RnaInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Rna(RnaVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(dummy_interval, edit),
                }),
            ));
        }

        // Try RNA-specific no product pattern (r.0)
        if let Ok((remaining, edit)) = parse_rna_no_product(input) {
            let dummy_pos = crate::hgvs::location::RnaPos {
                base: 1,
                offset: None,
                utr3: false,
            };
            let dummy_interval = RnaInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Rna(RnaVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(dummy_interval, edit),
                }),
            ));
        }

        // Try predicted RNA variant format: r.(interval edit)
        if let Ok((remaining, variant)) =
            parse_predicted_rna_variant(input, accession.clone(), gene_symbol.clone())
        {
            return Ok((remaining, variant));
        }

        let (input, interval) = parse_rna_interval(input)?;
        let (input, edit) = parse_na_edit(input)?;

        Ok((
            input,
            HgvsVariant::Rna(RnaVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, edit),
            }),
        ))
    }
}

/// Parse whole-RNA identity: r.=
fn parse_whole_rna_identity(input: &str) -> IResult<&str, crate::hgvs::edit::NaEdit> {
    map(tag("="), |_| {
        crate::hgvs::edit::NaEdit::whole_entity_identity()
    })
    .parse(input)
}

/// Parse whole-RNA unknown: r.?
fn parse_whole_rna_unknown(input: &str) -> IResult<&str, crate::hgvs::edit::NaEdit> {
    // Only match terminal ? (not ? followed by _ which indicates uncertain position)
    let (remaining, _) = tag("?").parse(input)?;
    if remaining.starts_with('_') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }
    Ok((remaining, crate::hgvs::edit::NaEdit::whole_entity_unknown()))
}

/// Parse RNA splice pattern: r.spl
fn parse_rna_splice(input: &str) -> IResult<&str, crate::hgvs::edit::NaEdit> {
    map(tag("spl"), |_| crate::hgvs::edit::NaEdit::Splice).parse(input)
}

/// Parse RNA no product pattern: r.0
fn parse_rna_no_product(input: &str) -> IResult<&str, crate::hgvs::edit::NaEdit> {
    // Match "0" but only as a standalone symbol (not followed by digits or other edit characters)
    let (remaining, _) = tag("0").parse(input)?;
    // Make sure it's not followed by more digits (would be a position)
    if remaining.chars().next().is_some_and(|c| c.is_ascii_digit()) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }
    Ok((remaining, crate::hgvs::edit::NaEdit::NoProduct))
}

/// Parse predicted RNA variant: r.(100_101insAUG)
fn parse_predicted_rna_variant(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    // Match opening paren
    let (input, _) = tag("(").parse(input)?;

    // Parse interval and edit inside parens
    let (input, interval) = parse_rna_interval(input)?;
    let (input, edit) = parse_na_edit(input)?;

    // Match closing paren
    let (input, _) = tag(")").parse(input)?;

    Ok((
        input,
        HgvsVariant::Rna(RnaVariant {
            accession,
            gene_symbol,
            loc_edit: LocEdit::new(interval, edit),
        }),
    ))
}

/// Parse position-based unknown phase shorthand for RNA: 100A>U(;)200del
fn parse_rna_position_unknown_phase(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    let mut variants = Vec::with_capacity(4);

    for part in input.split("(;)") {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }

        let (edit_remaining, interval) = parse_rna_interval(part)?;
        let (final_remaining, edit) = parse_na_edit(edit_remaining)?;

        if !final_remaining.trim().is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                final_remaining,
                nom::error::ErrorKind::Tag,
            )));
        }

        variants.push(HgvsVariant::Rna(RnaVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit: LocEdit::new(interval, edit),
        }));
    }

    if variants.len() < 2 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    Ok((
        "",
        HgvsVariant::Allele(AlleleVariant::new(variants, AllelePhase::Unknown)),
    ))
}

/// Parse RNA allele shorthand: [118_261del;118_373del], [100A>U(;)200del], or [100del];[0]
fn parse_rna_allele_shorthand(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    // Must start with [
    if !input.starts_with('[') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Check for trans allele pattern: [var];[var] or [var];[0] or [var];[?]
    if input.contains("];[") {
        return parse_rna_trans_allele_shorthand(input, accession, gene_symbol);
    }

    let close_bracket = input.rfind(']').ok_or_else(|| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Tag))
    })?;

    let content = &input[1..close_bracket];
    let remaining = &input[close_bracket + 1..];

    // Check for mixed phase: contains both ; and (;)
    let has_unknown_phase = content.contains("(;)");
    let has_cis_separator = {
        let temp = content.replace("(;)", "\x00");
        temp.contains(';')
    };

    let phase = if has_unknown_phase {
        AllelePhase::Unknown
    } else {
        AllelePhase::Cis
    };

    let mut variants = Vec::with_capacity(4);

    if has_unknown_phase && has_cis_separator {
        for unknown_group in content.split("(;)") {
            for part in unknown_group.split(';') {
                let part = part.trim();
                if part.is_empty() {
                    continue;
                }

                let (edit_remaining, interval) = parse_rna_interval(part)?;
                let (final_remaining, edit) = parse_na_edit(edit_remaining)?;

                if !final_remaining.trim().is_empty() {
                    return Err(nom::Err::Error(nom::error::Error::new(
                        final_remaining,
                        nom::error::ErrorKind::Tag,
                    )));
                }

                variants.push(HgvsVariant::Rna(RnaVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::new(interval, edit),
                }));
            }
        }
    } else {
        let separator = if has_unknown_phase { "(;)" } else { ";" };

        for part in content.split(separator) {
            let part = part.trim();
            if part.is_empty() {
                continue;
            }

            let (edit_remaining, interval) = parse_rna_interval(part)?;
            let (final_remaining, edit) = parse_na_edit(edit_remaining)?;

            if !final_remaining.trim().is_empty() {
                return Err(nom::Err::Error(nom::error::Error::new(
                    final_remaining,
                    nom::error::ErrorKind::Tag,
                )));
            }

            variants.push(HgvsVariant::Rna(RnaVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, edit),
            }));
        }
    }

    if variants.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    Ok((
        remaining,
        HgvsVariant::Allele(AlleleVariant::new(variants, phase)),
    ))
}

/// Parse RNA trans allele shorthand: [100del];[0], [100del];[?], [var1];[var2]
fn parse_rna_trans_allele_shorthand(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    let mut variants = Vec::with_capacity(2);
    let mut remaining = input;

    while !remaining.is_empty() {
        if !remaining.starts_with('[') {
            break;
        }

        let close_bracket = remaining.find(']').ok_or_else(|| {
            nom::Err::Error(nom::error::Error::new(
                remaining,
                nom::error::ErrorKind::Tag,
            ))
        })?;

        let content = &remaining[1..close_bracket].trim();

        let variant = if *content == "0" {
            HgvsVariant::NullAllele
        } else if *content == "?" {
            HgvsVariant::UnknownAllele
        } else {
            let (edit_remaining, interval) = parse_rna_interval(content)?;
            let (final_remaining, edit) = parse_na_edit(edit_remaining)?;

            if !final_remaining.trim().is_empty() {
                return Err(nom::Err::Error(nom::error::Error::new(
                    final_remaining,
                    nom::error::ErrorKind::Tag,
                )));
            }

            HgvsVariant::Rna(RnaVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, edit),
            })
        };

        variants.push(variant);
        remaining = &remaining[close_bracket + 1..];

        if remaining.starts_with(';') {
            remaining = &remaining[1..];
        }
    }

    if variants.len() < 2 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    Ok((
        remaining,
        HgvsVariant::Allele(AlleleVariant::new(variants, AllelePhase::Trans)),
    ))
}

/// Parse circular DNA variant (o.) - SVD-WG006
fn parse_circular_variant(
    accession: Accession,
    gene_symbol: Option<String>,
) -> impl FnMut(&str) -> IResult<&str, HgvsVariant> {
    move |input: &str| {
        let (input, _) = tag("o.").parse(input)?;
        let (input, interval) = parse_genome_interval(input)?;
        let (input, edit) = parse_na_edit(input)?;

        Ok((
            input,
            HgvsVariant::Circular(CircularVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(interval, edit),
            }),
        ))
    }
}

/// Parse a single variant (not an allele) from input
/// Returns the remaining input and the parsed variant
fn parse_single_variant(input: &str) -> Result<(&str, HgvsVariant), FerroError> {
    let trimmed = input.trim_start();

    // Parse accession
    let (remaining, accession) = parse_accession(trimmed).map_err(|e| FerroError::Parse {
        pos: 0,
        msg: format!("Failed to parse accession: {:?}", e),
        diagnostic: None,
    })?;

    // Parse gene symbol if present
    let (remaining, gene_symbol) = parse_gene_symbol(remaining).map_err(|e| FerroError::Parse {
        pos: trimmed.len() - remaining.len(),
        msg: format!("Failed to parse gene symbol: {:?}", e),
        diagnostic: None,
    })?;

    // Expect colon separator
    let remaining = remaining
        .strip_prefix(':')
        .ok_or_else(|| FerroError::Parse {
            pos: trimmed.len() - remaining.len(),
            msg: "Expected ':' after accession".to_string(),
            diagnostic: None,
        })?;

    // Dispatch based on type prefix
    let result = if remaining.starts_with('[') {
        // ClinVar-style allele: NM_xxx:[c.var1;c.var2] or NC_xxx:[g.var1;g.var2]
        return parse_clinvar_style_allele(remaining, &accession);
    } else if remaining.starts_with("c.") {
        parse_cds_variant(accession, gene_symbol).parse(remaining)
    } else if remaining.starts_with("p.") {
        parse_protein_variant(accession, gene_symbol).parse(remaining)
    } else if remaining.starts_with("g.") {
        parse_genome_variant(accession, gene_symbol).parse(remaining)
    } else if remaining.starts_with("n.") {
        parse_tx_variant(accession, gene_symbol).parse(remaining)
    } else if remaining.starts_with("r.") {
        parse_rna_variant(accession, gene_symbol).parse(remaining)
    } else if remaining.starts_with("m.") {
        parse_mt_variant(accession, gene_symbol).parse(remaining)
    } else if remaining.starts_with("o.") {
        parse_circular_variant(accession, gene_symbol).parse(remaining)
    } else {
        // Try to infer type from accession prefix and parse without type prefix
        let acc_prefix: &str = &accession.prefix;
        // For GenBank-style accessions like "AL513220", the entire name is the prefix,
        // so we also check if it starts with a known genomic prefix
        let is_genomic_prefix = acc_prefix == "NC"
            || acc_prefix == "NT"
            || acc_prefix == "NW"
            || acc_prefix == "AC"
            || acc_prefix == "NG"
            || acc_prefix.starts_with("AL")
            || acc_prefix.starts_with("BX")
            || acc_prefix.starts_with("CR")
            || acc_prefix.starts_with("CT")
            || acc_prefix.starts_with("CU");
        if is_genomic_prefix {
            // Genomic accession - try parsing as g. without the prefix
            let prefixed = format!("g.{}", remaining);
            match parse_genome_variant(accession.clone(), gene_symbol.clone()).parse(&prefixed) {
                Ok((rem, variant)) => {
                    // Adjust remaining - we consumed the same amount as the original remaining
                    let consumed = prefixed.len() - rem.len() - 2; // -2 for "g."
                    Ok((&remaining[consumed..], variant))
                }
                Err(_) => {
                    return Err(FerroError::Parse {
                        pos: trimmed.len() - remaining.len(),
                        msg: format!(
                            "Unknown variant type prefix: expected one of 'c.', 'g.', 'p.', 'n.', 'r.', 'm.', 'o.' but found '{}'",
                            remaining.chars().take(3).collect::<String>()
                        ),
                        diagnostic: None,
                    });
                }
            }
        } else if acc_prefix == "NM"
            || acc_prefix == "NR"
            || acc_prefix == "XM"
            || acc_prefix == "XR"
        {
            // Transcript accession - try parsing as c. without the prefix
            let prefixed = format!("c.{}", remaining);
            match parse_cds_variant(accession.clone(), gene_symbol.clone()).parse(&prefixed) {
                Ok((rem, variant)) => {
                    let consumed = prefixed.len() - rem.len() - 2; // -2 for "c."
                    Ok((&remaining[consumed..], variant))
                }
                Err(_) => {
                    return Err(FerroError::Parse {
                        pos: trimmed.len() - remaining.len(),
                        msg: format!(
                            "Unknown variant type prefix: expected one of 'c.', 'g.', 'p.', 'n.', 'r.', 'm.', 'o.' but found '{}'",
                            remaining.chars().take(3).collect::<String>()
                        ),
                        diagnostic: None,
                    });
                }
            }
        } else {
            return Err(FerroError::Parse {
                pos: trimmed.len() - remaining.len(),
                msg: format!(
                    "Unknown variant type prefix: expected one of 'c.', 'g.', 'p.', 'n.', 'r.', 'm.', 'o.' but found '{}'",
                    remaining.chars().take(3).collect::<String>()
                ),
                diagnostic: None,
            });
        }
    };

    match result {
        Ok((remaining, variant)) => Ok((remaining, variant)),
        Err(e) => Err(FerroError::Parse {
            pos: trimmed.len() - remaining.len(),
            msg: format!("Failed to parse variant: {:?}", e),
            diagnostic: None,
        }),
    }
}

/// Parse a ClinVar-style allele: NM_xxx:[c.var1;c.var2] or NC_xxx:[g.var1;g.var2]
/// This differs from the standard shorthand format NM_xxx:c.[var1;var2] by having
/// the variant type prefix (c., g., etc.) inside the brackets for each variant.
fn parse_clinvar_style_allele<'a>(
    input: &'a str,
    accession: &Accession,
) -> Result<(&'a str, HgvsVariant), FerroError> {
    // Must start with [
    if !input.starts_with('[') {
        return Err(FerroError::Parse {
            pos: 0,
            msg: "ClinVar-style allele must start with '['".to_string(),
            diagnostic: None,
        });
    }

    // Find the matching closing bracket
    let close_bracket = input.rfind(']').ok_or_else(|| FerroError::Parse {
        pos: 0,
        msg: "ClinVar-style allele must end with ']'".to_string(),
        diagnostic: None,
    })?;

    let content = &input[1..close_bracket];
    let remaining = &input[close_bracket + 1..];

    // Check for unknown phase marker (;)
    let has_unknown_phase = content.contains("(;)");

    // Split by either (;) or ; to get individual variant parts
    let parts: Vec<&str> = if has_unknown_phase {
        // Split on (;) first, then on ;
        content
            .split("(;)")
            .flat_map(|part| part.split(';'))
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .collect()
    } else {
        content
            .split(';')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .collect()
    };

    if parts.is_empty() {
        return Err(FerroError::Parse {
            pos: 1,
            msg: "ClinVar-style allele must contain at least one variant".to_string(),
            diagnostic: None,
        });
    }

    // Parse each variant with the shared accession prepended
    let mut variants = Vec::new();
    for part in parts {
        // Prepend shared accession to each variant part (e.g., "c.100A>G" -> "NM_xxx.1:c.100A>G")
        let full_variant = format!("{}:{}", accession, part);
        let (_, variant) = parse_single_variant(&full_variant)?;
        variants.push(variant);
    }

    let phase = if has_unknown_phase {
        AllelePhase::Unknown
    } else {
        AllelePhase::Cis
    };

    Ok((
        remaining,
        HgvsVariant::Allele(AlleleVariant::new(variants, phase)),
    ))
}

/// Parse a cis allele: [var1;var2;...]
fn parse_cis_allele(input: &str) -> Result<HgvsVariant, FerroError> {
    let input = input.trim();

    // Must start with [
    if !input.starts_with('[') {
        return Err(FerroError::Parse {
            pos: 0,
            msg: "Cis allele must start with '['".to_string(),
            diagnostic: None,
        });
    }

    // Find the matching ]
    let close_bracket = input.rfind(']').ok_or_else(|| FerroError::Parse {
        pos: input.len(),
        msg: "Missing closing ']' in cis allele".to_string(),
        diagnostic: None,
    })?;

    // Check nothing after the closing bracket
    let after_bracket = &input[close_bracket + 1..];
    if !after_bracket.trim().is_empty() {
        return Err(FerroError::Parse {
            pos: close_bracket + 1,
            msg: format!("Unexpected trailing characters: '{}'", after_bracket),
            diagnostic: None,
        });
    }

    // Extract content inside brackets
    let content = &input[1..close_bracket];

    // Split by semicolon and parse each variant
    let mut variants = Vec::with_capacity(4);
    for part in content.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let (remaining, variant) = parse_single_variant(part)?;
        if !remaining.trim().is_empty() {
            return Err(FerroError::Parse {
                pos: 0,
                msg: format!("Unexpected content after variant: '{}'", remaining),
                diagnostic: None,
            });
        }
        variants.push(variant);
    }

    if variants.is_empty() {
        return Err(FerroError::Parse {
            pos: 0,
            msg: "Empty allele".to_string(),
            diagnostic: None,
        });
    }

    Ok(HgvsVariant::Allele(AlleleVariant::new(
        variants,
        AllelePhase::Cis,
    )))
}

/// Parse a trans allele: [var1];[var2];...
/// Also supports special markers: [0] (null allele) and [?] (unknown allele)
fn parse_trans_allele(input: &str) -> Result<HgvsVariant, FerroError> {
    let input = input.trim();

    // Must start with [
    if !input.starts_with('[') {
        return Err(FerroError::Parse {
            pos: 0,
            msg: "Trans allele must start with '['".to_string(),
            diagnostic: None,
        });
    }

    // Split by ];[ pattern
    let mut variants = Vec::with_capacity(2);
    let mut remaining = input;

    while !remaining.is_empty() {
        // Must start with [
        if !remaining.starts_with('[') {
            return Err(FerroError::Parse {
                pos: input.len() - remaining.len(),
                msg: "Expected '[' in trans allele".to_string(),
                diagnostic: None,
            });
        }

        // Find closing bracket
        let close_bracket = remaining.find(']').ok_or_else(|| FerroError::Parse {
            pos: input.len() - remaining.len(),
            msg: "Missing closing ']' in trans allele".to_string(),
            diagnostic: None,
        })?;

        // Parse the variant inside the brackets
        let variant_str = &remaining[1..close_bracket].trim();

        // Check for special markers
        let variant = if *variant_str == "0" {
            HgvsVariant::NullAllele
        } else if *variant_str == "?" {
            HgvsVariant::UnknownAllele
        } else {
            let (var_remaining, var) = parse_single_variant(variant_str)?;
            if !var_remaining.trim().is_empty() {
                return Err(FerroError::Parse {
                    pos: 0,
                    msg: format!("Unexpected content after variant: '{}'", var_remaining),
                    diagnostic: None,
                });
            }
            var
        };
        variants.push(variant);

        // Move past the closing bracket
        remaining = &remaining[close_bracket + 1..];

        // Check for more variants (semicolon separator)
        if remaining.starts_with(';') {
            remaining = &remaining[1..];
        } else if !remaining.trim().is_empty() {
            return Err(FerroError::Parse {
                pos: input.len() - remaining.len(),
                msg: format!("Expected ';' or end of input, found: '{}'", remaining),
                diagnostic: None,
            });
        }
    }

    if variants.len() < 2 {
        return Err(FerroError::Parse {
            pos: 0,
            msg: "Trans allele requires at least two variants".to_string(),
            diagnostic: None,
        });
    }

    Ok(HgvsVariant::Allele(AlleleVariant::new(
        variants,
        AllelePhase::Trans,
    )))
}

/// Create an identity variant based on the type of the reference variant.
/// Used for mosaic/chimeric notation where "=" means reference allele.
fn create_identity_variant_from(reference: &HgvsVariant) -> Result<HgvsVariant, FerroError> {
    use crate::hgvs::edit::NaEdit;
    use crate::hgvs::location::{GenomePos, RnaPos, TxPos};

    match reference {
        HgvsVariant::Cds(cds) => {
            // Create c.= (whole-CDS identity)
            let dummy_pos = CdsPos {
                base: 1,
                offset: None,
                utr3: false,
            };
            let dummy_interval = CdsInterval::point(dummy_pos);
            Ok(HgvsVariant::Cds(CdsVariant {
                accession: cds.accession.clone(),
                gene_symbol: cds.gene_symbol.clone(),
                loc_edit: LocEdit::new(dummy_interval, NaEdit::whole_entity_identity()),
            }))
        }
        HgvsVariant::Genome(genome) => {
            // Create g.= (whole-genome identity)
            let dummy_interval = GenomeInterval::point(GenomePos::new(1));
            Ok(HgvsVariant::Genome(GenomeVariant {
                accession: genome.accession.clone(),
                gene_symbol: genome.gene_symbol.clone(),
                loc_edit: LocEdit::new(dummy_interval, NaEdit::whole_entity_identity()),
            }))
        }
        HgvsVariant::Tx(tx) => {
            // Create n.= (whole-transcript identity)
            let dummy_interval = TxInterval::point(TxPos::new(1));
            Ok(HgvsVariant::Tx(TxVariant {
                accession: tx.accession.clone(),
                gene_symbol: tx.gene_symbol.clone(),
                loc_edit: LocEdit::new(dummy_interval, NaEdit::whole_entity_identity()),
            }))
        }
        HgvsVariant::Rna(rna) => {
            // Create r.= (whole-RNA identity)
            let dummy_interval = RnaInterval::point(RnaPos::new(1));
            Ok(HgvsVariant::Rna(RnaVariant {
                accession: rna.accession.clone(),
                gene_symbol: rna.gene_symbol.clone(),
                loc_edit: LocEdit::new(dummy_interval, NaEdit::whole_entity_identity()),
            }))
        }
        HgvsVariant::Mt(mt) => {
            // Create m.= (whole-mitochondrial identity)
            let dummy_interval = GenomeInterval::point(GenomePos::new(1));
            Ok(HgvsVariant::Mt(MtVariant {
                accession: mt.accession.clone(),
                gene_symbol: mt.gene_symbol.clone(),
                loc_edit: LocEdit::new(dummy_interval, NaEdit::whole_entity_identity()),
            }))
        }
        HgvsVariant::Circular(circular) => {
            // Create o.= (whole-circular identity)
            let dummy_interval = GenomeInterval::point(GenomePos::new(1));
            Ok(HgvsVariant::Circular(CircularVariant {
                accession: circular.accession.clone(),
                gene_symbol: circular.gene_symbol.clone(),
                loc_edit: LocEdit::new(dummy_interval, NaEdit::whole_entity_identity()),
            }))
        }
        HgvsVariant::Protein(protein) => {
            // Create p.= (whole-protein identity)
            let dummy_interval = ProtInterval::point(ProtPos::new(AminoAcid::Met, 1));
            Ok(HgvsVariant::Protein(ProteinVariant {
                accession: protein.accession.clone(),
                gene_symbol: protein.gene_symbol.clone(),
                loc_edit: LocEdit::new(
                    dummy_interval,
                    crate::hgvs::edit::ProteinEdit::whole_protein_identity(),
                ),
            }))
        }
        _ => Err(FerroError::Parse {
            pos: 0,
            msg: "Cannot create identity variant for this variant type".to_string(),
            diagnostic: None,
        }),
    }
}

/// Parse a mosaic: var1/var2 (single slash)
/// Supports "=" as shorthand for reference allele (no change)
/// Supports accession inheritance: NM_000088.3:c.100A>G/c.200C>T
fn parse_mosaic_allele(input: &str) -> Result<HgvsVariant, FerroError> {
    let input = input.trim();

    // Split by single slash (but not double slash)
    let mut variants = Vec::with_capacity(2);
    let mut start = 0;
    let mut first_variant: Option<HgvsVariant> = None;
    let mut first_accession: Option<Accession> = None;
    let mut first_gene_symbol: Option<String> = None;

    while start < input.len() {
        // Find next single slash (not double)
        let mut found_slash = None;
        let bytes = input.as_bytes();
        for i in start..input.len() {
            if bytes[i] == b'/' {
                // Check if double slash
                if i + 1 < input.len() && bytes[i + 1] == b'/' {
                    // This is a chimeric separator, not mosaic
                    return Err(FerroError::Parse {
                        pos: i,
                        msg: "Found '//' (chimeric) but expected '/' (mosaic)".to_string(),
                        diagnostic: None,
                    });
                }
                found_slash = Some(i);
                break;
            }
        }

        let end = found_slash.unwrap_or(input.len());
        let variant_str = &input[start..end].trim();

        // Check for "=" shorthand meaning reference allele
        let variant = if *variant_str == "=" {
            // Need a reference variant to determine the type
            match &first_variant {
                Some(ref_var) => create_identity_variant_from(ref_var)?,
                None => {
                    return Err(FerroError::Parse {
                        pos: start,
                        msg: "Cannot use '=' as first variant in mosaic notation".to_string(),
                        diagnostic: None,
                    });
                }
            }
        } else {
            // Try parsing as a full variant first
            match parse_single_variant(variant_str) {
                Ok((remaining, var)) => {
                    if !remaining.trim().is_empty() {
                        return Err(FerroError::Parse {
                            pos: 0,
                            msg: format!("Unexpected content after variant: '{}'", remaining),
                            diagnostic: None,
                        });
                    }
                    var
                }
                Err(_) if first_accession.is_some() => {
                    // Try parsing as a variant with inherited accession
                    // Check if it starts with a type prefix
                    parse_variant_with_inherited_accession(
                        variant_str,
                        first_accession.as_ref().unwrap(),
                        first_gene_symbol.as_ref(),
                    )?
                }
                Err(e) => return Err(e),
            }
        };

        if first_variant.is_none() {
            first_variant = Some(variant.clone());
            // Extract accession and gene symbol from first variant
            match &variant {
                HgvsVariant::Cds(v) => {
                    first_accession = Some(v.accession.clone());
                    first_gene_symbol = v.gene_symbol.clone();
                }
                HgvsVariant::Genome(v) => {
                    first_accession = Some(v.accession.clone());
                    first_gene_symbol = v.gene_symbol.clone();
                }
                HgvsVariant::Tx(v) => {
                    first_accession = Some(v.accession.clone());
                    first_gene_symbol = v.gene_symbol.clone();
                }
                HgvsVariant::Rna(v) => {
                    first_accession = Some(v.accession.clone());
                    first_gene_symbol = v.gene_symbol.clone();
                }
                HgvsVariant::Protein(v) => {
                    first_accession = Some(v.accession.clone());
                    first_gene_symbol = v.gene_symbol.clone();
                }
                HgvsVariant::Mt(v) => {
                    first_accession = Some(v.accession.clone());
                    first_gene_symbol = v.gene_symbol.clone();
                }
                HgvsVariant::Circular(v) => {
                    first_accession = Some(v.accession.clone());
                    first_gene_symbol = v.gene_symbol.clone();
                }
                _ => {}
            }
        }
        variants.push(variant);

        if found_slash.is_some() {
            start = end + 1;
        } else {
            break;
        }
    }

    if variants.len() < 2 {
        return Err(FerroError::Parse {
            pos: 0,
            msg: "Mosaic allele requires at least two variants".to_string(),
            diagnostic: None,
        });
    }

    Ok(HgvsVariant::Allele(AlleleVariant::new(
        variants,
        AllelePhase::Mosaic,
    )))
}

/// Parse a variant that inherits its accession from a previous variant
/// Used for mosaic/chimeric notation like NM_000088.3:c.100A>G/c.200C>T
fn parse_variant_with_inherited_accession(
    input: &str,
    accession: &Accession,
    gene_symbol: Option<&String>,
) -> Result<HgvsVariant, FerroError> {
    let input = input.trim();

    // Check for type prefix and parse accordingly
    if let Some(rest) = input.strip_prefix("c.") {
        let (remaining, interval) = parse_cds_interval(rest).map_err(|e| FerroError::Parse {
            pos: 2,
            msg: format!("Failed to parse CDS interval: {:?}", e),
            diagnostic: None,
        })?;
        let (remaining, edit) = parse_na_edit(remaining).map_err(|e| FerroError::Parse {
            pos: input.len() - remaining.len(),
            msg: format!("Failed to parse edit: {:?}", e),
            diagnostic: None,
        })?;
        if !remaining.trim().is_empty() {
            return Err(FerroError::Parse {
                pos: input.len() - remaining.len(),
                msg: format!("Unexpected content: '{}'", remaining),
                diagnostic: None,
            });
        }
        Ok(HgvsVariant::Cds(CdsVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.cloned(),
            loc_edit: LocEdit::new(interval, edit),
        }))
    } else if let Some(rest) = input.strip_prefix("g.") {
        let (remaining, interval) = parse_genome_interval(rest).map_err(|e| FerroError::Parse {
            pos: 2,
            msg: format!("Failed to parse genome interval: {:?}", e),
            diagnostic: None,
        })?;
        let (remaining, edit) = parse_na_edit(remaining).map_err(|e| FerroError::Parse {
            pos: input.len() - remaining.len(),
            msg: format!("Failed to parse edit: {:?}", e),
            diagnostic: None,
        })?;
        if !remaining.trim().is_empty() {
            return Err(FerroError::Parse {
                pos: input.len() - remaining.len(),
                msg: format!("Unexpected content: '{}'", remaining),
                diagnostic: None,
            });
        }
        Ok(HgvsVariant::Genome(GenomeVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.cloned(),
            loc_edit: LocEdit::new(interval, edit),
        }))
    } else if let Some(rest) = input.strip_prefix("p.") {
        // Try whole-protein identity (p.= or p.(=))
        if let Ok((remaining, edit)) = parse_whole_protein_identity(rest) {
            let dummy_pos = ProtPos::new(AminoAcid::Met, 1);
            let dummy_interval = ProtInterval::point(dummy_pos);
            if !remaining.trim().is_empty() {
                return Err(FerroError::Parse {
                    pos: input.len() - remaining.len(),
                    msg: format!("Unexpected content: '{}'", remaining),
                    diagnostic: None,
                });
            }
            return Ok(HgvsVariant::Protein(ProteinVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.cloned(),
                loc_edit: LocEdit::new(dummy_interval, edit),
            }));
        }
        // Try whole-protein unknown (p.?)
        if let Ok((remaining, edit)) = parse_whole_protein_unknown(rest) {
            let dummy_pos = ProtPos::new(AminoAcid::Met, 1);
            let dummy_interval = ProtInterval::point(dummy_pos);
            if !remaining.trim().is_empty() {
                return Err(FerroError::Parse {
                    pos: input.len() - remaining.len(),
                    msg: format!("Unexpected content: '{}'", remaining),
                    diagnostic: None,
                });
            }
            return Ok(HgvsVariant::Protein(ProteinVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.cloned(),
                loc_edit: LocEdit::new(dummy_interval, edit),
            }));
        }
        // Try no protein (p.0)
        if let Ok((remaining, edit)) = parse_whole_protein_no_protein(rest) {
            let dummy_pos = ProtPos::new(AminoAcid::Met, 1);
            let dummy_interval = ProtInterval::point(dummy_pos);
            if !remaining.trim().is_empty() {
                return Err(FerroError::Parse {
                    pos: input.len() - remaining.len(),
                    msg: format!("Unexpected content: '{}'", remaining),
                    diagnostic: None,
                });
            }
            return Ok(HgvsVariant::Protein(ProteinVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.cloned(),
                loc_edit: LocEdit::new(dummy_interval, edit),
            }));
        }
        // Normal case: position + edit
        let (remaining, interval) = parse_prot_interval(rest).map_err(|e| FerroError::Parse {
            pos: 2,
            msg: format!("Failed to parse protein interval: {:?}", e),
            diagnostic: None,
        })?;
        let (remaining, edit) = parse_protein_edit(remaining).map_err(|e| FerroError::Parse {
            pos: input.len() - remaining.len(),
            msg: format!("Failed to parse protein edit: {:?}", e),
            diagnostic: None,
        })?;
        if !remaining.trim().is_empty() {
            return Err(FerroError::Parse {
                pos: input.len() - remaining.len(),
                msg: format!("Unexpected content: '{}'", remaining),
                diagnostic: None,
            });
        }
        Ok(HgvsVariant::Protein(ProteinVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.cloned(),
            loc_edit: LocEdit::new(interval, edit),
        }))
    } else if let Some(rest) = input.strip_prefix("n.") {
        let (remaining, interval) = parse_tx_interval(rest).map_err(|e| FerroError::Parse {
            pos: 2,
            msg: format!("Failed to parse transcript interval: {:?}", e),
            diagnostic: None,
        })?;
        let (remaining, edit) = parse_na_edit(remaining).map_err(|e| FerroError::Parse {
            pos: input.len() - remaining.len(),
            msg: format!("Failed to parse edit: {:?}", e),
            diagnostic: None,
        })?;
        if !remaining.trim().is_empty() {
            return Err(FerroError::Parse {
                pos: input.len() - remaining.len(),
                msg: format!("Unexpected content: '{}'", remaining),
                diagnostic: None,
            });
        }
        Ok(HgvsVariant::Tx(TxVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.cloned(),
            loc_edit: LocEdit::new(interval, edit),
        }))
    } else if let Some(rest) = input.strip_prefix("r.") {
        let (remaining, interval) = parse_rna_interval(rest).map_err(|e| FerroError::Parse {
            pos: 2,
            msg: format!("Failed to parse RNA interval: {:?}", e),
            diagnostic: None,
        })?;
        let (remaining, edit) = parse_na_edit(remaining).map_err(|e| FerroError::Parse {
            pos: input.len() - remaining.len(),
            msg: format!("Failed to parse edit: {:?}", e),
            diagnostic: None,
        })?;
        if !remaining.trim().is_empty() {
            return Err(FerroError::Parse {
                pos: input.len() - remaining.len(),
                msg: format!("Unexpected content: '{}'", remaining),
                diagnostic: None,
            });
        }
        Ok(HgvsVariant::Rna(RnaVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.cloned(),
            loc_edit: LocEdit::new(interval, edit),
        }))
    } else if let Some(rest) = input.strip_prefix("m.") {
        let (remaining, interval) = parse_genome_interval(rest).map_err(|e| FerroError::Parse {
            pos: 2,
            msg: format!("Failed to parse mitochondrial interval: {:?}", e),
            diagnostic: None,
        })?;
        let (remaining, edit) = parse_na_edit(remaining).map_err(|e| FerroError::Parse {
            pos: input.len() - remaining.len(),
            msg: format!("Failed to parse edit: {:?}", e),
            diagnostic: None,
        })?;
        if !remaining.trim().is_empty() {
            return Err(FerroError::Parse {
                pos: input.len() - remaining.len(),
                msg: format!("Unexpected content: '{}'", remaining),
                diagnostic: None,
            });
        }
        Ok(HgvsVariant::Mt(MtVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.cloned(),
            loc_edit: LocEdit::new(interval, edit),
        }))
    } else {
        Err(FerroError::Parse {
            pos: 0,
            msg: format!(
                "Unknown variant type prefix: expected one of 'c.', 'g.', 'p.', 'n.', 'r.', 'm.' but found '{}'",
                input.chars().take(3).collect::<String>()
            ),
            diagnostic: None,
        })
    }
}

/// Parse a chimeric: var1//var2 (double slash)
/// Supports "=" as shorthand for reference allele (no change)
fn parse_chimeric_allele(input: &str) -> Result<HgvsVariant, FerroError> {
    let input = input.trim();

    // Split by double slash
    let parts: Vec<&str> = input.split("//").collect();

    if parts.len() < 2 {
        return Err(FerroError::Parse {
            pos: 0,
            msg: "Chimeric allele requires '//' separator".to_string(),
            diagnostic: None,
        });
    }

    let mut variants = Vec::with_capacity(parts.len());
    let mut first_variant: Option<HgvsVariant> = None;

    for part in parts {
        let part = part.trim();
        if part.is_empty() {
            return Err(FerroError::Parse {
                pos: 0,
                msg: "Empty variant in chimeric allele".to_string(),
                diagnostic: None,
            });
        }

        // Check for "=" shorthand meaning reference allele
        let variant = if part == "=" {
            match &first_variant {
                Some(ref_var) => create_identity_variant_from(ref_var)?,
                None => {
                    return Err(FerroError::Parse {
                        pos: 0,
                        msg: "Cannot use '=' as first variant in chimeric notation".to_string(),
                        diagnostic: None,
                    });
                }
            }
        } else {
            let (remaining, var) = parse_single_variant(part)?;
            if !remaining.trim().is_empty() {
                return Err(FerroError::Parse {
                    pos: 0,
                    msg: format!("Unexpected content after variant: '{}'", remaining),
                    diagnostic: None,
                });
            }
            var
        };

        if first_variant.is_none() {
            first_variant = Some(variant.clone());
        }
        variants.push(variant);
    }

    Ok(HgvsVariant::Allele(AlleleVariant::new(
        variants,
        AllelePhase::Chimeric,
    )))
}

/// Parse an unknown phase allele: [var1(;)var2(;)...]
fn parse_unknown_phase_allele(input: &str) -> Result<HgvsVariant, FerroError> {
    let input = input.trim();

    // Must start with [ and end with ]
    if !input.starts_with('[') || !input.ends_with(']') {
        return Err(FerroError::Parse {
            pos: 0,
            msg: "Unknown phase allele must be wrapped in brackets".to_string(),
            diagnostic: None,
        });
    }

    // Extract content inside brackets
    let content = &input[1..input.len() - 1];

    // Split by (;) and parse each variant
    let mut variants = Vec::with_capacity(2);
    for part in content.split("(;)") {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let (remaining, variant) = parse_single_variant(part)?;
        if !remaining.trim().is_empty() {
            return Err(FerroError::Parse {
                pos: 0,
                msg: format!("Unexpected content after variant: '{}'", remaining),
                diagnostic: None,
            });
        }
        variants.push(variant);
    }

    if variants.is_empty() {
        return Err(FerroError::Parse {
            pos: 0,
            msg: "Empty allele".to_string(),
            diagnostic: None,
        });
    }

    Ok(HgvsVariant::Allele(AlleleVariant::new(
        variants,
        AllelePhase::Unknown,
    )))
}

/// Determine the type of allele from the input string
fn detect_allele_type(input: &str) -> Option<&'static str> {
    let input = input.trim();

    // Check for trans: [var];[var] pattern
    if input.starts_with('[') && input.contains("];[") {
        return Some("trans");
    }

    // Check for unknown phase: [var(;)var] pattern
    if input.starts_with('[') && input.ends_with(']') && input.contains("(;)") {
        return Some("unknown_phase");
    }

    // Check for cis: [var;var] pattern (single bracket wrapping all)
    if input.starts_with('[') && input.ends_with(']') && !input.contains("];[") {
        // Verify it contains ; inside brackets
        let inner = &input[1..input.len() - 1];
        if inner.contains(';') {
            return Some("cis");
        }
    }

    // Check for chimeric: var//var (double slash - check before single slash)
    if input.contains("//") {
        return Some("chimeric");
    }

    // Check for mosaic: var/var (single slash)
    // But not if it looks like a regular variant (no slash in accession part)
    if input.contains('/') {
        // Make sure the slash is between variants, not inside one
        // A simple heuristic: if there's a colon before the slash, it might be a mosaic
        if let Some(slash_pos) = input.find('/') {
            let before_slash = &input[..slash_pos];
            if before_slash.contains(':') {
                return Some("mosaic");
            }
        }
    }

    // Check for RNA fusion: accession:r.interval::accession:r.interval
    // The pattern has ":r." followed by "::" somewhere in the string
    if input.contains(":r.") && input.contains("::") {
        // Find :: and ensure it's between two r. parts
        if let Some(fusion_pos) = input.find("::") {
            let before = &input[..fusion_pos];
            let after = &input[fusion_pos + 2..];
            // Both parts should contain r. or the second should start with accession:r.
            if before.contains(":r.") && (after.contains(":r.") || after.starts_with("r.")) {
                return Some("rna_fusion");
            }
        }
    }

    None
}

/// Parse an RNA fusion variant - SVD-WG007
///
/// Format: `accession:r.interval::accession:r.interval`
/// Example: `NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924`
fn parse_rna_fusion(input: &str) -> Result<HgvsVariant, FerroError> {
    let input = input.trim();

    // Find the :: separator
    let fusion_pos = input.find("::").ok_or_else(|| FerroError::Parse {
        pos: 0,
        msg: "RNA fusion requires '::' separator".to_string(),
        diagnostic: None,
    })?;

    let five_prime_str = &input[..fusion_pos];
    let three_prime_str = &input[fusion_pos + 2..];

    // Parse the 5' partner
    let five_prime = parse_rna_fusion_breakpoint(five_prime_str)?;

    // Parse the 3' partner
    let three_prime = parse_rna_fusion_breakpoint(three_prime_str)?;

    Ok(HgvsVariant::RnaFusion(RnaFusionVariant::new(
        five_prime,
        three_prime,
    )))
}

/// Parse a single RNA fusion breakpoint (accession:r.interval)
fn parse_rna_fusion_breakpoint(input: &str) -> Result<RnaFusionBreakpoint, FerroError> {
    let input = input.trim();

    // Parse accession
    let (remaining, accession) = parse_accession(input).map_err(|e| FerroError::Parse {
        pos: 0,
        msg: format!("Failed to parse fusion accession: {:?}", e),
        diagnostic: None,
    })?;

    // Parse gene symbol if present
    let (remaining, gene_symbol) = parse_gene_symbol(remaining).map_err(|e| FerroError::Parse {
        pos: input.len() - remaining.len(),
        msg: format!("Failed to parse gene symbol: {:?}", e),
        diagnostic: None,
    })?;

    // Expect :r.
    let remaining = remaining
        .strip_prefix(":r.")
        .ok_or_else(|| FerroError::Parse {
            pos: input.len() - remaining.len(),
            msg: "Expected ':r.' in RNA fusion breakpoint".to_string(),
            diagnostic: None,
        })?;

    // Parse the RNA interval
    let (remaining, interval) = parse_rna_interval(remaining).map_err(|e| FerroError::Parse {
        pos: input.len() - remaining.len(),
        msg: format!("Failed to parse RNA interval: {:?}", e),
        diagnostic: None,
    })?;

    if !remaining.is_empty() {
        return Err(FerroError::Parse {
            pos: input.len() - remaining.len(),
            msg: format!("Unexpected content after RNA interval: '{}'", remaining),
            diagnostic: None,
        });
    }

    Ok(RnaFusionBreakpoint {
        accession,
        gene_symbol,
        interval,
    })
}

/// Parse a complete HGVS variant string
///
/// The parser is optimized with:
/// - Early type prefix checking to avoid unnecessary work
/// - Variant types ordered by frequency (c. > p. > g. > n. > r. > m.)
/// - Allele detection for compound variants
pub fn parse_variant(input: &str) -> Result<HgvsVariant, FerroError> {
    let input = input.trim();

    // Check for allele patterns and RNA fusion first
    if let Some(allele_type) = detect_allele_type(input) {
        return match allele_type {
            "cis" => parse_cis_allele(input),
            "trans" => parse_trans_allele(input),
            "mosaic" => parse_mosaic_allele(input),
            "chimeric" => parse_chimeric_allele(input),
            "unknown_phase" => parse_unknown_phase_allele(input),
            "rna_fusion" => parse_rna_fusion(input),
            _ => unreachable!(),
        };
    }

    // Parse as single variant
    let (remaining, variant) = parse_single_variant(input)?;

    if remaining.is_empty() {
        Ok(variant)
    } else {
        Err(FerroError::Parse {
            pos: input.len() - remaining.len(),
            msg: format!("Unexpected trailing characters: '{}'", remaining),
            diagnostic: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_genomic_substitution() {
        let variant = parse_variant("NC_000001.11:g.12345A>G").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(
            variant.accession().expect("Expected accession").full(),
            "NC_000001.11"
        );
        assert_eq!(format!("{}", variant), "NC_000001.11:g.12345A>G");
    }

    #[test]
    fn test_parse_genomic_deletion() {
        let variant = parse_variant("NC_000001.11:g.12345_12350del").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(format!("{}", variant), "NC_000001.11:g.12345_12350del");
    }

    #[test]
    fn test_parse_cds_substitution() {
        let variant = parse_variant("NM_000088.3:c.459A>G").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.459A>G");
    }

    #[test]
    fn test_parse_cds_deletion() {
        let variant = parse_variant("NM_000088.3:c.459del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
    }

    #[test]
    fn test_parse_cds_intronic() {
        let variant = parse_variant("NM_000088.3:c.459+5G>A").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.459+5G>A");
    }

    #[test]
    fn test_parse_cds_intronic_uncertain_offset() {
        // Uncertain positive offset (+?)
        let variant = parse_variant("NM_000088.3:c.100+?del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.100+?del");

        // Uncertain negative offset (-?)
        let variant = parse_variant("NM_000088.3:c.100-?del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.100-?del");

        // Range with uncertain offsets on both ends
        let variant = parse_variant("NM_000088.3:c.548-?_5193+?del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.548-?_5193+?del");

        // Complex uncertain position with range boundary
        let variant = parse_variant("NM_000088.3:c.(?_-1)_328+?del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.(?_-1)_328+?del");
    }

    #[test]
    fn test_parse_cds_inverted_intronic_range() {
        // Inverted intronic ranges where base positions appear "inverted" but are valid
        // These represent intronic regions spanning exon boundaries
        let variant = parse_variant("NM_000088.3:c.85-47_84+48del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.85-47_84+48del");

        let variant = parse_variant("NM_000088.3:c.86-46_85+47insA").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.86-46_85+47insA");

        // Delins with explicit deleted sequence (normalizes to delins)
        let variant = parse_variant("NM_000088.3:c.85-47_84+48delGCCAinsG").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.85-47_84+48delinsG");
    }

    #[test]
    fn test_parse_cds_utr() {
        let variant = parse_variant("NM_000088.3:c.-20G>A").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.-20G>A");

        let variant = parse_variant("NM_000088.3:c.*50G>A").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.*50G>A");
    }

    #[test]
    fn test_parse_whole_cds_unknown() {
        // Whole-CDS unknown: c.?
        let variant = parse_variant("NM_000088.3:c.?").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.?");
    }

    #[test]
    fn test_parse_transcript_variant() {
        let variant = parse_variant("NR_000001.1:n.100A>G").unwrap();
        assert!(matches!(variant, HgvsVariant::Tx(_)));
    }

    #[test]
    fn test_parse_protein_substitution() {
        let variant = parse_variant("NP_000079.2:p.Val600Glu").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
    }

    #[test]
    fn test_parse_protein_frameshift() {
        let variant = parse_variant("NP_000079.2:p.Lys23fs").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));

        // Frameshift with new amino acid
        let variant = parse_variant("NP_000079.2:p.Arg97ProfsTer23").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.Arg97ProfsTer23");
    }

    #[test]
    fn test_parse_insertion() {
        let variant = parse_variant("NM_000088.3:c.459_460insATG").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.459_460insATG");
    }

    #[test]
    fn test_parse_duplication() {
        let variant = parse_variant("NM_000088.3:c.459dup").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
    }

    #[test]
    fn test_parse_invalid() {
        let result = parse_variant("invalid");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_with_trailing() {
        let result = parse_variant("NC_000001.11:g.12345A>G extra");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_uncertain_genomic_deletion() {
        // Uncertain range deletion: (100_200)del
        let variant = parse_variant("NC_000001.11:g.(12345_12350)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(format!("{}", variant), "NC_000001.11:g.(12345)_(12350)del");
    }

    #[test]
    fn test_parse_uncertain_cds_deletion() {
        // Uncertain range deletion
        let variant = parse_variant("NM_000088.3:c.(100_200)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.(100)_(200)del");
    }

    #[test]
    fn test_parse_individual_uncertain_positions() {
        // Individual uncertain positions
        let variant = parse_variant("NC_000001.11:g.(100)_(200)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(format!("{}", variant), "NC_000001.11:g.(100)_(200)del");
    }

    #[test]
    fn test_parse_mixed_uncertainty() {
        // One certain, one uncertain
        let variant = parse_variant("NC_000001.11:g.100_(200)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(format!("{}", variant), "NC_000001.11:g.100_(200)del");
    }

    #[test]
    fn test_parse_rna_substitution() {
        let variant = parse_variant("NM_000088.3:r.100a>g").unwrap();
        assert!(matches!(variant, HgvsVariant::Rna(_)));
        // RNA uses lowercase nucleotides per HGVS spec
        assert_eq!(format!("{}", variant), "NM_000088.3:r.100a>g");
    }

    #[test]
    fn test_parse_rna_deletion() {
        let variant = parse_variant("NM_000088.3:r.100_200del").unwrap();
        assert!(matches!(variant, HgvsVariant::Rna(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:r.100_200del");
    }

    #[test]
    fn test_parse_rna_insertion() {
        let variant = parse_variant("NM_000088.3:r.100_101insaug").unwrap();
        assert!(matches!(variant, HgvsVariant::Rna(_)));
        // RNA uses lowercase nucleotides per HGVS spec
        assert_eq!(format!("{}", variant), "NM_000088.3:r.100_101insaug");
    }

    #[test]
    fn test_parse_rna_with_offset() {
        // Intronic RNA position
        let variant = parse_variant("NM_000088.3:r.100+5a>g").unwrap();
        assert!(matches!(variant, HgvsVariant::Rna(_)));
        // RNA uses lowercase nucleotides per HGVS spec
        assert_eq!(format!("{}", variant), "NM_000088.3:r.100+5a>g");
    }

    #[test]
    fn test_parse_rna_uncertain() {
        let variant = parse_variant("NM_000088.3:r.(100_200)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Rna(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:r.(100)_(200)del");
    }

    #[test]
    fn test_parse_rna_utr() {
        // 5' UTR position - RNA uses lowercase nucleotides per HGVS spec
        let variant = parse_variant("NM_000088.3:r.-14a>c").unwrap();
        assert!(matches!(variant, HgvsVariant::Rna(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:r.-14a>c");

        // 3' UTR position - RNA uses lowercase nucleotides per HGVS spec
        let variant = parse_variant("NM_000088.3:r.*41u>a").unwrap();
        assert!(matches!(variant, HgvsVariant::Rna(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:r.*41u>a");
    }

    #[test]
    fn test_parse_protein_identity() {
        // Whole protein identity (p.=)
        let variant = parse_variant("NP_000079.2:p.=").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.=");
    }

    #[test]
    fn test_parse_protein_identity_predicted() {
        // Predicted whole protein identity (p.(=))
        let variant = parse_variant("NP_000079.2:p.(=)").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.(=)");
    }

    #[test]
    fn test_parse_protein_position_identity() {
        // Position-specific identity (p.Val600=)
        let variant = parse_variant("NP_000079.2:p.Val600=").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.Val600=");
    }

    #[test]
    fn test_parse_protein_no_protein() {
        // No protein produced (p.0)
        let variant = parse_variant("NP_000079.2:p.0").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.0");
    }

    #[test]
    fn test_parse_protein_no_protein_predicted() {
        // Predicted no protein (p.0?)
        let variant = parse_variant("NP_000079.2:p.0?").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.0?");
    }

    #[test]
    fn test_parse_protein_start_codon_substitution() {
        // Start codon substitution (p.Met1Leu)
        let variant = parse_variant("NP_000079.2:p.Met1Leu").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.Met1Leu");
    }

    #[test]
    fn test_parse_protein_start_codon_deletion() {
        // Start codon deletion (p.Met1del)
        let variant = parse_variant("NP_000079.2:p.Met1del").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.Met1del");
    }

    #[test]
    fn test_parse_protein_start_codon_extension() {
        // Start codon N-terminal extension (p.Met1ext-5)
        let variant = parse_variant("NP_000079.2:p.Met1ext-5").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.Met1ext-5");
    }

    #[test]
    fn test_parse_protein_start_codon_uncertain() {
        // Uncertain effect on start codon (p.Met1?)
        let variant = parse_variant("NP_000079.2:p.Met1?").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.Met1?");
    }

    #[test]
    fn test_parse_protein_unusual_range() {
        // Unusual protein range with unknown amino acid (p.Met1_?4)
        // This is VEP-style notation where ?4 means "unknown amino acid at position 4"
        let variant = parse_variant("XP_005261260.1:p.Met1_?4").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        // The output normalizes ?4 to Xaa4 and adds implicit ? edit
        assert_eq!(format!("{}", variant), "XP_005261260.1:p.Met1_Xaa4?");

        // Another unusual range
        let variant = parse_variant("NP_000079.2:p.Gly100_?105").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.Gly100_Xaa105?");
    }

    #[test]
    fn test_parse_protein_uncertain_range_boundaries() {
        // Uncertain start boundary: (?_Met1)_Glu2del
        // Means somewhere from unknown to Met1, then definitely to Glu2
        let variant = parse_variant("NP_001005484.1:p.(?_Met1)_Glu2del").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_001005484.1:p.(?_Met1)_Glu2del");

        // Uncertain start boundary with known endpoint: (Met1_?)_Glu2del
        // Means somewhere from Met1 to unknown, then definitely to Glu2
        let variant = parse_variant("NP_001005484.1:p.(Met1_?)_Glu2del").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_001005484.1:p.(Met1_?)_Glu2del");

        // Uncertain end boundary: Met1_(Glu2_?)del
        // Means definitely Met1, then somewhere from Glu2 to unknown
        let variant = parse_variant("NP_001005484.1:p.Met1_(Glu2_?)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_001005484.1:p.Met1_(Glu2_?)del");

        // Both uncertain: (?_Met1)_(Glu2_?)del
        let variant = parse_variant("NP_001005484.1:p.(?_Met1)_(Glu2_?)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(
            format!("{}", variant),
            "NP_001005484.1:p.(?_Met1)_(Glu2_?)del"
        );
    }

    #[test]
    fn test_parse_protein_nested_uncertain_positions() {
        // Deeply nested uncertain positions: Met(?_1)_Glu2(?)
        // This is AminoAcid with position range followed by AminoAcid with uncertain position
        let variant = parse_variant("NP_000088.1:p.Met(?_1)_Glu2(?)").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        // Format output may vary based on internal representation

        // AminoAcid with just (?) - fully uncertain position
        let variant = parse_variant("NP_000088.1:p.Ter(?)").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));

        // Combined with range: His1811_Ter(?)
        let variant = parse_variant("NP_000088.1:p.His1811_Ter(?)").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));

        // AminoAcid with position range (?_N)
        let variant = parse_variant("NP_000088.1:p.Met(?_5)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));

        // AminoAcid with position range (N_?)
        let variant = parse_variant("NP_000088.1:p.Glu(3_?)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
    }

    #[test]
    fn test_parse_conversion() {
        // Gene conversion
        let variant = parse_variant("NM_000088.3:c.100_200conNM_000089.1:c.50_150").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(
            format!("{}", variant),
            "NM_000088.3:c.100_200conNM_000089.1:c.50_150"
        );
    }

    #[test]
    fn test_parse_genomic_conversion() {
        // Genomic conversion
        let variant = parse_variant("NC_000001.11:g.12345_12400conNC_000002.12:g.100_155").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(
            format!("{}", variant),
            "NC_000001.11:g.12345_12400conNC_000002.12:g.100_155"
        );
    }

    #[test]
    fn test_parse_position_based_conversion() {
        // Position-based conversion (same chromosome)
        let variant =
            parse_variant("NC_000017.11:g.42522624_42522669con42536337_42536382").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(
            format!("{}", variant),
            "NC_000017.11:g.42522624_42522669con42536337_42536382"
        );
    }

    #[test]
    fn test_parse_repeat_exact() {
        // Exact repeat count
        let variant = parse_variant("NM_000088.3:c.100CAG[12]").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.100CAG[12]");
    }

    #[test]
    fn test_parse_repeat_range() {
        // Range repeat count
        let variant = parse_variant("NM_000088.3:c.100CAG[10_15]").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.100CAG[10_15]");
    }

    #[test]
    fn test_parse_repeat_min_uncertain() {
        // Minimum with uncertain maximum
        let variant = parse_variant("NM_000088.3:c.100CAG[10_?]").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.100CAG[10_?]");
    }

    #[test]
    fn test_parse_repeat_max_uncertain() {
        // Maximum with uncertain minimum
        let variant = parse_variant("NM_000088.3:c.100CAG[?_20]").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.100CAG[?_20]");
    }

    #[test]
    fn test_parse_repeat_unknown() {
        // Fully unknown count
        let variant = parse_variant("NM_000088.3:c.100CAG[?]").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.100CAG[?]");
    }

    #[test]
    fn test_parse_unknown_end_position() {
        // Unknown end position: c.1_?del
        let variant = parse_variant("NM_000088.3:c.1_?del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.1_?del");
    }

    #[test]
    fn test_parse_unknown_start_position() {
        // Unknown start position: c.?_100del
        let variant = parse_variant("NM_000088.3:c.?_100del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.?_100del");
    }

    #[test]
    fn test_parse_genomic_unknown_position() {
        // Genomic with unknown positions
        let variant = parse_variant("NC_000001.11:g.100_?del").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(format!("{}", variant), "NC_000001.11:g.100_?del");

        let variant = parse_variant("NC_000001.11:g.?_200del").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(format!("{}", variant), "NC_000001.11:g.?_200del");
    }

    #[test]
    fn test_parse_transcript_unknown_position() {
        // Transcript with unknown positions
        let variant = parse_variant("NR_000001.1:n.100_?del").unwrap();
        assert!(matches!(variant, HgvsVariant::Tx(_)));
        assert_eq!(format!("{}", variant), "NR_000001.1:n.100_?del");
    }

    #[test]
    fn test_parse_whole_cds_identity() {
        // Whole-CDS identity: c.=
        let variant = parse_variant("NM_000088.3:c.=").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.=");
    }

    #[test]
    fn test_parse_position_specific_identity() {
        // Position-specific identity: c.100=
        let variant = parse_variant("NM_000088.3:c.100=").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.100=");
    }

    #[test]
    fn test_parse_whole_protein_unknown() {
        // Whole-protein unknown: p.?
        let variant = parse_variant("NP_000079.2:p.?").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.?");

        // Predicted whole-protein unknown: p.(?)
        let variant = parse_variant("NP_000079.2:p.(?)").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.(?)");
    }

    #[test]
    fn test_parse_position_specific_unknown() {
        // Position-specific unknown: p.Met1?
        let variant = parse_variant("NP_000079.2:p.Met1?").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.Met1?");
    }

    #[test]
    fn test_parse_predicted_protein_substitution() {
        // Predicted protein substitution: p.(Arg248Gln)
        let variant = parse_variant("NP_000079.2:p.(Arg248Gln)").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.(Arg248Gln)");
    }

    #[test]
    fn test_parse_predicted_protein_frameshift() {
        // Predicted protein frameshift: p.(Lys23fs)
        let variant = parse_variant("NP_000079.2:p.(Lys23fs)").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.(Lys23fs)");
    }

    #[test]
    fn test_parse_predicted_protein_deletion() {
        // Predicted protein deletion: p.(Val600del)
        let variant = parse_variant("NP_000079.2:p.(Val600del)").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        assert_eq!(format!("{}", variant), "NP_000079.2:p.(Val600del)");
    }

    #[test]
    fn test_parse_complex_cds_interval_exon_deletion() {
        // Complex interval with intronic breakpoints: c.(4185+1_4186-1)_(4357+1_4358-1)del
        let variant = parse_variant("NM_007294.4:c.(4185+1_4186-1)_(4357+1_4358-1)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(
            format!("{}", variant),
            "NM_007294.4:c.(4185+1_4186-1)_(4357+1_4358-1)del"
        );
    }

    #[test]
    fn test_parse_complex_cds_interval_with_unknown() {
        // Complex interval with unknown start: c.(?_-1)_(96+1_97-1)del
        let variant = parse_variant("NM_000546.6:c.(?_-1)_(96+1_97-1)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(
            format!("{}", variant),
            "NM_000546.6:c.(?_-1)_(96+1_97-1)del"
        );
    }

    #[test]
    fn test_parse_complex_genome_interval_uncertain_boundaries() {
        // Complex genomic interval with unknown boundaries: g.(?_43044294)_(43125364_?)del
        let variant = parse_variant("NC_000017.11:g.(?_43044294)_(43125364_?)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(
            format!("{}", variant),
            "NC_000017.11:g.(?_43044294)_(43125364_?)del"
        );
    }

    #[test]
    fn test_parse_simple_uncertain_boundary() {
        // Simple uncertain boundary: c.(4_6)_246del
        // The start position is uncertain (somewhere between 4 and 6)
        let variant = parse_variant("NM_000001.1:c.(4_6)_246del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000001.1:c.(4_6)_246del");

        // Also test with uncertain end position
        let variant = parse_variant("NM_000001.1:c.100_(200_250)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000001.1:c.100_(200_250)del");

        // Both positions uncertain
        let variant = parse_variant("NM_000001.1:c.(4_6)_(200_250)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000001.1:c.(4_6)_(200_250)del");
    }

    // Allele parsing tests

    #[test]
    fn test_parse_cis_allele() {
        // Cis allele: [var1;var2]
        let variant = parse_variant("[NM_000088.3:c.100A>G;NM_000088.3:c.200C>T]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Cis);
            assert_eq!(allele.variants.len(), 2);
        }
        assert_eq!(
            format!("{}", variant),
            "[NM_000088.3:c.100A>G;NM_000088.3:c.200C>T]"
        );
    }

    #[test]
    fn test_parse_trans_allele() {
        // Trans allele: [var1];[var2]
        let variant = parse_variant("[NM_000088.3:c.100A>G];[NM_000088.3:c.200C>T]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Trans);
            assert_eq!(allele.variants.len(), 2);
        }
        assert_eq!(
            format!("{}", variant),
            "[NM_000088.3:c.100A>G];[NM_000088.3:c.200C>T]"
        );
    }

    #[test]
    fn test_parse_mosaic_allele() {
        // Mosaic: var1/var2
        let variant = parse_variant("NM_000088.3:c.100A>G/NM_000088.3:c.200C>T").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Mosaic);
            assert_eq!(allele.variants.len(), 2);
        }
        assert_eq!(
            format!("{}", variant),
            "NM_000088.3:c.100A>G/NM_000088.3:c.200C>T"
        );
    }

    #[test]
    fn test_parse_mosaic_with_reference() {
        // Mosaic with reference allele shorthand: var/=
        let variant = parse_variant("NM_000088.3:c.123A>G/=").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Mosaic);
            assert_eq!(allele.variants.len(), 2);
            // First variant should be the mutation
            assert!(matches!(&allele.variants[0], HgvsVariant::Cds(_)));
            // Second variant should be identity (reference)
            if let HgvsVariant::Cds(cds) = &allele.variants[1] {
                assert!(cds
                    .loc_edit
                    .edit
                    .inner()
                    .map(|e| e.is_whole_entity_identity())
                    .unwrap_or(false));
            } else {
                panic!("Expected CDS variant");
            }
        }
        // Roundtrip: c.123A>G/c.=
        assert_eq!(
            format!("{}", variant),
            "NM_000088.3:c.123A>G/NM_000088.3:c.="
        );
    }

    #[test]
    fn test_parse_chimeric_with_reference() {
        // Chimeric with reference allele shorthand: var//=
        let variant = parse_variant("NM_000088.3:c.456C>T//=").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Chimeric);
            assert_eq!(allele.variants.len(), 2);
        }
        assert_eq!(
            format!("{}", variant),
            "NM_000088.3:c.456C>T//NM_000088.3:c.="
        );
    }

    #[test]
    fn test_parse_mosaic_reference_first_fails() {
        // Cannot use = as first variant
        let result = parse_variant("=/NM_000088.3:c.123A>G");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_chimeric_allele() {
        // Chimeric: var1//var2
        let variant = parse_variant("NM_000088.3:c.100A>G//NM_000088.3:c.200C>T").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Chimeric);
            assert_eq!(allele.variants.len(), 2);
        }
        assert_eq!(
            format!("{}", variant),
            "NM_000088.3:c.100A>G//NM_000088.3:c.200C>T"
        );
    }

    #[test]
    fn test_parse_unknown_phase_allele() {
        // Unknown phase allele: [var1(;)var2]
        let variant = parse_variant("[NM_000088.3:c.100A>G(;)NM_000088.3:c.200C>T]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Unknown);
            assert_eq!(allele.variants.len(), 2);
        }
        assert_eq!(
            format!("{}", variant),
            "[NM_000088.3:c.100A>G(;)NM_000088.3:c.200C>T]"
        );
    }

    #[test]
    fn test_parse_unknown_phase_allele_multiple() {
        // Unknown phase allele with 3 variants
        let variant =
            parse_variant("[NM_000088.3:c.100A>G(;)NM_000088.3:c.200C>T(;)NM_000088.3:c.300G>A]")
                .unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Unknown);
            assert_eq!(allele.variants.len(), 3);
        }
    }

    #[test]
    fn test_parse_cis_allele_shorthand() {
        // Cis allele shorthand: NM_000088.3:c.[145C>T;147C>G]
        let variant = parse_variant("NM_000088.3:c.[145C>T;147C>G]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Cis);
            assert_eq!(allele.variants.len(), 2);
            // Both variants should have the same accession
            for v in &allele.variants {
                assert_eq!(
                    v.accession().expect("Expected accession").full(),
                    "NM_000088.3"
                );
            }
        }
        // Output expands to full format
        assert_eq!(
            format!("{}", variant),
            "[NM_000088.3:c.145C>T;NM_000088.3:c.147C>G]"
        );
    }

    #[test]
    fn test_parse_unknown_phase_allele_shorthand() {
        // Unknown phase allele shorthand: NM_000088.3:c.[145C>T(;)147C>G]
        let variant = parse_variant("NM_000088.3:c.[145C>T(;)147C>G]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Unknown);
            assert_eq!(allele.variants.len(), 2);
        }
        assert_eq!(
            format!("{}", variant),
            "[NM_000088.3:c.145C>T(;)NM_000088.3:c.147C>G]"
        );
    }

    #[test]
    fn test_parse_mixed_phase_allele_shorthand() {
        // Mixed phase allele: some cis (;) and some unknown ((;))
        // Pattern: [123A>G;456C>T(;)789G>A] means 123 and 456 are cis, unknown phase with 789
        let variant = parse_variant("NM_000088.3:c.[123A>G;456C>T(;)789G>A]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            // Overall phase is Unknown when there's any unknown phase
            assert_eq!(allele.phase, AllelePhase::Unknown);
            assert_eq!(allele.variants.len(), 3);
        }
    }

    #[test]
    fn test_parse_mixed_phase_allele_shorthand_complex() {
        // More complex mixed phase: [A;B;C(;)D;E]
        let variant =
            parse_variant("NM_000088.3:c.[100A>G;200C>T;300G>A(;)400del;500dup]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Unknown);
            assert_eq!(allele.variants.len(), 5);
        }
    }

    #[test]
    fn test_parse_allele_with_interval_only_variant() {
        // Allele where first variant is interval-only (no edit) - implies identity/reference
        // ClinVar pattern: g.[interval;interval_with_edit]
        let variant = parse_variant("NC_000023.10:g.[100_200;300_400dup]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Cis);
            assert_eq!(allele.variants.len(), 2);
            // First variant should be identity
            if let HgvsVariant::Genome(v) = &allele.variants[0] {
                assert!(matches!(
                    v.loc_edit.edit,
                    crate::hgvs::uncertainty::Mu::Certain(
                        crate::hgvs::edit::NaEdit::Identity { .. }
                    )
                ));
            }
        }

        // Test with uncertain intervals (actual ClinVar pattern)
        let variant = parse_variant(
            "NC_000023.10:g.[(?_29619835)_(29843303_?);(?_30646799)_(30848980_?)dup]",
        )
        .unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.variants.len(), 2);
        }
    }

    #[test]
    fn test_parse_cis_allele_multiple_variants() {
        // Cis allele with 3 variants
        let variant =
            parse_variant("[NM_000088.3:c.100A>G;NM_000088.3:c.200C>T;NM_000088.3:c.300G>A]")
                .unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Cis);
            assert_eq!(allele.variants.len(), 3);
        }
    }

    #[test]
    fn test_parse_cis_allele_mixed_types() {
        // Cis allele with CDS and protein variants
        let variant = parse_variant("[NM_000088.3:c.100A>G;NP_000079.2:p.Val600Glu]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Cis);
            assert_eq!(allele.variants.len(), 2);
            assert!(matches!(allele.variants[0], HgvsVariant::Cds(_)));
            assert!(matches!(allele.variants[1], HgvsVariant::Protein(_)));
        }
    }

    // ClinVar-style allele tests (variant type prefix inside brackets)

    #[test]
    fn test_parse_cis_allele_clinvar_format() {
        // ClinVar uses format: NM_xxx:[c.var1;c.var2] instead of NM_xxx:c.[var1;var2]
        let variant = parse_variant("NM_006876.2:[c.1168A>G;c.1217C>T]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Cis);
            assert_eq!(allele.variants.len(), 2);
            // Both variants should be CDS variants with the shared accession
            assert!(
                matches!(&allele.variants[0], HgvsVariant::Cds(v) if v.accession.to_string() == "NM_006876.2")
            );
            assert!(
                matches!(&allele.variants[1], HgvsVariant::Cds(v) if v.accession.to_string() == "NM_006876.2")
            );
        }
    }

    #[test]
    fn test_parse_cis_allele_clinvar_format_genomic() {
        // ClinVar genomic format: NC_xxx:[g.var1;g.var2]
        let variant = parse_variant("NC_000017.11:[g.43094927del;g.43095845dup]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Cis);
            assert_eq!(allele.variants.len(), 2);
            assert!(matches!(&allele.variants[0], HgvsVariant::Genome(_)));
            assert!(matches!(&allele.variants[1], HgvsVariant::Genome(_)));
        }
    }

    #[test]
    fn test_parse_unknown_phase_allele_clinvar_format() {
        // ClinVar unknown phase format: NM_xxx:[c.var1(;)c.var2]
        let variant = parse_variant("NM_000088.3:[c.145C>T(;)c.147C>G]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Unknown);
            assert_eq!(allele.variants.len(), 2);
        }
    }

    // === SVD-WG006: Circular DNA (o.) tests ===

    #[test]
    fn test_parse_circular_substitution() {
        // Using NC_ (chromosome) for circular plasmid (common pattern)
        let variant = parse_variant("NC_001416.1:o.100A>G").unwrap();
        assert!(matches!(variant, HgvsVariant::Circular(_)));
        assert_eq!(format!("{}", variant), "NC_001416.1:o.100A>G");
    }

    #[test]
    fn test_parse_circular_duplication() {
        // Normal duplication on circular DNA
        let variant = parse_variant("NC_001416.1:o.100_200dup").unwrap();
        assert!(matches!(variant, HgvsVariant::Circular(_)));
        assert_eq!(format!("{}", variant), "NC_001416.1:o.100_200dup");
    }

    #[test]
    fn test_parse_circular_deletion() {
        // Deletion on circular DNA
        let variant = parse_variant("NC_001416.1:o.100_200del").unwrap();
        assert!(matches!(variant, HgvsVariant::Circular(_)));
        assert_eq!(format!("{}", variant), "NC_001416.1:o.100_200del");
    }

    #[test]
    fn test_parse_circular_insertion() {
        let variant = parse_variant("NC_001416.1:o.100_101insATG").unwrap();
        assert!(matches!(variant, HgvsVariant::Circular(_)));
        assert_eq!(format!("{}", variant), "NC_001416.1:o.100_101insATG");
    }

    // === SVD-WG007: RNA Fusion (::) tests ===

    #[test]
    fn test_parse_rna_fusion_basic() {
        // Basic fusion: NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924
        let variant = parse_variant("NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924").unwrap();
        assert!(matches!(variant, HgvsVariant::RnaFusion(_)));
        if let HgvsVariant::RnaFusion(fusion) = &variant {
            assert_eq!(fusion.five_prime.accession.full(), "NM_152263.2");
            assert_eq!(fusion.three_prime.accession.full(), "NM_002609.3");
        }
        assert_eq!(
            format!("{}", variant),
            "NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924"
        );
    }

    #[test]
    fn test_parse_rna_fusion_simple_positions() {
        let variant = parse_variant("NM_000546.5:r.1_200::NM_000245.3:r.100_500").unwrap();
        assert!(matches!(variant, HgvsVariant::RnaFusion(_)));
        assert_eq!(
            format!("{}", variant),
            "NM_000546.5:r.1_200::NM_000245.3:r.100_500"
        );
    }

    #[test]
    fn test_parse_rna_fusion_single_positions() {
        // Fusion with single position on each side
        let variant = parse_variant("NM_000001.1:r.100::NM_000002.1:r.200").unwrap();
        assert!(matches!(variant, HgvsVariant::RnaFusion(_)));
        assert_eq!(
            format!("{}", variant),
            "NM_000001.1:r.100::NM_000002.1:r.200"
        );
    }

    #[test]
    fn test_parse_rna_fusion_epcam_msh2() {
        // Real-world example: EPCAM::MSH2 fusion (Lynch syndrome)
        let variant = parse_variant("NM_002354.2:r.1_3469::NM_000251.2:r.3397_6871").unwrap();
        assert!(matches!(variant, HgvsVariant::RnaFusion(_)));
        if let HgvsVariant::RnaFusion(fusion) = &variant {
            assert_eq!(&*fusion.five_prime.accession.prefix, "NM");
            assert_eq!(&*fusion.three_prime.accession.prefix, "NM");
        }
    }

    // ============== Extended Pattern Coverage Tests ==============

    #[test]
    fn test_parse_nested_repeat_notation() {
        // Genotype notation: A[6][1] - 6 repeats on one allele, 1 on another
        let variant = parse_variant("NM_000088.3:c.100A[6][1]").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));

        // Three repeat counts
        let variant = parse_variant("NM_000088.3:c.100[4][5][6]").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
    }

    #[test]
    fn test_parse_complex_insertion_with_parts() {
        // Complex insertion with repeat
        let variant = parse_variant("NM_000088.3:c.100_101ins[A[10];T]").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));

        // Complex insertion with position range reference
        let variant = parse_variant("NM_000088.3:c.419_420ins[T;401_419]").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
    }

    #[test]
    fn test_parse_uniprot_protein_variant() {
        // 6-character UniProt accession
        let variant = parse_variant("P54802:p.Phe48Leu").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        if let HgvsVariant::Protein(prot) = &variant {
            assert_eq!(&*prot.accession.prefix, "P");
            assert_eq!(&*prot.accession.number, "54802");
        }

        // Another UniProt accession
        let variant = parse_variant("Q8TAM1:p.Arg34Pro").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
    }

    #[test]
    fn test_parse_uniprot_10char_accession() {
        // 10-character UniProt accession (first char is prefix, rest is number)
        let variant = parse_variant("A0A024R1R8:p.Val100Glu").unwrap();
        assert!(matches!(variant, HgvsVariant::Protein(_)));
        if let HgvsVariant::Protein(prot) = &variant {
            assert_eq!(&*prot.accession.prefix, "A");
            assert_eq!(&*prot.accession.number, "0A024R1R8");
            // Full accession should display correctly
            assert_eq!(prot.accession.full(), "A0A024R1R8");
        }
    }

    #[test]
    fn test_parse_inversion_with_length() {
        // Inversion with explicit length: inv3
        let variant = parse_variant("NM_000088.3:c.274_276inv3").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
    }

    #[test]
    fn test_parse_inversion_with_sequence() {
        // Inversion with explicit sequence: invATG
        let variant = parse_variant("NM_000088.3:c.274_276invATG").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
    }

    #[test]
    fn test_parse_uncertain_duplication_range() {
        // Uncertain duplication range: dup(731_741)
        let variant = parse_variant("NM_000088.3:c.722_742dup(731_741)").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));

        // Completely uncertain dup: dup?
        let variant = parse_variant("NM_000088.3:c.100dup?").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));

        // Uncertain dup with parens: dup(?)
        let variant = parse_variant("NM_000088.3:c.100dup(?)").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
    }

    #[test]
    fn test_parse_assembly_chromosome_notation() {
        // Assembly/chromosome notation
        let variant = parse_variant("GRCh38(chr1):g.100A>G").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        if let HgvsVariant::Genome(gen) = &variant {
            assert!(gen.accession.is_assembly_ref());
            assert_eq!(gen.accession.assembly.as_deref(), Some("GRCh38"));
            assert_eq!(gen.accession.chromosome.as_deref(), Some("chr1"));
        }

        // GRCh37 notation
        let variant = parse_variant("GRCh37(chrX):g.15000del").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
    }

    #[test]
    fn test_parse_circular_dna_notation() {
        // Circular DNA variant (o. prefix)
        let variant = parse_variant("NC_012920.1:o.100A>G").unwrap();
        assert!(matches!(variant, HgvsVariant::Circular(_)));
    }

    #[test]
    fn test_parse_complex_uncertain_positions() {
        // Already tested but verify roundtrip
        let variant = parse_variant("NM_000088.3:c.(?_-1)_328+?del").unwrap();
        assert_eq!(format!("{}", variant), "NM_000088.3:c.(?_-1)_328+?del");

        // Another complex pattern
        let variant = parse_variant("NM_000088.3:c.(?_1)_5074+?del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
    }

    #[test]
    fn test_parse_unknown_position_patterns() {
        // Unknown position with duplication: c.?dup
        let variant = parse_variant("NM_001412270.1:c.?dup").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_001412270.1:c.?dup");

        // Unknown position with offset at start: c.?-232_4484+?del
        let variant = parse_variant("NM_007294.3:c.?-232_4484+?del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_007294.3:c.?-232_4484+?del");

        // Same pattern with LRG accession
        let variant = parse_variant("LRG_292t1:c.?-232_4484+?del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "LRG_292t1:c.?-232_4484+?del");
    }

    #[test]
    fn test_parse_delins_external_reference() {
        // Delins with external GenBank reference
        let variant =
            parse_variant("NC_000008.11:g.86587460_86650711delins[KY923049.1:g.1_466]").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(
            format!("{}", variant),
            "NC_000008.11:g.86587460_86650711delins[KY923049.1:g.1_466]"
        );

        // CDS delins with external reference and intronic positions
        let variant =
            parse_variant("NM_019098.4:c.904-2824_1782-8208delins[KY923049.1:g.1_466]").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(
            format!("{}", variant),
            "NM_019098.4:c.904-2824_1782-8208delins[KY923049.1:g.1_466]"
        );
    }

    #[test]
    fn test_parse_complex_delins_array() {
        // Complex delins with sequence, positions, and inversions
        let variant = parse_variant("NC_000007.14:g.45043702_46521017delins[AGAAGGAAATTT;45310743_46521014;45043709_45310738inv]").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(format!("{}", variant), "NC_000007.14:g.45043702_46521017delins[AGAAGGAAATTT;45310743_46521014;45043709_45310738inv]");

        // Simple complex delins with just inversion
        let variant =
            parse_variant("NC_000016.9:g.78179358_78219143delins[78185355_78199419inv]").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(
            format!("{}", variant),
            "NC_000016.9:g.78179358_78219143delins[78185355_78199419inv]"
        );
    }

    #[test]
    fn test_parse_delins_sequence_repeat_range() {
        // Delins with sequence repeat range: delinsAAAGG[400_2000]
        let variant =
            parse_variant("NC_000004.12:g.39348425_39348479delinsAAAGG[400_2000]").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(
            format!("{}", variant),
            "NC_000004.12:g.39348425_39348479delinsAAAGG[400_2000]"
        );
    }

    #[test]
    fn test_parse_complex_insertion_with_external() {
        // Complex insertion with literal and external reference: ins[TCTT;KT192064.1:1_310]
        let variant = parse_variant("NM_017635.5:c.438_439ins[TCTT;KT192064.1:1_310]").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(
            format!("{}", variant),
            "NM_017635.5:c.438_439ins[TCTT;KT192064.1:1_310]"
        );
    }

    #[test]
    fn test_parse_telomere_qter_patterns() {
        // Telomere/qter notation - complex structural variants
        // NC_000009.12:g.12891379_qterdelins[T;NC_000020.11:g.47204889_qterinv]
        let variant =
            parse_variant("NC_000009.12:g.12891379_qterdelins[T;NC_000020.11:g.47204889_qterinv]")
                .unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(
            format!("{}", variant),
            "NC_000009.12:g.12891379_qterdelins[T;NC_000020.11:g.47204889_qterinv]"
        );

        // NC_000013.10:g.114819939_qterdelins[96729864_114814234inv;96735632_104289803]
        let variant = parse_variant(
            "NC_000013.10:g.114819939_qterdelins[96729864_114814234inv;96735632_104289803]",
        )
        .unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(
            format!("{}", variant),
            "NC_000013.10:g.114819939_qterdelins[96729864_114814234inv;96735632_104289803]"
        );
    }

    #[test]
    fn test_parse_predicted_substitution_in_parens() {
        // Predicted substitution in parens: c.(9740C>A)
        // This is a predicted variant where the whole edit is wrapped in parentheses
        let variant = parse_variant("NM_002016.2:c.(9740C>A)").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_002016.2:c.9740(C>A)");

        let variant = parse_variant("NM_006767.4:c.(742G>A)").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_006767.4:c.742(G>A)");
    }

    #[test]
    fn test_parse_uncertain_intron_exon_ranges() {
        // Uncertain intron/exon ranges - these are complex patterns
        // Pattern 1: c.(2727+1_2728-1)(2858+1_2859-1) - non-standard (no underscore between ranges)
        // Patterns 2 & 3 are position-only (no edit type) which is non-standard HGVS
    }

    #[test]
    fn test_parse_reverse_order_intervals() {
        // Reverse order intervals - where start > end (used for structural variants)
        // Feature 15: Reverse order intervals
        let variant =
            parse_variant("NC_000011.10:g.5238138_5153222insTATTT").expect("should parse");
        assert_eq!(
            variant.to_string(),
            "NC_000011.10:g.5238138_5153222insTATTT"
        );

        let variant = parse_variant("NC_000012.11:g.110593351_110576466dup").expect("should parse");
        assert_eq!(variant.to_string(), "NC_000012.11:g.110593351_110576466dup");

        // Also from Feature 16
        let variant = parse_variant("NC_000016.9:g.23634775_23621090dup").expect("should parse");
        assert_eq!(variant.to_string(), "NC_000016.9:g.23634775_23621090dup");
    }

    #[test]
    fn test_parse_position_only_cds() {
        // Feature 6: Position-Only CDS - positions without edit type
        let variant = parse_variant("NM_173651.4:c.5238_5240").expect("should parse");
        assert_eq!(variant.to_string(), "NM_173651.4:c.5238_5240");

        let variant = parse_variant("NM_134428.2:c.1486_1487").expect("should parse");
        assert_eq!(variant.to_string(), "NM_134428.2:c.1486_1487");

        // 3'UTR position-only
        let variant = parse_variant("NM_007375.3:c.*697").expect("should parse");
        assert_eq!(variant.to_string(), "NM_007375.3:c.*697");
    }

    #[test]
    fn test_parse_unknown_position_insertion() {
        // Unknown position insertion patterns
        // g.?_?ins(start_boundary)_(end_boundary)
        let _result = parse_variant("LRG_308:g.?_?ins(23632682_23625413)_(23625324_23619334)");

        let _result = parse_variant("NG_007406.1:g.?_?ins(23632682_23625413)_(23625324_23619334)");
    }

    #[test]
    fn test_parse_tx_downstream_positions() {
        // Feature 16: n.* downstream positions (after transcript end)
        let variant = parse_variant("NR_033294.1:n.*5C>G").expect("should parse");
        assert_eq!(variant.to_string(), "NR_033294.1:n.*5C>G");

        // With offset
        let variant = parse_variant("NR_033294.1:n.*5+10C>G").expect("should parse");
        assert_eq!(variant.to_string(), "NR_033294.1:n.*5+10C>G");
    }
}
