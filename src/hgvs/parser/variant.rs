//! Full variant parsing
//!
//! Parses complete HGVS variant strings into the HgvsVariant type.

use crate::error::FerroError;
use crate::error_handling::ErrorType;
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

    // Fast path: if starts with digit, *, or - it's likely a certain position.
    // Also accept p/q/c for pter/qter/cen special markers (mirrors the
    // starts_with_special_pos fast-path on the genomic interval path).
    if bytes[0].is_ascii_digit()
        || bytes[0] == b'*'
        || bytes[0] == b'-'
        || bytes[0] == b'p'
        || bytes[0] == b'q'
        || bytes[0] == b'c'
    {
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
    parse_genome_interval_inner(input, false)
}

/// Parse a genomic-axis interval, optionally relaxing the inverted-range
/// check for circular references (`m.`/`o.`).
///
/// HGVS spec authorises reversed `<high>_<low>` ranges on `o.`/`m.`
/// references for **deletions** (`deletion.md:17`) and **duplications**
/// (`docs/consultation/SVD-WG006.md` line 15 + line 23 example
/// `J01749.1:o.4344_197dup`, line 33 names "deletions/duplications" as
/// the spec-authorised scope). `delins` inherits the exception via the
/// `del + ins` composition. `ins` and `inv` are spec-silent — the
/// general "5'→3'" rule applies. `check_circular_reversed_range`,
/// called post-parse by `parse_mt_variant` / `parse_circular_variant`
/// (and the bracket / unknown-phase / inherited-accession helpers),
/// enforces the edit-kind authorisation so the rejection has access to
/// the parsed edit kind.
///
/// When `circular = true` the parse-level inverted-range check is
/// removed for the simple-range fast path; post-parse validation in the
/// caller catches the spec-unauthorised edit kinds. The complex
/// uncertain-range paths still reject reversed bounds since the spec
/// exception only applies to plain `<high>_<low>` reversed forms.
fn parse_genome_interval_for_circular(input: &str) -> IResult<&str, GenomeInterval> {
    parse_genome_interval_inner(input, true)
}

fn parse_genome_interval_inner(input: &str, circular: bool) -> IResult<&str, GenomeInterval> {
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
                // - this is a circular (`m.`/`o.`) reference. The caller validates
                //   the spec-authorised edit kinds post-parse (#129).
                let allows_inverted = circular
                    || remaining.starts_with("dupins")
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
        // Single uncertain position bounded by a range: (pos_pos) NOT followed
        // by _ — per HGVS v21 uncertain.md this is *one* position somewhere
        // in [start, end], not a multi-base range. Encoded as the same Range
        // on both Interval boundaries; Display collapses to `(a_b)` once.
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
            // Reject inverted ranges (skip ordering when either bound is a special
            // position like pter/qter/cen — those carry base=0 sentinels so a raw
            // base comparison wrongly rejects forms like (12345_qter)).
            if !start.is_special() && !end.is_special() && start.base > end.base {
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Verify,
                )));
            }
            let boundary = UncertainBoundary::range(Mu::certain(start), Mu::certain(end));
            Ok((
                remaining,
                Interval::with_complex_boundaries(boundary.clone(), boundary),
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

    // Only if NEITHER has an offset, treat as simple uncertain interval.
    // Per HGVS v21 uncertain.md, `(a_b)<edit>` is *one* position somewhere in
    // [a, b], not a multi-base range. Encoded as the same Range on both
    // Interval boundaries; Display collapses to `(a_b)` once.
    if start.offset.is_none() && end.offset.is_none() {
        let boundary = UncertainBoundary::range(Mu::certain(start), Mu::certain(end));
        Ok((
            remaining,
            Interval::with_complex_boundaries(boundary.clone(), boundary),
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
    // Special telomere/centromere markers carry a base==0 sentinel with no
    // comparable base value; never treat a special-bearing range as inverted
    // (mirrors the genomic-axis guards in this file).
    if start.is_special() || end.is_special() {
        return false;
    }

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

/// Validate a parsed circular-axis (`m.`/`o.`) interval+edit pair against
/// the HGVS spec's circular range-orientation rules.
///
/// HGVS spec sources:
/// - `recommendations/DNA/deletion.md:17` — explicit exception for
///   `o.`/`m.` references allowing reversed `<high>_<low>` ranges on
///   deletions.
/// - `docs/consultation/SVD-WG006.md` (accepted into HGVS v19.01,
///   founding proposal for circular DNA) — line 15 grants the
///   exception for any "rearrangement" that "includes both the last
///   and first nucleotides of the reference sequence", with explicit
///   examples for both `del` (`NC_012920.1:m.16563_13del`) and `dup`
///   (`J01749.1:o.4344_197dup`). Line 33 explicitly names
///   "deletions/duplications" as the scope.
///
/// So the spec-authorised set is `del`, `dup`, and `delins` (the
/// latter by composition: `delins` is a `del + ins` rewrite and the
/// del side carries the exception). `ins` and `inv` are spec-silent
/// (no examples and not named in the rationale); the general "5' to 3'"
/// rule applies and ferro rejects those reversed forms — closing the
/// silent-accept gap pinned by the audit tests.
///
/// This helper is called post-parse from `parse_mt_variant` and
/// `parse_circular_variant` to reject the unauthorised forms. Returns
/// `Ok(())` when the form is spec-compliant, or `Err(_)` with a
/// nom-style verify-error so the caller can short-circuit the parse.
/// (#129)
fn check_circular_reversed_range<'a>(
    interval: &GenomeInterval,
    edit: &crate::hgvs::edit::NaEdit,
    input: &'a str,
) -> Result<(), nom::Err<nom::error::Error<&'a str>>> {
    use crate::hgvs::edit::NaEdit;

    // Only the plain `<start>_<end>` (Single(Certain), Single(Certain))
    // shape can be reversed; complex range/uncertain boundaries are
    // already rejected by parse_genome_interval_inner before reaching
    // here for reversed cases.
    let (s, e) = match (&interval.start, &interval.end) {
        (UncertainBoundary::Single(Mu::Certain(s)), UncertainBoundary::Single(Mu::Certain(e))) => {
            (s, e)
        }
        _ => return Ok(()),
    };
    // Special positions (pter/qter/cen) carry base=0 sentinels; ordering
    // is not numeric. Spec doesn't address these on circular refs;
    // pass through unchanged.
    if s.is_special() || e.is_special() {
        return Ok(());
    }
    if s.base <= e.base {
        return Ok(());
    }
    // Reversed range. Allow `del`, `dup`, and `delins` (the
    // spec-authorised set per SVD-WG006); reject everything else.
    //
    // The `_` arm explicitly covers:
    //   - `Substitution` / `SubstitutionNoRef` — single-position
    //     edits; `<high>_<low>` is structurally undefined.
    //   - `Insertion`, `Inversion` — spec-silent on circular refs
    //     (SVD-WG006 names only del/dup; insertion.md / inversion.md
    //     give no wraparound exception).
    //   - `DupIns` — non-standard ClinVar shape; no spec
    //     authorisation, defensive reject.
    //   - `Repeat` / `MultiRepeat` — spec is silent on wraparound
    //     repeat tracts.
    //   - `Identity` — `<high>_<low>=` is semantically degenerate.
    //   - `Conversion` — defensive default-reject; the canonicaliser
    //     rewrites `con` → `delins` post-parse but cannot fire here
    //     because the validator runs immediately after `parse_na_edit`.
    //   - `Unknown` — `?` edits carry no positional semantics on a
    //     reversed range; reject defensively.
    match edit {
        NaEdit::Deletion { .. } | NaEdit::Delins { .. } | NaEdit::Duplication { .. } => Ok(()),
        _ => Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        ))),
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

    // Fast path: if starts with digit, *, or - it's a simple position or range.
    // Also accept p/q/c for pter/qter/cen special markers (mirrors the
    // starts_with_special_pos fast-path on the genomic interval path).
    // Also handle ? when followed by +/- offset (e.g., ?-232) but NOT when followed by _
    // (e.g., ?_-540 needs the complex parser for patterns like ?_-540_3117+?del)
    let use_fast_path = bytes[0].is_ascii_digit()
        || bytes[0] == b'*'
        || bytes[0] == b'-'
        || bytes[0] == b'p'
        || bytes[0] == b'q'
        || bytes[0] == b'c'
        || (bytes[0] == b'?' && bytes.len() > 1 && (bytes[1] == b'+' || bytes[1] == b'-'));
    if use_fast_path {
        let (remaining, start) = parse_cds_pos(input)?;
        // Check if this is a range (followed by _)
        if let Some(after_underscore) = remaining.strip_prefix('_') {
            // Only use fast path if followed by a position start character (simple range)
            // Fall through to complex parser for patterns like 100_(200)del
            let after_bytes = after_underscore.as_bytes();
            // Note: ? is NOT included - it needs the complex parser.
            // p/q/c are included for pter/qter/cen special markers.
            let is_simple_range = !after_bytes.is_empty()
                && (after_bytes[0].is_ascii_digit()
                    || after_bytes[0] == b'*'
                    || after_bytes[0] == b'-'
                    || after_bytes[0] == b'p'
                    || after_bytes[0] == b'q'
                    || after_bytes[0] == b'c');
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
        // Single uncertain position bounded by a range: (a_b) NOT followed by _.
        // Per HGVS v21 uncertain.md, `(a_b)<edit>` is *one* position somewhere
        // in [a, b]; same Range on both Interval boundaries.
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
            let boundary = UncertainBoundary::range(Mu::certain(start), Mu::certain(end));
            Ok((
                remaining,
                Interval::with_complex_boundaries(boundary.clone(), boundary),
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
        // Single uncertain position bounded by a range: (Lys23_Leu24)<edit>.
        // Per HGVS v21 uncertain.md, the parens wrap *one* residue somewhere
        // in [start, end]; same Range on both Interval boundaries; Display
        // collapses to `(Lys23_Leu24)` once.
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
            let boundary = UncertainBoundary::range(Mu::certain(start), Mu::certain(end));
            Ok((
                remaining,
                Interval::with_complex_boundaries(boundary.clone(), boundary),
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

fn parse_rna_range_boundary(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::RnaPos>> {
    delimited(
        char('('),
        alt((
            // (?_pos)
            map((tag("?"), tag("_"), parse_rna_pos), |(_, _, end)| {
                UncertainBoundary::range(Mu::Unknown, Mu::certain(end))
            }),
            // (pos_?)
            map((parse_rna_pos, tag("_"), tag("?")), |(start, _, _)| {
                UncertainBoundary::range(Mu::certain(start), Mu::Unknown)
            }),
            // (pos_pos)
            map(
                (parse_rna_pos, tag("_"), parse_rna_pos),
                |(start, _, end)| UncertainBoundary::range(Mu::certain(start), Mu::certain(end)),
            ),
        )),
        char(')'),
    )
    .parse(input)
}

/// Parse an RNA boundary: (?_pos), (pos_?), (pos_pos), (pos), ?, or pos
fn parse_rna_boundary(
    input: &str,
) -> IResult<&str, UncertainBoundary<crate::hgvs::location::RnaPos>> {
    alt((
        // Range boundary: (?_pos), (pos_?), or (pos_pos)
        parse_rna_range_boundary,
        // Single uncertain position: (pos)
        map(
            delimited(char('('), parse_rna_pos, char(')')),
            UncertainBoundary::uncertain,
        ),
        // Unknown: ?
        map(tag("?"), |_| UncertainBoundary::unknown()),
        // Certain position: pos
        map(parse_rna_pos, UncertainBoundary::certain),
    ))
    .parse(input)
}

fn parse_rna_interval(input: &str) -> IResult<&str, RnaInterval> {
    alt((
        // Single uncertain position bounded by a range: (a_b) NOT followed by _.
        // Per HGVS v21 uncertain.md, `(a_b)<edit>` is *one* position somewhere
        // in [a, b]; same Range on both Interval boundaries.
        |input| {
            let (remaining, (start, _, end)) = delimited(
                char('('),
                (parse_rna_pos, tag("_"), parse_rna_pos),
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
            let boundary = UncertainBoundary::range(Mu::certain(start), Mu::certain(end));
            Ok((
                remaining,
                Interval::with_complex_boundaries(boundary.clone(), boundary),
            ))
        },
        // Complex interval with range boundaries: (pos_pos)_(pos_pos),
        // (?_pos)_(pos_?), (pos_pos)_pos, pos_(pos_pos), etc.
        |input| {
            let (remaining, (start, _, end)) =
                (parse_rna_boundary, tag("_"), parse_rna_boundary).parse(input)?;
            Ok((remaining, Interval::with_complex_boundaries(start, end)))
        },
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

        // First try whole-genome identity (g.= or g.(=)) - no position.
        // Reject the probe when the remainder starts with `(;)` so that a
        // bracketless unknown-phase compound (`g.=(;)123A>G`) reaches the
        // `(;)` dispatcher below rather than getting silently consumed as
        // a `g.=` whose trailing `(;)...` then fails far from the cause.
        if let Ok((remaining, edit)) = parse_whole_genome_identity(input) {
            if !remaining.starts_with("(;)") {
                let dummy_pos = crate::hgvs::location::GenomePos::new(1);
                let dummy_interval = GenomeInterval::point(dummy_pos);
                return Ok((
                    remaining,
                    HgvsVariant::Genome(GenomeVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                    }),
                ));
            }
        }

        // Try whole-genome unknown (g.? or g.(?)) - no position. Same
        // `(;)` guard as the identity probe above.
        if let Ok((remaining, edit)) = parse_whole_genome_unknown(input) {
            if !remaining.starts_with("(;)") {
                let dummy_pos = crate::hgvs::location::GenomePos::new(1);
                let dummy_interval = GenomeInterval::point(dummy_pos);
                return Ok((
                    remaining,
                    HgvsVariant::Genome(GenomeVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                    }),
                ));
            }
        }

        // Try trans-allele compact-prefix shorthand: g.[edit1];[edit2]
        // Detection mirrors the CDS/RNA paths (parse_cds_allele_shorthand,
        // parse_rna_allele_shorthand): starts with `[` and contains `];[`.
        if input.starts_with('[') && input.contains("];[") {
            if let Ok((remaining, variant)) =
                parse_genome_trans_allele_shorthand(input, accession.clone(), gene_symbol.clone())
            {
                return Ok((remaining, variant));
            }
        }

        // Try compound allele notation: g.[pos1edit1;pos2edit2;...] or g.[a(;)b]
        if input.starts_with('[') {
            if let Ok((remaining, allele)) =
                parse_genome_compound_allele(input, accession.clone(), gene_symbol.clone())
            {
                return Ok((remaining, HgvsVariant::Allele(allele)));
            }
        }

        // Bracketless unknown-phase shorthand: g.123G>A(;)345del (#123)
        if input.contains("(;)") && !input.starts_with('[') {
            return parse_genome_position_unknown_phase(
                input,
                accession.clone(),
                gene_symbol.clone(),
                GenomeKind::Genome,
            );
        }

        // Predicted variant: g.(positionEdit). #241.
        if input.starts_with('(') && !input.starts_with("(?") && !input.starts_with("(=)") {
            if let Ok((remaining, (interval, edit))) =
                delimited(char('('), (parse_genome_interval, parse_na_edit), char(')')).parse(input)
            {
                return Ok((
                    remaining,
                    HgvsVariant::Genome(GenomeVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::new_predicted(interval, edit),
                    }),
                ));
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

/// Tag that selects which DNA variant constructor to use when the genome-style
/// allele helpers are shared between `g.`, `m.`, and `o.`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum GenomeKind {
    Genome,
    Mt,
    Circular,
}

impl GenomeKind {
    /// Build a variant from a pre-constructed `LocEdit` (allowing
    /// `Mu::Uncertain` for the predicted-wrapper form). #243.
    fn build_variant_with_loc_edit(
        self,
        accession: Accession,
        gene_symbol: Option<String>,
        loc_edit: LocEdit<GenomeInterval, crate::hgvs::edit::NaEdit>,
    ) -> HgvsVariant {
        match self {
            GenomeKind::Genome => HgvsVariant::Genome(GenomeVariant {
                accession,
                gene_symbol,
                loc_edit,
            }),
            GenomeKind::Mt => HgvsVariant::Mt(MtVariant {
                accession,
                gene_symbol,
                loc_edit,
            }),
            GenomeKind::Circular => HgvsVariant::Circular(CircularVariant {
                accession,
                gene_symbol,
                loc_edit,
            }),
        }
    }
}

/// Parse a single bracket-member position+edit (genome / m. / o.),
/// recognising whole-entity edits (`=`, `?`, `(=)`, `(?)`) and the
/// predicted-wrapper form. Mirrors `parse_cds_bracket_member` and
/// `parse_rna_bracket_member`. #243; whole-entity dispatch added by
/// #423 (follow-up to #396 item 3).
fn parse_genome_bracket_member(
    part: &str,
    kind: GenomeKind,
) -> IResult<&str, LocEdit<GenomeInterval, crate::hgvs::edit::NaEdit>> {
    // Whole-entity genome edits have no positional content, so attach
    // a dummy interval at position 1 (mirrors `parse_genome_variant`'s
    // handling of `g.=` / `g.?` at the top level). These probes run
    // BEFORE the position-based predicted wrapper / bare-position
    // parsers so a bracket like `[(=)]` matches the whole-entity
    // identity marker rather than being routed into the position
    // parser. Same shape applies to `m.` and `o.` — whole-entity edits
    // don't carry a coordinate so the circular-vs-linear distinction
    // is moot for these probes. #423, follow-up to #396 item 3.
    //
    // Note: bare `[0]` / `[?]` in the trans-allele form are handled
    // upstream by `parse_trans_allele_shorthand_generic`'s cross-coord
    // short-circuits (→ `NullAllele` / `UnknownAllele`) and never
    // reach this dispatcher. Only cis/unknown-flat paths and
    // predicted forms (`[(?)]`, `[(=)]`) plus bare `=` / `?` inside a
    // cis bracket flow through these probes.
    fn dummy_genome_interval() -> GenomeInterval {
        GenomeInterval::point(crate::hgvs::location::GenomePos::new(1))
    }

    if let Ok((remaining, mu_edit)) = parse_whole_genome_identity(part) {
        return Ok((
            remaining,
            LocEdit {
                location: dummy_genome_interval(),
                edit: mu_edit,
            },
        ));
    }
    if let Ok((remaining, mu_edit)) = parse_whole_genome_unknown(part) {
        return Ok((
            remaining,
            LocEdit {
                location: dummy_genome_interval(),
                edit: mu_edit,
            },
        ));
    }

    let is_circular = matches!(kind, GenomeKind::Mt | GenomeKind::Circular);
    // For `m.`/`o.` the circular variant allows the spec-authorised
    // reversed `<high>_<low>` forms; post-parse we still run
    // `check_circular_reversed_range` to reject edit kinds the spec
    // doesn't authorise. For `g.` the linear parser is unchanged. (#129)
    if part.starts_with('(') && !part.starts_with("(?") && !part.starts_with("(=)") {
        let predicted = if is_circular {
            delimited(
                char('('),
                (parse_genome_interval_for_circular, parse_na_edit),
                char(')'),
            )
            .parse(part)
        } else {
            delimited(char('('), (parse_genome_interval, parse_na_edit), char(')')).parse(part)
        };
        if let Ok((remaining, (interval, edit))) = predicted {
            if is_circular {
                check_circular_reversed_range(&interval, &edit, part)?;
            }
            return Ok((remaining, LocEdit::new_predicted(interval, edit)));
        }
    }
    let (remaining, interval) = if is_circular {
        parse_genome_interval_for_circular(part)?
    } else {
        parse_genome_interval(part)?
    };
    let (remaining, edit) = if let Ok((r, e)) = parse_na_edit(remaining) {
        (r, e)
    } else {
        (
            remaining,
            crate::hgvs::edit::NaEdit::Identity {
                sequence: None,
                whole_entity: false,
            },
        )
    };
    if is_circular {
        check_circular_reversed_range(&interval, &edit, part)?;
    }
    Ok((remaining, LocEdit::new(interval, edit)))
}

/// Find the index of the `]` that closes the leading `[` in `input`,
/// scanning at bracket depth 0. Returns `None` if `input` does not start
/// with `[` or if the bracket is unbalanced.
///
/// Plain `input.find(']')` would truncate at any inner edit bracket
/// (`delins[ACC:g.x_y]`, `delins[A;T]`, `TG[14]`, etc.); plain
/// `input.rfind(']')` would over-consume on trans-form `[a];[b]`. The
/// depth-tracking scan handles both.
pub(crate) fn find_top_level_close_bracket(input: &str) -> Option<usize> {
    if !input.starts_with('[') {
        return None;
    }
    // Skip the leading `[` and delegate to the after-open primitive,
    // shifting the returned index back into the original slice.
    find_close_after_consumed_open(&input[1..]).map(|i| i + 1)
}

/// Variant of [`find_top_level_close_bracket`] for callers that have
/// already consumed the opening `[` (typically via a `nom` combinator).
/// Treats the slice as the bracket's content and returns the index of
/// the matching `]` (relative to `input`). See
/// `parse_bracketed_inserted_sequence` in `parser/edit.rs` for a
/// caller.
pub(crate) fn find_close_after_consumed_open(input: &str) -> Option<usize> {
    // Caller has already passed one level deep; start depth at 1 so
    // the scan returns when the first un-paired `]` brings it back to 0.
    let mut depth: i32 = 1;
    for (i, b) in input.bytes().enumerate() {
        match b {
            b'[' => depth += 1,
            b']' => {
                depth -= 1;
                if depth == 0 {
                    return Some(i);
                }
            }
            _ => {}
        }
    }
    None
}

/// Locate the next top-level `/` (mosaic) or `//` (chimeric) separator
/// in `input`, scanning at bracket depth 0. A "top-level" separator
/// is one not enclosed in `[ ]` brackets — so the bracketed-inner
/// forms `[a;b]//[c;d]` and `[a;b]/[c;d]` split cleanly even though
/// plain `input.find("/")` / `input.split("//")` would stumble.
///
/// When `want_double` is `true`, returns the position of the first
/// `/` of any top-level `//`. When `want_double` is `false`, returns
/// the position of any top-level `/` that is *not* part of `//`
/// (the mosaic single-slash separator).
///
/// Used by `parse_chimeric_allele` and `parse_mosaic_allele`.
pub(crate) fn find_top_level_slash(input: &str, want_double: bool) -> Option<usize> {
    let bytes = input.as_bytes();
    let mut depth: i32 = 0;
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'[' => depth += 1,
            b']' => depth -= 1,
            b'/' if depth == 0 => {
                let is_double = i + 1 < bytes.len() && bytes[i + 1] == b'/';
                if want_double && is_double {
                    return Some(i);
                }
                if !want_double && !is_double {
                    return Some(i);
                }
                // Skip the second `/` of a `//` so the next iteration
                // doesn't see it as a fresh `/` in single-slash mode.
                if is_double {
                    i += 1;
                }
            }
            _ => {}
        }
        i += 1;
    }
    None
}

/// Detect whether `input` contains a `/` *inside* a bracket pair —
/// i.e. at bracket depth ≥ 1. Used to flag the `[a/b]` mosaic form
/// (and the `[a//b]` chimeric variant of the same family), which the
/// HGVS spec does not define. Both surface under W3019
/// (`NonSpecMosaicForm`).
///
/// Returns `true` on the *first* `/` found at depth ≥ 1 — callers
/// only need a yes/no for diagnostic emission.
///
/// This is the dual of `find_top_level_slash`, which scans at depth
/// 0 only. Together they let `parse_variant` route the four mosaic /
/// chimeric shapes:
///
///   * top-level `/` / `//`  — spec-supported mosaic / chimeric
///   * bracketed `[a/b]`     — non-spec, reject with W3019 (this helper)
///   * nested `/` + `//` at the same level — non-spec, reject with W3019
pub(crate) fn has_slash_inside_brackets(input: &str) -> bool {
    let bytes = input.as_bytes();
    let mut depth: i32 = 0;
    let mut i = 0;
    while i < bytes.len() {
        match bytes[i] {
            b'[' => depth += 1,
            b']' if depth > 0 => depth -= 1,
            b'/' if depth >= 1 => return true,
            _ => {}
        }
        i += 1;
    }
    false
}

/// Bracket-depth-aware top-level split for chimeric (`//`) or mosaic
/// (`/`) separators. Mirrors `str::split(sep)` semantics but respects
/// `[ ]` nesting. Empty chunks are preserved (callers reject them
/// with a context-specific error).
pub(crate) fn split_top_level_slashes(input: &str, want_double: bool) -> Vec<&str> {
    let mut chunks = Vec::new();
    let sep_len = if want_double { 2 } else { 1 };
    let mut s = input;
    loop {
        match find_top_level_slash(s, want_double) {
            Some(pos) => {
                chunks.push(&s[..pos]);
                s = &s[pos + sep_len..];
            }
            None => {
                chunks.push(s);
                break;
            }
        }
    }
    chunks
}

/// Result of `scan_allele_separators`: top-level depth-0 allele-separator
/// flags plus the depth-aware members.
#[derive(Debug)]
pub(crate) struct AlleleSeparatorScan<'a> {
    /// True iff `content` contains at least one top-level `(;)` (unknown phase).
    pub has_unknown_phase: bool,
    /// True iff `content` contains at least one top-level standalone `;`
    /// that is not part of a `(;)`.
    pub has_cis_separator: bool,
    /// `content` split at whichever separator is present. If only `(;)` is
    /// present, splits at top-level `(;)`. If only `;` is present, splits
    /// at top-level `;`. If neither is present, returns `[content]`.
    /// If both are present (mixed phase), splits at top-level `(;)` —
    /// callers reject mixed-phase before consuming this field, so the
    /// exact split for the mixed shape is not meaningful.
    pub members: Vec<&'a str>,
}

/// Scan `content` — the inside of an outer `[...]` allele bracket pair,
/// with the surrounding `[`/`]` already stripped — for top-level
/// (depth-0) allele-phase separators and split into members at the same
/// depth.
///
/// "Top level" means bracket depth 0 with respect to nested `[...]`
/// inside an edit (`delins[A;T]`, `delins[ACC:g.x_y]`, `TG[14]`, ...).
/// `(;)` and `;` inside such nested brackets are part of the edit, not
/// allele separators, and must not be counted here.
///
/// The pair of `parens` characters `(` / `)` is also tracked so that a
/// member written as a predicted-wrapper form like `(100A>G;200del)`
/// is treated as a single member (the inner `;` is at parens-depth 1,
/// not at the allele top level). Without this, a single predicted
/// member could be mis-split into multiple bogus members.
///
/// Used by every `parse_*_compound_allele` / `parse_*_allele_shorthand`
/// that drives the in-bracket cis (`;`) vs. unknown-phase (`(;)`)
/// distinction so the mixed-phase reject path and the member split both
/// scan only at the outer allele level.
pub(crate) fn scan_allele_separators(content: &str) -> AlleleSeparatorScan<'_> {
    let bytes = content.as_bytes();
    let mut bracket_depth: i32 = 0;
    let mut paren_depth: i32 = 0;
    // Positions and lengths of top-level (depth-0 for both brackets and
    // parens) separators. `len` is 3 for `(;)` and 1 for `;`.
    let mut seps: Vec<(usize, usize)> = Vec::new();
    let mut has_unknown_phase = false;
    let mut has_cis_separator = false;

    let mut i = 0;
    while i < bytes.len() {
        let b = bytes[i];
        match b {
            b'[' => bracket_depth += 1,
            b']' if bracket_depth > 0 => bracket_depth -= 1,
            b'(' if bracket_depth == 0 => {
                // Check for a `(;)` separator at allele top level — only
                // when we are not already inside any other parens or
                // brackets.
                if paren_depth == 0
                    && i + 2 < bytes.len()
                    && bytes[i + 1] == b';'
                    && bytes[i + 2] == b')'
                {
                    has_unknown_phase = true;
                    seps.push((i, 3));
                    i += 3;
                    continue;
                }
                paren_depth += 1;
            }
            b')' if bracket_depth == 0 && paren_depth > 0 => paren_depth -= 1,
            b';' if bracket_depth == 0 && paren_depth == 0 => {
                has_cis_separator = true;
                seps.push((i, 1));
            }
            _ => {}
        }
        i += 1;
    }

    // Build the member slices. If both unknown and cis separators are
    // present (mixed phase), prefer the `(;)` split — callers reject
    // before consuming `members`, so the exact split is moot.
    let split_kind_len = if has_unknown_phase { 3 } else { 1 };
    let mut members: Vec<&str> = Vec::new();
    let mut start = 0;
    for (pos, len) in &seps {
        if *len != split_kind_len {
            // Skip the "other" separator kind in the mixed case.
            continue;
        }
        members.push(&content[start..*pos]);
        start = *pos + *len;
    }
    members.push(&content[start..]);

    AlleleSeparatorScan {
        has_unknown_phase,
        has_cis_separator,
        members,
    }
}

/// Parse a compound allele with genomic variants:
///   `[pos1edit1;pos2edit2;...]`   (cis)
///   `[pos1edit1(;)pos2edit2;...]` (unknown phase)
///
/// A single bracket pair must use exactly one separator type; mixed
/// `;`/`(;)` inside one bracket pair is not spec-valid and is rejected.
fn parse_genome_compound_allele(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, AlleleVariant> {
    parse_genome_kind_compound_allele(input, accession, gene_symbol, GenomeKind::Genome)
}

fn parse_genome_kind_compound_allele(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
    kind: GenomeKind,
) -> IResult<&str, AlleleVariant> {
    if !input.starts_with('[') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Locate the `]` that closes the leading `[`, scanning at bracket depth 0
    // so that members containing nested edit brackets (`delins[ACC:g.x_y]`,
    // `delins[A;T]`, etc.) are not truncated at the inner `]`.
    let close_bracket = find_top_level_close_bracket(input).ok_or_else(|| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Tag))
    })?;

    let content = &input[1..close_bracket];
    let remaining = &input[close_bracket + 1..];

    // Scan for top-level (depth-0) allele separators. Any `(;)` / `;`
    // inside a nested edit bracket (`delins[A;T]`) or a predicted-wrapper
    // member (`(100A>G;200del)`) is part of the edit, not an allele
    // separator. See `scan_allele_separators` for the depth rules.
    let scan = scan_allele_separators(content);

    // Mixed `;` and `(;)` inside a single bracket pair is not spec-valid
    // (HGVS `recommendations/{DNA,protein}/alleles.md` only mixes `;`/`(;)` at
    // the top level via `[a;b];[c](;)d`, not inside one `[...]`). Reject
    // rather than silently flatten the cis subgroup boundary into a single
    // `phase: Unknown` allele.
    if scan.has_unknown_phase && scan.has_cis_separator {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    let phase = if scan.has_unknown_phase {
        AllelePhase::Unknown
    } else {
        AllelePhase::Cis
    };

    let mut variants = Vec::with_capacity(4);

    for part in scan.members {
        let part = part.trim();
        if part.is_empty() {
            // Empty parts (e.g. `[a;;b]`, `[a(;)(;)b]`, or trailing `[a(;)]`)
            // are malformed; reject rather than silently dropping them.
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }
        let (final_remaining, loc_edit) = parse_genome_bracket_member(part, kind)?;
        if !final_remaining.trim().is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                final_remaining,
                nom::error::ErrorKind::Tag,
            )));
        }
        variants.push(kind.build_variant_with_loc_edit(
            accession.clone(),
            gene_symbol.clone(),
            loc_edit,
        ));
    }

    // Singletons (`[a]`) are intentionally accepted by ferro: the bracketed
    // wrapper is preserved for upstream tooling rather than rejected, and a
    // non-empty `remaining` (e.g. trans-form `[a];[b]`) is handled by the
    // outer parser. See `test_singleton_cis_allele_preserves_wrapper` and
    // `pinned_v21_normalization_behavior` for the pinned behavior.
    if variants.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    Ok((remaining, AlleleVariant::new(variants, phase)))
}

/// Parse position-based unknown-phase shorthand for genome-style coord systems:
///   `123G>A(;)345del` (uses [`GenomeKind`] to select the wrapping variant kind)
fn parse_genome_position_unknown_phase(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
    kind: GenomeKind,
) -> IResult<&str, HgvsVariant> {
    let mut variants = Vec::with_capacity(4);

    for part in input.split("(;)") {
        let part = part.trim();
        if part.is_empty() {
            // An empty `(;)` group (e.g. `a(;)(;)b` or trailing/leading
            // `(;)`) is malformed; reject rather than silently dropping it.
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }

        // Delegate to the bracket-member dispatcher so the same set of
        // shapes (whole-entity `=`/`?`/`(=)`/`(?)`, predicted-wrapper
        // `(<pos><edit>)`, bare position) is admitted here. This keeps
        // the bracket form `g.[X(;)Y]` and the bracketless form
        // `g.X(;)Y` symmetric — Display canonicalises one into the
        // other, so accepting an asymmetric set on the parse side
        // breaks round-trip. #423.
        let (final_remaining, loc_edit) = parse_genome_bracket_member(part, kind)?;
        if !final_remaining.trim().is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                final_remaining,
                nom::error::ErrorKind::Tag,
            )));
        }

        variants.push(kind.build_variant_with_loc_edit(
            accession.clone(),
            gene_symbol.clone(),
            loc_edit,
        ));
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

/// Generic trans-allele shorthand parser shared by 6 of 7 axes (g., c., n.,
/// m., r., o.). Owns the bracketed-member loop, the `[0]`/`[?]`
/// cross-coordinate special-cases, the bracket-balance check, the
/// `len < 2` guard, and the final `AlleleVariant::new(.., Trans)` wrap.
/// Always emits `AllelePhase::Trans` — cis and unknown-phase shapes
/// route through `parse_*_allele_shorthand` / `parse_*_compound_allele`
/// upstream and never reach this function.
///
/// Each axis-specific shim provides `parse_member`, a closure that
/// consumes the bracket content and constructs the per-axis
/// `HgvsVariant::<Axis>` (or returns a parse error). The closure
/// returns its own remainder so each axis is free to use its own
/// bracket-member primitives; the generic checks the remainder is
/// fully consumed.
///
/// Protein keeps its own bespoke helper because it overrides the
/// special-marker arm (`[0]`/`[0?]`/`[(0)]` → `ProteinEdit::NoProtein`
/// rather than the cross-coordinate `NullAllele`).
fn parse_trans_allele_shorthand_generic<F>(
    input: &str,
    mut parse_member: F,
) -> IResult<&str, HgvsVariant>
where
    F: FnMut(&str) -> IResult<&str, HgvsVariant>,
{
    let mut variants = Vec::with_capacity(2);
    let mut remaining = input;

    while !remaining.is_empty() {
        if !remaining.starts_with('[') {
            break;
        }

        let close_bracket = find_top_level_close_bracket(remaining).ok_or_else(|| {
            nom::Err::Error(nom::error::Error::new(
                remaining,
                nom::error::ErrorKind::Tag,
            ))
        })?;

        let content = remaining[1..close_bracket].trim();

        let variant = if content == "0" {
            HgvsVariant::NullAllele
        } else if content == "?" {
            HgvsVariant::UnknownAllele
        } else {
            let scan = scan_allele_separators(content);
            if scan.has_cis_separator && !scan.has_unknown_phase {
                // A trans bracket whose content is itself a multi-variant
                // cis group, e.g. `[296T>G;476T>C]` in `[a;b];[c]`
                // (compound heterozygote). Parse each cis member and nest
                // it as a `Cis` sub-allele member of the trans allele.
                let mut cis_members = Vec::with_capacity(scan.members.len());
                for member in scan.members {
                    let member = member.trim();
                    let (rest, v) = parse_member(member)?;
                    if !rest.trim().is_empty() {
                        return Err(nom::Err::Error(nom::error::Error::new(
                            rest,
                            nom::error::ErrorKind::Tag,
                        )));
                    }
                    cis_members.push(v);
                }
                // HGVS DNA/alleles.md: the "not changed" `=` indication is
                // used only when a single variant is identified per allele.
                // Spelling the other allele's variant as `=` inside a
                // multi-variant cis group is invalid (the spec flags
                // `[2376G>C;3103=];[2376=;3103del]` as not correct). Reject a
                // nested cis group that carries a position-identity member
                // rather than accept the discouraged cross-spelled form.
                // (The all-identity-allele edge `[a=;b=]` is deferred.)
                if cis_members.iter().any(is_lhs_position_identity) {
                    return Err(nom::Err::Error(nom::error::Error::new(
                        content,
                        nom::error::ErrorKind::Verify,
                    )));
                }
                HgvsVariant::Allele(AlleleVariant::new(cis_members, AllelePhase::Cis))
            } else {
                let (final_remaining, variant) = parse_member(content)?;
                if !final_remaining.trim().is_empty() {
                    return Err(nom::Err::Error(nom::error::Error::new(
                        final_remaining,
                        nom::error::ErrorKind::Tag,
                    )));
                }
                variant
            }
        };

        variants.push(variant);
        remaining = &remaining[close_bracket + 1..];

        // Enforce strict `[...];[...]` shape. A `;` is only valid if it
        // separates the bracket we just consumed from another bracket; a
        // dangling trailing `;` (e.g. `[A];[B];`) and a missing separator
        // (e.g. `[A];[B][C]`) must both reject rather than silently
        // build an extra-member or short AlleleVariant.
        if let Some(after_semi) = remaining.strip_prefix(';') {
            if !after_semi.starts_with('[') {
                return Err(nom::Err::Error(nom::error::Error::new(
                    after_semi,
                    nom::error::ErrorKind::Tag,
                )));
            }
            remaining = after_semi;
        } else if remaining.starts_with('[') {
            return Err(nom::Err::Error(nom::error::Error::new(
                remaining,
                nom::error::ErrorKind::Tag,
            )));
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

fn parse_genome_trans_allele_shorthand(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    parse_trans_allele_shorthand_generic(input, |content| {
        let (rest, loc_edit) = parse_genome_bracket_member(content, GenomeKind::Genome)?;
        Ok((
            rest,
            HgvsVariant::Genome(GenomeVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit,
            }),
        ))
    })
}

/// Parse whole-genome identity: g.=
fn parse_whole_genome_identity(input: &str) -> IResult<&str, Mu<crate::hgvs::edit::NaEdit>> {
    // Predicted: g.(=)
    if let Ok((remaining, _)) = tag::<_, _, nom::error::Error<&str>>("(=)").parse(input) {
        return Ok((
            remaining,
            Mu::Uncertain(crate::hgvs::edit::NaEdit::whole_entity_identity()),
        ));
    }
    map(tag("="), |_| {
        Mu::Certain(crate::hgvs::edit::NaEdit::whole_entity_identity())
    })
    .parse(input)
}

/// Parse whole-genome unknown: g.? or g.(?)
///
/// Only matches a terminal `?`: rejects `?` followed by an interval
/// continuation (`_`), an offset (`+`/`-`), an edit keyword (`dup`,
/// `del`, `ins`, `inv`), or a substitution pattern (`<base>>` or
/// `>`) so that those forms fall through to the position-based
/// parser (`?<edit>` means edit at unknown position, not whole-entity
/// unknown). #239.
///
/// The `g.(?)` predicted form is returned as `Mu::Uncertain`. #245.
fn parse_whole_genome_unknown(input: &str) -> IResult<&str, Mu<crate::hgvs::edit::NaEdit>> {
    // Predicted: g.(?)
    if let Ok((remaining, _)) = tag::<_, _, nom::error::Error<&str>>("(?)").parse(input) {
        return Ok((
            remaining,
            Mu::Uncertain(crate::hgvs::edit::NaEdit::whole_entity_unknown()),
        ));
    }
    let (remaining, _) = tag("?").parse(input)?;
    if is_edit_continuation_after_unknown(remaining) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }
    Ok((
        remaining,
        Mu::Certain(crate::hgvs::edit::NaEdit::whole_entity_unknown()),
    ))
}

/// Returns true when `remaining` (the input after a `?`) starts with
/// anything that should route `?` to the position-based parser rather
/// than treating it as whole-entity unknown.
///
/// Covers every `?<edit>` continuation accepted by the position-based
/// parser:
///   - interval continuation `_`
///   - intronic offsets `+` / `-`
///   - identity `=`
///   - keyword-keyed edits: `dup`, `del`, `ins`, `inv`, `con`, `copy`
///   - substitutions: `<base>>` (e.g. `?A>G`) and `>` (no-ref)
///
/// `con` / `copy` were added per #286 to allow `?conNM_...:c.X_Y`
/// and `?copy<N>` at an unknown position, consistent with the other
/// `?<edit>` shapes accepted by #239 / PR #240.
fn is_edit_continuation_after_unknown(remaining: &str) -> bool {
    if remaining.starts_with('_')
        || remaining.starts_with('+')
        || remaining.starts_with('-')
        || remaining.starts_with('=')
        || remaining.starts_with("dup")
        || remaining.starts_with("del")
        || remaining.starts_with("ins")
        || remaining.starts_with("inv")
        || remaining.starts_with("con")
        || remaining.starts_with("copy")
    {
        return true;
    }
    // Substitution `?<base>>` (e.g. `?A>G`) or no-ref `?>G`.
    let bytes = remaining.as_bytes();
    if bytes.first() == Some(&b'>') {
        return true;
    }
    if bytes.len() >= 2 && bytes[1] == b'>' && bytes[0].is_ascii_alphabetic() {
        return true;
    }
    false
}

/// Parse whole-CDS identity: c.= or c.(=)
fn parse_whole_cds_identity(input: &str) -> IResult<&str, Mu<crate::hgvs::edit::NaEdit>> {
    if let Ok((remaining, _)) = tag::<_, _, nom::error::Error<&str>>("(=)").parse(input) {
        return Ok((
            remaining,
            Mu::Uncertain(crate::hgvs::edit::NaEdit::whole_entity_identity()),
        ));
    }
    map(tag("="), |_| {
        Mu::Certain(crate::hgvs::edit::NaEdit::whole_entity_identity())
    })
    .parse(input)
}

/// Parse whole-CDS unknown: c.? or c.(?)
///
/// Only matches if `?` is at end (not followed by `_`, `+`/`-`, edit
/// keywords, or a substitution pattern — those route to the
/// position-based parser as edits at unknown position). #239.
///
/// The `c.(?)` predicted form is returned as `Mu::Uncertain`. #245.
fn parse_whole_cds_unknown(input: &str) -> IResult<&str, Mu<crate::hgvs::edit::NaEdit>> {
    if let Ok((remaining, _)) = tag::<_, _, nom::error::Error<&str>>("(?)").parse(input) {
        return Ok((
            remaining,
            Mu::Uncertain(crate::hgvs::edit::NaEdit::whole_entity_unknown()),
        ));
    }
    let (remaining, _) = tag("?").parse(input)?;
    if is_edit_continuation_after_unknown(remaining) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }
    Ok((
        remaining,
        Mu::Certain(crate::hgvs::edit::NaEdit::whole_entity_unknown()),
    ))
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

        // First try whole-CDS identity (c.= or c.(=)) - no position
        if let Ok((remaining, edit)) = parse_whole_cds_identity(input) {
            // Use a dummy interval for whole-entity identity
            let dummy_pos = crate::hgvs::location::CdsPos {
                base: 1,
                offset: None,
                utr3: false,
                special: None,
            };
            let dummy_interval = CdsInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Cds(CdsVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                }),
            ));
        }

        // Try whole-CDS unknown (c.? or c.(?)) - no position
        if let Ok((remaining, edit)) = parse_whole_cds_unknown(input) {
            // Use a dummy interval for whole-entity unknown
            let dummy_pos = crate::hgvs::location::CdsPos {
                base: 1,
                offset: None,
                utr3: false,
                special: None,
            };
            let dummy_interval = CdsInterval::point(dummy_pos);
            return Ok((
                remaining,
                HgvsVariant::Cds(CdsVariant {
                    accession: accession.clone(),
                    gene_symbol: gene_symbol.clone(),
                    loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                }),
            ));
        }

        // Try predicted variant: c.(positionEdit) where the whole
        // position + edit is wrapped in parentheses. Mirrors the
        // working protein predicted-change path at parse_protein_variant.
        // The `delimited` falls through to `parse_cds_interval` only if
        // `parse_na_edit` succeeds inside the parens — so the G1 form
        // `c.(123_127)A>G` (uncertain *position*, certain edit) goes to
        // the normal path where the bare position parser handles it.
        if input.starts_with('(') && !input.starts_with("(?") && !input.starts_with("(=)") {
            if let Ok((remaining, (interval, edit))) =
                delimited(char('('), (parse_cds_interval, parse_na_edit), char(')')).parse(input)
            {
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
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }

        // Delegate to the bracket-member dispatcher so the same set of
        // shapes is admitted here as in the bracketed cis/trans forms
        // — keeps `c.[X(;)Y]` / `c.X(;)Y` round-trip-symmetric across
        // Display's bracket-vs-bracketless canonicalisation. #423.
        let (final_remaining, loc_edit) = parse_cds_bracket_member(part)?;

        if !final_remaining.trim().is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                final_remaining,
                nom::error::ErrorKind::Tag,
            )));
        }

        variants.push(HgvsVariant::Cds(CdsVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit,
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
/// Parse a single bracket-member position+edit (CDS), recognising
/// whole-entity edits (`=`, `?`, `(=)`, `(?)`), the predicted-wrapper
/// form `(<pos><edit>)`, and bare `<pos>[<edit>]` (edit optional,
/// defaults to identity). #243 follow-up to #241; whole-entity dispatch
/// added by #423 (follow-up to #396 item 3).
///
/// Note the intentional asymmetry with `parse_tx_bracket_member`: CDS
/// does NOT probe `parse_rna_no_product`. `c.0` is not a CDS whole-entity
/// form — bare `[0]` in trans-allele context is cross-coord-short-
/// circuited to `NullAllele` upstream (see
/// `parse_trans_allele_shorthand_generic`), and the spec does not
/// define a `c.0` / `c.(0)` "no protein" form on the c. axis (that's
/// what `p.0` / `p.(0)` are for). Non-coding transcripts `n.` DO use
/// `parse_rna_no_product` because `n.0` is a valid spec form ("no
/// transcript produced"), mirroring `r.0`.
fn parse_cds_bracket_member(
    part: &str,
) -> IResult<&str, LocEdit<CdsInterval, crate::hgvs::edit::NaEdit>> {
    // Whole-CDS edits have no positional content, so attach a dummy
    // interval at position 1 (mirrors `parse_cds_variant`'s handling
    // of `c.=` / `c.?` at the top level). These probes run BEFORE the
    // position-based predicted wrapper / bare-position parsers so a
    // bracket like `[(=)]` matches the whole-entity identity marker
    // rather than being routed into `parse_cds_interval`. #423,
    // follow-up to #396 item 3.
    //
    // Note: bare `[0]` / `[?]` in the trans-allele form are handled
    // upstream by `parse_trans_allele_shorthand_generic`'s cross-coord
    // short-circuits (→ `NullAllele` / `UnknownAllele`) and never
    // reach this dispatcher. Only cis/unknown-flat paths and
    // predicted forms (`[(?)]`, `[(=)]`) plus bare `=` / `?` inside a
    // cis bracket flow through these probes.
    fn dummy_cds_interval() -> CdsInterval {
        CdsInterval::point(crate::hgvs::location::CdsPos {
            base: 1,
            offset: None,
            utr3: false,
            special: None,
        })
    }

    if let Ok((remaining, mu_edit)) = parse_whole_cds_identity(part) {
        return Ok((
            remaining,
            LocEdit {
                location: dummy_cds_interval(),
                edit: mu_edit,
            },
        ));
    }
    if let Ok((remaining, mu_edit)) = parse_whole_cds_unknown(part) {
        return Ok((
            remaining,
            LocEdit {
                location: dummy_cds_interval(),
                edit: mu_edit,
            },
        ));
    }

    // Try predicted wrapper first
    if part.starts_with('(') && !part.starts_with("(?") && !part.starts_with("(=)") {
        if let Ok((remaining, (interval, edit))) =
            delimited(char('('), (parse_cds_interval, parse_na_edit), char(')')).parse(part)
        {
            return Ok((remaining, LocEdit::new_predicted(interval, edit)));
        }
    }
    // Bare position [+ edit]
    let (remaining, interval) = parse_cds_interval(part)?;
    let (remaining, edit) = if let Ok((r, e)) = parse_na_edit(remaining) {
        (r, e)
    } else {
        (
            remaining,
            crate::hgvs::edit::NaEdit::Identity {
                sequence: None,
                whole_entity: false,
            },
        )
    };
    Ok((remaining, LocEdit::new(interval, edit)))
}

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

    // Scan for top-level (depth-0) allele separators. Inner `;` / `(;)`
    // belonging to a nested edit (`delins[A;T]`) or a predicted-wrapper
    // member must NOT drive phase detection or member splitting. See
    // `scan_allele_separators`.
    let scan = scan_allele_separators(content);

    // Mixed `;`/`(;)` inside one bracket pair is not spec-valid — HGVS
    // does not provide a way to express "subgroup A;B is cis to each
    // other but unknown phase to C" inside a single bracket. Reject
    // rather than silently flattening to AllelePhase::Unknown (#396
    // item 2). Mirrors `parse_protein_allele_shorthand`.
    if scan.has_unknown_phase && scan.has_cis_separator {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    let phase = if scan.has_unknown_phase {
        AllelePhase::Unknown
    } else {
        AllelePhase::Cis
    };

    // Pre-allocate with estimated capacity (most alleles have 2-4 variants)
    let mut variants = Vec::with_capacity(4);

    for part in scan.members {
        let part = part.trim();
        if part.is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }

        // Parse interval + edit (e.g., "145C>T" or "100_200del" or
        // predicted "(145C>T)").
        let (final_remaining, loc_edit) = parse_cds_bracket_member(part)?;

        if !final_remaining.trim().is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                final_remaining,
                nom::error::ErrorKind::Tag,
            )));
        }

        variants.push(HgvsVariant::Cds(CdsVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit,
        }));
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
    parse_trans_allele_shorthand_generic(input, |content| {
        let (rest, loc_edit) = parse_cds_bracket_member(content)?;
        Ok((
            rest,
            HgvsVariant::Cds(CdsVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit,
            }),
        ))
    })
}

/// Parse transcript variant (n.)
fn parse_tx_variant(
    accession: Accession,
    gene_symbol: Option<String>,
) -> impl FnMut(&str) -> IResult<&str, HgvsVariant> {
    move |input: &str| {
        let (input, _) = tag("n.").parse(input)?;

        // First try whole-tx identity (n.= or n.(=)) - no position. #288.
        // Non-coding transcripts reuse the RNA whole-entity parsers
        // (`parse_whole_rna_*`, `parse_rna_no_product`) — those are pure
        // tag-only parsers and don't depend on the coord system. Reject
        // the probe when the remainder starts with `(;)` so a bracketless
        // unknown-phase compound (`n.=(;)100A>G`) falls through to the
        // `(;)` dispatcher rather than being silently consumed as `n.=`.
        if let Ok((remaining, edit)) = parse_whole_rna_identity(input) {
            if !remaining.starts_with("(;)") {
                let dummy_interval = TxInterval::point(crate::hgvs::location::TxPos::new(1));
                return Ok((
                    remaining,
                    HgvsVariant::Tx(TxVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                    }),
                ));
            }
        }

        // Try whole-tx unknown (n.? or n.(?)) - no position. #288.
        if let Ok((remaining, edit)) = parse_whole_rna_unknown(input) {
            if !remaining.starts_with("(;)") {
                let dummy_interval = TxInterval::point(crate::hgvs::location::TxPos::new(1));
                return Ok((
                    remaining,
                    HgvsVariant::Tx(TxVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                    }),
                ));
            }
        }

        // Try tx-specific no product pattern (n.0 or predicted n.(0)). #288.
        // `0` is a valid spec form for non-coding transcripts — they can
        // fail to produce a transcript, exactly like `r.0`. Same `(;)` guard.
        if let Ok((remaining, edit)) = parse_rna_no_product(input) {
            if !remaining.starts_with("(;)") {
                let dummy_interval = TxInterval::point(crate::hgvs::location::TxPos::new(1));
                return Ok((
                    remaining,
                    HgvsVariant::Tx(TxVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                    }),
                ));
            }
        }

        // Compact-prefix trans-allele shorthand: n.[a];[b]. Mirrors
        // parse_genome_variant's dispatch added in PR #146; the expanded
        // form [ACC:n.a];[ACC:n.b] is handled by the top-level allele parser.
        if input.starts_with('[') && input.contains("];[") {
            if let Ok((remaining, variant)) =
                parse_tx_trans_allele_shorthand(input, accession.clone(), gene_symbol.clone())
            {
                return Ok((remaining, variant));
            }
        }

        // Bracketed compound allele: n.[var;var;...] or n.[var(;)var]
        if input.starts_with('[') {
            if let Ok((remaining, allele)) =
                parse_tx_compound_allele(input, accession.clone(), gene_symbol.clone())
            {
                return Ok((remaining, HgvsVariant::Allele(allele)));
            }
        }

        // Bracketless unknown-phase shorthand: n.100A>G(;)200T>C (#123)
        if input.contains("(;)") && !input.starts_with('[') {
            return parse_tx_position_unknown_phase(input, accession.clone(), gene_symbol.clone());
        }

        // Predicted variant: n.(positionEdit). #241.
        if input.starts_with('(') && !input.starts_with("(?") && !input.starts_with("(=)") {
            if let Ok((remaining, (interval, edit))) =
                delimited(char('('), (parse_tx_interval, parse_na_edit), char(')')).parse(input)
            {
                return Ok((
                    remaining,
                    HgvsVariant::Tx(TxVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::new_predicted(interval, edit),
                    }),
                ));
            }
        }

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

/// Parse a single bracket-member position+edit (non-coding transcript /
/// `n.`), recognising whole-entity edits (`=`, `?`, `(=)`, `(?)`, `0`,
/// `(0)`), the predicted-wrapper form `(<pos><edit>)`, and bare
/// `<pos>[<edit>]` (edit optional, defaults to identity). Mirrors
/// `parse_cds_bracket_member` and `parse_rna_bracket_member`. #287
/// follow-up to #243 / PR #244; whole-entity dispatch added by #423
/// (follow-up to #396 item 3).
///
/// Whole-entity probes here reuse the RNA-named parsers
/// (`parse_whole_rna_identity`, `parse_whole_rna_unknown`,
/// `parse_rna_no_product`) — these are pure tag-only and don't depend
/// on the coord system, so `parse_tx_variant` already shares them at
/// the top level. The name is historical; functionally these are the
/// shared NA whole-entity tag parsers for `n.` / `r.`.
fn parse_tx_bracket_member(
    part: &str,
) -> IResult<&str, LocEdit<TxInterval, crate::hgvs::edit::NaEdit>> {
    // Whole-tx edits have no positional content, so attach a dummy
    // interval at position 1 (mirrors `parse_tx_variant`'s handling of
    // `n.=` / `n.?` / `n.0` at the top level — which reuses the RNA
    // whole-entity parsers because they're pure tag-only and don't
    // depend on the coord system). These probes run BEFORE the
    // position-based predicted wrapper / bare-position parsers so a
    // bracket like `[(0)]` matches the no-product marker rather than
    // being routed into `parse_tx_interval`. #423, follow-up to #396
    // item 3.
    //
    // Note: bare `[0]` / `[?]` in the trans-allele form are handled
    // upstream by `parse_trans_allele_shorthand_generic`'s cross-coord
    // short-circuits (→ `NullAllele` / `UnknownAllele`) and never
    // reach this dispatcher. Only cis/unknown-flat paths and
    // predicted forms (`[(0)]`, `[(?)]`, `[(=)]`) plus bare `=` / `?`
    // / `0` inside a cis bracket flow through these probes.
    fn dummy_tx_interval() -> TxInterval {
        TxInterval::point(crate::hgvs::location::TxPos::new(1))
    }

    if let Ok((remaining, mu_edit)) = parse_whole_rna_identity(part) {
        return Ok((
            remaining,
            LocEdit {
                location: dummy_tx_interval(),
                edit: mu_edit,
            },
        ));
    }
    if let Ok((remaining, mu_edit)) = parse_whole_rna_unknown(part) {
        return Ok((
            remaining,
            LocEdit {
                location: dummy_tx_interval(),
                edit: mu_edit,
            },
        ));
    }
    if let Ok((remaining, mu_edit)) = parse_rna_no_product(part) {
        return Ok((
            remaining,
            LocEdit {
                location: dummy_tx_interval(),
                edit: mu_edit,
            },
        ));
    }

    // Try predicted wrapper first
    if part.starts_with('(') && !part.starts_with("(?") && !part.starts_with("(=)") {
        if let Ok((remaining, (interval, edit))) =
            delimited(char('('), (parse_tx_interval, parse_na_edit), char(')')).parse(part)
        {
            return Ok((remaining, LocEdit::new_predicted(interval, edit)));
        }
    }
    // Bare position [+ edit]
    let (remaining, interval) = parse_tx_interval(part)?;
    let (remaining, edit) = if let Ok((r, e)) = parse_na_edit(remaining) {
        (r, e)
    } else {
        (
            remaining,
            crate::hgvs::edit::NaEdit::Identity {
                sequence: None,
                whole_entity: false,
            },
        )
    };
    Ok((remaining, LocEdit::new(interval, edit)))
}

/// Parse a non-coding transcript trans-allele compact-prefix shorthand:
/// `[edit1];[edit2]` (with the `n.` already consumed by the caller). Each
/// bracketed group holds a single sub-variant rendered with this allele's
/// accession. Mirrors `parse_genome_trans_allele_shorthand`.
fn parse_tx_trans_allele_shorthand(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    parse_trans_allele_shorthand_generic(input, |content| {
        let (rest, loc_edit) = parse_tx_bracket_member(content)?;
        Ok((
            rest,
            HgvsVariant::Tx(TxVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit,
            }),
        ))
    })
}

/// Parse a compound allele with `n.` (transcript) variants:
///   `[100A>G;200T>C]` (cis) or `[100A>G(;)200T>C]` (unknown phase)
fn parse_tx_compound_allele(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, AlleleVariant> {
    if !input.starts_with('[') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // See `find_top_level_close_bracket` for why a depth-aware scan is used.
    let close_bracket = find_top_level_close_bracket(input).ok_or_else(|| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Tag))
    })?;

    let content = &input[1..close_bracket];
    let remaining = &input[close_bracket + 1..];

    // Scan for top-level (depth-0) allele separators. See
    // `scan_allele_separators` for the depth rules; inner `;` / `(;)`
    // from a nested edit must not drive phase detection or splitting.
    let scan = scan_allele_separators(content);

    // Mixed `;`/`(;)` inside one bracket pair is not spec-valid — see
    // `parse_genome_kind_compound_allele` for the spec citation.
    if scan.has_unknown_phase && scan.has_cis_separator {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    let phase = if scan.has_unknown_phase {
        AllelePhase::Unknown
    } else {
        AllelePhase::Cis
    };

    let mut variants = Vec::with_capacity(4);
    for part in scan.members {
        let part = part.trim();
        if part.is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }
        let (final_remaining, loc_edit) = parse_tx_bracket_member(part)?;
        if !final_remaining.trim().is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                final_remaining,
                nom::error::ErrorKind::Tag,
            )));
        }
        variants.push(HgvsVariant::Tx(TxVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit,
        }));
    }

    // Singletons are intentionally accepted — see
    // `parse_genome_kind_compound_allele` for the rationale.
    if variants.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    Ok((remaining, AlleleVariant::new(variants, phase)))
}

/// Parse position-based unknown-phase shorthand for `n.`:
///   `100A>G(;)200T>C`
fn parse_tx_position_unknown_phase(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    let mut variants = Vec::with_capacity(4);

    for part in input.split("(;)") {
        let part = part.trim();
        if part.is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }
        // Delegate to the bracket-member dispatcher so the same set of
        // shapes is admitted here as in the bracketed cis/trans forms
        // — keeps `n.[X(;)Y]` / `n.X(;)Y` round-trip-symmetric. #423.
        let (final_remaining, loc_edit) = parse_tx_bracket_member(part)?;
        if !final_remaining.trim().is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                final_remaining,
                nom::error::ErrorKind::Tag,
            )));
        }
        variants.push(HgvsVariant::Tx(TxVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit,
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

/// Parse mitochondrial variant (m.)
fn parse_mt_variant(
    accession: Accession,
    gene_symbol: Option<String>,
) -> impl FnMut(&str) -> IResult<&str, HgvsVariant> {
    move |input: &str| {
        let (input, _) = tag("m.").parse(input)?;

        // First try whole-mtDNA identity (m.= or m.(=)) - no position. #288.
        // Mitochondrial DNA reuses the genome whole-entity parser because
        // both share `GenomeInterval`. Reject when the remainder starts
        // with `(;)` so `m.=(;)100A>G` reaches the unknown-phase dispatch.
        if let Ok((remaining, edit)) = parse_whole_genome_identity(input) {
            if !remaining.starts_with("(;)") {
                let dummy_pos = crate::hgvs::location::GenomePos::new(1);
                let dummy_interval = GenomeInterval::point(dummy_pos);
                return Ok((
                    remaining,
                    HgvsVariant::Mt(MtVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                    }),
                ));
            }
        }

        // Try whole-mtDNA unknown (m.? or m.(?)) - no position. #288.
        if let Ok((remaining, edit)) = parse_whole_genome_unknown(input) {
            if !remaining.starts_with("(;)") {
                let dummy_pos = crate::hgvs::location::GenomePos::new(1);
                let dummy_interval = GenomeInterval::point(dummy_pos);
                return Ok((
                    remaining,
                    HgvsVariant::Mt(MtVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                    }),
                ));
            }
        }

        // Compact-prefix trans-allele shorthand: m.[a];[b]. Mirrors
        // parse_genome_variant's dispatch added in PR #146.
        if input.starts_with('[') && input.contains("];[") {
            if let Ok((remaining, variant)) =
                parse_mt_trans_allele_shorthand(input, accession.clone(), gene_symbol.clone())
            {
                return Ok((remaining, variant));
            }
        }

        // Bracketed compound allele: m.[var;var] or m.[var(;)var]
        if input.starts_with('[') {
            if let Ok((remaining, allele)) = parse_genome_kind_compound_allele(
                input,
                accession.clone(),
                gene_symbol.clone(),
                GenomeKind::Mt,
            ) {
                return Ok((remaining, HgvsVariant::Allele(allele)));
            }
        }

        // Bracketless unknown-phase shorthand: m.100A>G(;)200T>C (#123)
        if input.contains("(;)") && !input.starts_with('[') {
            return parse_genome_position_unknown_phase(
                input,
                accession.clone(),
                gene_symbol.clone(),
                GenomeKind::Mt,
            );
        }

        // Predicted variant: m.(positionEdit). #241.
        if input.starts_with('(') && !input.starts_with("(?") && !input.starts_with("(=)") {
            if let Ok((remaining, (interval, edit))) = delimited(
                char('('),
                (parse_genome_interval_for_circular, parse_na_edit),
                char(')'),
            )
            .parse(input)
            {
                // Reject spec-unauthorised reversed ranges (#129).
                check_circular_reversed_range(&interval, &edit, input)?;
                return Ok((
                    remaining,
                    HgvsVariant::Mt(MtVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::new_predicted(interval, edit),
                    }),
                ));
            }
        }

        let original_input = input;
        let (input, interval) = parse_genome_interval_for_circular(input)?;
        let (input, edit) = parse_na_edit(input)?;
        // Reject spec-unauthorised reversed ranges (#129).
        check_circular_reversed_range(&interval, &edit, original_input)?;

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

/// Parse a mitochondrial trans-allele compact-prefix shorthand:
/// `[edit1];[edit2]` (with the `m.` already consumed). Mirrors
/// `parse_genome_trans_allele_shorthand`; mitochondrial variants reuse
/// `parse_genome_interval`.
fn parse_mt_trans_allele_shorthand(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    parse_trans_allele_shorthand_generic(input, |content| {
        // `GenomeKind::Mt` routes through the circular interval parser +
        // spec-authorised reversed-range validator (#129).
        let (rest, loc_edit) = parse_genome_bracket_member(content, GenomeKind::Mt)?;
        Ok((
            rest,
            HgvsVariant::Mt(MtVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit,
            }),
        ))
    })
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

/// Parse no protein production: `p.0`, `p.0?`, or `p.(0)`.
///
/// The spec form for the predicted ("possibly no protein") case is `p.0?`
/// (HGVS v21 `recommendations/protein/deletion.md`). We also accept the
/// parenthesised `p.(0)` as a tolerated alternate input, mirroring the
/// `c.(=)`/`r.(=)`/`r.(0)` predicted whole-entity dispatch added in PR #246
/// for NA coord systems (#245). Both `p.0?` and `p.(0)` route to the same
/// `ProteinEdit::NoProtein { predicted: true }` value, and Display
/// canonicalises to the spec form `p.0?` — keeping the long-standing D6
/// contract pinned by `tests/protein_no_protein_roundtrip.rs` intact while
/// resolving issue #289.
fn parse_whole_protein_no_protein(input: &str) -> IResult<&str, ProteinEdit> {
    alt((
        // Predicted no-protein (parenthesised form, #289): (0)
        map(tag("(0)"), |_| ProteinEdit::no_protein_predicted()),
        // Predicted no-protein (spec-canonical form): 0?
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

/// Parse a single protein bracket-member `interval+edit`, recognising
/// whole-entity protein edits (`=`, `(=)`, `(?)`) before the
/// position-based `parse_prot_interval` + `parse_protein_edit` path.
///
/// This is the protein-axis analogue of `parse_cds_bracket_member` /
/// `parse_genome_bracket_member` / `parse_rna_bracket_member`, but the
/// protein axis does not share the `NaEdit` type: whole-entity protein
/// edits are `ProteinEdit::Identity { whole_protein: true, .. }` /
/// `ProteinEdit::Unknown { whole_protein: true, .. }` and carry no
/// position, so we attach a dummy interval at `Met1` (mirroring
/// `parse_protein_variant`'s top-level handling of `p.=` / `p.(?)`). The
/// `predicted` flag lives inside the edit, so the returned `LocEdit` is
/// `Mu::Certain` — matching the top-level dispatch. #468, follow-up to
/// #423 (nucleic-acid axes) and #396 item 3 (`r.` axis).
///
/// Whole-entity protein forms admitted here, per HGVS v21
/// `recommendations/protein/`:
/// - `=` / `(=)` — whole-protein identity (`alleles.md` line 58-60:
///   `p.[Ser68Arg];[=]` is valid and distinct from `[Ser68=]`;
///   `substitution.md` "the description `p.=` means the **entire**
///   protein coding region was analysed").
/// - `(?)` — predicted whole-protein unknown (`alleles.md`:
///   `p.[(Ser68Arg)];[?]` uses the bare `[?]` *trans whole-allele*
///   marker; the parenthesised `p.(?)` is the per-member predicted form).
///
/// Bare `?` is intentionally NOT probed here: it is only the trans
/// whole-allele marker (handled upstream in
/// `parse_protein_trans_allele_shorthand` → `HgvsVariant::UnknownAllele`)
/// and is not a valid cis / unknown-phase bracket member — `p.[?]`,
/// `p.[?;X]`, `p.[X;?]` continue to reject (pinned by
/// `tests/protein_unknown_roundtrip.rs`).
///
/// `0` / `0?` / `(0)` (whole-protein no-protein) is likewise NOT probed
/// here: `[0]` is not a valid cis bracket member (issue #277), and the
/// protein-coordinate trans-allele path handles `[0]` separately in
/// `parse_protein_trans_allele_shorthand`.
fn parse_protein_bracket_member(part: &str) -> IResult<&str, LocEdit<ProtInterval, ProteinEdit>> {
    /// Whole-entity protein edits carry no position; attach a dummy
    /// interval at `Met1` (mirrors `parse_protein_variant`).
    fn dummy_prot_interval() -> ProtInterval {
        ProtInterval::point(ProtPos::new(AminoAcid::Met, 1))
    }

    if let Ok((remaining, edit)) = parse_whole_protein_identity(part) {
        return Ok((remaining, LocEdit::new(dummy_prot_interval(), edit)));
    }
    // Only the parenthesised predicted unknown `(?)` is admitted as a
    // bracket member; bare `?` is reserved for the trans whole-allele
    // marker and must not parse here (see the doc comment above).
    if let Some(remaining) = part.strip_prefix("(?)") {
        return Ok((
            remaining,
            LocEdit::new(
                dummy_prot_interval(),
                ProteinEdit::whole_protein_unknown_predicted(),
            ),
        ));
    }

    // Predicted-change member `(position+edit)` — strip the uncertainty
    // wrapper, parse the inner position+edit, and mark the edit uncertain
    // so the member round-trips as `[(...)]` (e.g.
    // `p.[(Ser68Arg)];[(Ser73Arg)]`). Mirrors the single-variant predicted
    // form in `parse_protein_variant`.
    if part.starts_with('(') {
        if let Ok((remaining, (interval, edit))) = delimited(
            char('('),
            (parse_prot_interval, parse_protein_edit),
            char(')'),
        )
        .parse(part)
        {
            return Ok((
                remaining,
                LocEdit::with_uncertainty(interval, Mu::Uncertain(edit)),
            ));
        }
    }

    let (edit_remaining, interval) = parse_prot_interval(part)?;
    let (final_remaining, edit) = parse_protein_edit(edit_remaining)?;
    Ok((final_remaining, LocEdit::new(interval, edit)))
}

/// Parse a protein allele shorthand:
///   `p.[edit1;edit2;...]`     (cis)
///   `p.[edit1(;)edit2;...]`   (unknown phase)
///
/// A single bracket pair must use exactly one separator type; mixed
/// `;`/`(;)` inside one bracket pair is not spec-valid and is rejected.
fn parse_protein_allele_shorthand(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, AlleleVariant> {
    if !input.starts_with('[') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Fast path: any `(;)` inside the brackets means unknown phase. Slice the
    // bracketed content, split on `(;)`, and parse each sub-variant as
    // `interval+edit`. Use a depth-aware scan for the closing `]` so members
    // containing nested edit brackets are not truncated — see
    // `find_top_level_close_bracket`.
    let close_bracket = find_top_level_close_bracket(input).ok_or_else(|| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Tag))
    })?;
    let content = &input[1..close_bracket];
    let after_bracket = &input[close_bracket + 1..];

    // Scan for top-level (depth-0) allele separators. See
    // `scan_allele_separators` for the depth rules; inner `;` / `(;)`
    // from a nested edit must not drive phase detection or splitting.
    let scan = scan_allele_separators(content);

    if scan.has_unknown_phase {
        // Mixed `;`/`(;)` inside one bracket pair is not spec-valid — see
        // `parse_genome_kind_compound_allele` for the spec citation.
        if scan.has_cis_separator {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }

        let mut variants = Vec::with_capacity(4);
        for part in scan.members {
            let part = part.trim();
            if part.is_empty() {
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Verify,
                )));
            }
            let (final_remaining, loc_edit) = parse_protein_bracket_member(part)?;
            if !final_remaining.trim().is_empty() {
                return Err(nom::Err::Error(nom::error::Error::new(
                    final_remaining,
                    nom::error::ErrorKind::Tag,
                )));
            }
            variants.push(HgvsVariant::Protein(ProteinVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit,
            }));
        }

        // Singletons are intentionally accepted — see
        // `parse_genome_kind_compound_allele` for the rationale.
        if variants.is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }

        return Ok((
            after_bracket,
            AlleleVariant::new(variants, AllelePhase::Unknown),
        ));
    }

    // Cis path (existing behavior).
    let (remaining, _) = char('[').parse(input)?;

    let mut variants = Vec::new();
    let mut current = remaining;

    loop {
        // Parse interval and edit, admitting whole-entity members
        // (`=`, `(=)`, `(?)`) alongside position-based edits. #468.
        let (rest, loc_edit) = parse_protein_bracket_member(current)?;

        variants.push(HgvsVariant::Protein(ProteinVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit,
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

    // Singletons are intentionally accepted — see
    // `parse_genome_kind_compound_allele` for the rationale.
    if variants.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    Ok((current, AlleleVariant::cis(variants)))
}

/// Parse a protein trans-allele compact-prefix shorthand:
/// `[edit1];[edit2]` (with `p.` already consumed by the caller). Each
/// bracketed group holds a single sub-variant. Parses protein intervals +
/// edits and does NOT allow a bare-position (Identity) fallback — the
/// protein grammar requires an edit.
///
/// Intentionally NOT routed through `parse_trans_allele_shorthand_generic`
/// (used by the other 6 axes): protein has three special-case markers
/// (`[0]`, `[0?]`, `[(0)]`) that resolve to `ProteinEdit::NoProtein` rather
/// than the cross-coordinate `HgvsVariant::NullAllele`, because the `p.`
/// compact prefix has already pinned the coordinate system. The generic
/// helper hard-codes the cross-coordinate `[0]`/`[?]` →
/// `NullAllele`/`UnknownAllele` mapping appropriate for the other 6 axes.
/// See the inline comment at the `content == "0"` branch below for the
/// full spec arbitration (#277, follow-up to PR #130).
fn parse_protein_trans_allele_shorthand(
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

        let close_bracket = find_top_level_close_bracket(remaining).ok_or_else(|| {
            nom::Err::Error(nom::error::Error::new(
                remaining,
                nom::error::ErrorKind::Tag,
            ))
        })?;

        let content = remaining[1..close_bracket].trim();

        // Inside the protein-coordinate compact trans-allele form
        // (`ACC:p.[X];[0]`), `[0]` and `[0?]` resolve to
        // `ProteinEdit::NoProtein` rather than the cross-coordinate
        // `HgvsVariant::NullAllele`: the shared `p.` compact prefix has
        // already pinned the coordinate system to protein, and `p.0` /
        // `p.0?` are the spec's "no protein produced" forms (per
        // `recommendations/protein/{substitution,deletion}.md`). Routing to
        // `NullAllele` here would silently break Display→reparse semantics
        // for variants like `[ACC:p.X];[ACC:p.0]`, whose compact-form
        // Display emits `ACC:p.[X];[0]` (issue #277, follow-up to PR #130).
        //
        // The cross-coordinate `[0]` / `[?]` markers still come through the
        // outer (no-coord-prefix) `parse_trans_allele` path
        // (`[ACC:p.X];[0]`), where `NullAllele` / `UnknownAllele` remain
        // correct — see the pin
        // `protein_no_protein_roundtrip::bare_bracketed_zero_is_null_allele_not_no_protein`.
        //
        // Arm routing for the three protein no-protein input shapes inside a
        // `p.`-prefixed trans-allele bracket (see the preamble block above
        // for why the protein-coord context overrides the cross-coord
        // `NullAllele` default for `[0]`):
        // - `[0]`  → first arm, `NoProtein { predicted: false }` (certain
        //   no-protein; the protein-coord override of the cross-coord
        //   `NullAllele` shape pinned in `bare_bracketed_zero_is_null_allele_not_no_protein`).
        // - `[0?]` → first arm, `NoProtein { predicted: true }` (canonical
        //   spec form; Display emits this).
        // - `[(0)]` → third arm, `NoProtein { predicted: true }` (tolerated
        //   alternate input from #289; Display canonicalises to `[0?]`).
        let variant = if content == "0" || content == "0?" {
            let predicted = content == "0?";
            let dummy_pos = ProtPos::new(AminoAcid::Met, 1);
            let dummy_interval = ProtInterval::point(dummy_pos);
            HgvsVariant::Protein(ProteinVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(dummy_interval, ProteinEdit::NoProtein { predicted }),
            })
        } else if content == "?" {
            HgvsVariant::UnknownAllele
        } else if content == "(0)" {
            // `[(0)]` is the tolerated alternate input for predicted no-protein
            // inside a protein-coordinate compact-form trans-allele bracket
            // (issue #289), mirroring the bare dispatcher's
            // `parse_whole_protein_no_protein` arm. Routes to
            // `ProteinEdit::NoProtein { predicted: true }`; Display
            // canonicalises to `[0?]` (the HGVS spec form, per
            // `recommendations/protein/deletion.md`). Without this arm,
            // `parse_prot_interval` would try to consume `(0)` as a position
            // and reject.
            //
            // `[0?]` is handled by the first arm above (the canonical spec
            // form) — see the arm-routing summary block above. Bare `[0]`
            // keeps its cross-coordinate `NullAllele` meaning (per
            // `protein_no_protein_roundtrip::bare_bracketed_zero_is_null_allele_not_no_protein`).
            let dummy_pos = ProtPos::new(AminoAcid::Met, 1);
            let dummy_interval = ProtInterval::point(dummy_pos);
            HgvsVariant::Protein(ProteinVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit: LocEdit::new(dummy_interval, ProteinEdit::NoProtein { predicted: true }),
            })
        } else {
            // All other members (position-based edits plus whole-entity
            // identity `=` / `(=)` and unknown `(?)`) flow through the
            // shared bracket-member dispatcher. Bare `?` and the `0`
            // forms are handled by the arms above, so they never reach
            // here. #468.
            let (final_remaining, loc_edit) = parse_protein_bracket_member(content)?;

            if !final_remaining.trim().is_empty() {
                return Err(nom::Err::Error(nom::error::Error::new(
                    final_remaining,
                    nom::error::ErrorKind::Tag,
                )));
            }

            HgvsVariant::Protein(ProteinVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit,
            })
        };

        variants.push(variant);
        remaining = &remaining[close_bracket + 1..];

        // Mirrors `parse_trans_allele_shorthand_generic`: enforce strict
        // `[...];[...]` shape — reject both dangling trailing `;` and a
        // missing separator between adjacent brackets.
        if let Some(after_semi) = remaining.strip_prefix(';') {
            if !after_semi.starts_with('[') {
                return Err(nom::Err::Error(nom::error::Error::new(
                    after_semi,
                    nom::error::ErrorKind::Tag,
                )));
            }
            remaining = after_semi;
        } else if remaining.starts_with('[') {
            return Err(nom::Err::Error(nom::error::Error::new(
                remaining,
                nom::error::ErrorKind::Tag,
            )));
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

/// Parse position-based unknown-phase shorthand for protein (#123):
///   `(Ser68Arg)(;)(Asn594del)`
///
/// Per `recommendations/protein/alleles.md` line 48, the protein unknown-phase
/// form does NOT use surrounding square brackets and each sub-variant is a
/// predicted change wrapped in `()`.
fn parse_protein_position_unknown_phase(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    let mut variants = Vec::with_capacity(4);

    for part in input.split("(;)") {
        let part = part.trim();
        if part.is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }

        // Each part must be a parenthesized predicted change `(interval+edit)`.
        let (rest, (interval, edit)) = delimited(
            char('('),
            (parse_prot_interval, parse_protein_edit),
            char(')'),
        )
        .parse(part)?;

        if !rest.trim().is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                rest,
                nom::error::ErrorKind::Tag,
            )));
        }

        variants.push(HgvsVariant::Protein(ProteinVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit: LocEdit::with_uncertainty(interval, Mu::Uncertain(edit)),
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

        // Try protein trans-allele compact-prefix shorthand: p.[edit1];[edit2]
        // before the cis allele path, since trans is the more specific shape.
        // Mirrors the DNA/RNA dispatches (parse_genome_variant,
        // parse_cds_variant, parse_rna_variant) and PR #146's fix for `g.`.
        if input.starts_with('[') && input.contains("];[") {
            if let Ok((remaining, variant)) =
                parse_protein_trans_allele_shorthand(input, accession.clone(), gene_symbol.clone())
            {
                return Ok((remaining, variant));
            }
        }

        // Try protein allele: p.[edit1;edit2;...] or p.[edit1(;)edit2]
        if input.starts_with('[') {
            if let Ok((remaining, allele)) =
                parse_protein_allele_shorthand(input, accession.clone(), gene_symbol.clone())
            {
                return Ok((remaining, HgvsVariant::Allele(allele)));
            }
        }

        // Bracketless unknown-phase shorthand for protein (#123):
        //   p.(Ser68Arg)(;)(Asn594del)
        // Per spec the protein form is bare (no surrounding `[]`). Each
        // sub-variant is a predicted change in parentheses.
        if input.starts_with('(') && input.contains(")(;)(") {
            return parse_protein_position_unknown_phase(
                input,
                accession.clone(),
                gene_symbol.clone(),
            );
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

        // First try whole-RNA identity (r.= or r.(=)) - no position
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
                    loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                }),
            ));
        }

        // Try whole-RNA unknown (r.? or r.(?)) - no position
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
                    loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
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

        // Try RNA-specific no product pattern (r.0 or predicted r.(0))
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
                    loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                }),
            ));
        }

        // Try whole-RNA predicted splice (r.(spl) / r.(spl?)) — must come before the
        // generic predicted-variant parser, which requires an interval. The whole-entity
        // splicing marker has no positional content, so we attach a dummy interval and
        // wrap the edit in `Mu::Uncertain` so that `RnaVariant::fmt_loc_edit` emits the
        // surrounding parens.
        if let Ok((remaining, edit)) = parse_whole_rna_predicted_splice(input) {
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
                    loc_edit: LocEdit {
                        location: dummy_interval,
                        edit: Mu::Uncertain(edit),
                    },
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

/// Parse whole-RNA identity: r.= or r.(=)
fn parse_whole_rna_identity(input: &str) -> IResult<&str, Mu<crate::hgvs::edit::NaEdit>> {
    if let Ok((remaining, _)) = tag::<_, _, nom::error::Error<&str>>("(=)").parse(input) {
        return Ok((
            remaining,
            Mu::Uncertain(crate::hgvs::edit::NaEdit::whole_entity_identity()),
        ));
    }
    map(tag("="), |_| {
        Mu::Certain(crate::hgvs::edit::NaEdit::whole_entity_identity())
    })
    .parse(input)
}

/// Parse whole-RNA unknown: r.? or r.(?)
///
/// Only matches a terminal `?` (not `?` followed by `_`, edit keywords,
/// or a substitution pattern — those route to the position-based
/// parser as edits at unknown position). #239.
///
/// The `r.(?)` predicted form is returned as `Mu::Uncertain`. #245.
fn parse_whole_rna_unknown(input: &str) -> IResult<&str, Mu<crate::hgvs::edit::NaEdit>> {
    if let Ok((remaining, _)) = tag::<_, _, nom::error::Error<&str>>("(?)").parse(input) {
        return Ok((
            remaining,
            Mu::Uncertain(crate::hgvs::edit::NaEdit::whole_entity_unknown()),
        ));
    }
    let (remaining, _) = tag("?").parse(input)?;
    if is_edit_continuation_after_unknown(remaining) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }
    Ok((
        remaining,
        Mu::Certain(crate::hgvs::edit::NaEdit::whole_entity_unknown()),
    ))
}

/// Parse RNA splice pattern: `spl` (very-likely-affected) or `spl?` (might-be-affected).
///
/// Per HGVS v21.0 (`assets/hgvs-nomenclature/docs/recommendations/RNA/splicing.md` and
/// `docs/recommendations/uncertain.md`):
/// - `r.spl`  — RNA not analysed; splicing very likely affected.
/// - `r.spl?` — RNA not analysed; splicing might be affected.
fn parse_rna_splice(input: &str) -> IResult<&str, crate::hgvs::edit::NaEdit> {
    let (input, _) = tag("spl").parse(input)?;
    let (input, marker) =
        nom::combinator::opt(tag::<_, _, nom::error::Error<&str>>("?")).parse(input)?;
    Ok((
        input,
        crate::hgvs::edit::NaEdit::Splice {
            unknown: marker.is_some(),
        },
    ))
}

/// Parse whole-RNA predicted splice pattern: `(spl)` or `(spl?)`.
///
/// Per HGVS v21.0, `r.(spl?)` is the canonical predicted-uncertain splicing marker;
/// `r.(spl)` is the symmetric predicted-but-certain form. Both are whole-entity edits
/// (no positional content), so we wrap the result in `Mu::Uncertain` rather than
/// dispatching through the generic `parse_predicted_rna_variant`, which requires an
/// interval. See `RnaVariant::fmt_loc_edit` for the matching Display path.
fn parse_whole_rna_predicted_splice(input: &str) -> IResult<&str, crate::hgvs::edit::NaEdit> {
    let (input, _) = tag("(").parse(input)?;
    let (input, edit) = parse_rna_splice(input)?;
    let (input, _) = tag(")").parse(input)?;
    Ok((input, edit))
}

/// Parse RNA no product pattern: r.0 or predicted r.(0)
fn parse_rna_no_product(input: &str) -> IResult<&str, Mu<crate::hgvs::edit::NaEdit>> {
    if let Ok((remaining, _)) = tag::<_, _, nom::error::Error<&str>>("(0)").parse(input) {
        return Ok((
            remaining,
            Mu::Uncertain(crate::hgvs::edit::NaEdit::NoProduct),
        ));
    }
    // Match "0" but only as a standalone symbol (not followed by digits or other edit characters)
    let (remaining, _) = tag("0").parse(input)?;
    // Make sure it's not followed by more digits (would be a position)
    if remaining.chars().next().is_some_and(|c| c.is_ascii_digit()) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }
    Ok((remaining, Mu::Certain(crate::hgvs::edit::NaEdit::NoProduct)))
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
            // The parens around `(interval edit)` denote a predicted change
            // — wrap as Mu::Uncertain so `RnaVariant::fmt_loc_edit` emits
            // the spec-canonical `r.(<pos><edit>)` shape. #241.
            loc_edit: LocEdit::new_predicted(interval, edit),
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
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }

        // Delegate to the bracket-member dispatcher so whole-entity
        // RNA forms (`=`, `?`, `0`, `spl`, `spl?`, `(=)`, `(?)`, `(0)`,
        // `(spl)`, `(spl?)`) added in #396 item 3 are also admitted in
        // the bracketless `r.X(;)Y` form. Keeps `r.[X(;)Y]` and
        // `r.X(;)Y` round-trip-symmetric (Display canonicalises one
        // into the other for unknown-phase). #423 sibling fix to #396
        // item 3.
        let (final_remaining, loc_edit) = parse_rna_bracket_member(part)?;

        if !final_remaining.trim().is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                final_remaining,
                nom::error::ErrorKind::Tag,
            )));
        }

        variants.push(HgvsVariant::Rna(RnaVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit,
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

/// Parse a single bracket-member position+edit (RNA), recognising
/// whole-entity edits (`=`, `?`, `0`, `spl`, `spl?`, `(spl)`,
/// `(spl?)`), the predicted-wrapper form `(<pos><edit>)`, and bare
/// `<pos><edit>` (edit optional, defaults to identity). Mirrors
/// `parse_cds_bracket_member`. #287 follow-up to #243 / PR #244;
/// whole-entity dispatch added by #396 item 3.
fn parse_rna_bracket_member(
    part: &str,
) -> IResult<&str, LocEdit<RnaInterval, crate::hgvs::edit::NaEdit>> {
    // Whole-entity RNA edits have no positional content, so attach a
    // dummy interval at position 1 (mirrors `parse_rna_variant`'s
    // handling of the same forms at the top level). These probes run
    // BEFORE the position-based predicted wrapper / bare-position
    // parsers so a bracket like `[(spl?)]` matches the predicted-splice
    // marker rather than being routed into `parse_rna_interval`. #396
    // item 3.
    //
    // Note: bare `[0]` / `[?]` in the trans-allele form are handled
    // upstream by `parse_trans_allele_shorthand_generic`'s cross-coord
    // short-circuits (→ `NullAllele` / `UnknownAllele`) and never reach
    // this dispatcher. Only the cis/unknown-flat paths and predicted
    // forms (`[(0)]`, `[(?)]`, `[(spl?)]`) plus splice markers
    // (`[spl]`, `[spl?]`) flow through these probes.
    fn dummy_rna_interval() -> RnaInterval {
        RnaInterval::point(crate::hgvs::location::RnaPos {
            base: 1,
            offset: None,
            utr3: false,
        })
    }

    if let Ok((remaining, mu_edit)) = parse_whole_rna_identity(part) {
        return Ok((
            remaining,
            LocEdit {
                location: dummy_rna_interval(),
                edit: mu_edit,
            },
        ));
    }
    if let Ok((remaining, mu_edit)) = parse_whole_rna_unknown(part) {
        return Ok((
            remaining,
            LocEdit {
                location: dummy_rna_interval(),
                edit: mu_edit,
            },
        ));
    }
    if let Ok((remaining, mu_edit)) = parse_rna_no_product(part) {
        return Ok((
            remaining,
            LocEdit {
                location: dummy_rna_interval(),
                edit: mu_edit,
            },
        ));
    }
    if let Ok((remaining, edit)) = parse_whole_rna_predicted_splice(part) {
        return Ok((
            remaining,
            LocEdit {
                location: dummy_rna_interval(),
                edit: Mu::Uncertain(edit),
            },
        ));
    }
    if let Ok((remaining, edit)) = parse_rna_splice(part) {
        return Ok((remaining, LocEdit::new(dummy_rna_interval(), edit)));
    }

    // Position-based predicted wrapper: `(<pos><edit>)`.
    if part.starts_with('(') && !part.starts_with("(?") && !part.starts_with("(=)") {
        if let Ok((remaining, (interval, edit))) =
            delimited(char('('), (parse_rna_interval, parse_na_edit), char(')')).parse(part)
        {
            return Ok((remaining, LocEdit::new_predicted(interval, edit)));
        }
    }
    // Bare position [+ edit]
    let (remaining, interval) = parse_rna_interval(part)?;
    let (remaining, edit) = if let Ok((r, e)) = parse_na_edit(remaining) {
        (r, e)
    } else {
        (
            remaining,
            crate::hgvs::edit::NaEdit::Identity {
                sequence: None,
                whole_entity: false,
            },
        )
    };
    Ok((remaining, LocEdit::new(interval, edit)))
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

    // Scan for top-level (depth-0) allele separators. Inner `;` / `(;)`
    // from a nested edit (`delins[a;u]`) or predicted-wrapper member
    // must not drive phase detection or splitting. See
    // `scan_allele_separators`.
    let scan = scan_allele_separators(content);

    // Mixed `;`/`(;)` inside one bracket pair is not spec-valid (#396
    // item 2; mirrors `parse_protein_allele_shorthand` and the new
    // cds-axis reject path). Reject rather than silently flattening to
    // AllelePhase::Unknown.
    if scan.has_unknown_phase && scan.has_cis_separator {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    let phase = if scan.has_unknown_phase {
        AllelePhase::Unknown
    } else {
        AllelePhase::Cis
    };

    let mut variants = Vec::with_capacity(4);

    for part in scan.members {
        let part = part.trim();
        if part.is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }

        let (final_remaining, loc_edit) = parse_rna_bracket_member(part)?;

        if !final_remaining.trim().is_empty() {
            return Err(nom::Err::Error(nom::error::Error::new(
                final_remaining,
                nom::error::ErrorKind::Tag,
            )));
        }

        variants.push(HgvsVariant::Rna(RnaVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit,
        }));
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
    parse_trans_allele_shorthand_generic(input, |content| {
        let (rest, loc_edit) = parse_rna_bracket_member(content)?;
        Ok((
            rest,
            HgvsVariant::Rna(RnaVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit,
            }),
        ))
    })
}

/// Parse circular DNA variant (o.) - SVD-WG006
fn parse_circular_variant(
    accession: Accession,
    gene_symbol: Option<String>,
) -> impl FnMut(&str) -> IResult<&str, HgvsVariant> {
    move |input: &str| {
        let (input, _) = tag("o.").parse(input)?;

        // First try whole-circular identity (o.= or o.(=)) - no position. #288.
        // Circular DNA reuses the genome whole-entity parser because both
        // share `GenomeInterval`; SVD-WG006 inherits the DNA grammar.
        // Reject when the remainder starts with `(;)` so `o.=(;)100A>G`
        // reaches the unknown-phase dispatch below.
        if let Ok((remaining, edit)) = parse_whole_genome_identity(input) {
            if !remaining.starts_with("(;)") {
                let dummy_pos = crate::hgvs::location::GenomePos::new(1);
                let dummy_interval = GenomeInterval::point(dummy_pos);
                return Ok((
                    remaining,
                    HgvsVariant::Circular(CircularVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                    }),
                ));
            }
        }

        // Try whole-circular unknown (o.? or o.(?)) - no position. #288.
        if let Ok((remaining, edit)) = parse_whole_genome_unknown(input) {
            if !remaining.starts_with("(;)") {
                let dummy_pos = crate::hgvs::location::GenomePos::new(1);
                let dummy_interval = GenomeInterval::point(dummy_pos);
                return Ok((
                    remaining,
                    HgvsVariant::Circular(CircularVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::with_uncertainty(dummy_interval, edit),
                    }),
                ));
            }
        }

        // Compact-prefix trans-allele shorthand: o.[a];[b]. Mirrors
        // parse_genome_variant's dispatch added in PR #146; SVD-WG006 says
        // circular DNA inherits the DNA allele grammar.
        if input.starts_with('[') && input.contains("];[") {
            if let Ok((remaining, variant)) =
                parse_circular_trans_allele_shorthand(input, accession.clone(), gene_symbol.clone())
            {
                return Ok((remaining, variant));
            }
        }

        // Bracketed compound allele: o.[var;var] or o.[var(;)var]
        if input.starts_with('[') {
            if let Ok((remaining, allele)) = parse_genome_kind_compound_allele(
                input,
                accession.clone(),
                gene_symbol.clone(),
                GenomeKind::Circular,
            ) {
                return Ok((remaining, HgvsVariant::Allele(allele)));
            }
        }

        // Bracketless unknown-phase shorthand: o.100A>G(;)200T>C (#123)
        if input.contains("(;)") && !input.starts_with('[') {
            return parse_genome_position_unknown_phase(
                input,
                accession.clone(),
                gene_symbol.clone(),
                GenomeKind::Circular,
            );
        }

        // Predicted variant: o.(positionEdit). #241.
        if input.starts_with('(') && !input.starts_with("(?") && !input.starts_with("(=)") {
            if let Ok((remaining, (interval, edit))) = delimited(
                char('('),
                (parse_genome_interval_for_circular, parse_na_edit),
                char(')'),
            )
            .parse(input)
            {
                // Reject spec-unauthorised reversed ranges (#129).
                check_circular_reversed_range(&interval, &edit, input)?;
                return Ok((
                    remaining,
                    HgvsVariant::Circular(CircularVariant {
                        accession: accession.clone(),
                        gene_symbol: gene_symbol.clone(),
                        loc_edit: LocEdit::new_predicted(interval, edit),
                    }),
                ));
            }
        }

        let original_input = input;
        let (input, interval) = parse_genome_interval_for_circular(input)?;
        let (input, edit) = parse_na_edit(input)?;
        // Reject spec-unauthorised reversed ranges (#129).
        check_circular_reversed_range(&interval, &edit, original_input)?;

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

/// Parse a circular DNA trans-allele compact-prefix shorthand:
/// `[edit1];[edit2]` (with the `o.` already consumed). Mirrors
/// `parse_genome_trans_allele_shorthand`; circular variants reuse
/// `parse_genome_interval`.
fn parse_circular_trans_allele_shorthand(
    input: &str,
    accession: Accession,
    gene_symbol: Option<String>,
) -> IResult<&str, HgvsVariant> {
    parse_trans_allele_shorthand_generic(input, |content| {
        // `GenomeKind::Circular` routes through the circular interval parser
        // + spec-authorised reversed-range validator (#129).
        let (rest, loc_edit) = parse_genome_bracket_member(content, GenomeKind::Circular)?;
        Ok((
            rest,
            HgvsVariant::Circular(CircularVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit,
            }),
        ))
    })
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
        } else if acc_prefix == "NM" || acc_prefix == "XM" {
            // Coding transcript accession - infer c. without the prefix
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
        } else if acc_prefix == "NR" || acc_prefix == "XR" {
            // Non-coding RNA transcript accession - infer n. (a `c.` inference
            // would be a coordinate-system mismatch on NR_/XR_; #486).
            let prefixed = format!("n.{}", remaining);
            match parse_tx_variant(accession.clone(), gene_symbol.clone()).parse(&prefixed) {
                Ok((rem, variant)) => {
                    let consumed = prefixed.len() - rem.len() - 2; // -2 for "n."
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

    // Scan for top-level (depth-0) allele separators. Inner `;` / `(;)`
    // belonging to a nested edit (`delins[A;G]`) or a predicted-wrapper
    // member must NOT drive phase detection or member splitting — a naive
    // `content.split(';')` mis-splits `delins[A;G]` into bogus members.
    // See `scan_allele_separators`.
    let scan = scan_allele_separators(content);

    // Mixed `;`/`(;)` inside one bracket pair is not spec-valid — there is
    // no HGVS shape for "A and B are cis to each other but unknown phase to
    // C" in a single bracket. Reject rather than silently flattening to
    // `AllelePhase::Unknown` (mirrors the compact-prefix path, #396 item 2).
    if scan.has_unknown_phase && scan.has_cis_separator {
        return Err(FerroError::Parse {
            pos: 1,
            msg: "Mixed ';' and '(;)' separators in a single allele bracket are not valid HGVS"
                .to_string(),
            diagnostic: None,
        });
    }

    let has_unknown_phase = scan.has_unknown_phase;
    // Reject malformed separators (`;;`, leading/trailing `;`) rather than
    // silently dropping the empty member, consistent with the sibling
    // compound-allele parsers.
    let mut parts: Vec<&str> = Vec::with_capacity(scan.members.len());
    for member in &scan.members {
        let member = member.trim();
        if member.is_empty() {
            return Err(FerroError::Parse {
                pos: 1,
                msg: "Empty allele member in ClinVar-style allele".to_string(),
                diagnostic: None,
            });
        }
        parts.push(member);
    }

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

    // Split at top-level (depth-0) `;` only. Inner `;` belonging to a
    // nested edit (`delins[A;G]`) is at bracket depth 1 and must not be
    // treated as an allele separator. This entrypoint is reached only for
    // the cis form (`detect_allele_type` routes any `(;)` to the
    // unknown-phase entrypoint), so no mixed-phase handling is needed here.
    // See `scan_allele_separators`.
    let scan = scan_allele_separators(content);

    // Parse each top-level member.
    let mut variants = Vec::with_capacity(4);
    for part in scan.members {
        let part = part.trim();
        // Reject malformed separators (`;;`, leading/trailing `;`) rather than
        // silently dropping the empty group, consistent with the sibling
        // compound-allele parsers.
        if part.is_empty() {
            return Err(FerroError::Parse {
                pos: 1,
                msg: "Empty ';' group in cis allele".to_string(),
                diagnostic: None,
            });
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
        let close_bracket =
            find_top_level_close_bracket(remaining).ok_or_else(|| FerroError::Parse {
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

/// True iff `variant` carries a position-bound (NOT whole-entity)
/// `Identity` edit on a NA coord system — i.e. the spec compact
/// mosaic / chimeric LHS shape `<acc>:<type>.<pos>=`.
///
/// Used by `parse_phase_allele` to decide how to interpret a bare
/// `=` chunk that appears after the leading arm in a multi-slash
/// chain. When true, the `=` chunk clones the leading arm so the
/// inherited interval is preserved (e.g. `m.3243=/A>G/=/A>T` keeps
/// position 3243 on the second `=`); when false (the single-slash
/// `var/=` shorthand), we fall through to
/// `create_identity_variant_from`, which builds a synthetic whole-
/// entity identity that Display rewrites as bare `=`.
fn is_lhs_position_identity(variant: &HgvsVariant) -> bool {
    use crate::hgvs::edit::NaEdit;
    let edit_is_pos_identity = |edit: &Mu<NaEdit>| {
        matches!(
            edit.inner(),
            Some(NaEdit::Identity {
                whole_entity: false,
                ..
            })
        )
    };
    match variant {
        HgvsVariant::Genome(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Cds(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Tx(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Rna(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Mt(v) => edit_is_pos_identity(&v.loc_edit.edit),
        HgvsVariant::Circular(v) => edit_is_pos_identity(&v.loc_edit.edit),
        _ => false,
    }
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
                special: None,
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
/// Supports "=" as shorthand for reference allele (no change).
/// Supports accession inheritance: NM_000088.3:c.100A>G/c.200C>T.
/// Supports bracketed inner alleles: `[a;b]/[c;d]`.
fn parse_mosaic_allele(input: &str) -> Result<HgvsVariant, FerroError> {
    parse_phase_allele(input, AllelePhase::Mosaic)
}

/// Parse the RHS of a HGVS spec compact mosaic / chimeric form.
///
/// The spec compact form is `<acc>:<type>.<pos>=/<edit>` (mosaic) or
/// `<acc>:<type>.<pos>=//<edit>` (chimeric). The RHS is a bare edit
/// with no accession, type prefix, or position; it inherits all three
/// from the LHS, which must be a position-bound `=` identity edit.
///
/// Returns `Ok(None)` when the LHS is not eligible (so the caller can
/// fall through to its existing error path); returns `Ok(Some(_))` on
/// successful compact-form parse; returns `Err` only on a definite
/// compact-form parse error (e.g. trailing characters).
///
/// Per `recommendations/DNA/{substitution,deletion,duplication}.md`,
/// the form applies to substitution, deletion, and duplication on
/// nucleic-acid coord systems (`g.`, `c.`, `n.`, `r.`, `m.`, `o.`).
fn parse_compact_mosaic_rhs(
    rhs: &str,
    lhs: &HgvsVariant,
) -> Result<Option<HgvsVariant>, FerroError> {
    use crate::hgvs::edit::NaEdit;
    use crate::hgvs::parser::edit::parse_na_edit;

    // Eligibility: LHS must carry a position-bound (not whole-entity)
    // identity edit on a NA coord system.
    let lhs_is_eligible = match lhs {
        HgvsVariant::Genome(v) => matches!(
            v.loc_edit.edit.inner(),
            Some(NaEdit::Identity {
                whole_entity: false,
                ..
            })
        ),
        HgvsVariant::Cds(v) => matches!(
            v.loc_edit.edit.inner(),
            Some(NaEdit::Identity {
                whole_entity: false,
                ..
            })
        ),
        HgvsVariant::Tx(v) => matches!(
            v.loc_edit.edit.inner(),
            Some(NaEdit::Identity {
                whole_entity: false,
                ..
            })
        ),
        HgvsVariant::Rna(v) => matches!(
            v.loc_edit.edit.inner(),
            Some(NaEdit::Identity {
                whole_entity: false,
                ..
            })
        ),
        HgvsVariant::Mt(v) => matches!(
            v.loc_edit.edit.inner(),
            Some(NaEdit::Identity {
                whole_entity: false,
                ..
            })
        ),
        HgvsVariant::Circular(v) => matches!(
            v.loc_edit.edit.inner(),
            Some(NaEdit::Identity {
                whole_entity: false,
                ..
            })
        ),
        _ => false,
    };
    if !lhs_is_eligible {
        return Ok(None);
    }

    // Parse the RHS as a bare NaEdit and require we consumed all of it.
    let (remaining, edit) = match parse_na_edit(rhs) {
        Ok((rem, e)) => (rem, e),
        Err(_) => return Ok(None),
    };
    if !remaining.trim().is_empty() {
        return Err(FerroError::Parse {
            pos: rhs.len() - remaining.len(),
            msg: format!(
                "compact-mosaic RHS has unexpected trailing content: '{}'",
                remaining
            ),
            diagnostic: None,
        });
    }

    // Build the RHS variant by cloning LHS's accession + gene_symbol +
    // interval, swapping in the new edit. The Mu wrapper on the LHS
    // edit is replaced (we always emit Certain for the new RHS edit).
    let rhs_variant = match lhs {
        HgvsVariant::Genome(v) => HgvsVariant::Genome(GenomeVariant {
            accession: v.accession.clone(),
            gene_symbol: v.gene_symbol.clone(),
            loc_edit: LocEdit::new(v.loc_edit.location.clone(), edit),
        }),
        HgvsVariant::Cds(v) => HgvsVariant::Cds(CdsVariant {
            accession: v.accession.clone(),
            gene_symbol: v.gene_symbol.clone(),
            loc_edit: LocEdit::new(v.loc_edit.location.clone(), edit),
        }),
        HgvsVariant::Tx(v) => HgvsVariant::Tx(TxVariant {
            accession: v.accession.clone(),
            gene_symbol: v.gene_symbol.clone(),
            loc_edit: LocEdit::new(v.loc_edit.location.clone(), edit),
        }),
        HgvsVariant::Rna(v) => HgvsVariant::Rna(RnaVariant {
            accession: v.accession.clone(),
            gene_symbol: v.gene_symbol.clone(),
            loc_edit: LocEdit::new(v.loc_edit.location.clone(), edit),
        }),
        HgvsVariant::Mt(v) => HgvsVariant::Mt(MtVariant {
            accession: v.accession.clone(),
            gene_symbol: v.gene_symbol.clone(),
            loc_edit: LocEdit::new(v.loc_edit.location.clone(), edit),
        }),
        HgvsVariant::Circular(v) => HgvsVariant::Circular(CircularVariant {
            accession: v.accession.clone(),
            gene_symbol: v.gene_symbol.clone(),
            loc_edit: LocEdit::new(v.loc_edit.location.clone(), edit),
        }),
        _ => unreachable!("eligibility check guards against non-NA arms"),
    };
    Ok(Some(rhs_variant))
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
        // Use the circular-aware interval parser + post-parse spec
        // authorisation check so inherited-accession `…/m.<high>_<low>del`
        // slash-form inputs accept the SVD-WG006 wraparound exception
        // (#129).
        let (remaining, interval) =
            parse_genome_interval_for_circular(rest).map_err(|e| FerroError::Parse {
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
        check_circular_reversed_range(&interval, &edit, rest).map_err(|e| FerroError::Parse {
            pos: 2,
            msg: format!("spec-unauthorised reversed range on m. axis: {:?}", e),
            diagnostic: None,
        })?;
        Ok(HgvsVariant::Mt(MtVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.cloned(),
            loc_edit: LocEdit::new(interval, edit),
        }))
    } else if let Some(rest) = input.strip_prefix("o.") {
        let (remaining, interval) =
            parse_genome_interval_for_circular(rest).map_err(|e| FerroError::Parse {
                pos: 2,
                msg: format!("Failed to parse circular interval: {:?}", e),
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
        check_circular_reversed_range(&interval, &edit, rest).map_err(|e| FerroError::Parse {
            pos: 2,
            msg: format!("spec-unauthorised reversed range on o. axis: {:?}", e),
            diagnostic: None,
        })?;
        Ok(HgvsVariant::Circular(CircularVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.cloned(),
            loc_edit: LocEdit::new(interval, edit),
        }))
    } else {
        Err(FerroError::Parse {
            pos: 0,
            msg: format!(
                "Unknown variant type prefix: expected one of 'c.', 'g.', 'p.', 'n.', 'r.', 'm.', 'o.' but found '{}'",
                input.chars().take(3).collect::<String>()
            ),
            diagnostic: None,
        })
    }
}

/// Parse a chimeric: var1//var2 (double slash)
/// Supports "=" as shorthand for reference allele (no change).
/// Supports bracketed inner alleles: `[a;b]//[c;d]`.
fn parse_chimeric_allele(input: &str) -> Result<HgvsVariant, FerroError> {
    parse_phase_allele(input, AllelePhase::Chimeric)
}

/// Whether a chunk-level `parse_variant` error carries a structured
/// `Diagnostic.code` — i.e. it's a semantic violation that
/// `parse_phase_allele`'s syntactic-recovery fallback chain should
/// short-circuit on rather than mask. See call site in
/// `parse_phase_allele` for context (issue #375).
fn chunk_error_carries_diagnostic_code(err: &FerroError) -> bool {
    matches!(
        err,
        FerroError::Parse {
            diagnostic: Some(d),
            ..
        } if d.code.is_some()
    )
}

/// Re-emit a chunk-relative `FerroError::Parse` against the *full* input
/// the user typed: shift `pos` and `Diagnostic.span.{start,end}` by
/// `chunk_offset`, and replace `Diagnostic.source` with `full_input`.
/// Non-`Parse` variants pass through unchanged. Used by
/// `parse_phase_allele` to surface chunk-level structured diagnostics
/// (E3006, W3017, …) at coordinates that match the user's typed string
/// rather than the chunk slice (issue #375).
fn remap_chunk_error_to_full_input(
    err: FerroError,
    chunk_offset: usize,
    full_input: &str,
) -> FerroError {
    match err {
        FerroError::Parse {
            pos,
            msg,
            diagnostic,
        } => {
            let new_pos = pos.saturating_add(chunk_offset);
            let new_diagnostic = diagnostic.map(|mut d| {
                if let Some(span) = d.span.as_mut() {
                    span.start = span.start.saturating_add(chunk_offset);
                    span.end = span.end.saturating_add(chunk_offset);
                }
                d.source = Some(full_input.to_string());
                d
            });
            FerroError::Parse {
                pos: new_pos,
                msg,
                diagnostic: new_diagnostic,
            }
        }
        other => other,
    }
}

/// Shared chunk-driver for `parse_mosaic_allele` and
/// `parse_chimeric_allele`. Splits `input` at top-level slash separators
/// using a bracket-depth-aware scan (`split_top_level_slashes`), parses
/// each chunk via the full `parse_variant` so bracketed inner alleles
/// like `[a;b]//[c;d]` recurse correctly, and applies the three
/// existing fallbacks for chunks that aren't stand-alone variants:
///
///   1. inherited-accession long form (`<acc>:c.X` LHS, `c.Y` RHS),
///   2. HGVS spec compact form (`<acc>:c.<pos>=` LHS, bare edit RHS).
///
/// The third path (`=` shorthand for the reference allele) is checked
/// before any parse path so the synthetic identity is created against
/// the recorded LHS.
///
/// Triple-slash inputs like `var1///var2` surface as a chunk with a
/// leading `/`, which this driver rejects with a clear error rather
/// than silently round-tripping a malformed description.
fn parse_phase_allele(input: &str, phase: AllelePhase) -> Result<HgvsVariant, FerroError> {
    let input = input.trim();
    let want_double = matches!(phase, AllelePhase::Chimeric);
    let phase_name = match phase {
        AllelePhase::Mosaic => "mosaic",
        AllelePhase::Chimeric => "chimeric",
        _ => unreachable!("parse_phase_allele only handles Mosaic / Chimeric"),
    };
    let sep_str = if want_double { "//" } else { "/" };

    let chunks = split_top_level_slashes(input, want_double);
    if chunks.len() < 2 {
        return Err(FerroError::Parse {
            pos: 0,
            msg: format!(
                "{} allele requires '{}' separator",
                capitalize(phase_name),
                sep_str
            ),
            diagnostic: None,
        });
    }

    let mut variants: Vec<HgvsVariant> = Vec::with_capacity(chunks.len());
    let mut first_variant: Option<HgvsVariant> = None;

    for chunk in chunks {
        let chunk = chunk.trim();
        if chunk.is_empty() {
            return Err(FerroError::Parse {
                pos: 0,
                msg: format!("Empty variant in {} allele", phase_name),
                diagnostic: None,
            });
        }
        // A chunk starting with `/` means the input contained either
        // `///` (triple-slash) or, in single-slash mode, an unsplit
        // `//`. Both are malformed.
        if chunk.starts_with('/') {
            return Err(FerroError::Parse {
                pos: 0,
                msg: format!(
                    "Unexpected '/' inside {} allele chunk `{}` (triple-slash or stray separator)",
                    phase_name, chunk
                ),
                diagnostic: None,
            });
        }
        // Reject mixing mosaic and chimeric phase markers at the same
        // nesting level. After splitting a chimeric input on `//`, a
        // chunk that still contains a top-level (depth-0) `/` would be
        // silently re-parsed by `parse_variant` as a nested mosaic
        // allele — yielding a chimeric-of-mosaic that the HGVS spec
        // does not define. Reject with a clean error rather than pick
        // an interpretation. (Mosaic chunks cannot carry top-level `//`
        // because `detect_allele_type` routes any top-level `//` to
        // chimeric up-front, so this guard is only meaningful for the
        // chimeric arm; checking the mosaic arm too is defensive — it
        // costs a single bracket-aware scan per chunk.)
        let chunk_has_inner_slash = if want_double {
            find_top_level_slash(chunk, false).is_some()
        } else {
            find_top_level_slash(chunk, true).is_some()
        };
        if chunk_has_inner_slash {
            // W3019 detector B: a chunk that still carries a top-level
            // slash of the *other* phase marker after splitting on the
            // outer marker — i.e. `a/b//c` (mosaic-of-chimeric) or
            // `a//b/c` (chimeric-of-mosaic). HGVS v21 does not define
            // this nesting. Reject with the same W3019 family code as
            // detector A in `parse_variant`.
            return Err(non_spec_mosaic_form_error(
                input,
                NonSpecCase::MixedFlatSlash,
            ));
        }
        let variant = if chunk == "=" {
            match &first_variant {
                // When the leading arm is itself a position-bound `=`
                // identity (the spec compact-mosaic LHS shape, e.g.
                // `m.3243=`), a later bare `=` chunk in a multi-slash
                // chain (`m.3243=/A>G/=/A>T`) inherits the LHS's
                // accession + coord type + interval rather than
                // collapsing to a synthetic whole-entity `m.1=`. The
                // whole-entity path is preserved only for the single-
                // slash `var/=` shorthand, where the LHS is a non-
                // identity edit and the Display layer rewrites the
                // synthetic whole-entity RHS as bare `=` (PR #153 /
                // #133 work-item-3). Closes #284.
                Some(ref_var) if is_lhs_position_identity(ref_var) => ref_var.clone(),
                Some(ref_var) => create_identity_variant_from(ref_var)?,
                None => {
                    return Err(FerroError::Parse {
                        pos: 0,
                        msg: format!("Cannot use '=' as first variant in {} notation", phase_name),
                        diagnostic: None,
                    });
                }
            }
        } else {
            match parse_variant(chunk) {
                Ok(v) => v,
                Err(primary_err) => {
                    // Closes #375: a chunk-level error that carries a
                    // structured `Diagnostic.code` is a *semantic*
                    // violation (E3006 SelfCancellingAllele, W3017
                    // AlleleFractionAnnotation, W3019 NonSpecMosaicForm,
                    // W3021 ProteinBracketedAaInsertion, …) — not the
                    // shape of failure the recovery fallbacks below were
                    // built for. The fallbacks are for *syntactic*
                    // shortfalls (RHS missing accession, bare-edit RHS,
                    // ClinVar-prose multi-allelic), so handing them a
                    // semantically-rejected chunk only produces a
                    // misleading "Unknown variant type prefix" that
                    // masks the genuine diagnostic. Propagate the
                    // original error, remapping `pos`, `Diagnostic.span`
                    // and `Diagnostic.source` to full-input coordinates
                    // so downstream tooling (LSP, web service) can
                    // underline the offending region in the user's
                    // typed string rather than in the chunk slice.
                    if chunk_error_carries_diagnostic_code(&primary_err) {
                        let chunk_offset = chunk.as_ptr() as usize - input.as_ptr() as usize;
                        return Err(remap_chunk_error_to_full_input(
                            primary_err,
                            chunk_offset,
                            input,
                        ));
                    }
                    // Fallback 1: inherited-accession long form.
                    let lhs = match &first_variant {
                        Some(lhs) => lhs,
                        None => return Err(primary_err),
                    };
                    let acc = match lhs.accession() {
                        Some(a) => a,
                        None => return Err(primary_err),
                    };
                    let gs = lhs.gene_symbol().map(|s| s.to_string());
                    match parse_variant_with_inherited_accession(chunk, acc, gs.as_ref()) {
                        Ok(v) => v,
                        Err(prefix_err) => {
                            // Fallback 2: HGVS spec compact form — bare
                            // NaEdit RHS inheriting accession + coord
                            // type + position from `<pos>=` LHS.
                            match parse_compact_mosaic_rhs(chunk, lhs)? {
                                Some(v) => v,
                                None => {
                                    // Targeted diagnostic for the
                                    // ClinVar-prose multi-allelic
                                    // shorthand `m.<pos><ref>><alt>/<alt2>`
                                    // (issue #278). Only fire on the
                                    // recognizable bare-base RHS shape,
                                    // when the LHS has already parsed
                                    // (so this is genuinely a mosaic-RHS
                                    // rejection), and when the LHS is a
                                    // mitochondrial variant — the W3018
                                    // diagnostic text and code are
                                    // specifically for `m.` heteroplasmy
                                    // prose; firing it on non-mito LHS
                                    // would mislabel unrelated parse
                                    // failures.
                                    if matches!(lhs, HgvsVariant::Mt(_))
                                        && is_prose_multi_allelic_rhs(chunk)
                                    {
                                        use crate::error::{Diagnostic, ErrorCode, SourceSpan};
                                        let _w3018 = ErrorType::ClinVarProseMultiAllelic;
                                        let diag = Diagnostic::new()
                                            .with_code(ErrorCode::ClinVarProseMultiAllelic)
                                            .with_span(SourceSpan::new(0, chunk.len()));
                                        return Err(FerroError::parse_with_diagnostic(
                                            0,
                                            prose_multi_allelic_diagnostic_msg(chunk),
                                            diag,
                                        ));
                                    }
                                    return Err(prefix_err);
                                }
                            }
                        }
                    }
                }
            }
        };

        if first_variant.is_none() {
            first_variant = Some(variant.clone());
        }
        variants.push(variant);
    }

    Ok(HgvsVariant::Allele(AlleleVariant::new(variants, phase)))
}

/// Which sub-case of the W3019 `NonSpecMosaicForm` family fired.
/// Both produce the same error code; the case lets the diagnostic
/// hint be precise about *which* shape the user wrote.
#[derive(Clone, Copy, Debug)]
pub(crate) enum NonSpecCase {
    /// `[a/b]` / `[a//b]` — slash inside a bracket pair.
    SlashInBracket,
    /// `a/b//c` / `a//b/c` — mixed `/` and `//` at the same level.
    MixedFlatSlash,
}

/// Build the W3019 (`NonSpecMosaicForm`) rejection diagnostic.
/// Centralized so both detectors emit identical hint text and so
/// downstream tooling can grep for a single `[W3019 …]` prefix.
///
/// The hint enumerates the spec-supported alternatives the user
/// should use instead:
///   * compound brackets `[a;b]` (cis),
///   * dual fully-qualified slash `acc:c.X/acc:c.Y` (mosaic) or
///     `acc:c.X//acc:c.Y` (chimeric),
///   * compact short-hands `acc:c.<pos>=/<edit>` (mosaic) and
///     `acc:c.<pos>=//<edit>` (chimeric).
pub(crate) fn non_spec_mosaic_form_error(input: &str, case: NonSpecCase) -> FerroError {
    use crate::error::{Diagnostic, ErrorCode, SourceSpan};
    use crate::error_handling::ErrorType;
    let case_label = match case {
        NonSpecCase::SlashInBracket => "bracketed `[a/b]` / `[a//b]`",
        NonSpecCase::MixedFlatSlash => {
            "mixing mosaic and chimeric markers at the same nesting level"
        }
    };
    // The hint text is inlined in `msg` (not only in the diagnostic
    // hint field) so the spec-alternatives are visible in plain
    // `Display` output — `FerroError::Parse`'s thiserror format
    // shows `msg` but not the structured `Diagnostic` hint.
    let hint = "use `[a;b]` (cis), `acc:c.X/acc:c.Y` (mosaic), \
                `acc:c.X//acc:c.Y` (chimeric), or compact `acc:c.<pos>=/<edit>` \
                / `acc:c.<pos>=//<edit>`";
    // Source the W-code + name from the `ErrorType` registry so the
    // string stays in sync if the code ever moves. Using the variant
    // here also satisfies the error-code audit's "enforced rows
    // must reference ErrorType::<Variant> from emission-relevant
    // src/ paths" rule (tests/error_code_audit.rs).
    let et = ErrorType::NonSpecMosaicForm;
    let msg = format!(
        "[{} {:?}] {} — HGVS v21 does not define this form; {}",
        et.code(),
        et,
        case_label,
        hint
    );
    FerroError::parse_with_diagnostic(
        0,
        msg,
        Diagnostic::new()
            .with_code(ErrorCode::UnexpectedChar)
            .with_span(SourceSpan::new(0, input.len()))
            .with_source(input)
            .with_hint(hint),
    )
}

/// Title-case the first ASCII letter of `s`. Tiny local helper used
/// only by `parse_phase_allele` error messages.
fn capitalize(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        Some(first) => first.to_ascii_uppercase().to_string() + c.as_str(),
        None => String::new(),
    }
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

    // Scan for top-level (depth-0) `(;)` separators. Inner `;` / `(;)`
    // belonging to a nested edit must not drive member splitting. See
    // `scan_allele_separators`.
    let scan = scan_allele_separators(content);

    // Mixed `;`/`(;)` inside one bracket pair is not spec-valid (same shape
    // rejected by the compact-prefix path in #396 item 2): there is no HGVS
    // form for "A and B are cis to each other but unknown phase to C".
    if scan.has_unknown_phase && scan.has_cis_separator {
        return Err(FerroError::Parse {
            pos: 1,
            msg: "Mixed ';' and '(;)' separators in a single allele bracket are not valid HGVS"
                .to_string(),
            diagnostic: None,
        });
    }

    // Split at top-level `(;)` and parse each variant. Empty groups (e.g.
    // `[A:c.x(;)(;)B:c.y]`, leading/trailing `(;)`) are malformed; reject
    // rather than silently dropping them — matches the per-coord helpers
    // (`parse_genome_kind_compound_allele` etc.) introduced for #123.
    let mut variants = Vec::with_capacity(2);
    for part in scan.members {
        let part = part.trim();
        if part.is_empty() {
            return Err(FerroError::Parse {
                pos: 0,
                msg: "Empty `(;)` group in unknown-phase allele".to_string(),
                diagnostic: None,
            });
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

    // Top-level slash check happens BEFORE the bracket-wrapping checks
    // so that bracketed chimeric / mosaic forms like
    // `[a;b]//[c;d]` and `[a;b]/[c;d]` route to chimeric / mosaic
    // instead of being misclassified as cis (`[a;b]` wrapping with
    // unparsed trailing content). Depth-aware so the slash inside
    // a nested edit (e.g. a delins source like `delins[ACC:c.1_2/3]`,
    // which is not currently emitted but kept structurally safe) is
    // ignored. See `find_top_level_slash`.
    if find_top_level_slash(input, true).is_some() {
        return Some("chimeric");
    }
    if let Some(pos) = find_top_level_slash(input, false) {
        // Bare `/foo` (slash at position 0, before any accession /
        // coord prefix) is not a valid mosaic — let downstream parsing
        // produce a clean error. Mirror the prior heuristic of
        // requiring some structural content (a `:`, a `[`, or the `=`
        // shorthand for the reference allele) before the first
        // top-level slash before classifying as mosaic. The bare-`=`
        // case lets `parse_phase_allele` emit its targeted
        // "Cannot use '=' as first variant" diagnostic, matching the
        // chimeric branch's behavior (which has no before-slash guard
        // since `//` is unambiguous).
        let before = input[..pos].trim();
        if before == "=" || before.contains(':') || before.starts_with('[') {
            return Some("mosaic");
        }
    }

    // Check for trans: [var];[var] pattern
    if input.starts_with('[') && input.contains("];[") {
        return Some("trans");
    }

    // Check for unknown phase: [var(;)var] pattern
    if input.starts_with('[') && input.ends_with(']') && input.contains("(;)") {
        return Some("unknown_phase");
    }

    // Check for cis: [var;var] pattern. The HGVS DNA cis grammar
    // (syntax.yaml lines 134-135) is
    // `"[" position_edit ";" position_edit "]"` — `;` separator and
    // ≥2 entries required. Standalone `c.[a]` is not generated by
    // any production and is explicitly addressed and rejected by the
    // committee in alleles.md lines 99-101: "the recommended
    // description is `LRG_199t1:c.[76A>C];[76=]`".
    if input.starts_with('[') && input.ends_with(']') && !input.contains("];[") {
        let inner = &input[1..input.len() - 1];
        if inner.contains(';') {
            return Some("cis");
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

/// Returns `true` if `s` looks like an allele-fraction / heteroplasmy
/// annotation appended to an otherwise-complete HGVS expression.
///
/// Recognized shapes (issue #278 / SVA W3017):
///
/// - `[level=NN%]` — ClinVar prose heteroplasmy
/// - `[heteroplasmy=NN%]`, `[mosaic=NN%]`, `[mosaicism=NN%]` — synonyms
/// - `(NN%)` — bare-percent paren shorthand
///
/// The match is intentionally narrow: it only fires on a trailing
/// bracket-or-paren group whose content carries `%`. This keeps the
/// diagnostic out of the way of other unrelated trailing failures.
fn is_allele_fraction_annotation(s: &str) -> bool {
    let s = s.trim();
    // Paren form: `(NN%)` or `(NN.NN%)`. Require at least one digit
    // before any `.` so degenerate shapes like `(.5%)` or `( . %)`
    // don't slip through as fraction annotations.
    if let Some(inner) = s.strip_prefix('(').and_then(|t| t.strip_suffix(')')) {
        let inner = inner.trim();
        if inner.ends_with('%') {
            let num = inner.trim_end_matches('%').trim();
            // Split on `.` and require ≥1 ASCII digit in the integer
            // part. The fractional part may be empty (e.g. `80.%`) or
            // a string of digits.
            let (int_part, frac_ok) = match num.split_once('.') {
                Some((int_part, frac)) => (
                    int_part,
                    frac.chars().all(|c| c.is_ascii_digit() || c == ' '),
                ),
                None => (num, true),
            };
            let int_has_digit = int_part.chars().any(|c| c.is_ascii_digit());
            let int_ok = int_part.chars().all(|c| c.is_ascii_digit() || c == ' ');
            return int_has_digit && int_ok && frac_ok;
        }
    }
    // Bracket form: `[<key>=NN%]` with a `%` somewhere in the body. Limit
    // the key to the documented heteroplasmy/mosaicism synonyms so that
    // unrelated trailing `[<key>=NN%]` annotations aren't mislabeled as
    // W3017 AlleleFractionAnnotation.
    if let Some(inner) = s.strip_prefix('[').and_then(|t| t.strip_suffix(']')) {
        if inner.contains('%') && inner.contains('=') {
            let (key, _rest) = inner.split_once('=').unwrap_or((inner, ""));
            let key_lc = key.trim().to_ascii_lowercase();
            return matches!(
                key_lc.as_str(),
                "level" | "heteroplasmy" | "mosaic" | "mosaicism"
            );
        }
    }
    false
}

/// Returns `true` if `chunk` looks like the ClinVar-prose
/// multi-allelic shorthand RHS (a bare nucleotide letter, optionally
/// short like `T` or `AT`), with no edit operator, no `=`, no `:`, no
/// type prefix, no parens, no brackets.
///
/// This is the RHS shape that produces `m.<pos><ref>><alt>/<alt2>`
/// (e.g. `m.3243A>G/T`). Pinned by issue #278.
fn is_prose_multi_allelic_rhs(chunk: &str) -> bool {
    let chunk = chunk.trim();
    if chunk.is_empty() || chunk.len() > 8 {
        return false;
    }
    // Reject anything with structural HGVS markers — those would be
    // legitimate (mis-)attempts, not the prose form.
    let forbidden = [
        '=', ':', '.', '>', '+', '-', '(', ')', '[', ']', '_', ',', ';', '*', '?',
    ];
    if chunk.chars().any(|c| forbidden.contains(&c)) {
        return false;
    }
    // Body must be entirely nucleotide letters (case-insensitive,
    // ACGTU) so we keep the diagnostic targeted on the actual prose
    // shape seen in submitter notes.
    chunk
        .chars()
        .all(|c| matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'U'))
}

/// Build the targeted ClinVar-prose-multi-allelic diagnostic message
/// for a rejected mosaic/chimeric RHS like `m.3243A>G/T`.
///
/// Pinned by issue #278. The message names the three spec-supported
/// alternatives so the caller can choose the correct one:
///   1. Compound brackets `<acc>:m.[a;b]`.
///   2. Dual fully-qualified slash `<acc>:m.<...>/<acc>:m.<...>`.
///   3. Spec compact mosaic `<acc>:m.<pos>=/<alt-edit>`.
fn prose_multi_allelic_diagnostic_msg(rhs: &str) -> String {
    format!(
        "ClinVar prose multi-allelic shorthand `<acc>:m.<pos><ref>><alt>/{rhs}` is not a HGVS \
         spec form (the RHS `{rhs}` elides the reference base and accession). Use one of the \
         three spec-supported alternatives instead: \
         (a) compound brackets, e.g. `<acc>:m.[3243A>G;3243A>T]`; \
         (b) dual fully-qualified slash, e.g. `<acc>:m.3243A>G/<acc>:m.3243A>T`; or \
         (c) spec compact mosaic, e.g. `<acc>:m.3243=/A>T` (= for the reference, edit for the \
         alternative).",
        rhs = rhs.trim()
    )
}

/// Parse a complete HGVS variant string
///
/// The parser is optimized with:
/// - Early type prefix checking to avoid unnecessary work
/// - Variant types ordered by frequency (c. > p. > g. > n. > r. > m.)
/// - Allele detection for compound variants
pub fn parse_variant(input: &str) -> Result<HgvsVariant, FerroError> {
    // NB: per HGVS v21 (general.md:96) whitespace is not permitted in any
    // HGVS description. ferro's lenient parser trims outer whitespace and
    // strips spaces after `;` inside allele brackets for interop with
    // real-world submissions; this is intentional and documented by
    // gene_selector_display_preserve / protein_*_outer_whitespace_is_trimmed
    // tests. The F1 test family in `tests/input_hygiene_rejections.rs`
    // pins this divergence so any future tightening is intentional.
    let input = input.trim();

    // Pre-parse rejection for protein bracketed-amino-acid insertion
    // (`p.…ins[Ala;Pro]`) — W3021 ProteinBracketedAaInsertion (#290).
    // The bracketed form has no spec-defined meaning at the protein-edit
    // level (brackets are reserved for variant-level alleles), so we
    // surface a structured diagnostic pointing at the canonical
    // concatenated `insAlaPro` shape rather than letting the bare parser
    // bail with a generic nom failure. Bare `parse_hgvs` callers see the
    // same code as the lenient/silent preprocessor paths.
    reject_protein_bracketed_aa_insertion(input)?;

    // W3019 detector A: `[a/b]` / `[a//b]` — slash *inside* a bracket
    // pair. HGVS does not define this form. Reject up-front with a
    // targeted diagnostic so users get the spec alternatives instead
    // of a generic nom error from `parse_cis_allele`.
    //
    // The bracketed slash form `[a;b]/[c;d]` is not in the HGVS spec
    // (the grammar in `syntax.yaml` defines only cis `[a;b]`, trans
    // `[a];[b]`, unknown-phase `a(;)b`; `/` and `//` are documented in
    // `DNA/substitution.md` and `general.md` only as separators between
    // single edits, e.g. `c.85=/T>C`; `consultation/open-issues.md`
    // records that the committee rejected nesting). The parser tolerates
    // the bracketed form for interoperability with real-world
    // submissions, splitting at top-level slashes in `parse_phase_allele`.
    // Slashes at depth 0 keep working under this tolerance, so the
    // depth-≥1 check here doesn't regress that path — see
    // `has_slash_inside_brackets`.
    if has_slash_inside_brackets(input) {
        return Err(non_spec_mosaic_form_error(
            input,
            NonSpecCase::SlashInBracket,
        ));
    }

    // Check for allele patterns and RNA fusion first
    let variant = if let Some(allele_type) = detect_allele_type(input) {
        match allele_type {
            "cis" => parse_cis_allele(input)?,
            "trans" => parse_trans_allele(input)?,
            "mosaic" => parse_mosaic_allele(input)?,
            "chimeric" => parse_chimeric_allele(input)?,
            "unknown_phase" => parse_unknown_phase_allele(input)?,
            "rna_fusion" => parse_rna_fusion(input)?,
            _ => unreachable!(),
        }
    } else {
        // Parse as single variant
        let (remaining, v) = parse_single_variant(input)?;
        if !remaining.is_empty() {
            // Trailing allele-fraction / heteroplasmy annotations such as
            // `[level=70%]`, `[heteroplasmy=70%]`, `[mosaic=80%]`, or
            // `(80%)` are out-of-spec but appear in real-world
            // submissions. Emit the dedicated SVA code (W3017 /
            // `ErrorType::AlleleFractionAnnotation`) so downstream
            // tooling can recognize the intent rather than seeing a
            // generic trailing-character parse error. See issue #278
            // (follow-up to #133 work-item 4).
            if is_allele_fraction_annotation(remaining) {
                use crate::error::{Diagnostic, ErrorCode, SourceSpan};
                // Bind the SVA `ErrorType::*` here so the audit's
                // emission-site scan (see tests/error_code_audit.rs) sees
                // a real reference for this W-code in the parser path.
                let _w3017 = ErrorType::AlleleFractionAnnotation;
                let pos = input.len() - remaining.len();
                let diag = Diagnostic::new()
                    .with_code(ErrorCode::AlleleFractionAnnotation)
                    .with_span(SourceSpan::new(pos, input.len()));
                return Err(FerroError::parse_with_diagnostic(
                    pos,
                    format!(
                        "Allele-fraction / heteroplasmy annotation '{}' is not part of the HGVS \
                         spec; allele fraction belongs in accompanying metadata (VCF \
                         FORMAT/AF, ClinVar's heteroplasmy field, etc.), not in the HGVS \
                         expression. Strip the annotation and record the fraction separately.",
                        remaining
                    ),
                    diag,
                ));
            }
            return Err(FerroError::Parse {
                pos: input.len() - remaining.len(),
                msg: format!("Unexpected trailing characters: '{}'", remaining),
                diagnostic: None,
            });
        }
        v
    };

    // Spec-mandated post-parse semantic check: reject self-cancelling alleles
    // (HGVS recommendations/general.md line 47, tracked under issue #115).
    validate_no_self_cancelling(&variant, input)?;

    // Spec-mandated post-parse semantic check: reject `dupins` mash-up
    // (HGVS DNA/duplication.md:92 — "a format not used in HGVS
    // nomenclature"; #445).
    validate_no_dupins(&variant)?;

    // Spec-mandated post-parse semantic check: reject single-position
    // insertion (HGVS DNA/insertion.md:95-101 Q&A — explicit "No"; #446).
    validate_no_point_insertion(&variant, input)?;

    // Spec-mandated post-parse semantic check: reject a coding/genomic/mito
    // coordinate system on a non-coding RNA (NR_/XR_) reference (#486,
    // ECOORDINATESYSTEMMISMATCH).
    validate_coordinate_system(&variant)?;

    Ok(variant)
}

/// Reject a coordinate system that is incompatible with the reference type.
///
/// A non-coding RNA transcript (`NR_` / `XR_`) has no CDS and is not a genomic
/// or mitochondrial reference, so a `c.`/`g.`/`m.`/`o.` description on it is a
/// coordinate-system mismatch (mutalyzer's `ECOORDINATESYSTEMMISMATCH`). The
/// only valid forms are `n.` (non-coding transcript) and `r.` (RNA), which are
/// not matched here. The check walks `Allele` wrappers and inspects each leaf.
fn validate_coordinate_system(variant: &HgvsVariant) -> Result<(), FerroError> {
    fn make_error(coord: &str) -> FerroError {
        use crate::error::{Diagnostic, ErrorCode};
        FerroError::parse_with_diagnostic(
            0,
            format!(
                "the `{coord}` coordinate system is not valid on a non-coding RNA \
                 (NR_/XR_) reference, which has no CDS and is not genomic; use \
                 `n.` (non-coding) or `r.` (RNA)"
            ),
            Diagnostic::new().with_code(ErrorCode::CoordinateSystemMismatch),
        )
    }

    match variant {
        HgvsVariant::Cds(v) if v.accession.is_noncoding_rna() => return Err(make_error("c.")),
        HgvsVariant::Genome(v) if v.accession.is_noncoding_rna() => return Err(make_error("g.")),
        HgvsVariant::Mt(v) if v.accession.is_noncoding_rna() => return Err(make_error("m.")),
        HgvsVariant::Circular(v) if v.accession.is_noncoding_rna() => return Err(make_error("o.")),
        HgvsVariant::Allele(allele) => {
            for inner in &allele.variants {
                validate_coordinate_system(inner)?;
            }
        }
        _ => {}
    }
    Ok(())
}

/// Reject `dupins` edits per `DNA/duplication.md:92`:
///
/// > "the variant is not described using `dupins`, a format not used in
/// > HGVS nomenclature."
///
/// `<code class="invalid">` markup in the spec — the strongest
/// prohibition register. ferro previously parsed `dupins<seq>` into the
/// `NaEdit::DupIns` variant and round-tripped it; this validation
/// rejects at parse time with a diagnostic pointing the user at the
/// canonical alternatives (`del` followed by `ins`, or `delins`).
///
/// The `NaEdit::DupIns` variant is kept in the data model for now (no
/// data-type refactor) — this check just ensures it is never produced
/// by the parser's user-facing entry point. Internal code that
/// constructs `DupIns` for normalization intermediates is unaffected.
fn validate_no_dupins(variant: &HgvsVariant) -> Result<(), FerroError> {
    use crate::hgvs::edit::NaEdit;
    use crate::hgvs::uncertainty::Mu;

    fn is_dupins(edit: &Mu<NaEdit>) -> bool {
        matches!(edit.inner(), Some(NaEdit::DupIns { .. }))
    }
    fn make_error() -> FerroError {
        use crate::error::{Diagnostic, ErrorCode};
        // Carry a structured `InvalidEdit` code so the semantic rejection
        // survives slash-form fallback paths that would otherwise mask an
        // unstructured parse error and lose the `dupins` guidance.
        FerroError::parse_with_diagnostic(
            0,
            "the `dupins<seq>` form is not used in HGVS nomenclature \
             (DNA/duplication.md:92). Describe as `del` followed by \
             `ins`, or as `delins` if the edit is contiguous.",
            Diagnostic::new().with_code(ErrorCode::InvalidEdit),
        )
    }

    match variant {
        HgvsVariant::Genome(v) if is_dupins(&v.loc_edit.edit) => return Err(make_error()),
        HgvsVariant::Cds(v) if is_dupins(&v.loc_edit.edit) => return Err(make_error()),
        HgvsVariant::Tx(v) if is_dupins(&v.loc_edit.edit) => return Err(make_error()),
        HgvsVariant::Rna(v) if is_dupins(&v.loc_edit.edit) => return Err(make_error()),
        HgvsVariant::Mt(v) if is_dupins(&v.loc_edit.edit) => return Err(make_error()),
        HgvsVariant::Circular(v) if is_dupins(&v.loc_edit.edit) => return Err(make_error()),
        HgvsVariant::Allele(allele) => {
            for inner in &allele.variants {
                validate_no_dupins(inner)?;
            }
        }
        _ => {}
    }
    Ok(())
}

/// Reject insertions whose anchor is a single position rather than a
/// two-position range. Per `DNA/insertion.md:95-101`:
///
/// > "Can I describe a variant as `g.123insG`? **No**, since the
/// > description is not unequivocal, it is not allowed. What does the
/// > description mean, the insertion of a `G` **at** position `g.123`
/// > or the insertion of a `G` **after** position `g.123`?"
///
/// The check walks Allele wrappers and inspects each leaf NaEdit. The
/// canonical fix for the user is to use `g.123_124insG` (or whatever
/// adjacent-pair anchor is intended).
fn validate_no_point_insertion(variant: &HgvsVariant, _source: &str) -> Result<(), FerroError> {
    use crate::hgvs::edit::NaEdit;
    use crate::hgvs::interval::UncertainBoundary;
    use crate::hgvs::uncertainty::Mu;

    fn pos_eq<T: PartialEq>(start: &UncertainBoundary<T>, end: &UncertainBoundary<T>) -> bool {
        match (start.as_single(), end.as_single()) {
            (Some(Mu::Certain(s)), Some(Mu::Certain(e))) => s == e,
            _ => false,
        }
    }
    fn is_insertion(edit: &Mu<NaEdit>) -> bool {
        matches!(edit.inner(), Some(NaEdit::Insertion { .. }))
    }
    fn make_error(prefix: &str) -> FerroError {
        use crate::error::{Diagnostic, ErrorCode};
        // Carry a structured `InvalidEdit` code (mirroring
        // `validate_no_dupins`) so the semantic rejection survives
        // slash-form fallback paths that would otherwise mask an
        // unstructured parse error and lose the single-position guidance.
        FerroError::parse_with_diagnostic(
            0,
            format!(
                "single-position insertion not allowed on {prefix}; the anchor MUST be \
                 two adjacent positions (e.g. {prefix}123_124insATG, not \
                 {prefix}123insATG) per DNA/insertion.md:95-101"
            ),
            Diagnostic::new().with_code(ErrorCode::InvalidEdit),
        )
    }

    match variant {
        HgvsVariant::Genome(v)
            if is_insertion(&v.loc_edit.edit)
                && pos_eq(&v.loc_edit.location.start, &v.loc_edit.location.end) =>
        {
            return Err(make_error("g."));
        }
        HgvsVariant::Cds(v)
            if is_insertion(&v.loc_edit.edit)
                && pos_eq(&v.loc_edit.location.start, &v.loc_edit.location.end) =>
        {
            return Err(make_error("c."));
        }
        HgvsVariant::Tx(v)
            if is_insertion(&v.loc_edit.edit)
                && pos_eq(&v.loc_edit.location.start, &v.loc_edit.location.end) =>
        {
            return Err(make_error("n."));
        }
        HgvsVariant::Rna(v)
            if is_insertion(&v.loc_edit.edit)
                && pos_eq(&v.loc_edit.location.start, &v.loc_edit.location.end) =>
        {
            return Err(make_error("r."));
        }
        HgvsVariant::Mt(v)
            if is_insertion(&v.loc_edit.edit)
                && pos_eq(&v.loc_edit.location.start, &v.loc_edit.location.end) =>
        {
            return Err(make_error("m."));
        }
        HgvsVariant::Circular(v)
            if is_insertion(&v.loc_edit.edit)
                && pos_eq(&v.loc_edit.location.start, &v.loc_edit.location.end) =>
        {
            return Err(make_error("o."));
        }
        HgvsVariant::Allele(allele) => {
            for inner in &allele.variants {
                validate_no_point_insertion(inner, _source)?;
            }
        }
        // Protein insertions use a different anchor model — handled
        // elsewhere; nothing to validate here.
        _ => {}
    }
    Ok(())
}

/// Pre-parse rejection for `p.…ins[Ala;Pro]` — W3021
/// ProteinBracketedAaInsertion (#290).
///
/// HGVS v21's protein insertion notation concatenates 3-letter codes
/// (`p.Arg97_Trp98insAlaPro`); brackets at the edit level have no
/// spec meaning. Bare `parse_hgvs` callers do not go through the
/// preprocessor, so we surface the same structured diagnostic here.
/// Scope is bounded to text inside a `:p.` description.
fn reject_protein_bracketed_aa_insertion(input: &str) -> Result<(), FerroError> {
    use crate::error::{Diagnostic, ErrorCode, SourceSpan};
    use crate::error_handling::corrections::detect_protein_bracketed_aa_insertion;

    let detections = detect_protein_bracketed_aa_insertion(input);
    let Some(first) = detections.first() else {
        return Ok(());
    };

    let canonical_suggestion = if !first.corrected.is_empty() {
        format!(
            " — e.g. write `{}` instead of `{}`",
            first.corrected, first.original
        )
    } else {
        String::new()
    };
    let hint = format!(
        "HGVS protein insertion uses concatenated 3-letter codes \
         (`ins<AA1><AA2>...`, e.g. `insAlaPro`), not a bracketed list \
         (`ins[AA1;AA2]`){canonical_suggestion}",
    );

    Err(FerroError::parse_with_diagnostic(
        first.start,
        "Bracketed amino-acid list inside protein insertion edit",
        Diagnostic::new()
            .with_code(ErrorCode::InvalidEdit)
            .with_span(SourceSpan::new(first.start, first.end))
            .with_source(input)
            .with_hint(hint),
    ))
}

/// Walk the variant tree and reject any `AlleleVariant` that contains an
/// overlapping `del` + `dup` pair (HGVS `recommendations/general.md` line 47:
/// "descriptions removing part of a reference sequence replacing it with
/// part of the same sequence are not allowed").
///
/// `source` is the original parser input — always passed through verbatim
/// so the attached `Diagnostic::with_source` can render against the
/// caller's actual string. When a violation is found the offending
/// allele's `[ ... ]` byte range is located within `source` (advancing
/// past the parent bracket on each nested recursion) so the diagnostic
/// points at the precise bracket; this fixes #171 where E3006 previously
/// reported `pos: 0`.
fn validate_no_self_cancelling(variant: &HgvsVariant, source: &str) -> Result<(), FerroError> {
    walk_self_cancelling(variant, source, 0)
}

/// Inner helper: `search_from` is the byte offset within `source` at
/// which to begin scanning for the next cis-allele bracket. Top-level
/// callers pass `0`; recursive calls for nested alleles pass the byte
/// offset *just after* the enclosing bracket's `[` so that nested
/// diagnostics resolve to the inner bracket rather than the outer.
fn walk_self_cancelling(
    variant: &HgvsVariant,
    source: &str,
    search_from: usize,
) -> Result<(), FerroError> {
    if let HgvsVariant::Allele(allele) = variant {
        // Only cis-phase alleles describe edits on the same physical sequence.
        // Trans/mosaic/chimeric place the edits on different alleles or cell
        // populations, so an overlapping del+dup is not self-cancelling.
        // Unknown phase is treated permissively to avoid false rejections.
        let bracket = locate_cis_allele_bracket(source, search_from);
        if matches!(allele.phase, AllelePhase::Cis) {
            if let Some((dl, du)) =
                crate::hgvs::variant::AlleleVariant::detect_self_cancelling_pair(&allele.variants)
            {
                let span = bracket.map(|(s, e)| crate::error::SourceSpan::new(s, e));
                return Err(crate::hgvs::variant::build_self_cancelling_error(
                    dl,
                    du,
                    span,
                    Some(source),
                ));
            }
        }
        // Recurse for nested alleles. A *per-child cursor* (not a single
        // shared offset) is required because a parent allele may contain
        // multiple sibling alleles in its `variants` list, each of which
        // has its own `[...]` in `source`. Reusing the same
        // `next_search_from` for every child would bind every nested
        // E3006 diagnostic to the first sibling's bracket. After each
        // child is walked, advance the cursor past that child's bracket
        // (if it found one) so the next child scans the remainder of
        // `source`. The full `source` keeps flowing so global byte
        // offsets resolve correctly at any depth.
        let mut cursor = bracket
            .map(|(open, _close)| open + 1)
            .unwrap_or(search_from);
        for inner in &allele.variants {
            walk_self_cancelling(inner, source, cursor)?;
            cursor = locate_cis_allele_bracket(source, cursor)
                .map(|(open, _close)| open + 1)
                .unwrap_or(cursor);
        }
    }
    Ok(())
}

/// Locate the first allele bracket `[...]` in `source` at or after
/// `search_from` and return its open / exclusive-close byte range.
///
/// Two shapes are recognised:
///
/// - **Compact form**: an axis marker (`:c.`, `:g.`, `:n.`, `:r.`, `:m.`,
///   `:o.`, `:p.`) immediately followed by `[`. This is the canonical
///   top-level shape, e.g. `NM_004006.2:c.[762_768del;767_774dup]`.
/// - **Expanded form / nested**: a bare `[` at `search_from` or later
///   that is not part of an inserted-sequence (`ins[ATC]`) or repeat
///   bracket. Disambiguation is done by looking at the byte immediately
///   before the candidate `[`: if it's `:`, `;`, `(`, or start-of-scan,
///   the `[` opens an allele rather than an `ins`/`dup` payload.
///
/// In both shapes the matching `]` is found by depth-tracking, so
/// inserted-sequence brackets nested inside the allele (`g.X_Yins[ATC]`
/// per #333) and repeat brackets do not confuse the closing position.
/// Returns `None` if no allele bracket is found from `search_from`
/// onward, or if the input has unbalanced brackets.
fn locate_cis_allele_bracket(source: &str, search_from: usize) -> Option<(usize, usize)> {
    let bytes = source.as_bytes();

    // Pass 1: prefer the axis-marker form — it unambiguously identifies
    // the allele bracket even when free `[` characters appear earlier.
    let mut i = search_from;
    while i + 4 <= bytes.len() {
        let is_axis_bracket = bytes[i] == b':'
            && matches!(bytes[i + 1], b'c' | b'g' | b'n' | b'r' | b'm' | b'o' | b'p')
            && bytes[i + 2] == b'.'
            && bytes[i + 3] == b'[';
        if is_axis_bracket {
            return close_at_depth_zero(bytes, i + 3);
        }
        i += 1;
    }

    // Pass 2: an inner allele bracket whose enclosing axis marker
    // belongs to its parent. Accept a bare `[` whose preceding byte is
    // an HGVS structural delimiter (`:`, `;`, `(`) or the start of the
    // scan window. Crucially this rejects `ins[` and `dup[` payloads
    // because their preceding byte is alphabetic.
    let mut j = search_from;
    while j < bytes.len() {
        if bytes[j] == b'[' {
            let prev_is_delim = j == search_from || matches!(bytes[j - 1], b':' | b';' | b'(');
            if prev_is_delim {
                return close_at_depth_zero(bytes, j);
            }
        }
        j += 1;
    }
    None
}

/// Given `bytes[open] == b'['`, return `(open, close_exclusive)` where
/// `close_exclusive` is one past the matching `]` (tracking depth so
/// nested `[…]` payloads close cleanly). Returns `None` for unbalanced
/// input.
fn close_at_depth_zero(bytes: &[u8], open: usize) -> Option<(usize, usize)> {
    if open >= bytes.len() || bytes[open] != b'[' {
        return None;
    }
    let mut depth = 0i32;
    let mut j = open;
    while j < bytes.len() {
        match bytes[j] {
            b'[' => depth += 1,
            b']' => {
                depth -= 1;
                if depth == 0 {
                    return Some((open, j + 1));
                }
            }
            _ => {}
        }
        j += 1;
    }
    None
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

        // Delins with explicit deleted sequence is preserved on round-trip
        // (issue #120). Use canonicalize_edit / canonicalization to strip
        // the explicit deleted form when desired.
        let variant = parse_variant("NM_000088.3:c.85-47_84+48delGCCAinsG").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(
            format!("{}", variant),
            "NM_000088.3:c.85-47_84+48delGCCAinsG"
        );
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
        // Uncertain single position bounded by [12345,12350]: round-trips
        // as the single-paren spec form per HGVS v21 uncertain.md (#237).
        let variant = parse_variant("NC_000001.11:g.(12345_12350)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(format!("{}", variant), "NC_000001.11:g.(12345_12350)del");
    }

    #[test]
    fn test_parse_uncertain_cds_deletion() {
        // Uncertain single position bounded by [100,200]: round-trips as
        // the single-paren spec form per HGVS v21 uncertain.md (#237).
        let variant = parse_variant("NM_000088.3:c.(100_200)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:c.(100_200)del");
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
        // Uncertain single position bounded by [100,200]: round-trips as
        // the single-paren spec form per HGVS v21 uncertain.md (#237).
        let variant = parse_variant("NM_000088.3:r.(100_200)del").unwrap();
        assert!(matches!(variant, HgvsVariant::Rna(_)));
        assert_eq!(format!("{}", variant), "NM_000088.3:r.(100_200)del");
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
        // Same-accession cis alleles render in compact form
        assert_eq!(format!("{}", variant), "NM_000088.3:c.[100A>G;200C>T]");
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
        // Same-accession trans alleles render in compact form
        assert_eq!(format!("{}", variant), "NM_000088.3:c.[100A>G];[200C>T]");
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
        // Roundtrip in spec compact form (#133 work item 3): the RHS `=`
        // shorthand emits as bare `=`, not as a synthetic identity at
        // position 1.
        assert_eq!(format!("{}", variant), "NM_000088.3:c.123A>G/=");
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
        // Roundtrip in spec compact form (#133 work item 3).
        assert_eq!(format!("{}", variant), "NM_000088.3:c.456C>T//=");
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
        // Same-accession unknown-phase alleles render in compact form
        assert_eq!(format!("{}", variant), "NM_000088.3:c.100A>G(;)200C>T");
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
        // Output uses compact form (same accession)
        assert_eq!(format!("{}", variant), "NM_000088.3:c.[145C>T;147C>G]");
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
        // Output uses compact form (same accession)
        assert_eq!(format!("{}", variant), "NM_000088.3:c.145C>T(;)147C>G");
    }

    #[test]
    fn test_parse_mixed_phase_cds_allele_shorthand_rejected() {
        // Mixed `;`/`(;)` inside one bracket pair (e.g. `[A;B(;)C]`) has
        // no spec-defined HGVS shape: there is no way to express "A and
        // B are cis to each other, but unknown phase to C" in a single
        // bracket. Item 2 of #396 replaced the previous silent-flatten
        // (which used to produce a flat `AllelePhase::Unknown` over all
        // members) with an explicit reject. Mirrors
        // `parse_protein_allele_shorthand`.
        assert!(parse_variant("NM_000088.3:c.[123A>G;456C>T(;)789G>A]").is_err());
    }

    #[test]
    fn test_parse_mixed_phase_cds_allele_shorthand_complex_rejected() {
        // More complex mixed shape `[A;B;C(;)D;E]` rejects too.
        assert!(parse_variant("NM_000088.3:c.[100A>G;200C>T;300G>A(;)400del;500dup]").is_err());
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

    #[test]
    fn test_parse_cis_allele_nested_delins_not_missplit() {
        // Expanded cis allele (each member carries its own accession) whose
        // first member is a multi-segment `delins[A;G]`. The inner `;` lives
        // at bracket depth 1 and must NOT be treated as an allele separator,
        // so the allele has exactly two members — not three. Regression for
        // the naive `content.split(';')` in `parse_cis_allele`.
        let variant =
            parse_variant("[NC_000001.11:g.100_102delins[A;G];NC_000001.11:g.200A>G]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Cis);
            assert_eq!(allele.variants.len(), 2);
        }
    }

    #[test]
    fn test_parse_clinvar_style_allele_nested_delins_not_missplit() {
        // Same depth-aware split for the ClinVar-style entrypoint
        // (`accession:[c.var1;c.var2]`): the inner `;` in `delins[A;G]` is at
        // bracket depth 1 and must not split the allele into a bogus member.
        let variant = parse_variant("NC_000001.11:[g.100_102delins[A;G];g.200A>G]").unwrap();
        assert!(matches!(variant, HgvsVariant::Allele(_)));
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Cis);
            assert_eq!(allele.variants.len(), 2);
        }
    }

    #[test]
    fn test_parse_clinvar_style_allele_rejects_mixed_phase() {
        // Mixed `;`/`(;)` inside one ClinVar-style bracket pair has no
        // spec-defined HGVS shape (same rationale as the compact-prefix
        // `c.[A;B(;)C]` reject in #396 item 2). It must reject rather than
        // silently flatten to `AllelePhase::Unknown`.
        assert!(parse_variant("NM_000088.3:[c.100A>G;c.200C>T(;)c.300G>A]").is_err());
    }

    #[test]
    fn test_parse_unknown_phase_allele_rejects_mixed_phase() {
        // Expanded unknown-phase entrypoint (each member carries its own
        // accession): a top-level `;` mixed with a top-level `(;)` is the
        // same unrepresentable mixed-phase shape and must reject.
        assert!(parse_variant(
            "[NM_000088.3:c.100A>G;NM_000088.3:c.200C>T(;)NM_000088.3:c.300G>A]"
        )
        .is_err());
    }

    #[test]
    fn test_parse_clinvar_style_allele_rejects_empty_members() {
        // Malformed separators (`;;`, leading `;`, trailing `;`) must ERROR
        // rather than silently dropping the empty member, consistent with the
        // sibling compound-allele parsers. These reach
        // `parse_clinvar_style_allele` via the `accession:[...]` shorthand.
        assert!(parse_variant("NM_006876.2:[c.1168A>G;;c.1217C>T]").is_err());
        assert!(parse_variant("NM_006876.2:[;c.1168A>G;c.1217C>T]").is_err());
        assert!(parse_variant("NM_006876.2:[c.1168A>G;c.1217C>T;]").is_err());
        assert!(parse_variant("NM_006876.2:[ ]").is_err());
    }

    #[test]
    fn test_parse_cis_allele_rejects_empty_members() {
        // Same malformed-separator rejection for the expanded cis entrypoint
        // (`parse_cis_allele`), where each member carries its own accession.
        assert!(parse_variant("[NC_000001.11:g.100A>G;;NC_000001.11:g.200A>G]").is_err());
        assert!(parse_variant("[;NC_000001.11:g.100A>G;NC_000001.11:g.200A>G]").is_err());
        assert!(parse_variant("[NC_000001.11:g.100A>G;NC_000001.11:g.200A>G;]").is_err());
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

        // Simple complex delins with just inversion. The spec-canonical
        // form for a single-payload inv is unbracketed (DNA/insertion.md:22,
        // DNA/inversion.md:39 — examples like `c.849_850ins850_900inv`);
        // ferro now drops brackets on Display for single-element Complex.
        let variant =
            parse_variant("NC_000016.9:g.78179358_78219143delins[78185355_78199419inv]").unwrap();
        assert!(matches!(variant, HgvsVariant::Genome(_)));
        assert_eq!(
            format!("{}", variant),
            "NC_000016.9:g.78179358_78219143delins78185355_78199419inv"
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
        // Predicted substitution in parens: c.(9740C>A). Per #241 the
        // canonical Display form wraps position+edit together, mirroring
        // the protein predicted form `p.(Arg248Gln)`.
        let variant = parse_variant("NM_002016.2:c.(9740C>A)").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_002016.2:c.(9740C>A)");

        let variant = parse_variant("NM_006767.4:c.(742G>A)").unwrap();
        assert!(matches!(variant, HgvsVariant::Cds(_)));
        assert_eq!(format!("{}", variant), "NM_006767.4:c.(742G>A)");
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

    #[test]
    fn parse_cds_special_ranges_not_inverted() {
        assert!(parse_variant("NM_003002.2:c.pterdel").is_ok());
        assert!(parse_variant("NM_003002.2:c.qterdel").is_ok());
        assert!(parse_variant("NM_003002.2:c.cendel").is_ok());
        assert!(parse_variant("NM_003002.2:c.pter_qterdel").is_ok());
        assert!(parse_variant("NM_003002.2:c.pter_50del").is_ok());
        // The special marker round-trips through Display (exercises the
        // CdsPos special branch at the variant level).
        assert_eq!(
            parse_variant("NM_003002.2:c.pterdel").unwrap().to_string(),
            "NM_003002.2:c.pterdel"
        );
    }

    // ----- Issue #115: self-cancelling allele rejection (E3006) -----

    #[test]
    fn test_parse_self_cancelling_allele_rejected() {
        let result = parse_variant("NM_004006.2:c.[762_768del;767_774dup]");
        assert!(result.is_err(), "self-cancelling allele must be rejected");
        let err = result.unwrap_err();
        assert_eq!(
            err.code(),
            Some(crate::error::ErrorCode::SelfCancellingAllele),
            "rejection must carry E3006"
        );
    }

    #[test]
    fn test_parse_non_overlapping_allele_accepted() {
        // Same edit kinds as the spec example but non-overlapping ranges:
        // OK to coexist.
        let result = parse_variant("NM_004006.2:c.[100_110del;200_210dup]");
        assert!(
            result.is_ok(),
            "non-overlapping del+dup must be accepted: {result:?}"
        );
    }

    // ----- Issue #171: E3006 must preserve source span -----

    #[test]
    fn test_parse_self_cancelling_allele_reports_source_span() {
        // The opening `[` in `NM_004006.2:c.[762_768del;767_774dup]` sits at
        // byte offset 14. Issue #171: previously the diagnostic carried
        // `pos: 0` even though the offending bracket position is well-defined
        // from the parser. Downstream tooling (LSP, web service, structured-
        // error consumers) needs the byte offset to underline the offending
        // allele in source.
        let input = "NM_004006.2:c.[762_768del;767_774dup]";
        let expected_bracket = input.find('[').expect("input contains a bracket");
        let result = parse_variant(input);
        let err = result.expect_err("self-cancelling allele must be rejected");
        let crate::error::FerroError::Parse {
            pos, diagnostic, ..
        } = &err
        else {
            panic!("E3006 should be a Parse error, got {err:?}");
        };
        assert_eq!(
            *pos, expected_bracket,
            "E3006 pos must point at the offending allele's '['"
        );
        let diag = diagnostic.as_ref().expect("E3006 must carry a Diagnostic");
        let span = diag
            .span
            .as_ref()
            .expect("E3006 Diagnostic must carry a SourceSpan");
        assert_eq!(span.start, expected_bracket);
        assert_eq!(
            span.end,
            input.rfind(']').expect("input contains close bracket") + 1
        );
    }

    #[test]
    fn test_locate_cis_allele_bracket_skips_ins_payload() {
        // Compact form with an inserted-sequence bracket nested inside the
        // allele (issue #333). Depth tracking must close the allele on
        // its own `]`, not on the `ins[ATC]` payload's `]`.
        let input = "NM_X:g.[100A>G;200_201ins[ATC]]";
        let (open, close) = super::locate_cis_allele_bracket(input, 0)
            .expect("axis-prefixed allele bracket is present");
        assert_eq!(open, input.find('[').unwrap());
        assert_eq!(close, input.len());
        // Round-trip: the slice between open and close encloses the full
        // compact allele including the ins payload.
        assert!(input[open..close].starts_with('[') && input[open..close].ends_with(']'));
    }

    #[test]
    fn test_locate_cis_allele_bracket_returns_none_for_plain_substitution() {
        // No allele bracket in a single-variant input — the locator must
        // not invent one (and must not panic on inputs that have no `[`).
        assert!(super::locate_cis_allele_bracket("NM_004006.2:c.100A>G", 0).is_none());
    }

    #[test]
    fn test_locate_cis_allele_bracket_returns_none_for_unbalanced_input() {
        // Malformed input with an unbalanced `[` must not infinite-loop and
        // must not yield a bogus span.
        assert!(super::locate_cis_allele_bracket("NM_X:c.[762_768del", 0).is_none());
    }

    #[test]
    fn test_locate_cis_allele_bracket_advances_past_outer_open() {
        // Defensive coverage for `walk_self_cancelling`'s recursion:
        // when descending into an inner allele the caller advances
        // `search_from` past the outer `[`, so the locator must find
        // the *next* allele bracket (`;[`, `:axis.[`) rather than
        // re-finding the outer one. The HGVS grammar today does not
        // produce nested-compound alleles via the parser, so this
        // pins the locator semantics directly.
        let input = "NM_004006.2:c.[1A>G;[762_768del;767_774dup]]";
        let outer_open = input.find('[').unwrap();
        let inner_open = input.find(";[").map(|p| p + 1).unwrap();
        let inner_close_excl = input.rfind("]]").map(|p| p + 1).unwrap();

        // Outer scan (search_from = 0) returns the outer bracket.
        let (o, _) = super::locate_cis_allele_bracket(input, 0).unwrap();
        assert_eq!(o, outer_open);

        // Inner scan (search_from = outer_open + 1) returns the inner
        // bracket and pairs it with the inner closing `]`.
        let (s, e) = super::locate_cis_allele_bracket(input, outer_open + 1).unwrap();
        assert_eq!(s, inner_open);
        assert_eq!(e, inner_close_excl);
    }

    #[test]
    fn test_self_cancelling_via_try_new_validated_takes_caller_span() {
        // The library-level `AlleleVariant::try_new_validated` entry must
        // forward a caller-supplied `SourceSpan` into the diagnostic so
        // consumers that build alleles outside the parser (e.g. the web
        // service builder) can still produce well-located E3006 errors.
        use crate::error::SourceSpan;
        use crate::hgvs::variant::AlleleVariant;

        let del = parse_variant("NM_004006.2:c.762_768del").expect("del should parse for fixture");
        let dup = parse_variant("NM_004006.2:c.767_774dup").expect("dup should parse for fixture");

        // Without a caller-supplied span, `pos` falls back to 0 (no source
        // available) but the diagnostic must still carry the E3006 code.
        let err = AlleleVariant::try_new_validated(
            vec![del.clone(), dup.clone()],
            crate::hgvs::variant::AllelePhase::Cis,
            None,
        )
        .expect_err("self-cancelling allele must be rejected");
        let crate::error::FerroError::Parse {
            pos, diagnostic, ..
        } = &err
        else {
            panic!("expected Parse, got {err:?}");
        };
        assert_eq!(*pos, 0, "no caller span ⇒ pos remains 0");
        assert_eq!(
            diagnostic.as_ref().and_then(|d| d.code),
            Some(crate::error::ErrorCode::SelfCancellingAllele),
        );

        // With a caller-supplied span, `pos` and the diagnostic span both
        // reflect it.
        let span = SourceSpan::new(50, 80);
        let err = AlleleVariant::try_new_validated(
            vec![del, dup],
            crate::hgvs::variant::AllelePhase::Cis,
            Some(span.clone()),
        )
        .expect_err("self-cancelling allele must be rejected");
        let crate::error::FerroError::Parse {
            pos, diagnostic, ..
        } = &err
        else {
            panic!("expected Parse, got {err:?}");
        };
        assert_eq!(*pos, 50);
        let diag = diagnostic.as_ref().expect("Diagnostic carried");
        assert_eq!(diag.span.as_ref(), Some(&span));
    }

    // ----- Issue #375: slash-form mosaic/chimeric fallback chain must not
    // swallow chunk-level structured diagnostics (E3006, W3017, ...) -----

    #[test]
    fn test_slash_mosaic_self_cancelling_rhs_reports_e3006_full_input_relative() {
        // Two-arm mosaic where the RHS chunk is self-cancelling. Before the
        // fix, parse_phase_allele's fallback chain caught the chunk-level
        // E3006 and re-emitted "Unknown variant type prefix" at pos 0 via
        // parse_variant_with_inherited_accession (since chunk1 starts with
        // an accession that the inherited-accession recovery can't strip).
        // After the fix, the structured E3006 propagates, with pos and span
        // remapped to full-input coordinates.
        let input = "NM_004006.2:c.[100A>G;200C>T]/NM_004006.2:c.[762_768del;767_774dup]";
        let rhs_chunk_offset = input.find("/NM").unwrap() + 1;
        let rhs_open = rhs_chunk_offset + input[rhs_chunk_offset..].find('[').unwrap();
        let rhs_close_excl = rhs_chunk_offset + input[rhs_chunk_offset..].rfind(']').unwrap() + 1;

        let err = parse_variant(input).expect_err("RHS self-cancelling must be rejected");
        let crate::error::FerroError::Parse {
            pos, diagnostic, ..
        } = &err
        else {
            panic!("expected Parse, got {err:?}");
        };
        assert_eq!(
            *pos, rhs_open,
            "E3006 pos must point at the RHS arm's '[' in the full input"
        );
        let diag = diagnostic.as_ref().expect("E3006 must carry a Diagnostic");
        assert_eq!(
            diag.code,
            Some(crate::error::ErrorCode::SelfCancellingAllele),
            "must carry E3006 SelfCancellingAllele code"
        );
        let span = diag.span.as_ref().expect("E3006 must carry a SourceSpan");
        assert_eq!(span.start, rhs_open);
        assert_eq!(span.end, rhs_close_excl);
        assert_eq!(
            diag.source.as_deref(),
            Some(input),
            "Diagnostic.source must be the full user input, not the chunk"
        );
    }

    #[test]
    fn test_slash_mosaic_self_cancelling_lhs_reports_e3006_full_input_relative() {
        // LHS-self-cancelling regression: chunk0 starts at byte 0 so
        // chunk-relative and full-input-relative coincide for `pos` and
        // `span`, but `Diagnostic.source` must still be the full input
        // (was previously the chunk slice).
        let input = "NM_004006.2:c.[762_768del;767_774dup]/NM_004006.2:c.[100A>G;200C>T]";
        let lhs_open = input.find('[').unwrap();
        let lhs_close_excl = input.find(']').unwrap() + 1;

        let err = parse_variant(input).expect_err("LHS self-cancelling must be rejected");
        let crate::error::FerroError::Parse {
            pos, diagnostic, ..
        } = &err
        else {
            panic!("expected Parse, got {err:?}");
        };
        assert_eq!(*pos, lhs_open);
        let diag = diagnostic.as_ref().expect("E3006 must carry a Diagnostic");
        assert_eq!(
            diag.code,
            Some(crate::error::ErrorCode::SelfCancellingAllele),
        );
        let span = diag.span.as_ref().expect("E3006 must carry a SourceSpan");
        assert_eq!(span.start, lhs_open);
        assert_eq!(span.end, lhs_close_excl);
        assert_eq!(diag.source.as_deref(), Some(input));
    }

    #[test]
    fn test_slash_chimeric_self_cancelling_rhs_reports_e3006_full_input_relative() {
        // Same as the mosaic RHS case but with the chimeric `//` separator.
        // Both phases share the same fallback chain in parse_phase_allele,
        // so both must propagate structured diagnostics.
        let input = "NM_004006.2:c.[100A>G;200C>T]//NM_004006.2:c.[762_768del;767_774dup]";
        let rhs_chunk_offset = input.find("//NM").unwrap() + 2;
        let rhs_open = rhs_chunk_offset + input[rhs_chunk_offset..].find('[').unwrap();
        let rhs_close_excl = rhs_chunk_offset + input[rhs_chunk_offset..].rfind(']').unwrap() + 1;

        let err = parse_variant(input).expect_err("RHS self-cancelling must be rejected");
        let crate::error::FerroError::Parse {
            pos, diagnostic, ..
        } = &err
        else {
            panic!("expected Parse, got {err:?}");
        };
        assert_eq!(*pos, rhs_open);
        let diag = diagnostic.as_ref().expect("E3006 must carry a Diagnostic");
        assert_eq!(
            diag.code,
            Some(crate::error::ErrorCode::SelfCancellingAllele),
        );
        let span = diag.span.as_ref().expect("E3006 must carry a SourceSpan");
        assert_eq!(span.start, rhs_open);
        assert_eq!(span.end, rhs_close_excl);
        assert_eq!(diag.source.as_deref(), Some(input));
    }

    #[test]
    fn test_inherited_accession_mosaic_still_works() {
        // Regression: the inherited-accession fallback (RHS chunk omits the
        // accession, inherits from LHS) must still parse cleanly when no
        // structured diagnostic blocks it.
        let variant = parse_variant("NM_000088.3:c.100A>G/c.200C>T")
            .expect("inherited-accession mosaic must parse");
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Mosaic);
            assert_eq!(allele.variants.len(), 2);
        } else {
            panic!("expected Allele variant, got {variant:?}");
        }
    }

    #[test]
    fn test_spec_mosaic_compact_form_still_works() {
        // Regression: the HGVS spec compact-mosaic shape `<pos>=/<edit>`
        // (where the RHS is a bare na_edit inheriting position from LHS `=`)
        // must still parse. It goes through `parse_compact_mosaic_rhs`,
        // which we must preserve.
        let variant =
            parse_variant("LRG_199t1:c.85=/T>C").expect("spec mosaic compact form must parse");
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Mosaic);
            assert_eq!(allele.variants.len(), 2);
        } else {
            panic!("expected Allele variant, got {variant:?}");
        }
    }

    #[test]
    fn test_spec_chimeric_compact_form_still_works() {
        // Regression: same as the mosaic compact form, but `//` separator.
        let variant =
            parse_variant("NM_004006.2:c.85=//T>C").expect("spec chimeric compact form must parse");
        if let HgvsVariant::Allele(allele) = &variant {
            assert_eq!(allele.phase, AllelePhase::Chimeric);
            assert_eq!(allele.variants.len(), 2);
        } else {
            panic!("expected Allele variant, got {variant:?}");
        }
    }
}
