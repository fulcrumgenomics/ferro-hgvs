//! Fast-path parsers for common HGVS patterns
//!
//! These specialized parsers bypass nom's combinator overhead for the most common
//! patterns, providing significant speedups for typical real-world workloads.
//!
//! # Performance
//!
//! The fast-path parser provides a **~1.7x end-to-end speedup** on the full ClinVar
//! corpus. It handles ~72% of real ClinVar variants (RefSeq / Ensembl / LRG /
//! Assembly `g.`/`c.` substitutions and plain deletions/duplications);
//! everything else falls back to the generic parser.
//! The per-pattern micro-benchmark speedup (fast path alone vs. generic parser) is
//! higher — roughly 45-58% faster for the patterns listed below — but the corpus-level
//! figure accounts for the fallback fraction and the fast-path overhead on
//! non-eligible inputs.
//!
//! | Pattern | Micro-benchmark speedup |
//! |---------|------------------------|
//! | `NC_000001.11:g.12345A>G` (RefSeq genomic) | ~50% faster |
//! | `NM_000088.3:c.459A>G` (RefSeq coding) | ~57% faster |
//! | `ENST00000357033.8:c.100A>G` (Ensembl) | ~54% faster |
//! | `GRCh38(chr1):g.12345A>G` (Assembly) | ~47% faster |
//!
//! # Tradeoffs
//!
//! The fast-path adds a small overhead (~3-6%) for patterns it cannot optimize:
//! - Non-coding RNA: `n.100A>G`
//! - RNA variants: `r.100a>g`
//! - Non-coding RNA accessions: NR_/XR_ (any coordinate system)
//! - Insertions / delins / inversions / repeats
//!
//! This overhead comes from the quick-rejection checks needed to identify
//! fast-path candidates. For workloads dominated by substitutions, the net
//! effect is strongly positive.
//!
//! # Supported Fast-Paths
//!
//! Substitutions (position + single base change), including intronic
//! (`c.100+5A>G`) and UTR (`c.*100A>G`, `c.-50A>G`) coding positions, and plain
//! deletions / duplications (a single position or `a_b` range followed by a bare
//! `del` / `dup`, with no trailing sequence or length) for:
//! - **RefSeq**: NC_, NM_, NG_, XM_ accessions with `g.` or `c.` variants
//!   (`NR_`/`XR_` non-coding RNA transcripts are deferred to the generic
//!   parser — their only valid coordinate systems are `n.`/`r.`)
//! - **Ensembl**: ENST, ENSG accessions with `g.` or `c.` variants
//! - **LRG**: LRG_ accessions with `g.` or `c.` variants
//! - **Assembly**: GRCh37/38, hg19/38 notation with `g.` variants
//!
//! Deletion/duplication forms that carry a deleted sequence or explicit length
//! (`delA`, `del3`, `dupAG`), insertions, delins, inversions, and repeats are
//! left to the generic parser, which owns their canonicalization.
//!
//! # When to Use
//!
//! Call [`super::parse_hgvs`] for all new code — the fast path is now the
//! default, so callers get the performance benefit automatically without
//! selecting a separate entry point.
//!
//! [`super::parse_hgvs_fast`] is retained for backward compatibility only; it
//! delegates unconditionally to [`super::parse_hgvs`] and offers no additional
//! performance over calling [`super::parse_hgvs`] directly.

use crate::hgvs::edit::{Base, NaEdit};
use crate::hgvs::interval::{CdsInterval, GenomeInterval};
use crate::hgvs::location::{CdsPos, GenomePos};
use crate::hgvs::variant::{Accession, CdsVariant, GenomeVariant, HgvsVariant, LocEdit};

/// Lookup table for valid IUPAC base characters
/// Supports A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N
#[inline]
const fn is_iupac_base(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'U'
            | b'R'
            | b'Y'
            | b'S'
            | b'W'
            | b'K'
            | b'M'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'N'
    )
}

/// Scan consecutive ASCII digits and compute their numeric value
#[inline]
fn scan_digits(bytes: &[u8], start: usize) -> (u64, usize) {
    let mut end = start;
    let mut value = 0u64;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        value = value * 10 + (bytes[end] - b'0') as u64;
        end += 1;
    }
    (value, end)
}

/// True for a base byte the fast path may parse in a DNA substitution. This is
/// [`is_iupac_base`] minus the RNA base `U`: HGVS forbids `U` in a DNA (`g.` /
/// `c.`) description (#486, `ENODNA`), so a `U` base must defer to the generic
/// parser, which rejects it. (Lowercase bases already fail `is_iupac_base`.)
#[inline]
const fn is_fast_dna_base(b: u8) -> bool {
    is_iupac_base(b) && b != b'U'
}

/// Scan an HGVS *position*: a run of ASCII digits that is canonical (no leading
/// zero) and fits in `u64`. Returns `(value, end)`, or `None` if the position
/// is absent, has a leading zero (`0`, `00`, `007`, `01` — all rejected by the
/// generic parser), or overflows `u64`; the fast path defers to the generic
/// parser on `None`. Distinct from [`scan_digits`], which is also used for
/// accession numbers where leading zeros are valid (e.g. `NM_000088`).
#[inline]
fn scan_position(bytes: &[u8], start: usize) -> Option<(u64, usize)> {
    if start >= bytes.len() || !bytes[start].is_ascii_digit() || bytes[start] == b'0' {
        return None;
    }
    let mut end = start;
    let mut value: u64 = 0;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        value = value
            .checked_mul(10)?
            .checked_add((bytes[end] - b'0') as u64)?;
        end += 1;
    }
    Some((value, end))
}

/// Result of fast-path parsing attempt
#[allow(clippy::large_enum_variant)]
pub enum FastPathResult {
    /// Successfully parsed with fast path
    Success(HgvsVariant),
    /// Pattern not recognized - fall back to generic parser
    Fallback,
}

/// The trailing edit kind a fast-path candidate carries.
///
/// Only the simplest, unambiguous forms are recognized; everything else defers
/// to the generic parser. All three recognized tails are exactly three bytes
/// (`X>Y`, `del`, `dup`), so the downstream intronic/UTR exclusion checks that
/// key off `len - 3` apply uniformly.
#[derive(Clone, Copy)]
enum FastEdit {
    Substitution,
    Deletion,
    Duplication,
}

/// Classify the trailing edit of a fast-path candidate, or `None` to defer.
///
/// - Substitution: `…<base>><base>`.
/// - Plain deletion / duplication: `…<digit>del` / `…<digit>dup`. The required
///   preceding digit means a position immediately precedes the keyword, so only
///   the bare forms (`100del`, `100_200dup`) match. Forms that carry a trailing
///   sequence or length (`delA`, `del3`, `dupAG`) end in a base or digit — not
///   `l`/`p` — and `delins`/`dupins` end in a base, so all of them fall through
///   to the generic parser, which owns their canonicalization.
#[inline]
fn classify_fast_edit(bytes: &[u8], len: usize) -> Option<FastEdit> {
    if len >= 3
        && bytes[len - 2] == b'>'
        && is_fast_dna_base(bytes[len - 1])
        && is_fast_dna_base(bytes[len - 3])
    {
        return Some(FastEdit::Substitution);
    }
    if len >= 4 && bytes[len - 4].is_ascii_digit() {
        match &bytes[len - 3..] {
            b"del" => return Some(FastEdit::Deletion),
            b"dup" => return Some(FastEdit::Duplication),
            _ => {}
        }
    }
    None
}

/// Attempt to parse using fast-path for known patterns
///
/// Returns `FastPathResult::Success` if the pattern was recognized and parsed,
/// or `FastPathResult::Fallback` if the generic parser should be used.
#[inline]
pub fn try_fast_path(input: &str) -> FastPathResult {
    let bytes = input.as_bytes();
    let len = bytes.len();

    if len < 10 {
        // Minimum: "X_1.1:g.1A>G" = 12 chars
        return FastPathResult::Fallback;
    }

    // Classify the trailing edit up front. This is the quick-rejection gate:
    // anything that is not a recognized substitution / plain del / plain dup
    // tail defers immediately, avoiding parse work for the patterns the fast
    // path does not handle.
    let edit_kind = match classify_fast_edit(bytes, len) {
        Some(k) => k,
        None => return FastPathResult::Fallback,
    };

    // Quick exclusion: find colon and check variant type + complex patterns
    // Combined check to minimize overhead
    if let Some(colon_pos) = memchr::memchr(b':', bytes) {
        if colon_pos + 2 < len {
            let type_char = bytes[colon_pos + 1];
            let dot_char = bytes[colon_pos + 2];

            // Exclude :n. (non-coding RNA) and :r. (RNA) - not supported
            if (type_char == b'n' || type_char == b'r') && dot_char == b'.' {
                return FastPathResult::Fallback;
            }

            // For :c. del/dup, exclude UTR (`*`/`-`) and intronic (`+`/`-`)
            // coordinates: the del/dup builders parse only plain in-CDS
            // positions. Substitutions handle UTR/intronic positions themselves
            // (see `scan_cds_sub_position`), so they are not excluded here.
            if type_char == b'c'
                && dot_char == b'.'
                && colon_pos + 3 < len
                && !matches!(edit_kind, FastEdit::Substitution)
            {
                let pos_start = bytes[colon_pos + 3];
                // c.*100 (UTR3) or c.-50 (UTR5/negative position)
                if pos_start == b'*' || pos_start == b'-' {
                    return FastPathResult::Fallback;
                }

                // Check for intronic offset using memchr (faster than loop).
                // Look for + or - in the region between the position start and
                // the trailing keyword. Guard the bounds: for a very short `c.`
                // body `colon_pos + 4` can exceed `len - 3`, and an unchecked
                // slice would panic instead of falling back. When the region is
                // empty there is no room for an intronic offset anyway, so skip.
                let region_start = colon_pos + 4;
                let region_end = len - 3;
                if region_start < region_end {
                    let search_region = &bytes[region_start..region_end];
                    if memchr::memchr2(b'+', b'-', search_region).is_some() {
                        return FastPathResult::Fallback;
                    }
                }
            }
        }
    }

    // Dispatch based on first character
    match bytes[0] {
        // RefSeq accessions: NC_, NM_, NP_, NG_, XM_, XP_
        // (NR_/XR_ are deferred to the generic parser before try_refseq_fast_path)
        b'N' | b'X' => {
            if bytes.len() > 2 && bytes[2] == b'_' {
                // `NR_`/`XR_` are non-coding RNA transcripts (the only RefSeq
                // nucleotide prefixes whose second letter is `R`). Per HGVS
                // (refseq.md: the prefix declares the reference *type*; `n.`
                // is "based on a transcript not coding for a protein"), their
                // only valid coordinate systems are `n.`/`r.`. The fast path
                // parses only `g.`/`c.` substitutions, so parsing one here
                // would mint a `Genome`/`Cds` variant — silently accepting a
                // coordinate-system/reference-type mismatch that the generic
                // parser rejects (#486, ECOORDINATESYSTEMMISMATCH). Defer to
                // `parse_variant` so it owns that spec rule and the two entry
                // points stay observationally identical.
                if bytes[1] == b'R' {
                    return FastPathResult::Fallback;
                }
                return try_refseq_fast_path(input, bytes, edit_kind);
            }
            FastPathResult::Fallback
        }
        // Ensembl accessions: ENST, ENSG, ENSP
        b'E' => {
            if bytes.len() >= 4 && bytes[1] == b'N' && bytes[2] == b'S' {
                return try_ensembl_fast_path(input, bytes, edit_kind);
            }
            FastPathResult::Fallback
        }
        // LRG accessions
        b'L' => {
            if bytes.len() >= 4 && bytes[1] == b'R' && bytes[2] == b'G' && bytes[3] == b'_' {
                return try_lrg_fast_path(input, bytes, edit_kind);
            }
            FastPathResult::Fallback
        }
        // Assembly notation: GRCh37, GRCh38
        b'G' => {
            if bytes.len() >= 6 && bytes[1] == b'R' && bytes[2] == b'C' && bytes[3] == b'h' {
                return try_assembly_fast_path(input, bytes, edit_kind);
            }
            FastPathResult::Fallback
        }
        // Assembly notation: hg18, hg19, hg38
        b'h' => {
            if bytes.len() >= 4 && bytes[1] == b'g' {
                return try_assembly_fast_path(input, bytes, edit_kind);
            }
            FastPathResult::Fallback
        }
        _ => FastPathResult::Fallback,
    }
}

/// Fast-path for RefSeq accessions (NC_, NM_, NP_, etc.)
///
/// Pattern: `[NX][A-Z]_DIGITS.VERSION:TYPE.POSITION[EDIT]`
#[inline]
fn try_refseq_fast_path(input: &str, bytes: &[u8], edit_kind: FastEdit) -> FastPathResult {
    // Parse prefix: two letters
    let prefix_end = 2;
    if !bytes[0].is_ascii_uppercase() || !bytes[1].is_ascii_uppercase() {
        return FastPathResult::Fallback;
    }

    // Skip underscore
    if bytes[2] != b'_' {
        return FastPathResult::Fallback;
    }

    // Parse accession number (digits)
    let (_number_value, number_end) = scan_digits(bytes, 3);
    if number_end == 3 {
        return FastPathResult::Fallback; // No digits found
    }
    let number_str = &input[3..number_end];

    // Parse optional version
    let (version, version_end) = if number_end < bytes.len() && bytes[number_end] == b'.' {
        let (v, ve) = scan_digits(bytes, number_end + 1);
        if ve == number_end + 1 {
            return FastPathResult::Fallback; // Dot but no version
        }
        (Some(v as u32), ve)
    } else {
        (None, number_end)
    };

    // Expect colon
    if version_end >= bytes.len() || bytes[version_end] != b':' {
        return FastPathResult::Fallback;
    }

    // Get type prefix (c., g., p., n., etc.)
    let type_start = version_end + 1;
    if type_start + 2 > bytes.len() || bytes[type_start + 1] != b'.' {
        return FastPathResult::Fallback;
    }

    let prefix = &input[0..prefix_end];
    let accession = Accession::with_style(
        prefix, number_str, version, false, // RefSeq uses underscore style
    );

    // Dispatch based on type
    let type_char = bytes[type_start];
    let edit_start = type_start + 2;

    match type_char {
        b'g' => dispatch_genome(input, bytes, edit_start, accession, edit_kind),
        b'c' => dispatch_cds(input, bytes, edit_start, accession, edit_kind),
        b'p' => FastPathResult::Fallback, // Protein is more complex
        _ => FastPathResult::Fallback,
    }
}

/// Fast-path for Ensembl accessions (ENST, ENSG, ENSP)
///
/// Pattern: `ENS[TGPSR]DIGITS.VERSION:TYPE.POSITION[EDIT]`
#[inline]
fn try_ensembl_fast_path(input: &str, bytes: &[u8], edit_kind: FastEdit) -> FastPathResult {
    // Check ENS prefix
    if bytes.len() < 15 || bytes[0] != b'E' || bytes[1] != b'N' || bytes[2] != b'S' {
        return FastPathResult::Fallback;
    }

    // Check type character
    let type_char = bytes[3];
    if !matches!(type_char, b'T' | b'G' | b'P' | b'E' | b'R') {
        return FastPathResult::Fallback;
    }

    // Parse digits (Ensembl IDs have 11-15 digits to accommodate various formats)
    let (_number_value, number_end) = scan_digits(bytes, 4);
    let digit_count = number_end - 4;
    if number_end == 4 || !(11..=15).contains(&digit_count) {
        return FastPathResult::Fallback;
    }
    let number_str = &input[4..number_end];

    // Parse optional version
    let (version, version_end) = if number_end < bytes.len() && bytes[number_end] == b'.' {
        let (v, ve) = scan_digits(bytes, number_end + 1);
        if ve == number_end + 1 {
            return FastPathResult::Fallback;
        }
        (Some(v as u32), ve)
    } else {
        (None, number_end)
    };

    // Expect colon
    if version_end >= bytes.len() || bytes[version_end] != b':' {
        return FastPathResult::Fallback;
    }

    // Get variant type prefix
    let var_type_start = version_end + 1;
    if var_type_start + 2 > bytes.len() || bytes[var_type_start + 1] != b'.' {
        return FastPathResult::Fallback;
    }

    let prefix = &input[0..4]; // ENST, ENSG, etc.
    let accession = Accession::with_style(
        prefix, number_str, version, true, // Ensembl style (no underscore)
    );

    let var_type_char = bytes[var_type_start];
    let edit_start = var_type_start + 2;

    match var_type_char {
        b'g' => dispatch_genome(input, bytes, edit_start, accession, edit_kind),
        b'c' => dispatch_cds(input, bytes, edit_start, accession, edit_kind),
        b'p' => FastPathResult::Fallback,
        _ => FastPathResult::Fallback,
    }
}

/// Fast-path for LRG accessions
///
/// Pattern: `LRG_DIGITS:TYPE.POSITION[EDIT]`
#[inline]
fn try_lrg_fast_path(input: &str, bytes: &[u8], edit_kind: FastEdit) -> FastPathResult {
    // Check LRG_ prefix
    if bytes.len() < 8
        || bytes[0] != b'L'
        || bytes[1] != b'R'
        || bytes[2] != b'G'
        || bytes[3] != b'_'
    {
        return FastPathResult::Fallback;
    }

    // Parse LRG number
    let (_number_value, number_end) = scan_digits(bytes, 4);
    if number_end == 4 {
        return FastPathResult::Fallback;
    }
    let number_str = &input[4..number_end];

    // Check for transcript suffix (t1, p1) or direct colon
    let (full_number, version_end) =
        if number_end < bytes.len() && (bytes[number_end] == b't' || bytes[number_end] == b'p') {
            // LRG with transcript: LRG_123t1
            let (_tx_num, tx_end) = scan_digits(bytes, number_end + 1);
            if tx_end == number_end + 1 {
                return FastPathResult::Fallback;
            }
            (&input[4..tx_end], tx_end)
        } else {
            (number_str, number_end)
        };

    // Expect colon
    if version_end >= bytes.len() || bytes[version_end] != b':' {
        return FastPathResult::Fallback;
    }

    // Get variant type prefix
    let type_start = version_end + 1;
    if type_start + 2 > bytes.len() || bytes[type_start + 1] != b'.' {
        return FastPathResult::Fallback;
    }

    let accession = Accession::with_style(
        "LRG",
        full_number,
        None, // LRG doesn't use versions
        false,
    );

    let type_char = bytes[type_start];
    let edit_start = type_start + 2;

    match type_char {
        b'g' => dispatch_genome(input, bytes, edit_start, accession, edit_kind),
        b'c' => dispatch_cds(input, bytes, edit_start, accession, edit_kind),
        b'p' => FastPathResult::Fallback,
        _ => FastPathResult::Fallback,
    }
}

/// Fast-path for assembly notation (GRCh37, GRCh38, hg19, hg38)
///
/// Pattern: `ASSEMBLY(CHROM):TYPE.POSITION[EDIT]`
#[inline]
fn try_assembly_fast_path(input: &str, bytes: &[u8], edit_kind: FastEdit) -> FastPathResult {
    // Find opening paren
    let paren_pos = bytes.iter().position(|&b| b == b'(');
    if paren_pos.is_none() {
        return FastPathResult::Fallback;
    }
    let paren_pos = paren_pos.unwrap();

    // Validate assembly name
    let assembly = &input[0..paren_pos];
    if !matches!(
        assembly,
        "GRCh37" | "GRCh38" | "hg19" | "hg38" | "hg18" | "GRCh36"
    ) {
        return FastPathResult::Fallback;
    }

    // Find closing paren
    let close_paren = bytes[paren_pos + 1..].iter().position(|&b| b == b')');
    if close_paren.is_none() {
        return FastPathResult::Fallback;
    }
    let close_paren = paren_pos + 1 + close_paren.unwrap();

    let chromosome = &input[paren_pos + 1..close_paren];

    // Expect colon after closing paren
    if close_paren + 1 >= bytes.len() || bytes[close_paren + 1] != b':' {
        return FastPathResult::Fallback;
    }

    // Get variant type prefix
    let type_start = close_paren + 2;
    if type_start + 2 > bytes.len() || bytes[type_start + 1] != b'.' {
        return FastPathResult::Fallback;
    }

    let accession = Accession::from_assembly(assembly, chromosome);

    let type_char = bytes[type_start];
    let edit_start = type_start + 2;

    match type_char {
        b'g' => dispatch_genome(input, bytes, edit_start, accession, edit_kind),
        _ => FastPathResult::Fallback, // Assembly is typically genomic only
    }
}

/// Dispatch a genomic (`g.`) edit to the matching fast-path builder.
#[inline]
fn dispatch_genome(
    input: &str,
    bytes: &[u8],
    edit_start: usize,
    accession: Accession,
    edit_kind: FastEdit,
) -> FastPathResult {
    match edit_kind {
        FastEdit::Substitution => {
            try_parse_genome_substitution(input, bytes, edit_start, accession)
        }
        FastEdit::Deletion => try_parse_genome_del_dup(bytes, edit_start, accession, false),
        FastEdit::Duplication => try_parse_genome_del_dup(bytes, edit_start, accession, true),
    }
}

/// Dispatch a coding (`c.`) edit to the matching fast-path builder.
#[inline]
fn dispatch_cds(
    input: &str,
    bytes: &[u8],
    edit_start: usize,
    accession: Accession,
    edit_kind: FastEdit,
) -> FastPathResult {
    match edit_kind {
        FastEdit::Substitution => try_parse_cds_substitution(input, bytes, edit_start, accession),
        FastEdit::Deletion => try_parse_cds_del_dup(bytes, edit_start, accession, false),
        FastEdit::Duplication => try_parse_cds_del_dup(bytes, edit_start, accession, true),
    }
}

/// Scan a fast-path genomic interval at `start`: either a single canonical
/// position or a `a_b` range. Returns the interval and the index just past the
/// position text, or `None` to defer. Defers on a non-canonical/overflowing
/// position and on an inverted range (`a > b`): the generic parser rejects an
/// inverted *deletion* and accepts an inverted *duplication*, so rather than
/// replicate that split we cede both inverted cases to it (they are rare and
/// `Fallback` stays observationally identical either way).
#[inline]
fn scan_genome_interval(bytes: &[u8], start: usize) -> Option<(GenomeInterval, usize)> {
    let (p1, e1) = scan_position(bytes, start)?;
    if e1 < bytes.len() && bytes[e1] == b'_' {
        let (p2, e2) = scan_position(bytes, e1 + 1)?;
        if p1 > p2 {
            return None;
        }
        Some((
            GenomeInterval::new(GenomePos::new(p1), GenomePos::new(p2)),
            e2,
        ))
    } else {
        Some((GenomeInterval::point(GenomePos::new(p1)), e1))
    }
}

/// Scan a fast-path coding interval at `start`. Like [`scan_genome_interval`]
/// but builds `CdsPos` and defers any position beyond `i64::MAX` (which would
/// wrap). UTR (`*`/`-`) and intronic (`+`/`-`) coordinates are already excluded
/// by the dispatcher's caller, so only plain in-CDS positions reach here.
#[inline]
fn scan_cds_interval(bytes: &[u8], start: usize) -> Option<(CdsInterval, usize)> {
    let (p1, e1) = scan_position(bytes, start)?;
    if p1 > i64::MAX as u64 {
        return None;
    }
    if e1 < bytes.len() && bytes[e1] == b'_' {
        let (p2, e2) = scan_position(bytes, e1 + 1)?;
        if p2 > i64::MAX as u64 || p1 > p2 {
            return None;
        }
        Some((
            CdsInterval::new(CdsPos::new(p1 as i64), CdsPos::new(p2 as i64)),
            e2,
        ))
    } else {
        Some((CdsInterval::point(CdsPos::new(p1 as i64)), e1))
    }
}

/// Scan a fast-path CDS *substitution* position at `start`: a plain in-CDS
/// position, a 5'UTR (`-N`) or 3'UTR (`*N`) position, each with an optional
/// intronic offset (`+M`/`-M`). Returns the [`CdsPos`] and the index just past
/// it, or `None` to defer.
///
/// Defers on anything it does not reproduce exactly: leading zeros / `0` /
/// `u64`-or-`i64` overflow (via [`scan_position`]), uncertain offsets (`+?` /
/// `-?`), and `?` / `(` / special markers (which start with a non-`*`/`-`/digit
/// byte and so fail [`scan_position`]). The generic parser owns those — and
/// because the caller falls back on `None`, deferring stays observationally
/// identical.
#[inline]
fn scan_cds_sub_position(bytes: &[u8], start: usize) -> Option<(CdsPos, usize)> {
    if start >= bytes.len() {
        return None;
    }
    // `scan_position` yields a canonical, non-zero `u64`; cap it at `i64::MAX`
    // since `CdsPos::base`/`offset` are `i64`.
    let to_i64 = |v: u64| (v <= i64::MAX as u64).then_some(v as i64);

    let (base, utr3, after_base) = match bytes[start] {
        b'*' => {
            let (v, j) = scan_position(bytes, start + 1)?;
            (to_i64(v)?, true, j)
        }
        b'-' => {
            let (v, j) = scan_position(bytes, start + 1)?;
            (-to_i64(v)?, false, j)
        }
        _ => {
            let (v, j) = scan_position(bytes, start)?;
            (to_i64(v)?, false, j)
        }
    };

    let mut i = after_base;
    let offset = if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        let neg = bytes[i] == b'-';
        // Uncertain offset (`+?` / `-?`) defers to the generic parser.
        if i + 1 < bytes.len() && bytes[i + 1] == b'?' {
            return None;
        }
        let (v, j) = scan_position(bytes, i + 1)?;
        let mag = to_i64(v)?;
        i = j;
        Some(if neg { -mag } else { mag })
    } else {
        None
    };

    Some((
        CdsPos {
            base,
            offset,
            utr3,
            special: None,
        },
        i,
    ))
}

/// Build a plain genomic deletion or duplication (`g.NNNdel`, `g.A_Bdup`).
///
/// Only the bare forms (no trailing sequence or length) reach here — the edit
/// keyword must immediately end the string. `sequence`/`length`/`uncertain_extent`
/// are `None`, matching the generic parser's output for these inputs.
#[inline]
fn try_parse_genome_del_dup(
    bytes: &[u8],
    edit_start: usize,
    accession: Accession,
    is_dup: bool,
) -> FastPathResult {
    let (interval, after) = match scan_genome_interval(bytes, edit_start) {
        Some(x) => x,
        None => return FastPathResult::Fallback,
    };
    // The classified 3-char `del`/`dup` keyword must immediately end the string;
    // anything between the position and it (e.g. an extra coordinate) defers.
    if after + 3 != bytes.len() {
        return FastPathResult::Fallback;
    }
    let edit = if is_dup {
        NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: None,
        }
    } else {
        NaEdit::Deletion {
            sequence: None,
            length: None,
        }
    };
    FastPathResult::Success(HgvsVariant::Genome(GenomeVariant {
        accession,
        gene_symbol: None,
        loc_edit: LocEdit::new(interval, edit),
    }))
}

/// Build a plain coding deletion or duplication (`c.NNNdel`, `c.A_Bdup`).
#[inline]
fn try_parse_cds_del_dup(
    bytes: &[u8],
    edit_start: usize,
    accession: Accession,
    is_dup: bool,
) -> FastPathResult {
    let (interval, after) = match scan_cds_interval(bytes, edit_start) {
        Some(x) => x,
        None => return FastPathResult::Fallback,
    };
    if after + 3 != bytes.len() {
        return FastPathResult::Fallback;
    }
    let edit = if is_dup {
        NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: None,
        }
    } else {
        NaEdit::Deletion {
            sequence: None,
            length: None,
        }
    };
    FastPathResult::Success(HgvsVariant::Cds(CdsVariant {
        accession,
        gene_symbol: None,
        loc_edit: LocEdit::new(interval, edit),
    }))
}

/// Try to parse a simple genomic substitution (g.NNNA>G)
///
/// This handles the most common case: a single position with a simple substitution.
#[inline]
fn try_parse_genome_substitution(
    _input: &str,
    bytes: &[u8],
    edit_start: usize,
    accession: Accession,
) -> FastPathResult {
    // Parse a canonical position (no leading zero, fits u64); a non-canonical
    // or overflowing position defers to the generic parser, which rejects it.
    let (position, pos_end) = match scan_position(bytes, edit_start) {
        Some(p) => p,
        None => return FastPathResult::Fallback,
    };

    // Check for simple substitution pattern: A>G at end of string
    // pos_end should point to ref base, pos_end+1 to '>', pos_end+2 to alt base
    if pos_end + 3 != bytes.len() {
        return FastPathResult::Fallback; // Not a simple substitution at end
    }

    if !is_fast_dna_base(bytes[pos_end])
        || bytes[pos_end + 1] != b'>'
        || !is_fast_dna_base(bytes[pos_end + 2])
    {
        return FastPathResult::Fallback;
    }

    let reference = Base::from_char(bytes[pos_end] as char).unwrap();
    let alternative = Base::from_char(bytes[pos_end + 2] as char).unwrap();

    let pos = GenomePos::new(position);
    let interval = GenomeInterval::point(pos);
    let edit = NaEdit::Substitution {
        reference,
        alternative,
    };

    FastPathResult::Success(HgvsVariant::Genome(GenomeVariant {
        accession,
        gene_symbol: None,
        loc_edit: LocEdit::new(interval, edit),
    }))
}

/// Try to parse a CDS substitution: a plain (`c.459A>G`), intronic
/// (`c.100+5A>G`, `c.100-5A>G`), or UTR (`c.*100A>G`, `c.-50A>G`) position
/// followed by a single base change. The position is parsed by
/// [`scan_cds_sub_position`], which defers any form it cannot reproduce exactly.
#[inline]
fn try_parse_cds_substitution(
    _input: &str,
    bytes: &[u8],
    edit_start: usize,
    accession: Accession,
) -> FastPathResult {
    // Parse the CDS position: plain in-CDS, 5'UTR (`-N`), or 3'UTR (`*N`), each
    // with an optional intronic offset. Anything else defers.
    let (pos, pos_end) = match scan_cds_sub_position(bytes, edit_start) {
        Some(p) => p,
        None => return FastPathResult::Fallback,
    };

    // Check for simple substitution pattern: A>G at end of string
    if pos_end + 3 != bytes.len() {
        return FastPathResult::Fallback;
    }

    if !is_fast_dna_base(bytes[pos_end])
        || bytes[pos_end + 1] != b'>'
        || !is_fast_dna_base(bytes[pos_end + 2])
    {
        return FastPathResult::Fallback;
    }

    let reference = Base::from_char(bytes[pos_end] as char).unwrap();
    let alternative = Base::from_char(bytes[pos_end + 2] as char).unwrap();

    let interval = CdsInterval::point(pos);
    let edit = NaEdit::Substitution {
        reference,
        alternative,
    };

    FastPathResult::Success(HgvsVariant::Cds(CdsVariant {
        accession,
        gene_symbol: None,
        loc_edit: LocEdit::new(interval, edit),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_refseq_genomic_substitution() {
        match try_fast_path("NC_000001.11:g.12345A>G") {
            FastPathResult::Success(variant) => {
                assert!(matches!(variant, HgvsVariant::Genome(_)));
                if let HgvsVariant::Genome(g) = variant {
                    assert_eq!(&*g.accession.prefix, "NC");
                    assert_eq!(&*g.accession.number, "000001");
                    assert_eq!(g.accession.version, Some(11));
                }
            }
            FastPathResult::Fallback => panic!("Expected success for RefSeq genomic"),
        }
    }

    #[test]
    fn test_refseq_cds_substitution() {
        match try_fast_path("NM_000088.3:c.459A>G") {
            FastPathResult::Success(variant) => {
                assert!(matches!(variant, HgvsVariant::Cds(_)));
                if let HgvsVariant::Cds(c) = variant {
                    assert_eq!(&*c.accession.prefix, "NM");
                    assert_eq!(&*c.accession.number, "000088");
                    assert_eq!(c.accession.version, Some(3));
                }
            }
            FastPathResult::Fallback => panic!("Expected success for RefSeq CDS"),
        }
    }

    #[test]
    fn test_ensembl_genomic_substitution() {
        match try_fast_path("ENSG00000141510.5:g.12345A>G") {
            FastPathResult::Success(variant) => {
                assert!(matches!(variant, HgvsVariant::Genome(_)));
                if let HgvsVariant::Genome(g) = variant {
                    assert_eq!(&*g.accession.prefix, "ENSG");
                    assert!(g.accession.ensembl_style);
                }
            }
            FastPathResult::Fallback => panic!("Expected success for Ensembl genomic"),
        }
    }

    #[test]
    fn test_ensembl_cds_substitution() {
        match try_fast_path("ENST00000012345.1:c.100A>G") {
            FastPathResult::Success(variant) => {
                assert!(matches!(variant, HgvsVariant::Cds(_)));
            }
            FastPathResult::Fallback => panic!("Expected success for Ensembl CDS"),
        }
    }

    #[test]
    fn test_lrg_substitution() {
        match try_fast_path("LRG_1:g.12345A>G") {
            FastPathResult::Success(variant) => {
                assert!(matches!(variant, HgvsVariant::Genome(_)));
            }
            FastPathResult::Fallback => panic!("Expected success for LRG"),
        }
    }

    #[test]
    fn test_assembly_substitution() {
        match try_fast_path("GRCh38(chr1):g.12345A>G") {
            FastPathResult::Success(variant) => {
                assert!(matches!(variant, HgvsVariant::Genome(_)));
            }
            FastPathResult::Fallback => panic!("Expected success for assembly"),
        }
    }

    #[test]
    fn test_hg_assembly_substitution() {
        match try_fast_path("hg38(chr1):g.12345A>G") {
            FastPathResult::Success(variant) => {
                assert!(matches!(variant, HgvsVariant::Genome(_)));
            }
            FastPathResult::Fallback => panic!("Expected success for hg assembly"),
        }
    }

    #[test]
    fn test_fallback_for_complex_patterns() {
        // Uncertain intronic offset still defers (plain intronic/UTR subs are
        // now fast-pathed; `+?`/`-?` are not)
        assert!(matches!(
            try_fast_path("NM_000088.3:c.100+?A>G"),
            FastPathResult::Fallback
        ));

        // Unknown / parenthesized positions still defer
        assert!(matches!(
            try_fast_path("NM_000088.3:c.?A>G"),
            FastPathResult::Fallback
        ));

        // Insertion (plain del/dup are now fast-pathed, but ins still defers —
        // its inserted sequence is parsed by the generic parser)
        assert!(matches!(
            try_fast_path("NC_000001.11:g.12345_12346insA"),
            FastPathResult::Fallback
        ));

        // Deletion-insertion (delins) still defers
        assert!(matches!(
            try_fast_path("NC_000001.11:g.12345_12346delinsAC"),
            FastPathResult::Fallback
        ));

        // Deletion with an explicit trailing sequence/length still defers
        assert!(matches!(
            try_fast_path("NC_000001.11:g.12345delA"),
            FastPathResult::Fallback
        ));

        // Protein
        assert!(matches!(
            try_fast_path("NP_000079.2:p.Arg100Gly"),
            FastPathResult::Fallback
        ));
    }

    #[test]
    fn test_fast_path_intronic_and_utr_substitution() {
        use crate::hgvs::parser::variant::parse_variant;
        // Intronic and UTR CDS substitutions are fast-pathed; each must equal
        // the generic parser's output exactly.
        let cases = [
            "NM_000088.3:c.100+5A>G",  // intronic, downstream
            "NM_000088.3:c.100-5A>G",  // intronic, upstream
            "NM_000088.3:c.*100A>G",   // 3' UTR
            "NM_000088.3:c.-50A>G",    // 5' UTR
            "NM_000088.3:c.*100+5A>G", // 3' UTR with intronic offset
            "NM_000088.3:c.-50-5A>G",  // 5' UTR with intronic offset
            "ENST00000357033.8:c.100+1G>A",
        ];
        for input in cases {
            let fast = match try_fast_path(input) {
                FastPathResult::Success(v) => v,
                FastPathResult::Fallback => panic!("expected fast-path Success for {input:?}"),
            };
            let generic = parse_variant(input).expect("generic should parse");
            assert_eq!(fast, generic, "fast-path/generic mismatch for {input:?}");
        }
    }

    #[test]
    fn test_fast_path_plain_deletion_and_duplication() {
        use crate::hgvs::parser::variant::parse_variant;
        // Plain deletions/duplications (no trailing sequence or length) are
        // fast-pathed; each must equal the generic parser's output exactly.
        let cases = [
            "NC_000001.11:g.12345del",
            "NC_000001.11:g.12345_12346del",
            "NC_000001.11:g.12345dup",
            "NC_000001.11:g.12345_12346dup",
            "NM_000088.3:c.459del",
            "NM_000088.3:c.459_460del",
            "NM_000088.3:c.459dup",
            "NM_000088.3:c.459_460dup",
            "ENST00000357033.8:c.100del",
            "LRG_1:g.5000del",
            "GRCh38(chr1):g.12345del",
        ];
        for input in cases {
            let fast = match try_fast_path(input) {
                FastPathResult::Success(v) => v,
                FastPathResult::Fallback => panic!("expected fast-path Success for {input:?}"),
            };
            let generic = parse_variant(input).expect("generic should parse");
            assert_eq!(fast, generic, "fast-path/generic mismatch for {input:?}");
        }
    }

    #[test]
    fn test_fallback_for_unknown_patterns() {
        // Unknown accession type
        assert!(matches!(
            try_fast_path("ABC_12345.1:g.100A>G"),
            FastPathResult::Fallback
        ));

        // Too short
        assert!(matches!(
            try_fast_path("N:g.1A>G"),
            FastPathResult::Fallback
        ));
    }
}
