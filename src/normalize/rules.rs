//! HGVS normalization rules
//!
//! Additional rules for HGVS-compliant normalization beyond shuffling.
//!
//! # Coordinate Systems
//!
//! This module uses two coordinate systems:
//!
//! | Context | Basis | Type |
//! |---------|-------|------|
//! | Array indexing | 0-based | `usize` |
//! | HGVS output positions | 1-based | `u64` |
//!
//! Functions that take `pos: u64` parameters use 0-based positions for internal
//! array indexing. Return values marked with "1-indexed" comments use HGVS-style
//! 1-based positions.
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use std::str::FromStr;

use crate::coords::{hgvs_pos_to_index, index_to_hgvs_pos};
use crate::error::FerroError;
use crate::hgvs::edit::{InsertedPart, InsertedSequence, NaEdit, Sequence};
use crate::reference::ReferenceProvider;
use crate::sequence::complement_base;

/// Check if an edit type needs normalization
///
/// "Needs normalization" gates entry into `normalize_na_edit`, which performs
/// both reference-allele validation AND 3'-shifting. Most listed kinds need the
/// shift; `MultiRepeat` and `Substitution` are exceptions — both are flagged
/// `true` for the *validation* half only, not the shift. A compound
/// multi-unit tract has no canonical shuffle target, so `normalize_na_edit`'s
/// `MultiRepeat` arm passes it through unchanged after validating
/// (`validate_multirepeat_tract`, #953). Substitutions are never shuffled
/// either way; routing them through validates the stated reference base
/// against the loaded reference (#1052) and returns the edit unchanged (see
/// the `Substitution` arm in `normalize_na_edit`'s `alt_seq` match). If the
/// `MultiRepeat` arm ever gains a real canonicalization, this flag already
/// routes it correctly; until then the two are reconciled and documented, not
/// silently contradictory.
pub fn needs_normalization(edit: &NaEdit) -> bool {
    if matches!(
        edit,
        NaEdit::Deletion { .. }
            | NaEdit::Insertion { .. }
            | NaEdit::Duplication { .. }
            | NaEdit::Delins { .. }
            | NaEdit::Inversion { .. }
            | NaEdit::Repeat { .. }
            | NaEdit::MultiRepeat { .. }
    ) {
        return true;
    }
    // Substitutions are routed through `normalize_na_edit` for the *validation*
    // half only (reference-base check), like `MultiRepeat` — never for shuffling.
    // A `ref == alt` sub additionally rewrites to identity (`=`); a real sub is
    // returned unchanged after validation (see the `Substitution` arm in
    // `normalize_na_edit`'s `alt_seq` match). Intronic subs are held OUT of this
    // path by the per-axis guards (they cannot be validated without genomic
    // projection and must not error) — see #1052.
    matches!(edit, NaEdit::Substitution { .. })
}

/// A substitution whose stated reference and alternative bases differ (i.e. not
/// the degenerate `ref == alt` identity case). Used by the per-axis intronic
/// guards to hold real substitutions out of the genomic intronic projection.
pub fn is_real_substitution(edit: &NaEdit) -> bool {
    matches!(
        edit,
        NaEdit::Substitution { reference, alternative } if reference != alternative
    )
}

/// A real substitution (see [`is_real_substitution`]) whose `LocEdit` wraps
/// it in `Mu::Uncertain` (HGVS `(...)`, e.g. `c.(10A>G)`).
///
/// Used by every axis normalizer (`normalize_genome`, `normalize_mt`,
/// `normalize_cds`, `normalize_tx`, `normalize_rna`) as a single shared
/// definition for the pre-validation carve-out: an uncertain/predicted-
/// wrapped substitution must stay a silent pass-through in every mode,
/// regardless of whether the stated ref matches the reference — validating
/// it would emit a new `RefSeqMismatch` (lenient) or a hard
/// `ReferenceMismatch` rejection (strict) on input that normalized silently
/// before #1052. `Mu::inner()` unwraps both `Certain` and `Uncertain`, so
/// the unwrapped `edit` alone can't distinguish them; the caller must pass
/// the `Mu` wrapper (typically `variant.loc_edit.edit`) alongside it.
///
/// This also protects against a second, independent hazard: every axis's
/// post-validation result construction (`LocEdit::new(...)`) unconditionally
/// emits `Mu::Certain`, silently dropping an `Uncertain` wrapper on
/// reconstruction. Short-circuiting here — before that reconstruction runs —
/// avoids relying on a fix to that (separate, pre-existing) defect.
pub fn is_uncertain_real_substitution(edit: &NaEdit, edit_mu: &crate::hgvs::Mu<NaEdit>) -> bool {
    is_real_substitution(edit) && edit_mu.is_uncertain()
}

/// The largest chance-match probability at which an inverted-copy match is
/// taken as evidence rather than coincidence.
///
/// Derived rather than chosen: at `1/N` over a representation-stability corpus
/// of `N` rows, the *expected number of false positives over that corpus* is
/// below one. A bare length cutoff would be a number with no such account behind
/// it, and the spec supplies none — `DNA/insertion.md:22` scopes the copy-range
/// form to "larger insert sequences" and defines "larger" nowhere.
///
/// **Read it as an order of magnitude — about `1e-5` — not as a row count.**
/// The literal denominator is a *retired* figure:
/// `examples/dump_normalized_corpus.rs` records 85,642 as one entry in a
/// history, superseded by 86,398, and says in terms that it deliberately does
/// not name a current entry because "the current total is whatever the generator
/// prints today". Two more reasons the denominator will not bear close reading:
/// it counts the protein axis, while only nucleotide rows can carry an inserted
/// revcomp payload and only `g.` rows reach this code; and the account treats
/// the corpus as `N` independent trials while **two** candidate spans are tested
/// per insertion junction, so the budget is optimistic by at least 2×. None of
/// that moves the order of magnitude, which is all this constant claims.
///
/// Under even base composition this behaves like a floor of nine bases
/// (`0.25^9 = 3.8e-6`, `0.25^8 = 1.5e-5`), but it is deliberately not written as
/// one: in a skewed or homopolymeric window a short match is much likelier than
/// `4^-n`, and a flat length rule would wave exactly those through.
///
/// The **effect** is what is guarded, not the derivation — see
/// `tests/it/recommended_form_pins.rs::the_coincidence_floor_sits_between_eight_and_nine_bases`.
/// Before that pin existed, weakening this constant 85× left all 41 tests in the
/// three dedicated corpora green.
const MAX_CHANCE_MATCH_PROBABILITY: f64 = 1.0 / 85_642.0;

/// Half the width of the window whose base composition
/// [`chance_match_probability`] keys on, in bases either side of the insertion
/// junction.
///
/// The composition estimate must be a property of the **reference**, not of the
/// caller. `ref_seq` is the normalization window — `[start − w, end + w)` for
/// `NormalizeConfig::window_size`, a public field, then grown by
/// `normalize_in_grown_window` until the 3'-shift settles — so keying on the
/// composition of the whole slice made the emitted description depend on a knob
/// that exists to bound the shuffle and has no statistical meaning. On a contig
/// with an AT→GC boundary, one variant rendered `ins<literal>` at
/// `window_size = 80` and `ins<range>inv` at `window_size = 120`.
///
/// 50 rather than the default `window_size` of 100 so that the fixed window sits
/// comfortably *inside* the supplied slice at the default setting and at every
/// larger one, including after the growth loop enlarges it. Near a contig end
/// the window still clamps, but it clamps to where the reference ends rather
/// than to how the caller was configured.
const COMPOSITION_HALF_WIDTH: usize = 50;

/// Probability that an unrelated payload would reverse-complement-match `span`
/// by chance, under the base composition of `window`.
///
/// For each base of the span, an unrelated payload matches at that position if
/// it happens to carry that base's complement, which it does with the frequency
/// of that complement in the surrounding sequence. Positions are treated as
/// independent, which is the conservative direction here: real sequence is
/// correlated, so the true chance of a match is *higher* than this product and
/// the test is if anything too permissive rather than too strict.
///
/// Counts are Laplace-smoothed so a base absent from the window yields a small
/// probability rather than zero — without it, a match on a base the window never
/// shows would read as infinitely significant.
fn chance_match_probability(window: &[u8], span: &[u8]) -> f64 {
    fn index(base: u8) -> Option<usize> {
        match base.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        }
    }
    fn complement(base: u8) -> u8 {
        match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other,
        }
    }

    let mut counts = [0usize; 4];
    let mut total = 0usize;
    for &base in window {
        if let Some(i) = index(base) {
            counts[i] += 1;
            total += 1;
        }
    }
    if total == 0 {
        return 1.0;
    }

    let mut probability = 1.0f64;
    for &base in span {
        let Some(i) = index(complement(base)) else {
            return 1.0;
        };
        probability *= (counts[i] as f64 + 1.0) / (total as f64 + 4.0);
    }
    probability
}

/// The reference span whose reverse complement is `payload`, when that span
/// immediately abuts the insertion junction on its **5' side** and the match is
/// unlikely enough to be evidence of an inverted duplication.
///
/// The inverted sibling of [`insertion_is_duplication`]. `DNA/inversion.md:69`
/// gives the spelling — `g.234_235ins123_234inv` — and five Notes restate it,
/// each rejecting `dupinv`.
///
/// Three conditions, only the last of which is ours:
///
/// 1. **More than one nucleotide.** `DNA/inversion.md` defines an inversion as a
///    change where "**more than one nucleotide** replacing the original sequence
///    are the reverse complement of the original sequence", so a single base is
///    outside the concept entirely.
/// 2. **Not a plain duplication.** A payload that equals its neighbour *and* its
///    neighbour's reverse complement is self-complementary, which makes it a
///    duplication — and `general.md:55` ranks duplication (4) above insertion
///    (5), so it belongs to [`insertion_is_duplication`], not here. Declining
///    keeps `#179`'s self-complementary cases unchanged.
/// 3. **Unlikely by chance**, per [`MAX_CHANCE_MATCH_PROBABILITY`].
///
/// `gap` is the 1-based position 5' of the junction, in the same frame `ref_seq`
/// is indexed in; the returned span is 1-based inclusive in that same frame, and
/// the caller converts it to axis coordinates exactly as it does the anchor.
/// Keeping that conversion in one place is deliberate — a payload that travelled
/// a different path from its own anchor is how a range insert ends up naming
/// bases far from the ones it means.
pub(crate) fn inverted_adjacent_copy_span(
    ref_seq: &[u8],
    gap: u64,
    payload: &[u8],
) -> Option<(u64, u64)> {
    let len = payload.len();
    // `DNA/inversion.md`: an inversion is a change where "more than one
    // nucleotide" replacing the original are its reverse complement, so a single
    // base is outside the concept.
    if len < 2 {
        return None;
    }
    let gap_idx = usize::try_from(gap).ok()?;
    if gap_idx > ref_seq.len() {
        return None;
    }

    // `:69` gives two spellings, "depending on whether the inverted copy is 5'
    // or 3' of the original copy". The 5'-abutting source (`g.234_235ins
    // 123_234inv`, the span ENDING at the junction) is tried first, and a match
    // there wins, because the spec supplies no tie-break and preferring one
    // deterministically is what keeps the answer independent of which side was
    // examined first.
    //
    // Both sides match exactly when the `2 * len` bases straddling the junction
    // are a TANDEM REPEAT of period `len` — candidate 0 is `A = ref[gap-len..gap]`
    // and candidate 1 is `B = ref[gap..gap+len]`, and both reverse-complement-
    // match one payload iff `revcomp(A) == revcomp(B)` iff `A == B`. Nothing
    // palindromic is required, and the difference is operational: tandem repeats
    // of 9–30 bases are ubiquitous in real sequence while palindromes of that
    // length are not, so this tie-break is exercised far more often than a
    // self-complementarity reading would suggest. Pinned by
    // `a_tandem_repeat_makes_both_candidates_match_and_the_five_prime_span_wins`.
    let candidates = [
        gap_idx.checked_sub(len).map(|start| (start, gap_idx)),
        gap_idx
            .checked_add(len)
            .filter(|end| *end <= ref_seq.len())
            .map(|end| (gap_idx, end)),
    ];

    for (start_idx, end_idx) in candidates.into_iter().flatten() {
        let span_bytes = &ref_seq[start_idx..end_idx];
        let Ok(span) = std::str::from_utf8(span_bytes) else {
            continue;
        };
        if !crate::sequence::reverse_complement(span)
            .as_bytes()
            .eq_ignore_ascii_case(payload)
        {
            continue;
        }
        // A payload equal to its neighbour AND to that neighbour's reverse
        // complement is self-complementary, which makes this a duplication —
        // `general.md:55` ranks duplication (4) above insertion (5), so it
        // belongs to `insertion_is_duplication`. Declining keeps #179's
        // self-complementary cases unchanged.
        if payload.eq_ignore_ascii_case(span_bytes) {
            continue;
        }
        // Key the coincidence floor on a fixed window centred on the junction,
        // not on `ref_seq` — see `COMPOSITION_HALF_WIDTH`.
        let composition_start = gap_idx.saturating_sub(COMPOSITION_HALF_WIDTH);
        let composition_end = gap_idx
            .saturating_add(COMPOSITION_HALF_WIDTH)
            .min(ref_seq.len());
        let composition = &ref_seq[composition_start..composition_end];
        if chance_match_probability(composition, span_bytes) > MAX_CHANCE_MATCH_PROBABILITY {
            continue;
        }
        return Some((start_idx as u64 + 1, end_idx as u64));
    }
    None
}

/// Check if an insertion should be represented as a duplication
///
/// In HGVS, if an inserted sequence is identical to the sequence
/// immediately 5' (before) or 3' (after) of the insertion point,
/// the variant should be represented as a duplication.
///
/// IMPORTANT: This only applies to INSERTIONS, not deletions.
/// Deletions always stay as deletions - they just shift 3'.
pub fn insertion_is_duplication(ref_seq: &[u8], pos: u64, inserted_seq: &[u8]) -> bool {
    let ins_len = inserted_seq.len();
    let pos_idx = pos as usize;

    // Guard against positions beyond the reference sequence
    if ref_seq.is_empty() || pos_idx > ref_seq.len() {
        return false;
    }

    // Check if sequence before the insertion point matches the inserted sequence
    // (this would be a duplication of the preceding sequence)
    //
    // Case-insensitive on both sides (#1318): `ref_seq` may be soft-masked while
    // `inserted_seq` comes from the author's upper-case description, and a raw
    // byte test then reports "not a duplication" for one that is — the same
    // class of miss the delins path carried, measured on `acgtacgt` + `ACGT`.
    if pos_idx >= ins_len && pos_idx <= ref_seq.len() {
        let before_start = pos_idx - ins_len;
        if ref_seq[before_start..pos_idx].eq_ignore_ascii_case(inserted_seq) {
            return true;
        }
    }

    // Also check if sequence after matches (for 5' duplications)
    if pos_idx + ins_len <= ref_seq.len()
        && ref_seq[pos_idx..pos_idx + ins_len].eq_ignore_ascii_case(inserted_seq)
    {
        return true;
    }

    false
}

/// In-phase correctness gate for insertion->dup / insertion->repeat conversion
/// (issue #882).
///
/// Returns `true` iff inserting `inserted_seq` at `ins_pos` (0-based index of
/// the base immediately 5' of the insertion point; the insertion sits between
/// `ins_pos` and `ins_pos + 1`) yields exactly the same full sequence as the
/// tandem edit that replaces `ref_seq[region_start..region_end_excl]` with
/// `unit` repeated `total_count` times.
///
/// `region_start`, `region_end_excl`, `unit`, and `total_count` are the
/// POST-alignment values the caller is about to emit (the dup/repeat's own
/// window and unit), so this validates the literal emitted edit — never the raw
/// pre-selection tract. A conversion is spec-valid only when the tandem form
/// decodes to the identical sequence as the input insertion.
///
/// Total: never panics; returns `false` on any out-of-range or degenerate
/// argument.
fn tandem_edit_preserves_insertion(
    ref_seq: &[u8],
    ins_pos: usize,
    inserted_seq: &[u8],
    region_start: usize,
    region_end_excl: usize,
    unit: &[u8],
    total_count: u64,
) -> bool {
    // Bounds guards. `ref_seq[..=ins_pos]` is `ref_seq[..ins_pos + 1]`, which
    // panics when `ins_pos == ref_seq.len()`, so guard `ins_pos >= len`.
    if unit.is_empty()
        || ins_pos >= ref_seq.len()
        || region_start > region_end_excl
        || region_end_excl > ref_seq.len()
    {
        return false;
    }

    // Sequence produced by applying the insertion.
    let mut from_insertion = Vec::with_capacity(ref_seq.len() + inserted_seq.len());
    from_insertion.extend_from_slice(&ref_seq[..=ins_pos]);
    from_insertion.extend_from_slice(inserted_seq);
    from_insertion.extend_from_slice(&ref_seq[ins_pos + 1..]);

    // Sequence produced by the emitted tandem edit: replace the region with
    // `unit` repeated `total_count` times. Compute the capacity with checked
    // arithmetic so a degenerate `total_count` can't truncate (32-bit `usize`)
    // or overflow into a panic; per the contract, bail out `false` if it would.
    let Ok(total_count_usize) = usize::try_from(total_count) else {
        return false;
    };
    let capacity = unit
        .len()
        .checked_mul(total_count_usize)
        .and_then(|body| body.checked_add(region_start))
        .and_then(|head| head.checked_add(ref_seq.len().saturating_sub(region_end_excl)));
    let Some(capacity) = capacity else {
        return false;
    };
    let mut from_tandem = Vec::with_capacity(capacity);
    from_tandem.extend_from_slice(&ref_seq[..region_start]);
    for _ in 0..total_count {
        from_tandem.extend_from_slice(unit);
    }
    from_tandem.extend_from_slice(&ref_seq[region_end_excl..]);

    from_insertion == from_tandem
}

/// Determine the canonical representation of an indel
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CanonicalForm {
    /// Standard deletion
    Deletion,
    /// Duplication (del of a repeated sequence)
    Duplication,
    /// Deletion-insertion
    Delins,
    /// Standard insertion
    Insertion,
    /// Repeat notation (e.g., A[9] for homopolymer)
    Repeat {
        /// The repeated base
        base: u8,
        /// Total count after the edit
        count: u64,
        /// Start position (1-based HGVS, inclusive) of the repeat region
        start: u64,
        /// End position (1-based HGVS, inclusive) of the repeat region
        end: u64,
    },
}

/// Result of analyzing a potential repeat region
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepeatAnalysis {
    /// Whether this is a homopolymer (single-base repeat)
    pub is_homopolymer: bool,
    /// The repeated base (for homopolymer)
    pub base: Option<u8>,
    /// Start position (0-based) of the repeat in reference
    pub ref_start: usize,
    /// End position (0-based, exclusive) of the repeat in reference
    pub ref_end: usize,
    /// Count of repeats in reference
    pub ref_count: u64,
}

/// Analyze a position to find the extent of any homopolymer repeat
///
/// Given a position in the reference, finds the full extent of any
/// single-base repeat that includes this position.
///
/// Returns None if the position is not in a repeat (count < 2).
pub fn find_homopolymer_at(ref_seq: &[u8], pos: usize) -> Option<RepeatAnalysis> {
    if pos >= ref_seq.len() {
        return None;
    }

    let base = ref_seq[pos];

    // Find start of the repeat (scan left)
    let mut start = pos;
    while start > 0 && ref_seq[start - 1] == base {
        start -= 1;
    }

    // Find end of the repeat (scan right)
    let mut end = pos + 1;
    while end < ref_seq.len() && ref_seq[end] == base {
        end += 1;
    }

    let count = (end - start) as u64;

    // Only consider it a repeat if count >= 2
    if count < 2 {
        return None;
    }

    Some(RepeatAnalysis {
        is_homopolymer: true,
        base: Some(base),
        ref_start: start,
        ref_end: end,
        ref_count: count,
    })
}

/// Check if an insertion in a tandem repeat should become repeat notation.
///
/// Per HGVS spec: when an insertion adds 2 or more copies of a tandem
/// repeat unit already present adjacent to the insertion point, the
/// variant is expressed as repeat notation `unit[N+k]`. A single-copy
/// addition is left as a duplication (handled by `insertion_is_duplication`).
///
/// `is_coding` enables the spec's codon-frame exception (repeated.md): in
/// c. context, repeat notation requires `unit_len % 3 == 0`. When the gate
/// blocks the rewrite this returns `None` and the caller falls back to a
/// literal `ins` description. A second gate, `tandem_edit_preserves_insertion`
/// (#882), likewise returns `None` when the emitted repeat is out of phase —
/// i.e. it would decode to a different sequence than the input insertion — so
/// an out-of-phase insertion stays a literal `ins` rather than a spurious
/// repeat.
///
/// Returns `Some((first_byte, total_count, ref_start_1based,
/// ref_end_1based, unit_bytes))` when repeat notation applies, `None`
/// otherwise. `first_byte` is preserved for backwards-compat with the
/// homopolymer caller; `unit_bytes` is the full tandem unit.
///
/// # Direction-awareness (audit per #357)
///
/// Unlike `insertion_to_duplication` (which picks an explicit
/// direction-aligned slot), `insertion_to_repeat` returns the canonical
/// `unit[N+k]` window spanning the full reference tract. The window covers the
/// entire tract regardless of where shuffle would have placed the original
/// insertion; only the tract's *phase* is direction-dependent, delegated to
/// [`normalize_repeat`] (see below). For a homopolymer or a tract bounded by
/// non-unit bases the 3'- and 5'-most phases coincide, so the displayed range is
/// still direction-invariant there (#357,
/// `insertion_to_repeat_output_is_direction_invariant_on_n_axis`); an ambiguous
/// alternating tract, by contrast, is phased to the direction-appropriate end of
/// the tract (#8), which is what makes the emitted repeat a fixed point under
/// re-normalization.
///
/// # The codon gate is answered for the tract, not the input (#1210)
///
/// `is_coding` is the caller's verdict for the *input* edit's span. Where
/// `cds_span` is supplied it takes precedence, and the gate is re-answered for
/// the tandem tract this function discovers — the span the emitted repeat
/// actually occupies. The two differ exactly when the edit shifts across
/// `cds_end`, which is where the non-idempotency lived. Pass `None` (genomic,
/// `n.`, or any caller with no CDS context) to keep `is_coding` verbatim.
///
/// This function owns only the insertion-specific work: rotation-aware
/// discovery of the adjacent tandem tract (`find_tandem_extent`) and the #882
/// preservation gate. The canonicalization arithmetic — the 3'-phase
/// (`three_prime_align_tract`), the copy count, and the c.-codon gate — is
/// **delegated to [`normalize_repeat`]** (#866), making that function the single
/// source of truth for ins→repeat canonicalization (the #852 non-idempotency
/// lived in exactly this seam). The two keep separate tract-*finders*
/// (`find_tandem_extent` here, [`count_tandem_repeats`] there) because their
/// anchoring contracts genuinely differ; the cross-check test
/// `insertion_to_repeat_agrees_with_normalize_repeat_phase` pins their agreement.
pub fn insertion_to_repeat(
    ref_seq: &[u8],
    pos: u64,
    inserted_seq: &[u8],
    is_coding: bool,
    cds_span: Option<(u64, u64)>,
    direction: crate::normalize::config::ShuffleDirection,
) -> Option<(u8, u64, u64, u64, Vec<u8>)> {
    if inserted_seq.is_empty() {
        return None;
    }

    let base_unit = smallest_repeat_unit(inserted_seq);
    let added_copies = (inserted_seq.len() / base_unit.len()) as u64;
    if added_copies < 2 {
        // Single-copy addition is a duplication, not repeat notation.
        return None;
    }

    // #1846: the same early exit `insertion_to_duplication` carries, for the
    // same reason. This loop rotates and rebuilds a `base_unit.len()`-byte `Vec`
    // once per base of the unit; when the unit outgrows `ref_seq` every
    // rotation's `find_tandem_extent` returns `None`, so `best` is already
    // `None` and the Θ(unit²) work is guaranteed futile. An aperiodic payload
    // returns above (`added_copies == 1` routes to the dup path), so this fires
    // only for a payload that really is 2+ copies of an oversized unit — e.g.
    // two copies of a multi-megabase unit — which the `added_copies` guard does
    // not bound. See the fuller note on `insertion_to_duplication`.
    if base_unit.len() > ref_seq.len() {
        return None;
    }

    // The c.-context codon gate (repeated.md: unit_len % 3 == 0) is applied
    // inside `normalize_repeat` below, not here — see the delegation note.

    // The inserted sequence may be written as any cyclic rotation of the
    // reference repeat unit (e.g. `insCAGCAG` against a `GCA` tract). Try every
    // rotation of `base_unit` and pick the one that yields the largest
    // contiguous tandem run anchored at the insertion point.
    //
    // Each rotation is resolved to the repeat it would actually emit — and
    // gated — before it can win (#1355). Scoring first and gating only the
    // winner let a rotation that abuts a longer but out-of-phase tract discard a
    // valid in-phase one and take the whole result down with it, which is
    // byte-for-byte the defect #1349 corrected in `insertion_to_duplication`.
    // The two functions are a deliberate symmetric pair, so they now share the
    // shape as well: resolve per candidate, `continue` on rejection.
    let u_len = base_unit.len();
    let mut best: Option<RepeatCandidate> = None;
    for r in 0..u_len {
        let mut rotated = Vec::with_capacity(u_len);
        rotated.extend_from_slice(&base_unit[r..]);
        rotated.extend_from_slice(&base_unit[..r]);
        let Some((ref_start, ref_count)) = find_tandem_extent(ref_seq, pos as usize, &rotated)
        else {
            continue;
        };
        if ref_count == 0 {
            continue;
        }
        // Ties still fall to iteration order, as before: strictly-greater keeps
        // the first rotation found. Adding a direction-aware tie-break here
        // would be a separate behaviour change from the ordering fix
        // (`insertion_to_duplication` grew one for #403; this function has never
        // had one).
        if best
            .as_ref()
            .is_some_and(|current| ref_count <= current.ref_count)
        {
            // Cannot become `best` whatever the gates say, and resolving it
            // costs a `normalize_repeat` call plus two `ref_seq`-sized buffers.
            continue;
        }
        let Some(candidate) = resolve_repeat_rotation(
            ref_seq,
            pos,
            inserted_seq,
            ref_start,
            ref_count,
            rotated,
            added_copies,
            is_coding,
            cds_span,
            direction,
        ) else {
            continue;
        };
        best = Some(candidate);
    }

    let best = best?;
    Some((
        best.rotated_unit[0],
        best.total_count,
        index_to_hgvs_pos(best.start_idx),
        best.emitted_end_hgvs,
        best.rotated_unit,
    ))
}

/// One rotation's resolved repeat emission, after both gates have accepted it.
struct RepeatCandidate {
    /// Copies of the rotated unit in the matched reference tract — the score.
    ref_count: u64,
    /// 0-based start of the emitted repeat window.
    start_idx: usize,
    /// 1-based inclusive HGVS end of the emitted repeat window.
    emitted_end_hgvs: u64,
    /// The phase-aligned unit `normalize_repeat` settled on, which may differ
    /// from the rotation fed in.
    rotated_unit: Vec<u8>,
    /// Total copies the emitted repeat asserts.
    total_count: u64,
}

/// Resolve one cyclic rotation to the repeat it would emit, or `None` if either
/// gate rejects it.
///
/// Split out of [`insertion_to_repeat`] so every rotation is judged on the values
/// the caller would actually receive (#1355), mirroring [`align_dup_slot`]'s role
/// in `insertion_to_duplication` (#1349). The difference between the two is why
/// this one is bigger: a duplication slot is pure arithmetic, whereas a repeat
/// has to go through `normalize_repeat`, which owns the window, the count and the
/// codon decision.
///
/// # Both rejections are a `continue` at the call site, and both mean `None`
///
/// - The **#882 phase gate** (`tandem_edit_preserves_insertion`): a rotation that
///   fails it is not a repeat of *this* insertion at all, so it must not compete,
///   and must not shadow a lower-scoring rotation that is. This is the #1355 fix.
/// - The **#1210 codon gate**, reached as `RepeatNormResult::Insertion`: in a `c.`
///   context a non-codon-length unit may not be written as repeat notation, so
///   the caller emits the literal `ins` instead.
///
/// The codon gate previously aborted the whole function (`return None`) because
/// it was only ever reached by the winner. Inside a selection loop that would be
/// wrong in a new way — a losing rotation could abort a search that had a valid
/// answer — so it became a `continue` too. **This is observably different only
/// when the highest-scoring rotation is codon-gated while a lower-scoring one is
/// not**, and the two answers there are "literal `ins`" versus "a repeat over the
/// other tract". `continue` is the right one, for the same reason the phase gate
/// is: all rotations of a unit share its length, so the codon gate can differ
/// between them only through `gate_is_coding` — that is, through which *tract*
/// each rotation found. The gate is a statement about a tract, not about the
/// insertion, and a tract that clears both gates describes the same edit
/// correctly.
///
/// The FAIL-set ledgers this could have moved are the manifest-gated conformance
/// axes — `axis_normalized`, `axis_genomic`, `axis_errors` and both idempotency
/// axes. They were re-run against a real prepared reference when this landed; the
/// measurement lives in the commit message rather than here, since a doc comment
/// asserting "the suite is green" rots the moment anything else changes.
#[allow(clippy::too_many_arguments)]
fn resolve_repeat_rotation(
    ref_seq: &[u8],
    pos: u64,
    inserted_seq: &[u8],
    ref_start: usize,
    ref_count: u64,
    unit: Vec<u8>,
    added_copies: u64,
    is_coding: bool,
    cds_span: Option<(u64, u64)>,
    direction: crate::normalize::config::ShuffleDirection,
) -> Option<RepeatCandidate> {
    let ref_end_excl = ref_start + ref_count as usize * unit.len();

    // #1210: answer the codon gate for the span the REPEAT will occupy, not for
    // the span the input insertion occupied.
    //
    // `is_coding` is precomputed by the caller for the *input* edit. That is the
    // wrong question here: the discovered tandem tract is a different span, and a
    // CDS-resident input can 3'-shift clean out of the CDS, so the two disagree
    // precisely at the `cds_end` boundary. `r.15delinsgaa` on a transcript whose
    // CDS ends at `r.15` reduces to an insertion that lands in the 3'UTR at
    // `r.*1_*2`; gating it on the input's CDS-residency suppressed the repeat on
    // pass 1, while pass 2 — whose input was by then UTR-resident — emitted
    // `r.*1a[3]`, so normalization was not idempotent.
    //
    // `cds_span` exists for exactly this ("lets a site that decides about a
    // DIFFERENT span ask about that span instead", #1185); this is the second
    // site to need it. `None` keeps the caller's precomputed verdict, so callers
    // with no CDS context (genomic / `n.`) and the unit tests below are
    // unaffected.
    //
    // Per rotation since #1355, which is also more correct than it was: each
    // rotation finds its own tract, so each gets its own residency answer instead
    // of the winner's.
    //
    // #1389: the window to ask about is the one the repeat is **emitted** over,
    // which is not the tract `find_tandem_extent` discovered. `normalize_repeat`
    // owns a direction-aware phase alignment and can return a shifted window, and
    // the two disagree exactly at a CDS boundary: for `CTGTGTT` with
    // `cds_span = (3, 6)`, the raw `TG` tract is 1-based `2..=5` and so reads as
    // non-resident (`2 < 3`), switching the gate off — while the emitted window is
    // `3..=6`, wholly CDS-resident, and got written as `GT[5]` with a two-base
    // unit that `repeated.md`'s codon rule forbids. The 5' direction of the same
    // fixture emits `TG` over `2..=5`, which genuinely is not resident, so before
    // this the two directions disagreed about whether the repeat was admissible.
    //
    // So probe for the emitted window with the gate lifted, then delegate with the
    // gate derived from it. Two calls rather than re-implementing the codon test
    // here, which keeps `normalize_repeat` the single source of truth for the
    // codon decision (#866) — the probe is only ever used to learn *which window*
    // to ask about. This is sound because `is_coding` enters `normalize_repeat` at
    // exactly one point, `codon_blocks_repeat`, which routes the result variant
    // only: the window, unit and counts are all computed before it and do not
    // depend on it. `insertion_to_repeat`'s `ref_count <= best.ref_count` skip
    // means only score-improving rotations pay for the second call.
    //
    // Frames line up: `normalize_repeat` returns 1-based inclusive HGVS positions
    // and `cds_span` is 1-based inclusive in the same transcript frame, so the
    // comparison is direct. With no `cds_span` there is no window question to ask
    // — the caller's precomputed verdict stands and no probe is issued, so callers
    // with no CDS context (genomic / `n.`) and the unit tests below are
    // unaffected.
    let target = ref_count + added_copies;
    let end_pos_incl = ref_end_excl - 1;
    let gate_is_coding = match cds_span {
        Some((cds_start, cds_end)) => {
            let probe = normalize_repeat(
                ref_seq,
                ref_start,
                end_pos_incl,
                &unit,
                target,
                // Gate lifted: this call answers "which window?", not "may it be
                // repeat notation?". Passing the real flag would let a gated
                // result hide the window the gate has to be asked about.
                false,
                direction,
            );
            match probe {
                RepeatNormResult::Repeat { start, end, .. } => start >= cds_start && end <= cds_end,
                // Unreachable for the same reason the delegation below asserts:
                // `added_copies >= 2`. Fall back to the caller's verdict rather
                // than panicking here, so the assertion stays in one place.
                _ => is_coding,
            }
        }
        None => is_coding,
    };

    // #866: delegate the canonical-window + count + codon-gate decision to
    // `normalize_repeat` (the single source of truth for ins→repeat
    // canonicalization). Re-express the discovered tract as the well-formed
    // explicit repeat range `unit[ref_count]` over `[ref_start, ref_end_excl)` and
    // expand it by the inserted copies. Because `find_tandem_extent` and
    // `count_tandem_repeats` walk the tract identically, `normalize_repeat`
    // re-discovers the same tract, absorbs no flanks, and applies the SAME
    // direction-aware phase alignment + codon gate — so the emitted
    // window/unit/count agree with `normalize_repeat` on the same reference tract
    // in the same direction. `added_copies >= 2` ⇒ the target is
    // `>= ref_count + 2`, so only `Repeat` (or the codon-gated `Insertion` →
    // `None`) is reachable; `Deletion`/`Duplication`/`Unchanged` cannot occur.
    //
    // `target` and `end_pos_incl` are bound above, where the #1389 probe needs
    // them; this call is the authoritative one.
    let (start_idx, emitted_end_hgvs, rotated_unit, total_count) = match normalize_repeat(
        ref_seq,
        ref_start,
        end_pos_incl,
        &unit,
        target,
        gate_is_coding,
        direction,
    ) {
        RepeatNormResult::Repeat {
            start,
            end,
            sequence,
            count,
        } => (hgvs_pos_to_index(start), end, sequence, count),
        // Codon gate (c. + non-codon unit) reroutes the expansion to `Insertion`;
        // the caller then emits the literal `ins` form. Matches the pre-#866 gate.
        RepeatNormResult::Insertion { .. } => return None,
        // `added_copies >= 2` ⇒ `target = ref_count + added_copies >= ref_count + 2`,
        // so a contraction / one-copy-expansion / identity result cannot arise here.
        // Assert it loudly rather than silently returning `None`, so a future change
        // that relaxes the `added_copies >= 2` guard is caught at test time.
        RepeatNormResult::Deletion { .. }
        | RepeatNormResult::Duplication { .. }
        | RepeatNormResult::Unchanged => unreachable!(
            "insertion_to_repeat delegates with added_copies >= 2, so target >= \
                 ref_count + 2; only Repeat or the codon-gated Insertion are reachable"
        ),
    };
    // `emitted_end_hgvs` is a 1-based INCLUSIVE HGVS position, numerically equal to
    // the 0-based EXCLUSIVE aligned end; recover the exclusive end for the #882 gate
    // with no double correction.
    let end_excl = hgvs_pos_to_index(emitted_end_hgvs) + 1;

    // #882: only emit the repeat when it decodes to the same sequence as the
    // input insertion.
    if !tandem_edit_preserves_insertion(
        ref_seq,
        pos as usize,
        inserted_seq,
        start_idx,
        end_excl,
        &rotated_unit,
        total_count,
    ) {
        return None;
    }

    Some(RepeatCandidate {
        ref_count,
        start_idx,
        emitted_end_hgvs,
        rotated_unit,
        total_count,
    })
}

/// Description of how a single-copy insertion can be re-expressed as a
/// 3'-rule canonical duplication.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct InsToDupResult {
    /// The phase-aligned tandem unit chosen by rotation iteration. May
    /// differ from the `inserted_seq` argument when the alt was spelled
    /// as a non-zero cyclic rotation of the reference unit (issue #132).
    pub unit: Vec<u8>,
    /// 1-based HGVS start of the direction-aligned unit-aligned duplication
    /// position inside the rotation-aligned reference tract (inclusive).
    /// Resolves to the most-3' slot when the caller passed
    /// `ShuffleDirection::ThreePrime`, the most-5' slot under `FivePrime`
    /// (closes #340 subgroup (i)).
    pub start: u64,
    /// 1-based HGVS end of the direction-aligned unit-aligned duplication
    /// position inside the rotation-aligned reference tract (inclusive).
    /// Always equals `start + unit.len() - 1`.
    pub end: u64,
    /// Number of full reference copies of `unit` in the chosen tract. The
    /// caller uses this to decide whether the tract-aligned dup position
    /// should win over a post-shuffle simple-dup candidate: a multi-copy
    /// tract (`ref_count >= 2`) has a meaningful phase that must be
    /// preserved; a single-copy tract is just one unit, and a longer
    /// post-shuffle partial-match extension is more 3' (issue #180).
    pub ref_count: u64,
}

/// Compute the canonical 3'-rule duplication position for a single-copy
/// insertion whose alt is a (possibly rotated) tandem-repeat unit adjacent
/// to a reference tract.
///
/// # Why this helper exists (the `insertion_to_repeat` asymmetry)
///
/// The 3' shuffle (`shuffle.rs`) walks one base at a time. When the inserted
/// alt is a non-zero cyclic rotation of an adjacent reference unit, the very
/// first-base check fails (`alt[0] != ref[ins_point]`) and shuffle never
/// starts. `insertion_to_repeat` already iterates rotations to handle the
/// 2+-copy case; this helper does the same for the 1-copy (duplication) case,
/// restoring symmetry between the single- and multi-copy paths.
///
/// # Algorithm
///
/// 1. Compute `unit = smallest_repeat_unit(inserted_seq)`. Reject when
///    `inserted_seq.len() / unit.len() != 1` (multi-copy is repeat
///    notation's job, handled by `insertion_to_repeat`).
/// 2. For each rotation `r in 0..unit.len()`, build the rotated unit and
///    call `find_tandem_extent` to locate the maximal tandem run abutting
///    the insertion point.
/// 3. Resolve each rotation all the way to the duplication it would emit
///    (step 5), and hold it to the `tandem_edit_preserves_insertion` phase
///    gate (#882) — a rotation that would decode to a different sequence than
///    the input is not a candidate at all.
/// 4. Among the rotations that clear that gate, pick the one yielding the
///    largest tandem run (ties broken by direction, below). The gate is
///    evaluated only for a rotation that would otherwise take the lead, which
///    selects the same winner while skipping the gate's `ref_seq`-sized
///    reconstructions for rotations that have already lost.
/// 5. Return the unit-aligned duplication slot inside that tract whose
///    position matches the requested shuffle direction:
///    - `ThreePrime`: `dup_start_idx = tract.end - unit.len()` — the
///      most-3' unit (current default; pre-#340 behavior).
///    - `FivePrime`: `dup_start_idx = ref_start` — the most-5' unit, so
///      a homopolymer insertion under 5'-direction shuffle reports the
///      leftmost-canonical dup rather than the 3'-most one (closes #340
///      subgroup (i)).
///
///    Both end at `dup_start_idx + unit.len() - 1`, converted to 1-based
///    HGVS positions.
///
/// Returns `None` when no rotation matches an adjacent tract (true
/// non-tandem insertion), **and also** when every matching rotation's emitted
/// duplication is out of phase — i.e. it would decode to a different sequence
/// than the input insertion, so it is not a tandem duplication per
/// `duplication.md` (the `tandem_edit_preserves_insertion` gate, #882). In
/// both cases the caller leaves it as a plain `ins`.
///
/// The phase gate is applied per rotation *during* selection, not once to the
/// winner (#1349). Scoring first and gating afterwards let a high-scoring but
/// out-of-phase rotation shadow a lower-scoring in-phase one and take the whole
/// result down with it: for `insTCC` at `GTACGT|CCCCTCCTCCTCCTCCC`, the alt's
/// own phase abuts a one-copy `TCC` tract (the correct dup) while the `CCT`
/// rotation abuts a two-copy tract further 3' that is out of phase with the
/// insertion. `CCT` won on count, failed the gate, and the correct candidate was
/// never reconsidered — costing the 5' path its `ins → dup` promotion and
/// leaving a non-idempotent boundary `delins` behind.
pub(crate) fn insertion_to_duplication(
    ref_seq: &[u8],
    pos: u64,
    inserted_seq: &[u8],
    direction: crate::normalize::config::ShuffleDirection,
) -> Option<InsToDupResult> {
    use crate::normalize::config::ShuffleDirection;

    if inserted_seq.is_empty() {
        return None;
    }

    let base_unit = smallest_repeat_unit(inserted_seq);
    let added_copies = inserted_seq.len() / base_unit.len();
    if added_copies != 1 {
        // 2+-copy insertions are handled by `insertion_to_repeat`.
        return None;
    }

    // #1846: bail before the rotation scan when the unit cannot fit the window.
    // Every rotation below has length `base_unit.len()`, and `find_tandem_extent`
    // skips every anchor once `anchor + unit.len() > ref_seq.len()` — so a unit
    // longer than `ref_seq` matches no anchor in any rotation, and `best` is
    // already determined to stay `None`. Without this exit the loop still runs
    // `base_unit.len()` iterations, each allocating and filling a
    // `base_unit.len()`-byte `Vec`, i.e. Θ(L²) bytes copied to reach that fixed
    // `None`. A spec-legal `ins[ACC:g.A_B]` expands to a payload of the named
    // range's length (a 38-char cross-reference named 10.8 MB), so an unbounded
    // scan here never terminates on such input; the early exit makes it O(L).
    // `ref_seq` is the fetched shuffle window (default 100 bases, grown by
    // `normalize_in_grown_window`), so this fires only far above any window the
    // normalizer actually uses and cannot change which inputs convert to a dup.
    if base_unit.len() > ref_seq.len() {
        return None;
    }

    /// One rotation resolved all the way to the duplication it would emit, and
    /// already past the phase gate — i.e. a candidate the caller could receive.
    struct DupCandidate {
        /// Emitted unit, post-alignment (the aligner may re-rotate it).
        unit: Vec<u8>,
        /// 0-based start of the emitted duplication slot.
        dup_start_idx: usize,
        /// 0-based end of the emitted duplication slot (inclusive).
        dup_end_idx: usize,
        /// Start of the matched reference tract, used only for the tie-break.
        ref_start: usize,
        /// Copies of `unit` in the matched tract — the primary score.
        ref_count: u64,
    }

    let u_len = base_unit.len();
    let mut best: Option<DupCandidate> = None;
    for r in 0..u_len {
        let mut rotated = Vec::with_capacity(u_len);
        rotated.extend_from_slice(&base_unit[r..]);
        rotated.extend_from_slice(&base_unit[..r]);
        let Some((ref_start, ref_count)) = find_tandem_extent(ref_seq, pos as usize, &rotated)
        else {
            continue;
        };
        if ref_count == 0 {
            continue;
        }

        // Resolve this rotation to the duplication it would emit, so the phase
        // gate below judges the same values the caller would receive.
        let (unit, dup_start_idx) =
            align_dup_slot(ref_seq, ref_start, ref_count, rotated, direction);
        let dup_end_idx = dup_start_idx + unit.len() - 1;

        // Selection rules:
        // - Higher `ref_count` wins outright (preserves the multi-copy
        //   tract phase per #132 / #180).
        // - On ties, the rotation whose anchor is direction-canonical
        //   wins: most-5' (smallest `ref_start`) for `FivePrime`,
        //   most-3' (largest `ref_start`, which orders the same as
        //   `tract_end` since all rotations share `u_len`) for
        //   `ThreePrime`. Without this tie-break, the `ref_count == 1`
        //   case — an isolated single tandem unit adjacent to the
        //   insertion point — falls to first-found iteration order and
        //   disagrees with biocommons on the anchor position. Closes
        //   #403.
        let is_better = match best.as_ref() {
            None => true,
            Some(current) => match ref_count.cmp(&current.ref_count) {
                std::cmp::Ordering::Greater => true,
                std::cmp::Ordering::Less => false,
                std::cmp::Ordering::Equal => match direction {
                    ShuffleDirection::FivePrime => ref_start < current.ref_start,
                    ShuffleDirection::ThreePrime => ref_start > current.ref_start,
                },
            },
        };
        if !is_better {
            // Cannot become `best` whatever the gate says, and the gate rebuilds
            // two `ref_seq`-sized buffers per call — so don't pay for it.
            continue;
        }

        // #882 phase gate, applied per candidate (#1349): only emit a dup when it
        // decodes to the same sequence as the input insertion. A rotation that
        // fails here is not a duplication of this insertion at all, so it must
        // not compete — and must not shadow a lower-scoring rotation that is.
        // Gating a rotation that already lost the comparison above would change
        // nothing, since a failing rotation never becomes `best` either way.
        if !tandem_edit_preserves_insertion(
            ref_seq,
            pos as usize,
            inserted_seq,
            dup_start_idx,
            dup_end_idx + 1,
            &unit,
            2,
        ) {
            continue;
        }

        best = Some(DupCandidate {
            unit,
            dup_start_idx,
            dup_end_idx,
            ref_start,
            ref_count,
        });
    }

    let best = best?;

    Some(InsToDupResult {
        unit: best.unit,
        start: index_to_hgvs_pos(best.dup_start_idx),
        end: index_to_hgvs_pos(best.dup_end_idx),
        ref_count: best.ref_count,
    })
}

/// Resolve a matched tandem tract to the duplication slot to emit for `direction`,
/// returning the (possibly re-rotated) unit and the 0-based start of the slot.
///
/// Split out of `insertion_to_duplication` so each rotation can be resolved —
/// and phase-checked — independently during selection (#1349).
fn align_dup_slot(
    ref_seq: &[u8],
    ref_start: usize,
    ref_count: u64,
    rotated: Vec<u8>,
    direction: crate::normalize::config::ShuffleDirection,
) -> (Vec<u8>, usize) {
    use crate::normalize::config::ShuffleDirection;

    let mut unit = rotated;
    let tract_end_idx = ref_start + (ref_count as usize) * unit.len();
    // `find_tandem_extent` returns `ref_start` aligned on a `unit`-length
    // anchor probe, so picking `ref_start` for the 5' end never lands the
    // dup slot mid-unit even when the inserted alt was a non-zero
    // rotation of the canonical reference unit.
    let dup_start_idx = match direction {
        ShuffleDirection::ThreePrime => {
            // Push to the *true* 3'-most slot, crossing any off-phase tract
            // extension. `tract_end_idx - unit.len()` only 3'-aligns within the
            // abutting rotation's phase; when the run continues one base further
            // in the opposite phase (e.g. inserting `TG` at the 5' edge of a
            // `TGTGTGTGT` run — the abutting `TG`/`GT` tract ends one base short
            // of the run's true 3' end), stopping there violates the 3'rule.
            // `three_prime_align_tract` walks the remaining flank while
            // `ref[end] == rotated[0]` — exactly the duplication slide rule
            // `ref[s] == ref[s+L]` — rotating the unit to the run's 3' phase,
            // mirroring the repeat path (#852). HGVS DNA duplication.md: "the
            // most 3' position possible … is arbitrarily assigned to have been
            // changed (3'rule)." Verified against mutalyzer (#864).
            let (_, aligned_end, aligned_unit) =
                three_prime_align_tract(ref_seq, ref_start, tract_end_idx, &unit);
            unit = aligned_unit;
            aligned_end - unit.len()
        }
        ShuffleDirection::FivePrime => {
            // Mirror of the 3' branch: push to the *true* 5'-most slot, crossing
            // any off-phase tract extension. `ref_start` only 5'-aligns within the
            // abutting rotation's phase; when the run continues one base further 5'
            // in the opposite phase (e.g. inserting `TG` at the 3' edge of a
            // `GTGTGTG` run — the abutting `TG` tract starts one base short of the
            // run's true 5' end), stopping there under-shifts by one and is not a
            // fixed point (the emitted `dup` re-shuffles one base further 5').
            // `five_prime_align_tract` walks the remaining 5' flank while
            // `ref[start-1] == rotated.last()` — the mirror duplication slide rule
            // — rotating the unit to the run's 5' phase (#8's dup-path sibling).
            let (aligned_start, _, aligned_unit) =
                five_prime_align_tract(ref_seq, ref_start, tract_end_idx, &unit);
            unit = aligned_unit;
            aligned_start
        }
    };
    (unit, dup_start_idx)
}

/// Smallest unit `U` such that `seq = U * k` for some integer k ≥ 1.
pub(crate) fn smallest_repeat_unit(seq: &[u8]) -> &[u8] {
    let n = seq.len();
    for u in 1..=n {
        if !n.is_multiple_of(u) {
            continue;
        }
        let unit = &seq[..u];
        if seq.chunks_exact(u).all(|c| c == unit) {
            return unit;
        }
    }
    seq // fallback (unreachable for n>0)
}

/// Find the maximal tandem run of `unit` near position `pos` in
/// `ref_seq`. `pos` is the 0-based index of the base immediately 5' of
/// the insertion point — i.e., the insertion is between `pos` and
/// `pos+1`. The tract may not be unit-aligned with `pos`, so we probe
/// anchor offsets in a window around `pos` and pick the run that
/// abuts/contains the insertion point.
///
/// Returns `(ref_start_0based, count)` of the chosen tandem run.
fn find_tandem_extent(ref_seq: &[u8], pos: usize, unit: &[u8]) -> Option<(usize, u64)> {
    let u = unit.len();
    if u == 0 {
        return None;
    }

    // Insertion point in "between-bases" coords: between pos and pos+1.
    let ins_point = pos + 1;

    // Probe candidate anchor offsets within [ins_point - u, ins_point].
    // For each anchor `a`, compute the maximal tandem run that includes
    // anchor `a` (walk left and right in steps of `u`). Accept runs that
    // abut or span the insertion point. Among accepted runs, pick the
    // largest count.
    let lo = ins_point.saturating_sub(u);
    let hi = ins_point.min(ref_seq.len());

    let mut best: Option<(usize, u64)> = None;

    for anchor in lo..=hi {
        if anchor + u > ref_seq.len() {
            continue;
        }
        if &ref_seq[anchor..anchor + u] != unit {
            continue;
        }

        // Anchor identifies one unit-aligned occurrence; extend to the full tract.
        // The `?`-via-let-else is defensive: we just verified the anchor matches,
        // so `extend_tandem_tract` shouldn't return None here.
        let Some(TandemTract {
            start,
            end,
            ref_count: count,
        }) = extend_tandem_tract(ref_seq, anchor..anchor + u, unit)
        else {
            continue;
        };

        // Require the run [start, end) to abut/span ins_point: the
        // insertion point must lie within [start, end].
        if ins_point < start || ins_point > end {
            continue;
        }

        match best {
            None => best = Some((start, count)),
            Some((_, bc)) if count > bc => best = Some((start, count)),
            _ => {}
        }
    }

    best
}

/// Information about a maximal tandem repeat tract in a reference sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct TandemTract {
    /// 0-based, inclusive start of the tract.
    pub start: usize,
    /// 0-based, exclusive end of the tract.
    pub end: usize,
    /// Number of full `unit`-length copies in the tract.
    pub ref_count: u64,
}

/// Walk left from `anchor.start` and right from `anchor.end` extending
/// unit-by-unit while `ref_seq` matches `unit`. The anchor span must lie on
/// `unit`-length boundaries — i.e., `anchor.len() % unit.len() == 0` — and
/// `ref_seq[anchor]` must be `unit` repeated some non-negative number of times.
///
/// Returns the maximal tract enclosing the anchor, or `None` if the anchor
/// itself isn't unit-periodic (defensive — `deletion_to_repeat` already checks
/// before calling, but we don't trust future callers).
pub(crate) fn extend_tandem_tract(
    ref_seq: &[u8],
    anchor: std::ops::Range<usize>,
    unit: &[u8],
) -> Option<TandemTract> {
    let u = unit.len();
    if u == 0 || anchor.start > anchor.end || anchor.end > ref_seq.len() {
        return None;
    }
    if !(anchor.end - anchor.start).is_multiple_of(u) {
        return None;
    }
    // Defensive: verify the anchor itself is unit-periodic.
    if !ref_seq[anchor.start..anchor.end]
        .chunks_exact(u)
        .all(|chunk| chunk == unit)
    {
        return None;
    }

    // Walk left in steps of u.
    let mut start = anchor.start;
    while start >= u && &ref_seq[start - u..start] == unit {
        start -= u;
    }
    // Walk right in steps of u.
    let mut end = anchor.end;
    while end + u <= ref_seq.len() && &ref_seq[end..end + u] == unit {
        end += u;
    }

    let ref_count = ((end - start) / u) as u64;
    Some(TandemTract {
        start,
        end,
        ref_count,
    })
}

/// Description of how a deletion can be re-expressed as repeat notation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct DelToRepeatResult {
    /// The phase-aligned tandem unit (length 1 for homopolymer, longer for
    /// multi-base tandem).
    pub unit: Vec<u8>,
    /// Post-deletion unit count (`N - k`).
    pub count: u64,
    /// 1-based HGVS start of the canonical reference tract (inclusive).
    pub start: u64,
    /// 1-based HGVS end of the canonical reference tract (inclusive).
    pub end: u64,
}

/// Compute the canonical `unit[N-k]` description for a post-3'-shift deletion,
/// or return `None` if the deletion does not meet B2 trigger conditions.
///
/// # Shuffle phase-alignment lemma (why no rotation iteration)
///
/// The 3' shuffle (`shuffle.rs`) extends `del_end` while
/// `ref[del_start] == ref[del_end]`, terminating exactly when the period
/// breaks. For a deletion sitting inside a tandem tract with period `p`, the
/// post-shift `del_slice = ref_seq[del_start..del_end]` is therefore
/// unit-aligned with the surrounding tract, and `smallest_repeat_unit(del_slice)`
/// equals the tract's natural unit. So `extend_tandem_tract` walking from the
/// post-shift `del_start` finds the full tract on the r=0 phase — no rotation
/// iteration is required here.
///
/// # Trigger conditions
///
/// All five must hold; otherwise returns `None` (caller emits plain `del`):
///
/// 1. `del_len % unit_len == 0` (length is a multiple of the unit's period).
/// 2. The reference tract enclosing the deletion has `ref_count >= 2`.
/// 3. `k = del_len / unit_len >= 2` (single-unit removal stays as `del`).
/// 4. `post_count = ref_count - k >= 1` (full tract removal stays as `del`).
/// 5. The deletion lies entirely within the tract (handled by `extend_tandem_tract`'s
///    defensive periodicity check on the anchor).
///
/// The bounds check rejects 0-length deletions (`del_start == del_end`).
/// Callers passing insertion-point-shaped zero-width ranges should use
/// `extend_tandem_tract` directly instead.
///
/// # The codon gate is answered for the tract, not the input (#1270)
///
/// A sixth condition applies on the `c.`/`r.` axes and is deliberately not in
/// the list above, because it is about the *emitted* description rather than the
/// trigger: `repeated.md` requires `unit_len % 3 == 0` for repeat notation in a
/// coding context.
///
/// `is_coding` is the caller's verdict for the span the *input* deletion
/// occupied. Where `cds_span` is supplied it takes precedence, and the gate is
/// re-answered for the tandem tract this function discovers — the span the
/// emitted repeat actually occupies. The two differ exactly when the tract runs
/// past `cds_end`, whose 3'UTR part the carve-out exempts. Pass `None`
/// (genomic, `n.`, or any caller with no CDS context) to keep `is_coding`
/// verbatim.
///
/// This mirrors [`insertion_to_repeat`], which took the same argument for #1185
/// and #1210; before #1270 only that path re-asked, so a `T` run spanning
/// `cds_end` gave `c.5_6insTT` -> `c.4_*5T[16]` while the matching `c.5_6del`
/// gave a plain `c.*4_*5del`. The gate is therefore evaluated *after*
/// `extend_tandem_tract` rather than before it, since the tract is what it now
/// asks about.
pub(crate) fn deletion_to_repeat(
    ref_seq: &[u8],
    del_start: usize,
    del_end: usize,
    is_coding: bool,
    cds_span: Option<(u64, u64)>,
) -> Option<DelToRepeatResult> {
    if del_start >= del_end || del_end > ref_seq.len() {
        return None;
    }
    let del_slice = &ref_seq[del_start..del_end];
    let unit_slice = smallest_repeat_unit(del_slice);
    let p = unit_slice.len();
    if p == 0 || !(del_end - del_start).is_multiple_of(p) {
        return None;
    }

    let tract = extend_tandem_tract(ref_seq, del_start..del_end, unit_slice)?;
    if tract.ref_count < 2 {
        return None;
    }

    // HGVS spec (repeated.md): in c. context, repeat notation requires
    // unit_len % 3 == 0. Falls back to plain del when the gate blocks.
    //
    // Answered for the **tract**, not for the input span — the #1185/#1210
    // argument that `insertion_to_repeat` above already applies, and the gate
    // is therefore checked after `extend_tandem_tract` rather than before it.
    // `is_coding` is the caller's verdict for the span the *input* deletion
    // occupied, and a CDS-resident deletion inside a tract that runs past
    // `cds_end` sits on a tract the carve-out exempts. Measured on a `T` run
    // spanning `cds_end`: `c.5_6insTT` produced `c.4_*5T[16]` while the
    // matching `c.5_6del` produced a plain `c.*4_*5del`, because only the
    // insertion path re-asked (#1270).
    //
    // Frames line up as they do there: `tract.start` / `tract.end` are 0-based
    // indices into the transcript-frame `ref_seq` with `end` exclusive, and
    // `cds_span` is 1-based inclusive in that same frame. `None` keeps the
    // caller's verdict verbatim, so genomic / `n.` callers and the unit tests
    // below are unaffected.
    let gate_is_coding = match cds_span {
        Some((cds_start, cds_end)) => {
            let tract_start = tract.start as u64 + 1;
            let tract_end = tract.end as u64;
            tract_start >= cds_start && tract_end <= cds_end
        }
        None => is_coding,
    };
    if gate_is_coding && !p.is_multiple_of(3) {
        return None;
    }

    let k = ((del_end - del_start) / p) as u64;
    if k < 2 {
        return None;
    }

    let post_count = tract.ref_count - k;
    if post_count == 0 {
        return None;
    }

    Some(DelToRepeatResult {
        unit: unit_slice.to_vec(),
        count: post_count,
        start: index_to_hgvs_pos(tract.start),
        end: index_to_hgvs_pos(tract.end - 1),
    })
}

/// True when `last` is the complement of `first`, i.e. inverting the pair
/// changes nothing.
///
/// Both sides fold case, and an unmodelled byte on either side answers `false`
/// — both properties come from [`complement_base`], which is now the crate's
/// single byte-level complement (#1318). This used to work around a local
/// helper whose fallback returned unmodelled bytes unchanged; with that helper
/// gone the wrapper is a thin, honest comparison.
fn is_complementary_pair(first: u8, last: u8) -> bool {
    // `complement_base` never answers with an unmodelled byte, so if `last` is
    // one the comparison is false — which is the wanted behaviour, and the
    // reason this does not need a separate "is `last` a base?" test.
    complement_base(first) == Some(last.to_ascii_uppercase())
}

/// True when inverting `base` on its own changes nothing, i.e. it is one of the
/// self-complementary symbols `W` (A/T), `S` (C/G) and `N`.
///
/// A byte outside the modelled alphabet answers **`false`**, not `true`: there
/// is no complement to report, so "inverting it changes nothing" is not a claim
/// this can make. The previous helper returned such bytes unchanged, which made
/// every unmodelled symbol look self-complementary — the #1249 trap.
fn is_self_complementary(base: u8) -> bool {
    is_complementary_pair(base, base)
}

/// The `reference > alternative` pair that describes inverting a single base:
/// the base itself and its complement.
///
/// Returns `None` when the base is its own complement (nothing is substituted)
/// or when either side is not a typed [`crate::hgvs::edit::Base`], so callers
/// never fabricate a substitution they cannot spell.
pub fn complementary_substitution(
    base: u8,
) -> Option<(crate::hgvs::edit::Base, crate::hgvs::edit::Base)> {
    use crate::hgvs::edit::Base;

    if is_self_complementary(base) {
        return None;
    }
    Some((
        Base::from_char(base as char)?,
        Base::from_char(complement_base(base)? as char)?,
    ))
}

/// Shorten an inversion by removing outer bases that cancel out
///
/// When the first base of an inversion is complementary to the last base,
/// inverting them produces no net change, so they can be excluded.
/// This is repeated until no more cancellation is possible.
///
/// Returns the shortened (start, end) positions (0-indexed), or `None` if the
/// inversion reduces to identity.
///
/// A residue of **two or more** bases is still an inversion. A residue of
/// exactly **one** base is not identity: that base is replaced by its own
/// complement, which `DNA/inversion.md:16` says to describe as a substitution
/// ("a one-nucleotide inversion should be described as a substitution"). Such a
/// residue is returned as a one-base range so the caller can emit that
/// substitution; only a self-complementary base (`W`, `S`, `N`) collapses to
/// identity. Conflating the one- and zero-base cases silently dropped real
/// variants — see #1249.
///
/// Note for external callers: a span that once yielded `None` may now yield a
/// one-base range.
pub fn shorten_inversion(ref_seq: &[u8], start: usize, end: usize) -> Option<(usize, usize)> {
    if start >= end || end > ref_seq.len() {
        return None;
    }

    let mut s = start;
    let mut e = end;

    // Keep shrinking while the outer *pair* is complementary. The bound is
    // `s + 1 < e`, not `s < e`: at `e == s + 1` the span is a single base and
    // `ref_seq[s]` and `ref_seq[e - 1]` are the same byte, so a
    // self-complementary centre would cancel against itself and drive `e` below
    // `s`. The lone-centre case is classified after the loop instead.
    while s + 1 < e {
        let first = ref_seq[s];
        let last = ref_seq[e - 1]; // e is exclusive

        // Check if first base is complement of last base
        if is_complementary_pair(first, last) {
            s += 1;
            e -= 1;
        } else {
            break;
        }
    }

    // Every base cancelled against its partner, so the inversion is a genuine
    // no-op. This is the only span that shrinks to nothing, and it requires an
    // even length.
    if e == s {
        return None;
    }

    // An odd-length span always leaves a centre base. Inverting one base
    // replaces it with its complement, so it is identity only when the base is
    // its own complement.
    if e == s + 1 && is_self_complementary(ref_seq[s]) {
        return None;
    }

    Some((s, e))
}

/// Canonical form for a delins per HGVS edit-type priority (sub > del > inv > dup > ins).
///
/// All position fields are 0-indexed. Sub/del/ins/inv/keep-as-delins variants
/// carry positions taken from the *trimmed* delins (after greedy shared-affix
/// elimination on both ends); `Identity` and `Duplication` reflect the full
/// input range because shared-affix trimming would destroy those classifications.
/// The caller converts each 0-indexed position to HGVS 1-based for emission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DelinsCanonical {
    /// Inserted sequence equals the deleted reference. Emit as identity (`=`)
    /// at the input range.
    Identity,
    /// 1-base substitution at 0-indexed `position` (in `ref_seq`). Emit `ref>alt`.
    Substitution {
        position: usize,
        reference: crate::hgvs::edit::Base,
        alternative: crate::hgvs::edit::Base,
    },
    /// Pure deletion (insert consumed entirely by shared-affix trimming) over
    /// half-open 0-indexed `[start, end)`.
    Deletion { start: usize, end: usize },
    /// Pure insertion (deleted reference consumed entirely by shared-affix
    /// trimming). `after_index` is the 0-indexed position of the base
    /// immediately AFTER the insertion point — equivalently, the 1-based HGVS
    /// position of the base immediately BEFORE. HGVS form:
    /// `c.{after_index}_{after_index + 1}ins{sequence}`.
    Insertion {
        after_index: usize,
        sequence: Vec<u8>,
    },
    /// N>=2 delins whose insert is the reverse complement of the deleted
    /// reference (possibly after shared-affix trimming), with the
    /// complementary-outer-bases shortening rule applied. Half-open 0-indexed
    /// `[start, end)`; by construction `end > start + 1`.
    Inversion { start: usize, end: usize },
    /// N -> 2N delins where insert is the (full-range) deleted sequence
    /// repeated twice. Range fields are the full input range, not a trimmed
    /// sub-range — see the doc on `canonicalize_delins` for why duplication
    /// must be detected before trimming.
    Duplication { start: usize, end: usize },
    /// None of the higher-priority forms apply. The (possibly trimmed) delins
    /// occupies half-open 0-indexed `[start, end)` with `sequence` as the
    /// inserted bases. If no trimming was possible these match the input.
    KeepAsDelins {
        start: usize,
        end: usize,
        sequence: Vec<u8>,
    },
}

/// Per-position sub-edit kind emitted by `decompose_delins` when a delins
/// span decomposes into the spec-canonical edit-priority form (`inv`
/// sub-span recognition from issue #160, plus the sub-only branch from
/// issue #165 / tracking issue #81 item A10).
///
/// All position fields are 0-indexed offsets into `ref_seq` (the same input
/// passed to `decompose_delins`), matching the convention used by
/// `DelinsCanonical` so callers can convert with the same `index_to_hgvs_pos`
/// helper.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DelinsSubedit {
    /// 1-base substitution at `position`. Renders as `g.<position+1>R>A`.
    Substitution {
        position: usize,
        reference: crate::hgvs::edit::Base,
        alternative: crate::hgvs::edit::Base,
    },
    /// Inversion over half-open `[start, end)`; by construction `end > start + 1`.
    Inversion { start: usize, end: usize },
    /// 1-base identity at `position` — interior position unchanged from
    /// reference (e.g. codon-frame merge synthesized middle base, issue
    /// #79). The `base` field carries the unchanged reference byte so
    /// callers can reconstruct the codon-frame triplet's alt sequence
    /// when preserving the spec's `general.md:35-38` exception during
    /// decomposition. Always an IUPAC `Base` variant — `decompose_delins`
    /// abandons the whole decomposition if it sees a non-IUPAC byte at
    /// an identity position. Renders as `g.<position+1>=` when emitted
    /// directly.
    IdentityAt {
        position: usize,
        base: crate::hgvs::edit::Base,
    },
}

/// Decompose a delins into a sequence of canonical sub-edits when a maximal
/// contiguous mismatch run is *entirely* an inversion, or the span contains
/// at least one position whose alt byte equals the ref byte (an interior
/// identity). Returns `Some(edits)` only when the resulting decomposition has
/// more than one element AND at least one element is an `Inversion` or an
/// `IdentityAt`.
///
/// Implements:
/// - Issue #160 (item A2 + A10 inv-branch of tracking issue #81), as
///   corrected by issue #1034: an `Inversion` is emitted only when a whole
///   maximal contiguous mismatch run is a reverse complement, per the HGVS
///   edit-priority rule applied to the
///   *maximal contiguous run* (`DNA/inversion.md`). A reverse-complement
///   sub-run of a longer contiguous change is NOT carved out — that change
///   stays a single `delins`.
/// - Issue #165 (item A10 sub-only branch of #81): an interior position
///   matching the reference (an `IdentityAt`) is the spec's signal that
///   the surrounding mismatches are independent variants under
///   `general.md:34` (variants separated by ≥1 unchanged nt are
///   described individually) and per the `substitution > delins`
///   priority. The caller decides whether to emit such triplets as
///   separate variants or to preserve the spec's codon-frame exception
///   (`general.md:35-38`) by re-grouping eligible `[Sub; Identity; Sub]`
///   patterns; this function reports the decomposition either way.
///
/// Returns `None` for:
/// - Inputs whose post-scan emission contains only adjacent
///   `Substitution` elements (no `Inversion`, no `IdentityAt`). Such a
///   decomposition would just re-merge into the input delins under the
///   adjacency rule (`substitution.md`: consecutive nucleotide changes
///   are described as `delins`, issue #182), so returning `None` lets
///   the caller short-circuit and avoid pointless churn.
/// - Inputs where the entire span is a single inversion (already handled
///   by the full-span check in `canonicalize_delins`).
/// - Unequal-length delins (`alt.len() != ref.len()`).
/// - Length-1 inputs (substitution range, never a delins shape).
/// - Inputs containing a non-IUPAC byte at any position (substitution
///   or identity); abandoning the decomposition keeps the original
///   delins form so a subsequent round-trip lands on the same string.
/// - Out-of-bounds or empty ranges (`start >= end` or
///   `end > ref_seq.len()`).
///
/// Algorithm (left-to-right scan over maximal contiguous mismatch runs):
/// - At an identity position (`alt[i] == ref[i]`): emit `IdentityAt(i,
///   ref[i])`, advance `i += 1`. The byte is recorded so a downstream
///   codon-frame preservation pass can reconstruct a 3-base delins's alt
///   sequence. An identity also bounds the maximal contiguous run — variants
///   separated by ≥1 unchanged nt are independent (`general.md:34`).
/// - Otherwise collect the maximal contiguous run `[rs, re)` of positions
///   where `alt != ref`, then:
///   - If the **whole** run is an inversion (`re - rs >= 2` and
///     `alt[rs..re] == revcomp(ref[rs..re])`): emit `Inversion(rs, re)`.
///   - Else: emit a `Substitution` for each position in the run.
///
/// A reverse-complement **sub-run** is never carved out of a longer
/// contiguous change: only a run that is *entirely* a reverse complement is
/// typed as `inv`, per `DNA/inversion.md` and the issue #1034 fix. The
/// surrounding `Substitution`s re-merge into a single `delins` under the
/// caller's adjacency rule (`substitution.md` / issue #182).
///
/// `shorten_inversion`'s complementary-outer-pair peeling is not applied
/// here: within a pure mismatch run a peelable outer pair would imply
/// `alt == ref` at the boundary (an identity, not a mismatch), so the run
/// bounds are already minimal for a full-run inversion.
pub fn decompose_delins(
    ref_seq: &[u8],
    start: usize,
    end: usize,
    inserted_seq: &[u8],
) -> Option<Vec<DelinsSubedit>> {
    use crate::hgvs::edit::Base;

    // Bounds + length precondition. Sub-span decomposition is only meaningful
    // when alt and ref have the same length and the span has at least 2 bases.
    if start >= end || end > ref_seq.len() {
        return None;
    }
    let n = end - start;
    if inserted_seq.len() != n || n < 2 {
        return None;
    }
    let deleted = &ref_seq[start..end];

    let mut emitted: Vec<DelinsSubedit> = Vec::new();
    let mut has_inv = false;
    let mut has_identity = false;
    let mut i = 0;
    while i < n {
        if deleted[i].eq_ignore_ascii_case(&inserted_seq[i]) {
            // Identity position: record the unchanged ref byte so the caller
            // can rebuild a codon-frame triplet's alt sequence (issue #165),
            // and let it bound the maximal contiguous run. Abandon the whole
            // decomposition on a non-IUPAC byte, mirroring the substitution
            // branch below: an identity at a non-IUPAC byte cannot be
            // re-rendered as a 3-base delins alt without silently coercing the
            // unknown byte to `N`, which would diverge from the next
            // round-trip's input.
            let b = Base::from_char(deleted[i] as char)?;
            emitted.push(DelinsSubedit::IdentityAt {
                position: start + i,
                base: b,
            });
            has_identity = true;
            i += 1;
            continue;
        }

        // Maximal contiguous run of mismatches `[run_start, run_end)`.
        //
        // Case-insensitive for the same reason as the identity test above, and
        // it must be: this is the same rule written twice, and folding case in
        // only one of them is worse than folding it in neither. On a
        // soft-masked reference the run would swallow positions that the test
        // above calls identities, so `has_identity` never became true and the
        // final gate returned `None` for a span whose uppercase twin decomposes
        // — defeating the fix above for exactly the input it targets (#1318).
        let run_end_differs = |k: usize| !deleted[k].eq_ignore_ascii_case(&inserted_seq[k]);
        let run_start = i;
        let mut run_end = i + 1;
        while run_end < n && run_end_differs(run_end) {
            run_end += 1;
        }

        // Only type the run as `inv` when the ENTIRE run is a reverse
        // complement (issue #1034). A reverse-complement sub-run is never
        // carved out of a longer contiguous change — the surrounding
        // substitutions re-merge into a single `delins` under the caller's
        // adjacency rule.
        if run_end - run_start >= 2
            && is_revcomp(
                &deleted[run_start..run_end],
                &inserted_seq[run_start..run_end],
            )
        {
            emitted.push(DelinsSubedit::Inversion {
                start: start + run_start,
                end: start + run_end,
            });
            has_inv = true;
        } else {
            // Emit a per-position substitution for each mismatch. Bases must
            // be IUPAC; a non-IUPAC byte cannot be expressed as `Base`, so
            // abandon the whole decomposition (the caller keeps the delins
            // as-is). Mirrors the `canonicalize_delins` substitution branch.
            for k in run_start..run_end {
                let r = Base::from_char(deleted[k] as char)?;
                let a = Base::from_char(inserted_seq[k] as char)?;
                emitted.push(DelinsSubedit::Substitution {
                    position: start + k,
                    reference: r,
                    alternative: a,
                });
            }
        }
        i = run_end;
    }

    // Trigger: commit only if the decomposition has > 1 element AND it
    // contains at least one `Inversion` or `IdentityAt`. The former is the
    // issue #160 path (`inv > delins` priority). The latter is the issue
    // #165 / item A10 path (`sub > delins` priority, with `IdentityAt`
    // marking the gap between non-adjacent subs that must split per
    // `general.md:34`).
    //
    // A pure-`Substitution` emission of length > 1 has no interior gap and
    // would be re-merged back into the input delins by the caller's
    // adjacency rule (`substitution.md` / issue #182), so returning `None`
    // there lets the caller short-circuit the round-trip.
    //
    // A single-`Inversion` result is the same shape that the existing
    // full-span check in `canonicalize_delins` already produces, so there
    // is no point splitting for it either.
    if emitted.len() >= 2 && (has_inv || has_identity) {
        Some(emitted)
    } else {
        None
    }
}

/// Classify a delins into its HGVS canonical form.
///
/// Implements item A2 of issue #81 alongside the previously separate
/// identity / substitution / duplication checks, plus shared-affix trimming
/// so a delins that is structurally a smaller edit (sub / del / ins / inv)
/// surrounded by identical context canonicalizes to that smaller form per
/// the HGVS minimal-form rule. Priority follows the HGVS edit-priority
/// recommendation (substitution > deletion > inversion > duplication >
/// insertion); `Identity` short-circuits ahead of all of them because an
/// identity is never a real edit.
///
/// Algorithm:
/// 1. Bounds / degenerate input -> `KeepAsDelins` (untrimmed).
/// 2. Full-range `Identity` (`insert == deleted`).
/// 3. Full-range `Duplication` (`insert == deleted + deleted`). Must precede
///    trimming, since greedy shared-affix trimming on a duplication would
///    consume the entire deleted range and falsely reclassify it as a pure
///    insertion of the duplicated unit.
/// 4. Greedy shared-affix trimming on both ends.
/// 5. Reclassify the trimmed range as `Insertion` (deleted empty),
///    `Deletion` (insert empty), `Substitution` (1->1), `Inversion`
///    (N>=2 revcomp + outer-pair shortening), or `KeepAsDelins`.
///
/// Spec references (`assets/hgvs-nomenclature/docs/recommendations/DNA/`):
/// - `inversion.md`: an inversion requires "more than one nucleotide" and is
///   defined as the inserted sequence being the reverse complement of the
///   deleted reference; a one-nucleotide complement is a substitution.
/// - `delins.md`: a delins is "not a substitution, inversion or conversion".
pub fn canonicalize_delins(
    ref_seq: &[u8],
    start: usize,
    end: usize,
    inserted_seq: &[u8],
) -> DelinsCanonical {
    use crate::hgvs::edit::Base;

    if start >= end || end > ref_seq.len() || inserted_seq.is_empty() {
        return DelinsCanonical::KeepAsDelins {
            start,
            end,
            sequence: inserted_seq.to_vec(),
        };
    }

    let deleted = &ref_seq[start..end];

    // 1. Identity (insert == deleted reference). Caught at full range; trimming
    //    would otherwise consume the entire range and lose this classification.
    //
    //    Case-insensitive (#1318): on a soft-masked reference `deleted` arrives
    //    lower-case while `inserted_seq` comes from the author's upper-case
    //    description, so a raw byte test missed an identity that is one. That
    //    is not merely a mis-classification — the range then fell through to
    //    the trimming below, which now consumes it entirely, and the pure-
    //    insertion branch's `debug_assert!("Identity case caught above")` is
    //    exactly the claim this step has to make true.
    if deleted.eq_ignore_ascii_case(inserted_seq) {
        return DelinsCanonical::Identity;
    }

    // 2. Duplication (N -> 2N, full-range insert == deleted+deleted). Must be
    //    checked before trimming: a true duplication has identical prefixes
    //    (insert[..N] == deleted == ref[start..end]) so greedy trimming would
    //    eat the entire deleted range and downgrade dup to ins, violating the
    //    sub > del > inv > dup > ins priority.
    if inserted_seq.len() == 2 * deleted.len() {
        let (first_half, second_half) = inserted_seq.split_at(deleted.len());
        // Case-insensitive for the same reason as step 1 (#1318): a
        // soft-masked reference must classify as a duplication exactly when the
        // upper-case one does.
        if first_half.eq_ignore_ascii_case(deleted) && second_half.eq_ignore_ascii_case(deleted) {
            return DelinsCanonical::Duplication { start, end };
        }
    }

    // 3. Trim shared affixes on both ends. After this, sub/del/ins/inv are
    //    all classified on the trimmed range, so equivalent edits canonicalize
    //    identically regardless of how much identical context the input
    //    carried.
    let (k_prefix, l_suffix) = shared_affix_lengths(deleted, inserted_seq);
    let trim_start = start + k_prefix;
    let trim_end = end - l_suffix;
    let trim_insert = &inserted_seq[k_prefix..inserted_seq.len() - l_suffix];

    // 4a. Pure insertion (trim consumed the entire deleted range).
    if trim_start == trim_end {
        // Identity (insert == deleted) is the only way greedy trim could leave
        // both sides empty; that case returned at step 1.
        debug_assert!(!trim_insert.is_empty(), "Identity case caught above");
        return DelinsCanonical::Insertion {
            after_index: trim_start,
            sequence: trim_insert.to_vec(),
        };
    }

    // 4b. Pure deletion (trim consumed the entire inserted sequence).
    if trim_insert.is_empty() {
        return DelinsCanonical::Deletion {
            start: trim_start,
            end: trim_end,
        };
    }

    let trim_deleted = &ref_seq[trim_start..trim_end];

    // 4c. 1-base substitution at the trimmed position. Falls through on
    //     non-IUPAC bytes (Base::from_char returns None) because we cannot
    //     express the substitution without a typed Base — better to keep the
    //     variant as a delins than fabricate one. The trim_deleted.len() >= 2
    //     gates on the inversion / dup checks below also exclude this path,
    //     so non-IUPAC 1->1 inputs end up at KeepAsDelins as intended.
    if trim_deleted.len() == 1 && trim_insert.len() == 1 {
        if let (Some(reference), Some(alternative)) = (
            Base::from_char(trim_deleted[0] as char),
            Base::from_char(trim_insert[0] as char),
        ) {
            return DelinsCanonical::Substitution {
                position: trim_start,
                reference,
                alternative,
            };
        }
    }

    // 4d. Inversion at trimmed range (revcomp + outer-pair shortening).
    //     Shared-affix trimming can reveal an inversion that the full-range
    //     check missed, e.g. `ACGAGT -> ACTCGT`: full-range revcomp does not
    //     match (revcomp(ACGAGT) = ACTCGT — does match here), but cases like
    //     a non-palindromic shared prefix only become revcomp after trimming.
    if trim_deleted.len() >= 2
        && trim_insert.len() == trim_deleted.len()
        && is_revcomp(trim_deleted, trim_insert)
    {
        // Invariant: outer-pair shortening of a true revcomp cannot collapse
        // to identity here. A full collapse would mean trim_deleted is its own
        // reverse complement, i.e. trim_deleted == trim_insert. But then
        // greedy shared-affix trimming would have consumed the entire range,
        // leaving trim_start == trim_end (handled by the Insertion / Identity
        // branches above).
        let (s, e) = shorten_inversion(ref_seq, trim_start, trim_end).expect(
            "revcomp delins cannot collapse to identity under shortening; \
             that case is handled by the Identity / Insertion branches above",
        );
        debug_assert!(e > s + 1, "Inversion interval must contain >=2 bases");
        return DelinsCanonical::Inversion { start: s, end: e };
    }

    // 5. Nothing reduced to a higher form. Emit minimal (trimmed) delins.
    DelinsCanonical::KeepAsDelins {
        start: trim_start,
        end: trim_end,
        sequence: trim_insert.to_vec(),
    }
}

/// Compute greedy shared-prefix and shared-suffix lengths between `deleted`
/// and `inserted`, with the constraint that the prefix and suffix together
/// cannot consume more than `min(deleted.len(), inserted.len())` bytes from
/// either side (so the two regions never overlap on either string).
fn shared_affix_lengths(deleted: &[u8], inserted: &[u8]) -> (usize, usize) {
    // Case-insensitive on both sides (#1318): on a soft-masked reference
    // `deleted` arrives lower-case while `inserted` comes from the author's
    // upper-case description, so a raw byte test found no shared affix for
    // sequences that are identical. That mattered beyond a missed trim — it left
    // a palindromic delins untrimmed on the way into the inversion branch, whose
    // `expect` assumes trimming has already handled the collapse-to-identity
    // case, and the mismatch panicked.
    let same = |a: u8, b: u8| a.eq_ignore_ascii_case(&b);
    let max_total = deleted.len().min(inserted.len());
    let mut k = 0;
    while k < max_total && same(deleted[k], inserted[k]) {
        k += 1;
    }
    let mut l = 0;
    while k + l < max_total
        && same(
            deleted[deleted.len() - 1 - l],
            inserted[inserted.len() - 1 - l],
        )
    {
        l += 1;
    }
    (k, l)
}

/// Bytewise reverse-complement equality check.
///
/// Returns true iff `inserted` is the reverse complement of `deleted`, both
/// of equal length, comparing case-insensitively on **both** sides.
///
/// Allocation-free. Folding case matters here for the same reason it does in
/// [`is_complementary_pair`]: on a soft-masked reference `deleted` arrives
/// lower-case while `inserted` comes from the author's upper-case description
/// (#1318).
pub(crate) fn is_revcomp(deleted: &[u8], inserted: &[u8]) -> bool {
    deleted.len() == inserted.len()
        && deleted
            .iter()
            .rev()
            .zip(inserted.iter())
            .all(|(d, i)| complement_base(*d) == Some(i.to_ascii_uppercase()))
}

/// Check if a duplication in a homopolymer should become repeat notation
///
/// When duplicating bases within a homopolymer repeat,
/// the result should be expressed as repeat notation.
/// Result of duplication to repeat conversion
#[derive(Debug, Clone)]
pub enum DupToRepeatResult {
    /// Single-base homopolymer repeat
    Homopolymer {
        base: u8,
        count: u64,
        start: u64, // 1-based (HGVS), use index_to_hgvs_pos() for conversion
        end: u64,   // 1-based (HGVS), inclusive
    },
    /// Multi-base tandem repeat
    TandemRepeat {
        unit: Vec<u8>,
        count: u64,
        start: u64, // 1-based (HGVS), use index_to_hgvs_pos() for conversion
        end: u64,   // 1-based (HGVS), inclusive
    },
    /// Codon-frame gate triggered: structural conditions for repeat
    /// notation are met, but `unit_len % 3 != 0` in c. context, so the
    /// canonical form is a literal `ins` of the duplicated sequence at
    /// the 3' flanking position. `start` and `end` are the two flanking
    /// 1-based HGVS positions (last tract base and the next base).
    GatedInsertion {
        start: u64,
        end: u64,
        sequence: Vec<u8>,
    },
}

pub fn duplication_to_repeat(
    ref_seq: &[u8],
    start: u64,
    end: u64,
    is_coding: bool,
) -> Option<DupToRepeatResult> {
    let start_idx = start as usize;
    let end_idx = end as usize;

    if start_idx >= ref_seq.len() || end_idx > ref_seq.len() || start_idx >= end_idx {
        return None;
    }

    let dup_seq = &ref_seq[start_idx..end_idx];
    if dup_seq.is_empty() {
        return None;
    }

    let dup_len = dup_seq.len();

    // Check if all duplicated bases are the same (homopolymer)
    // IMPORTANT: Only convert to repeat notation when duplicating 2+ bases.
    // Single-base duplications stay as simple dups (e.g., c.5266dup stays as dup, not C[4])
    let first = dup_seq[0];
    if dup_len >= 2 && dup_seq.iter().all(|&b| b == first) {
        // Find the full homopolymer extent
        if let Some(analysis) = find_homopolymer_at(ref_seq, start_idx) {
            if analysis.base == Some(first) {
                // Codon-frame gate (repeated.md): homopolymer unit_len=1 is
                // never codon-aligned, so c. context blocks repeat notation.
                // Emit GatedInsertion so the caller renders as `ins<dup_seq>`
                // at the 3' flanking position of the reference tract.
                if is_coding {
                    let last_tract_idx = analysis.ref_start + analysis.ref_count as usize - 1;
                    return Some(DupToRepeatResult::GatedInsertion {
                        start: index_to_hgvs_pos(last_tract_idx),
                        end: index_to_hgvs_pos(last_tract_idx) + 1,
                        sequence: dup_seq.to_vec(),
                    });
                }
                let total_count = analysis.ref_count + dup_len as u64;
                return Some(DupToRepeatResult::Homopolymer {
                    base: first,
                    count: total_count,
                    start: index_to_hgvs_pos(analysis.ref_start),
                    end: index_to_hgvs_pos(analysis.ref_start + analysis.ref_count as usize - 1),
                });
            }
        }
    }

    // Check for multi-base tandem repeat
    // The duplicated sequence could be multiple copies of a repeat unit
    // IMPORTANT: Only convert to repeat notation when duplicating MULTIPLE copies
    // of the repeat unit. Single-copy duplications stay as simple dups.
    // Try different unit lengths from 1 to half the dup length
    for unit_len in 1..=dup_len / 2 {
        if !dup_len.is_multiple_of(unit_len) {
            continue;
        }

        let unit = &dup_seq[0..unit_len];
        let copies_in_dup = dup_len / unit_len;

        // Only convert to repeat if duplicating 2+ copies of the unit
        if copies_in_dup < 2 {
            continue;
        }

        // Check if dup_seq is made of repeated copies of unit
        let is_repeat = (0..copies_in_dup).all(|i| {
            let chunk = &dup_seq[i * unit_len..(i + 1) * unit_len];
            chunk == unit
        });

        if !is_repeat {
            continue;
        }

        // Found a repeat unit. Now find the full extent in the reference.
        if let Some((ref_count, rep_start, rep_end)) =
            count_tandem_repeats(ref_seq, start_idx, unit)
        {
            // Codon-frame gate (repeated.md): in c., repeat notation requires
            // unit_len % 3 == 0. Structural conditions are met but the gate
            // forces a literal `ins<dup_seq>` at the 3' tract flanking
            // position instead.
            if is_coding && !unit_len.is_multiple_of(3) {
                let last_tract_idx = rep_end - 1;
                return Some(DupToRepeatResult::GatedInsertion {
                    start: index_to_hgvs_pos(last_tract_idx),
                    end: index_to_hgvs_pos(last_tract_idx) + 1,
                    sequence: dup_seq.to_vec(),
                });
            }
            // Total count = reference count + duplicated copies
            let total_count = ref_count + copies_in_dup as u64;
            // rep_end is exclusive (0-based), so last position is rep_end - 1
            return Some(DupToRepeatResult::TandemRepeat {
                unit: unit.to_vec(),
                count: total_count,
                start: index_to_hgvs_pos(rep_start),
                end: index_to_hgvs_pos(rep_end - 1),
            });
        }
    }

    // Note: We do NOT convert single-copy tandem repeat duplications to repeat notation
    // e.g., c.360_362dupGCA stays as dup, not GCA[8]

    None
}

/// Find tandem repeat at a position for a given repeat unit
///
/// Given a position and a repeat unit sequence, finds how many times that
/// sequence is repeated at that position in the reference.
///
/// Returns (count, start, end) where:
/// - count is the number of repeats found
/// - start is the 0-indexed start of the repeat region
/// - end is the 0-indexed exclusive end
///
/// # Anchoring contract (do NOT replace with `extend_tandem_tract`)
///
/// This finder scans **backward** from `pos` in unit-length steps, so it accepts
/// a tract that merely *ends at or before* `pos` — `pos` itself need not match
/// `unit`. `extend_tandem_tract` (used by `find_tandem_extent`) instead
/// *requires* the anchor span `ref_seq[pos..pos+unit_len]` to equal `unit`. The
/// two are therefore **not interchangeable** (#866): e.g.
/// `count_tandem_repeats(b"ATATXX", 4, b"AT") == Some((2, 0, 4))` (the `AT`
/// tract ends at index 4) while `extend_tandem_tract(b"ATATXX", 4..6, b"AT")
/// == None` (`ref[4..6] == "XX"`). `normalize_repeat` relies on this
/// boundary-anchor behavior, so delegating here would silently drop valid
/// tracts. The genuinely shared piece across the ins→repeat / repeat paths is
/// the 3'-phase step `three_prime_align_tract` (#852), not the tract finder.
///
/// # Case
///
/// Both scans compare with [`slice::eq_ignore_ascii_case`], because `ref_seq`
/// comes from the reference and `repeat_unit` comes from the description. A
/// reference FASTA is routinely soft-masked and a description is conventionally
/// written in upper case, so a raw byte comparison finds no run at all on a
/// masked tract — and, since this returns coordinates rather than a verdict,
/// that surfaces as a *different canonical form* rather than as an error
/// (#1491). Four of fifteen repeat inputs measured changed answer with case
/// alone: `g.259_267CAG[2]` lowered to `g.265_267del` on an upper-case
/// reference and stayed an unlowered repeat on the soft-masked twin, and three
/// single-position anchors widened to their tract on one and not the other.
///
/// This is the convention the rest of this module already follows — the dup,
/// delins-trimming and inversion-typing comparisons all fold — and folding
/// costs nothing on an unmasked reference, where `eq_ignore_ascii_case` and
/// `==` agree by definition. Note also that only *coordinates* leave this
/// function, so folding here cannot put a lower-case base into an emitted
/// sequence.
///
/// **Case is all this folds — it is not alphabet normalization.**
/// `eq_ignore_ascii_case` maps `u` to `U`, never `U` to `T`, so an `r.` tract
/// whose reference is spelled with uracil is still invisible here. The SPDI
/// call site handles that separately by running the fetched window through
/// `apply_alphabet` before it ever calls this function (#1452,
/// `spdi::convert::fetch_normalized_reference_bases`), which folds case *and*
/// rewrites `U` to `T`. The two are complementary rather than redundant, and
/// the overlap is easy to misread as duplication: with this fold in place,
/// reverting the SPDI-side normalization leaves every soft-masking test green
/// and fails exactly one — `a_uracil_spelled_tract_is_found_on_the_rna_axis`.
/// That single test is what pins the difference, so do not "simplify" either
/// site away on the grounds that the other covers it.
pub fn count_tandem_repeats(
    ref_seq: &[u8],
    pos: usize,
    repeat_unit: &[u8],
) -> Option<(u64, usize, usize)> {
    if repeat_unit.is_empty() || pos >= ref_seq.len() {
        return None;
    }

    let unit_len = repeat_unit.len();

    // Check if the repeat unit matches at the position
    if pos + unit_len > ref_seq.len() {
        return None;
    }

    // Try to find repeats starting at or before this position
    // First, scan backwards to find the start of the repeat region
    let mut start = pos;
    while start >= unit_len {
        let candidate = &ref_seq[start - unit_len..start];
        if candidate.eq_ignore_ascii_case(repeat_unit) {
            start -= unit_len;
        } else {
            break;
        }
    }

    // Now count forward from the start
    let mut end = start;
    let mut count = 0u64;
    while end + unit_len <= ref_seq.len() {
        if ref_seq[end..end + unit_len].eq_ignore_ascii_case(repeat_unit) {
            count += 1;
            end += unit_len;
        } else {
            break;
        }
    }

    if count >= 1 {
        Some((count, start, end))
    } else {
        None
    }
}

/// Result of repeat normalization
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RepeatNormResult {
    /// Convert to deletion (specified count < reference count)
    Deletion {
        start: u64, // 1-based (HGVS), inclusive start of region to delete
        end: u64,   // 1-based (HGVS), inclusive end of region to delete
    },
    /// Convert to duplication (specified count = reference count + 1)
    Duplication {
        start: u64, // 1-based (HGVS), inclusive start of duplicated region
        end: u64,   // 1-based (HGVS), inclusive end of duplicated region
        sequence: Vec<u8>,
    },
    /// Convert to insertion. Used when the codon-frame gate blocks repeat
    /// notation in c. context for an expansion of >=2 unit copies. `start`
    /// and `end` are the two flanking 1-based HGVS positions; `sequence` is
    /// the literal inserted sequence (canonical unit repeated `specified -
    /// ref_count` times).
    Insertion {
        start: u64,
        end: u64,
        sequence: Vec<u8>,
    },
    /// Keep as repeat notation with canonical position
    Repeat {
        start: u64, // 1-based (HGVS), inclusive start
        end: u64,   // 1-based (HGVS), inclusive end
        sequence: Vec<u8>,
        count: u64,
    },
    /// No change needed
    Unchanged,
}

/// Apply the HGVS 3'-rule to a tandem-repeat tract, returning the most-3' phase.
///
/// `ref_start`/`ref_end` are 0-based with `ref_end` **exclusive**:
/// `ref_seq[ref_start..ref_end]` is a whole number of copies of `unit`. While the
/// base immediately past the tract equals the rotating unit's head, the window
/// slides one base 3' and the unit rotates left by one — a pure re-phasing that
/// preserves the width (copy count). Returns the rotated
/// `(start, end_exclusive, rotated_unit)`.
///
/// Shared by [`normalize_repeat`] and [`insertion_to_repeat`] so their
/// ins->repeat canonicalization agrees in a single pass (idempotency; #852, #866).
fn three_prime_align_tract(
    ref_seq: &[u8],
    ref_start: usize,
    ref_end: usize,
    unit: &[u8],
) -> (usize, usize, Vec<u8>) {
    let mut start = ref_start;
    let mut end = ref_end;
    let mut rotated = unit.to_vec();
    while end < ref_seq.len() && !rotated.is_empty() && ref_seq[end] == rotated[0] {
        rotated.rotate_left(1);
        start += 1;
        end += 1;
    }
    (start, end, rotated)
}

/// 5' mirror of [`three_prime_align_tract`] for `ShuffleDirection::FivePrime`.
///
/// While the base immediately 5' of the tract equals the rotating unit's *tail*,
/// the window slides one base 5' and the unit rotates right by one — the
/// direction-mirrored re-phasing that preserves the copy count. Returns the
/// rotated `(start, end_exclusive, rotated_unit)`.
///
/// This is what keeps a `FivePrime` repeat idempotent AND confluent over an
/// ambiguous alternating tract: `AG[4]` anchored at the tract's 3' phase and
/// `GA[4]` anchored at its 5' phase both re-anchor to the single 5'-most form
/// (#8, the tandem-repeat mirror of the #403 direction tie-break). `ThreePrime`
/// keeps [`three_prime_align_tract`], so 3' output is unchanged.
fn five_prime_align_tract(
    ref_seq: &[u8],
    ref_start: usize,
    ref_end: usize,
    unit: &[u8],
) -> (usize, usize, Vec<u8>) {
    let mut start = ref_start;
    let mut end = ref_end;
    let mut rotated = unit.to_vec();
    while start > 0
        && rotated
            .last()
            .is_some_and(|&tail| ref_seq[start - 1] == tail)
    {
        rotated.rotate_right(1);
        start -= 1;
        end -= 1;
    }
    (start, end, rotated)
}

/// Normalize a repeat variant
///
/// Re-base an anchored repeat's count from the tract the caller spelled onto
/// the tract the rotation search named — or decline the re-phase.
///
/// `[N]` counts copies of the unit the caller **spelled**, against that unit's
/// own maximal tract; the emitted description counts against the winning
/// rotation's tract, so the difference has to be carried or reference copies
/// are silently created or deleted (#1618). `N + winner - literal` leaves the
/// net base change unchanged, and a pure re-phase (`DNA/repeated.md:97`'s 11
/// `TG` -> 11 `GT`) has the two counts equal and so is untouched.
///
/// `None` means **decline** — the honest answer for two shapes, neither of
/// which any count on the winner's tract can express:
///
/// * **The tracts are disjoint.** The description's window has moved to a
///   different stretch of reference, so re-basing conserves the net length and
///   nothing else. On `AATGTGTGGTAA`, `g.265TG[6]` re-based onto the 1-copy
///   `GT` tract is `g.265_266GT[4]`, which denotes `TGTGTGGTGTGTGT` where the
///   input denotes `TGTGTGTGTGTGGT` — the same 14 bases in a different order,
///   as `spdi::compare_denoted_sequences` reports.
/// * **The re-based count would underflow.** `g.265TG[1]` on that same
///   reference is `1 + 1 - 3`. The `saturating_sub` this replaces clamped it to
///   `0` and routed to a two-base deletion, where the input denotes **four**
///   bases fewer — a silent loss of two bases, in the arm the re-basing rule
///   itself added. (The saturation's original comment justified it by the
///   disjoint case, where `literal == winner` and nothing ever saturates.)
///
/// The caller answers a decline by keeping the spelled unit's own tract, which
/// is what `[N]` was written against, so the output denotes exactly the input's
/// bases.
///
/// A `None` literal tract takes **no** adjustment: the spelled unit tiles
/// nowhere at or 5' of the anchor, so there are no copies to re-base against,
/// and `hgvs_to_spdi` reads such an input as denoting no sequence at all — no
/// baseline exists to conserve. `the_none_literal_arm_takes_no_adjustment` pins
/// the asymmetry that leaves, as a limit rather than as a fix.
fn rebase_count_onto_tract(
    specified_count: u64,
    literal_tract: Option<(u64, usize, usize)>,
    winner: (u64, usize, usize),
) -> Option<u64> {
    let Some((literal_count, literal_start, literal_end)) = literal_tract else {
        return Some(specified_count);
    };
    let (winner_count, winner_start, winner_end) = winner;
    // Half-open spans, so touching-but-not-overlapping is disjoint.
    if literal_end <= winner_start || winner_end <= literal_start {
        return None;
    }
    specified_count
        .checked_add(winner_count)?
        .checked_sub(literal_count)
}

/// Given a repeat notation like CAT[1], determines the appropriate
/// normalized representation by comparing to the reference.
///
/// Arguments:
/// - ref_seq: The reference sequence (0-indexed)
/// - pos: The 0-indexed start of the input repeat range
/// - end_pos: The 0-indexed inclusive end of the input repeat range. For a
///   single-position anchor (`c.4CAT[1]`, where only the repeat start is named)
///   this equals `pos`; for an explicit range (`g.3_4GT[5]`, `c.*16_*18T[8]`) it
///   is the range's last base. See the flank-absorption note below.
/// - repeat_unit: The repeated sequence (e.g., b"CAT")
/// - specified_count: The count specified in the variant
///
/// Returns the normalized representation
///
/// Shares the 3'-phase step (`three_prime_align_tract`) with
/// [`insertion_to_repeat`]; kept in lockstep by the cross-check test
/// `insertion_to_repeat_agrees_with_normalize_repeat_phase`. Uses
/// [`count_tandem_repeats`] (+ an off-phase rotation/offset search) as its
/// tract-finder rather than `find_tandem_extent` because it anchors at an
/// existing position (possibly mid- or end-of-tract), not an insertion point —
/// see #866 and the note on [`count_tandem_repeats`].
pub fn normalize_repeat(
    ref_seq: &[u8],
    pos: usize,
    end_pos: usize,
    repeat_unit: &[u8],
    specified_count: u64,
    is_coding: bool,
    direction: crate::normalize::config::ShuffleDirection,
) -> RepeatNormResult {
    // Match `count_tandem_repeats`'s pre-refactor contract: an empty unit
    // is meaningless and falls through to `Unchanged`. Without this guard,
    // `smallest_repeat_unit(b"")` returns `b""` and the period division below
    // panics.
    if repeat_unit.is_empty() {
        return RepeatNormResult::Unchanged;
    }

    // Canonicalize to the smallest period so callers spelling a non-minimal
    // unit (e.g. `ATAT[1]` over an `AT[4]` tract) reach the right branch.
    // Without this, a contraction misses B2 because k is computed in the
    // caller-spelled unit (1 ATAT removed) instead of the canonical unit
    // (2 ATs removed).
    let canonical_unit = smallest_repeat_unit(repeat_unit);
    let copies_per_input_unit = (repeat_unit.len() / canonical_unit.len()) as u64;
    let mut specified_count = specified_count * copies_per_input_unit;

    // Count how many times the canonical unit appears in the reference.
    //
    // For a single-position anchor (`end_pos == pos`) the unit can land off the
    // tract's phase boundary — a reverse-complemented projected repeat puts its
    // single anchor at the genomic 3' end of the tract, often mid-unit (#852).
    // Examples vs the `GT`/`GC` tracts: `g.4913TG[5]` (anchor reads `GT`, not the
    // queried `TG`) and `g.4920GC[5]` (anchor is the tract's last base, so a
    // boundary-anchored count undercounts). `count_tandem_repeats` requires the
    // unit to tile *at* the anchor, so the literal lookup misses or undercounts.
    // Search unit rotations AND small anchor offsets (`pos - a`) for the maximal
    // tandem tract that CONTAINS `pos`, mirroring `insertion_to_repeat`'s rotation
    // handling; `three_prime_align_tract` below then fixes the final phase. Bounded
    // (`unit_len^2` probes). An explicit range is spelled unit-aligned by the
    // caller, so it keeps the literal unit (no search) — only length>=2 units on a
    // single anchor can be off-phase.
    //
    // The search is a MAXIMIZATION, not a fallback (it used to run only when the
    // literal lookup returned `None`). A literal match that succeeds can still be
    // the *shorter* phase: for `…ACTGAAGTGTTACTTT…` anchored at the tract's last
    // base, literal `TG` tiles 1 copy while the `GT` rotation tiles 2 copies
    // across the same anchor. Taking the literal hit unconditionally made
    // normalization non-idempotent — pass 1 emitted a window covering only the
    // short tract, and pass 2 read that as an explicit range, re-discovered the
    // real (longer) tract and absorbed the extra copies into the count, so window
    // AND count both grew (`r.-124ug[14]` → `r.-125_-124gu[14]` →
    // `r.-127_-124gu[15]`). Seed `best` with the literal result and let a strictly
    // longer spanning rotation displace it; ties keep the literal, so the only
    // behavior change is preferring a genuinely larger tract.
    let mut working_unit = canonical_unit.to_vec();
    let literal_tract = count_tandem_repeats(ref_seq, pos, &working_unit);
    // Copies of the caller's own unit at the anchor, kept so the count can be
    // re-based onto whatever tract the maximization below actually names.
    //
    // `[N]` counts copies of the unit the caller SPELLED, against that unit's
    // own tract; re-phasing to a rotation that tiles a different tract changes
    // what the count counts, so the difference has to be carried or reference
    // copies are silently created or deleted (#1618).
    //
    // `DNA/repeated.md:97-99` is the spec's own worked case: re-spelling 11 `TG`
    // copies as `GT` slides the tract one base 3', swallowing the first `T` of
    // the neighbouring run, and the spec decrements that run `T[7]` -> `T[6]` so
    // the total is conserved. Same principle, applied to a re-phase that changes
    // the tract's size rather than to a pure slide.
    let tract = match literal_tract {
        _ if end_pos == pos && working_unit.len() >= 2 => {
            // Hold the literal seed to the SAME span test the rotations below
            // must pass. `count_tandem_repeats` can return a tract lying entirely
            // 5' of the anchor — e.g. `ref="ACACTT", pos=4, unit="AC"` back-scans
            // to `start=0` and counts forward to `end=4`, so `pos < e` is false.
            // Seeding `best` with such a tract would let a non-spanning literal
            // out-count, and thereby suppress, a rotation that genuinely spans the
            // anchor.
            let mut best: Option<(u64, usize, usize, Vec<u8>)> = literal_tract
                .filter(|&(_, s, e)| s <= pos && pos < e)
                .map(|(c, s, e)| (c, s, e, canonical_unit.to_vec()));
            for a in 0..canonical_unit.len() {
                let Some(anchor) = pos.checked_sub(a) else {
                    continue;
                };
                for r in 0..canonical_unit.len() {
                    let mut rotated = canonical_unit.to_vec();
                    rotated.rotate_left(r);
                    if let Some((c, s, e)) = count_tandem_repeats(ref_seq, anchor, &rotated) {
                        // Only accept a tract that actually spans the anchor, and
                        // keep the longest such tract.
                        if s <= pos && pos < e && best.as_ref().is_none_or(|(bc, ..)| c > *bc) {
                            best = Some((c, s, e, rotated));
                        }
                    }
                }
            }
            match best {
                Some((c, s, e, u)) => {
                    // Re-base the count onto the tract the maximization named,
                    // or decline the re-phase when no count on that tract can
                    // denote the input's bases. The literal tract is the
                    // baseline even when it does not span the anchor — that is
                    // the #1618 shape exactly (`GTGT` at g.259..262 is 2 copies
                    // of `GT` but only 1 of `TG`, whose tract ends AT the anchor
                    // and so was filtered out of `best` above).
                    //
                    // **Copies, not bases.** This was first written as a
                    // per-edge base measurement,
                    // `(ls.saturating_sub(s) + e.saturating_sub(le)) / u.len()`,
                    // which reads the two edges independently and so cannot see
                    // a tract that NARROWS. Where the literal tract spans the
                    // anchor the two agree exactly — `best` is seeded with it,
                    // so a winner must be strictly longer, and back/forward-scan
                    // maximality bounds each edge's disagreement to under one
                    // unit, which the truncating division then absorbs
                    // (exhaustively checked: 0 disagreements over ~590k
                    // (sequence, anchor, unit) triples). Where it does NOT span
                    // the anchor there is no seed, so the winner may hold FEWER
                    // copies than the literal tract, and the base form reports a
                    // gain for a loss: on `AATGTGTGGTAA`, `g.265TG[6]` re-phases
                    // a 3-copy `TG` tract onto a 1-copy `GT` one and the base
                    // form added a copy instead of removing two.
                    //
                    // **Where the literal tract is `None`, nothing is re-based**
                    // — see the note on `count_tandem_repeats` above and
                    // `the_none_literal_arm_takes_no_adjustment`, which pins the
                    // resulting asymmetry as a limit. The spelled unit tiles
                    // nowhere at or 5' of the anchor, so there are no copies to
                    // count the caller's `[N]` against and `hgvs_to_spdi` reads
                    // such an input as denoting no sequence at all: there is no
                    // baseline to conserve, which is why inventing one here
                    // would be a ruling rather than a rescue.
                    match rebase_count_onto_tract(specified_count, literal_tract, (c, s, e)) {
                        Some(rebased) => {
                            specified_count = rebased;
                            working_unit = u;
                            Some((c, s, e))
                        }
                        // The re-phase cannot carry the count. Keep the spelled
                        // unit's own tract, which is what `[N]` was written
                        // against, so the output denotes exactly the input's
                        // bases.
                        None => literal_tract,
                    }
                }
                None => None,
            }
        }
        // Explicit range, or a 1-base unit: the caller spelled it unit-aligned,
        // so the literal lookup stands as-is.
        other => other,
    };
    let Some((ref_count, mut ref_start, mut ref_end)) = tract else {
        return RepeatNormResult::Unchanged;
    };
    let canonical_unit: &[u8] = &working_unit;

    // Flank absorption for an under-specified explicit range. `[N]` counts the
    // alt copies of the *stated* reference units; reference repeat copies inside
    // the maximal tract but outside the stated range are untouched by the edit
    // and remain in the variant, so they must be added to the canonical count
    // (which re-anchors to the full tract `count_tandem_repeats` found). This is
    // mutalyzer-verified real-world behavior: `LRG_303:g.3_4GT[5]` ->
    // `g.1_4GT[6]` (the 5' copy g.1_2 absorbs, 5+1=6); `c.*16_*18T[8]` ->
    // `c.*16_*19T[9]` (the 3' copy *19 absorbs). A single-position anchor
    // (`end_pos == pos`, e.g. `c.4CAT[1]`) names only the repeat start and means
    // "the whole repeat -> N copies", so it absorbs nothing. A well-formed
    // explicit range (stated span == true tract) also absorbs nothing (the
    // flanks below are zero).
    if end_pos > pos {
        let unit_len = canonical_unit.len();
        let three_prime_bases = ref_end.saturating_sub(end_pos + 1);
        let five_prime_bases = pos.saturating_sub(ref_start);
        let absorbed = ((three_prime_bases + five_prime_bases) / unit_len) as u64;
        specified_count += absorbed;
    }

    // Direction-aware unit rotation. For `ThreePrime`, repeated.md L44 ("applying
    // the 3'rule, the repeat has to be described as an AGC repeat") — slide to the
    // 3'-most phase via the shared helper, so `insertion_to_repeat` reaches the
    // identical phase (idempotency; #852). For `FivePrime`, slide to the mirror
    // 5'-most phase so a `dup`→`unit[N]` emission (which the shuffle stage already
    // 5'-anchors) is a fixed point rather than drifting 3' (#8). `ref_end` stays
    // exclusive (the later `index_to_hgvs_pos(ref_end - 1)` is unchanged).
    use crate::normalize::config::ShuffleDirection;
    let (ns, ne, rotated_unit) = match direction {
        ShuffleDirection::ThreePrime => {
            three_prime_align_tract(ref_seq, ref_start, ref_end, canonical_unit)
        }
        ShuffleDirection::FivePrime => {
            five_prime_align_tract(ref_seq, ref_start, ref_end, canonical_unit)
        }
    };
    ref_start = ns;
    ref_end = ne;
    let canonical_unit: &[u8] = &rotated_unit;

    let unit_len = canonical_unit.len() as u64;
    // HGVS spec (repeated.md): in c. context, repeat notation requires
    // unit_len % 3 == 0. When the gate blocks, contraction-with-survivors
    // routes to Deletion and expansion-of->=2-copies routes to Insertion.
    let codon_blocks_repeat = is_coding && !canonical_unit.len().is_multiple_of(3);

    if specified_count < ref_count {
        let k = ref_count - specified_count;
        if k >= 2 && specified_count >= 1 && !codon_blocks_repeat {
            // B2 (symmetric with A7): >=2 unit reduction with surviving units → repeat
            RepeatNormResult::Repeat {
                start: index_to_hgvs_pos(ref_start),
                end: index_to_hgvs_pos(ref_end - 1),
                sequence: canonical_unit.to_vec(),
                count: specified_count,
            }
        } else {
            // 1-unit reduction, full tract removal, or codon-frame-gated → deletion
            // (HGVS prioritization: deletion outranks unranked repeat[0])
            //
            // Removing any `k` copies of a pure tandem unit yields the same
            // haplotype, so the direction rule picks which copies are named: the
            // 3'-most under `ThreePrime`, the 5'-most under `FivePrime`. Naming
            // the 3'-most copies unconditionally made the emitted `del` a
            // non-fixed-point under 5' (it re-shuffled to the 5'-most slot).
            let del_len = (k as usize) * unit_len as usize;
            let (del_start_idx, del_end_idx) = match direction {
                ShuffleDirection::ThreePrime => (ref_end - del_len, ref_end - 1),
                ShuffleDirection::FivePrime => (ref_start, ref_start + del_len - 1),
            };
            RepeatNormResult::Deletion {
                start: index_to_hgvs_pos(del_start_idx),
                end: index_to_hgvs_pos(del_end_idx),
            }
        }
    } else if specified_count == ref_count + 1 {
        // Convert to duplication - we're adding exactly one copy.
        // dup is always permitted; the spec exception only forbids `[N]`.
        //
        // Duplicating any one copy of a pure tandem unit yields the same
        // haplotype, so the direction rule picks which copy is named: the last
        // (3'-most) under `ThreePrime`, the first (5'-most) under `FivePrime`.
        // Naming the last copy unconditionally made the emitted `dup` a
        // non-fixed-point under 5' (it re-shuffled to the 5'-most slot).
        // `ref_end` is exclusive, so the tract's last position is `ref_end - 1`.
        let (dup_start_idx, dup_end_idx) = match direction {
            ShuffleDirection::ThreePrime => (ref_end - canonical_unit.len(), ref_end - 1),
            ShuffleDirection::FivePrime => (ref_start, ref_start + canonical_unit.len() - 1),
        };
        RepeatNormResult::Duplication {
            start: index_to_hgvs_pos(dup_start_idx),
            end: index_to_hgvs_pos(dup_end_idx),
            sequence: canonical_unit.to_vec(),
        }
    } else if specified_count == ref_count {
        // Same as reference - this is identity (no change)
        RepeatNormResult::Unchanged
    } else if codon_blocks_repeat {
        // Expansion of >=2 copies in c. with non-codon-aligned unit: spec
        // mandates `ins<literal>` form (e.g., c.1741_1742insTATATATA), not
        // `[N]`.
        //
        // Adding the copies at either end of a pure tandem tract yields the same
        // haplotype, so the direction rule picks the flank the insertion is named
        // against: the 3' flank (between `ref_end - 1` and `ref_end`) under
        // `ThreePrime`, the 5' flank (between `ref_start - 1` and `ref_start`)
        // under `FivePrime`. Naming the 3' flank unconditionally made the emitted
        // `ins` a non-fixed-point under 5' (it re-shuffled to the 5' flank).
        //
        // A tract flush against the sequence start (`ref_start == 0`) has no 5'
        // flank base to name, so it keeps the 3' flank; that is the transcript-
        // start boundary case the CDS-start clamp in `normalize_na_edit` owns.
        let added_copies = specified_count - ref_count;
        let mut inserted = Vec::with_capacity((added_copies as usize) * canonical_unit.len());
        for _ in 0..added_copies {
            inserted.extend_from_slice(canonical_unit);
        }
        let flank_left = match direction {
            ShuffleDirection::FivePrime if ref_start > 0 => index_to_hgvs_pos(ref_start - 1),
            _ => index_to_hgvs_pos(ref_end - 1),
        };
        let flank_right = flank_left + 1;
        RepeatNormResult::Insertion {
            start: flank_left,
            end: flank_right,
            sequence: inserted,
        }
    } else {
        // Default expansion (>=2 copies, not gated): repeat notation.
        // The repeat region describes the REFERENCE tract (per HGVS spec)
        // ref_end is exclusive, so last position is ref_end - 1
        RepeatNormResult::Repeat {
            start: index_to_hgvs_pos(ref_start),
            end: index_to_hgvs_pos(ref_end - 1),
            sequence: canonical_unit.to_vec(),
            count: specified_count,
        }
    }
}

/// Get the canonical form for an edit
///
/// HGVS rules:
/// - Deletions ALWAYS stay as deletions (just shift 3')
/// - Insertions become duplications if they match adjacent sequence
/// - Duplications stay as duplications
/// - Delins stays as delins
pub fn get_canonical_form(edit: &NaEdit, ref_seq: &[u8], start: u64, _end: u64) -> CanonicalForm {
    use crate::hgvs::edit::InsertedSequence;

    match edit {
        NaEdit::Deletion { .. } => {
            // Deletions always stay as deletions - never convert to dup
            CanonicalForm::Deletion
        }
        NaEdit::Insertion { sequence } => {
            // Insertions may become duplications if they match adjacent sequence
            // Only check Literal sequences (others like Count, Range don't have actual bases)
            if let InsertedSequence::Literal(seq) = sequence {
                // Convert Sequence (Vec<Base>) to bytes
                let seq_bytes: Vec<u8> = seq.bases().iter().map(|b| *b as u8).collect();
                if insertion_is_duplication(ref_seq, start, &seq_bytes) {
                    return CanonicalForm::Duplication;
                }
            }
            CanonicalForm::Insertion
        }
        NaEdit::Delins { .. } => CanonicalForm::Delins,
        NaEdit::Duplication { .. } => CanonicalForm::Duplication,
        _ => CanonicalForm::Deletion, // Default fallback
    }
}

/// Apply minimal notation rules to an edit without requiring reference data.
///
/// This function applies HGVS minimal notation rules:
/// - Deletions: remove explicit sequence and length (del12 → del, delATG → del)
/// - Delins: remove explicit deleted sequence from delins notation
///   (delATGinsGGG → delinsGGG)
/// - Duplications: remove explicit sequence and length (dupATG → dup, dup12 → dup)
///
/// This is useful for canonicalizing variants even when reference sequence
/// is not available for full 3'/5' normalization.
pub fn canonicalize_edit(edit: &NaEdit) -> NaEdit {
    match edit {
        NaEdit::Deletion { .. } => NaEdit::Deletion {
            sequence: None,
            length: None,
        },
        NaEdit::Duplication {
            uncertain_extent, ..
        } => NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: uncertain_extent.clone(),
        },
        NaEdit::Delins { sequence, .. } => {
            // Keep only the inserted sequence; strip the explicit deleted
            // sequence/length per the HGVS spec recommendation that the
            // short form `delinsXXX` is preferred.
            NaEdit::Delins {
                sequence: sequence.clone(),
                deleted: None,
                deleted_length: None,
                substitution_reference: None,
            }
        }
        // Inversion: §HGVS v21.0 DNA/inversion.md — "the recommendation
        // is not to describe the inverted nucleotide sequence." Strip
        // both `sequence` and `length` so `g.100_105invATGCC` /
        // `g.100_105inv5` collapse to the canonical `g.100_105inv`.
        // Closes-after: #352.
        NaEdit::Inversion { .. } => NaEdit::Inversion {
            sequence: None,
            length: None,
        },
        // A4: a substitution where ref == alt (e.g. `c.100A>A`) is degenerate;
        // the HGVS v21 spec calls the form "not allowed" (recommendations/DNA/
        // other.md) and gives `c.100=` as the canonical alternative. The rule
        // is purely syntactic on the edit's stated bases, so it lives here in
        // the no-reference canonicalization path alongside the other minimal-
        // notation rewrites — it must fire regardless of provider availability.
        NaEdit::Substitution {
            reference,
            alternative,
        } if reference == alternative => NaEdit::position_identity(),
        // Other edits pass through unchanged
        _ => edit.clone(),
    }
}

/// Apply minimal notation to a variant without reference data.
/// Returns true if the edit was modified.
pub fn should_canonicalize(edit: &NaEdit) -> bool {
    match edit {
        NaEdit::Deletion { sequence, length } => sequence.is_some() || length.is_some(),
        NaEdit::Duplication {
            sequence, length, ..
        } => sequence.is_some() || length.is_some(),
        // Companion to the Delins arm in `canonicalize_edit`: a top-level
        // `<a>_<b>delXinsY` (or `<a>_<b>del<N>insY`) is redundant per
        // recommendations/DNA/delins.md — the spec says `<a>_<b>delinsY`
        // is the canonical form. Strip the explicit deleted bases/length
        // regardless of provider availability. Closes #338.
        NaEdit::Delins {
            deleted,
            deleted_length,
            ..
        } => deleted.is_some() || deleted_length.is_some(),
        // Companion to the Inversion arm in `canonicalize_edit`: strip
        // redundant `sequence` / `length` per §DNA/inversion.md. Closes-
        // after: #352.
        NaEdit::Inversion { sequence, length } => sequence.is_some() || length.is_some(),
        // Companion to the A4 arm in `canonicalize_edit`: route degenerate
        // substitutions through the no-reference canonicalize path so the
        // rewrite fires even when the provider is empty.
        NaEdit::Substitution {
            reference,
            alternative,
        } => reference == alternative,
        _ => false,
    }
}

/// Canonicalize a `con` (conversion) edit to its SVD-WG009 `delins` equivalent.
///
/// Per HGVS Nomenclature DNA delins recommendations and Community
/// Consultation SVD-WG009 (accepted Nov 2020), the `con` form is "no longer
/// used" and should be described as `delins`. This helper performs the
/// purely-syntactic rewrite:
///
/// - Same-reference position-only source (e.g. `42536337_42536382`):
///   emits `Delins{PositionRange{start, end}}`, displayed as
///   `delins42536337_42536382`.
/// - Cross-reference source (anything else, e.g. `NM_000089.1:c.789_1011`):
///   emits `Delins{Reference(source)}`, displayed as
///   `delins[NM_000089.1:c.789_1011]`. The square brackets and source
///   reference-type prefix are required by SVD-WG009.
///
/// Returns `None` for any non-`Conversion` edit so callers can use it as
/// an early-return probe.
///
/// This is a pure transformation; it does not require reference data and
/// does not validate the source. Validation is delegated to the existing
/// `delins`-source parser on any subsequent round-trip.
pub fn canonicalize_conversion_to_delins(edit: &NaEdit) -> Option<NaEdit> {
    use crate::hgvs::edit::InsertedSequence;

    let source = match edit {
        NaEdit::Conversion { source } => source,
        _ => return None,
    };

    // Try same-reference position-only form: `<digits>_<digits>`.
    // HGVS positions are 1-based and ordered, so only emit a PositionRange
    // when both endpoints are >= 1 and start <= end. Anything else falls
    // through to Reference so we don't fabricate a structurally-invalid
    // delins range (e.g. `0_0`, reversed `10_2`).
    if let Some((s, e)) = source.split_once('_') {
        if !s.is_empty()
            && !e.is_empty()
            && s.bytes().all(|b| b.is_ascii_digit())
            && e.bytes().all(|b| b.is_ascii_digit())
        {
            if let (Ok(start), Ok(end)) = (s.parse::<u64>(), e.parse::<u64>()) {
                if start >= 1 && end >= 1 && start <= end {
                    return Some(NaEdit::Delins {
                        sequence: InsertedSequence::PositionRange { start, end },
                        deleted: None,
                        deleted_length: None,
                        substitution_reference: None,
                    });
                }
            }
        }
    }

    // Cross-reference source (or anything that isn't a clean integer pair):
    // wrap in `Reference`, which displays as `[<source>]`.
    //
    // `InsertedSequence::Reference` carries the *inner* payload (e.g.
    // `NC_000022.11:g.100_200`) and its `Display` impl wraps the value
    // in a single `[...]` pair. The parser for `con<source>` captures
    // the entire token after `con` verbatim, so a bracketed input like
    // `con[NC_...:g.X_Y]` lands here with `source == "[NC_...:g.X_Y]"`.
    // Storing that as-is in `Reference` would yield `delins[[NC_...]]`
    // on round-trip (doubled brackets) and leak the same doubling into
    // any downstream `format!("ins[{}]", reference)` error message.
    // Strip a single matching outer `[...]` pair so the stored payload
    // is exactly the unbracketed inner reference, matching what the
    // `delins[...]` parser path stores for the same shape (issue #333).
    let inner = strip_outer_brackets(source);
    Some(NaEdit::Delins {
        sequence: InsertedSequence::Reference(inner.to_string()),
        deleted: None,
        deleted_length: None,
        substitution_reference: None,
    })
}

/// Map the RNA base `U`/`u` to thymine (`T`) throughout an `r.` edit's literal
/// base sequences, so the edit can be compared against the DNA-stored reference
/// during normalization (#736).
///
/// The parser keeps `u` as a distinct [`Base::U`] (`b'U'`), separate from
/// [`Base::T`] (`b'T'`). `normalize_rna` compares the edit's inserted / delins
/// bases against the transcript reference, which is stored as DNA (`T`), so a
/// `Base::U` never matches and insertions / delins silently fail to canonicalize
/// or 3'-shift. Re-expressing the edit in DNA before normalization fixes this;
/// `RnaVariant`'s `Display` (`to_rna_string`) renders `T`→`u` again on output,
/// so the rewrite is display-neutral.
///
/// Returns `Some` with the rewritten edit iff at least one `U` was present, and
/// `None` for an already-DNA edit — mirroring
/// [`canonicalize_conversion_to_delins`]'s `None`-means-no-op convention so
/// callers can use it as an early-return probe. The walked shapes are a superset
/// of the literal-base carriers checked by the parser's `validate_no_u_in_dna`
/// (substitution, substitution-no-ref, deletion, duplication, insertion, delins,
/// repeat, multi-repeat), additionally covering the inversion, breakpoint-insertion,
/// and dupins shapes that `validate_no_u_in_dna` deliberately skips; every other
/// shape carries no literal base alphabet and is returned unchanged.
pub fn rna_uracil_to_thymine(edit: &NaEdit) -> Option<NaEdit> {
    use crate::hgvs::edit::{Base, RepeatUnit};

    fn map_base(b: Base) -> Base {
        if matches!(b, Base::U) {
            Base::T
        } else {
            b
        }
    }
    fn map_seq(s: &Sequence) -> Sequence {
        Sequence::new(s.bases().iter().copied().map(map_base).collect())
    }
    fn map_opt_seq(s: &Option<Sequence>) -> Option<Sequence> {
        s.as_ref().map(map_seq)
    }
    fn map_part(p: &InsertedPart) -> InsertedPart {
        match p {
            InsertedPart::Literal(s) => InsertedPart::Literal(map_seq(s)),
            InsertedPart::Repeat { base, count } => InsertedPart::Repeat {
                base: map_base(*base),
                count: count.clone(),
            },
            other => other.clone(),
        }
    }
    fn map_ins(ins: &InsertedSequence) -> InsertedSequence {
        match ins {
            InsertedSequence::Literal(s) => InsertedSequence::Literal(map_seq(s)),
            InsertedSequence::Repeat { base, count } => InsertedSequence::Repeat {
                base: map_base(*base),
                count: count.clone(),
            },
            InsertedSequence::SequenceRepeat { sequence, count } => {
                InsertedSequence::SequenceRepeat {
                    sequence: map_seq(sequence),
                    count: count.clone(),
                }
            }
            InsertedSequence::Complex(parts) => {
                InsertedSequence::Complex(parts.iter().map(map_part).collect())
            }
            other => other.clone(),
        }
    }

    let mapped = match edit {
        NaEdit::Substitution {
            reference,
            alternative,
        } => NaEdit::Substitution {
            reference: map_base(*reference),
            alternative: map_base(*alternative),
        },
        NaEdit::SubstitutionNoRef { alternative } => NaEdit::SubstitutionNoRef {
            alternative: map_base(*alternative),
        },
        NaEdit::Deletion { sequence, length } => NaEdit::Deletion {
            sequence: map_opt_seq(sequence),
            length: *length,
        },
        NaEdit::Duplication {
            sequence,
            length,
            uncertain_extent,
        } => NaEdit::Duplication {
            sequence: map_opt_seq(sequence),
            length: *length,
            uncertain_extent: uncertain_extent.clone(),
        },
        NaEdit::Insertion { sequence } => NaEdit::Insertion {
            sequence: map_ins(sequence),
        },
        NaEdit::BreakpointInsertion { sequence } => NaEdit::BreakpointInsertion {
            sequence: map_ins(sequence),
        },
        NaEdit::Delins {
            sequence,
            deleted,
            deleted_length,
            substitution_reference,
        } => NaEdit::Delins {
            sequence: map_ins(sequence),
            deleted: map_opt_seq(deleted),
            deleted_length: *deleted_length,
            substitution_reference: map_opt_seq(substitution_reference),
        },
        NaEdit::DupIns { sequence } => NaEdit::DupIns {
            sequence: map_ins(sequence),
        },
        NaEdit::Inversion { sequence, length } => NaEdit::Inversion {
            sequence: map_opt_seq(sequence),
            length: *length,
        },
        NaEdit::Repeat {
            sequence,
            count,
            additional_counts,
            trailing,
        } => NaEdit::Repeat {
            sequence: map_opt_seq(sequence),
            count: count.clone(),
            additional_counts: additional_counts.clone(),
            trailing: map_opt_seq(trailing),
        },
        NaEdit::MultiRepeat { units } => NaEdit::MultiRepeat {
            units: units
                .iter()
                .map(|u| RepeatUnit {
                    sequence: map_seq(&u.sequence),
                    count: u.count.clone(),
                })
                .collect(),
        },
        // Every other shape carries no literal base alphabet.
        other => other.clone(),
    };

    (mapped != *edit).then_some(mapped)
}

/// If `s` is wrapped in a single matching outer `[...]` pair — i.e.
/// the leading `[` at index 0 closes at the trailing `]` at
/// `s.len() - 1` under standard bracket-depth bookkeeping — return the
/// substring between the brackets. Otherwise return `s` unchanged.
///
/// The depth check guards against degenerate sources like `[A][B]`
/// (two adjacent groups) where stripping the first and last byte would
/// yield `A][B`, an invalid payload. It also leaves alone strings that
/// happen to start with `[` but never close back to depth-0 at the end
/// (e.g. `[A]B`).
fn strip_outer_brackets(s: &str) -> &str {
    let bytes = s.as_bytes();
    if bytes.len() < 2 || bytes[0] != b'[' || bytes[bytes.len() - 1] != b']' {
        return s;
    }
    let mut depth: usize = 0;
    for (i, &b) in bytes.iter().enumerate() {
        match b {
            b'[' => depth += 1,
            b']' => {
                if depth == 0 {
                    // Unbalanced: an inner `]` precedes any `[`.
                    return s;
                }
                depth -= 1;
                if depth == 0 && i != bytes.len() - 1 {
                    // Outer `[` closed before the end — the trailing
                    // `]` is a separate group, not the matching pair.
                    return s;
                }
            }
            _ => {}
        }
    }
    if depth == 0 {
        &s[1..s.len() - 1]
    } else {
        s
    }
}

/// Coordinate system context for `ins[...]` payload expansion.
///
/// Determines how inner position ranges (e.g. the `100_120` in
/// `g.X_Yins[100_120]`) are interpreted when looked up against the
/// reference. `g.` and `m.` ranges are direct genomic positions on the
/// outer accession; `n.` ranges are direct transcript positions; `c.`
/// ranges are CDS-relative and must be translated to transcript
/// coordinates via the transcript's `cds_start`.
///
/// This parameter is a deliberate addition beyond the
/// `(seq, accession, provider)` triple the helper signatures take
/// elsewhere — without it the helper cannot disambiguate
/// `c.X_Yins[100_120]` (CDS) from `n.X_Yins[100_120]` (transcript) for
/// the same accession. See issue #333.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InsCoordKind {
    /// Genomic (`g.`, `m.`, `o.`) or non-coding transcript (`n.`):
    /// positions are 1-based positions on the outer accession.
    ///
    /// **Note on `n.`:** non-coding transcript positions are 1-based
    /// transcript coordinates (not CDS-relative — there is no CDS to be
    /// relative to). The helper fetches bytes directly without any
    /// offset translation. The locked decision "`c./n.` range coords
    /// interpret as CDS" applies only to `c.`; `n.` is naturally
    /// transcript-coordinate by definition of the coord system.
    Direct,
    /// Coding DNA (`c.`): positions are 1-based CDS-relative; the helper
    /// translates to transcript coordinates via `cds_start` looked up
    /// from the provider.
    Cds,
    /// RNA (`r.`): position numbering follows the associated DNA reference —
    /// CDS-relative (== `c.`) for a coding transcript, transcript-relative
    /// (== `n.`) for a non-coding one. The coding/non-coding split needs the
    /// provider, so this is resolved to `Cds`/`Direct` at fetch time in
    /// `fetch_position_range_bases` via `Transcript::is_coding()`. Produced only
    /// by `parse_cross_reference` for an `r.` cross-reference payload. #773.
    Rna,
}

/// Fetch a 1-based inclusive byte range `[start_1based..=end_1based]`
/// from the provider, treating the bytes as the literal HGVS reference
/// sequence at those positions.
///
/// For `kind == InsCoordKind::Cds` the helper translates the input
/// positions to transcript positions via the transcript's `cds_start`
/// before fetching. Negative / `c.0` / UTR-3 CDS coordinates are not
/// supported here — those don't arise in well-formed `ins[A_B]` payloads
/// (the HGVS spec example `c.849_850ins858_895` uses positive CDS
/// integers only) and would route through `FerroError::Unsupported` at
/// the parser level.
fn fetch_position_range_bases<P: ReferenceProvider>(
    accession: &str,
    start_1based: u64,
    end_1based: u64,
    kind: InsCoordKind,
    provider: &P,
) -> Result<String, FerroError> {
    if start_1based == 0 || end_1based < start_1based {
        return Err(FerroError::InvalidCoordinates {
            msg: format!(
                "ins[{}_{}] range must satisfy 1 <= start <= end",
                start_1based, end_1based
            ),
        });
    }

    // Resolve an `r.` axis to its effective DNA-numbering frame. r. follows the
    // associated DNA numbering: CDS-relative (== `c.`) for a coding transcript,
    // transcript-relative (== `n.`) for a non-coding one. Deciding this needs
    // the provider, so it happens here rather than in `parse_cross_reference`.
    // We reuse the existing `Cds`/`Direct` arithmetic below so the cross-ref
    // path cannot drift from `normalize_rna`. No `U->T` step is applied: the
    // provider stores transcript bases as DNA (T), so fetched bases are already
    // DNA — translating them would be a no-op (the `g./c./n.` paths splice
    // provider output directly for the same reason). #773.
    let kind = match kind {
        InsCoordKind::Rna => {
            if provider.get_transcript(accession)?.is_coding() {
                InsCoordKind::Cds
            } else {
                InsCoordKind::Direct
            }
        }
        other => other,
    };

    // Translate CDS → transcript coordinates if needed. The helper only
    // accepts positive CDS positions; UTR-3 (`*N`) and 5' UTR (`-N`) are
    // not expressible in the `Complex(vec![PositionRange{..}])` shape
    // we land on from the parser (those carry sigils that the simple
    // u64 PositionRange cannot represent), so we conservatively reject
    // any CDS lookup whose translated tx position would underflow.
    let (tx_start, tx_end) = match kind {
        InsCoordKind::Direct => (start_1based, end_1based),
        InsCoordKind::Cds => {
            let transcript = provider.get_transcript(accession)?;
            let cds_start = transcript
                .cds_start
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!(
                        "transcript {} has no CDS start; cannot expand c. ins[A_B] payload",
                        accession
                    ),
                })?;
            // `cds_start` is documented as a 1-based transcript position
            // (`reference::transcript::Transcript::cds_start`). A zero
            // value violates that invariant and would let `c.1` translate
            // to `tx_start = 0`, then underflow on the `tx_start - 1`
            // 0-based conversion at the `get_sequence` call below. Reject
            // it explicitly with the same `ConversionError` shape as the
            // missing-CDS branch above.
            if cds_start == 0 {
                return Err(FerroError::ConversionError {
                    msg: format!(
                        "transcript {} has invalid CDS start 0 (expected 1-based); \
                             cannot expand c. ins[A_B] payload",
                        accession
                    ),
                });
            }
            // c.N (positive) maps to tx position cds_start + N - 1.
            // Guard against u64 overflow in the addition. The subtraction
            // is safe because `start_1based` and `end_1based` were already
            // validated as `>= 1` by the caller.
            let tx_start = cds_start.checked_add(start_1based - 1).ok_or_else(|| {
                FerroError::ConversionError {
                    msg: format!(
                        "c. ins[{}_{}] translation overflowed transcript coordinates \
                             on accession {}",
                        start_1based, end_1based, accession
                    ),
                }
            })?;
            let tx_end = cds_start.checked_add(end_1based - 1).ok_or_else(|| {
                FerroError::ConversionError {
                    msg: format!(
                        "c. ins[{}_{}] translation overflowed transcript coordinates \
                             on accession {}",
                        start_1based, end_1based, accession
                    ),
                }
            })?;
            (tx_start, tx_end)
        }
        // `Rna` is always resolved to `Cds` or `Direct` by the block above
        // before this match is reached — this arm is unreachable in practice.
        InsCoordKind::Rna => unreachable!("Rna kind must be resolved before this match"),
    };

    // `provider.get_sequence(id, s, e)` follows the convention used
    // throughout the normalizer: returns the slice `seq[s..e]` in
    // 0-based half-open indexing. HGVS 1-based position `p` lives at
    // index `p - 1`. So fetching HGVS positions `tx_start..=tx_end`
    // (inclusive) requires `(tx_start - 1, tx_end)`.
    provider.get_sequence(accession, tx_start - 1, tx_end)
}

/// A parsed cross-reference payload: the inner accession, how its positions
/// are interpreted, and the position range.
struct CrossReferenceParse {
    /// The payload's accession string (e.g. `NM_000088.3` or the compound
    /// `NG_012337.3(NM_012459.2)`), as written before the `:`.
    accession: String,
    /// How the positions are translated when fetching bases (`Cds` translates
    /// via the foreign transcript's `cds_start`; `Direct` fetches directly;
    /// `Rna` is resolved coding-aware to `Cds`/`Direct` in
    /// `fetch_position_range_bases` via `Transcript::is_coding()`).
    kind: InsCoordKind,
    /// One-based start position.
    start: u64,
    /// One-based end position (equal to `start` for a single-position payload).
    end: u64,
    /// Whether the payload range carried a trailing `inv`, meaning the fetched
    /// bases are inserted reverse-complemented. Only ever set for a two-part
    /// `A_Binv` range — a single position has no orientation to invert, matching
    /// the same-reference parse paths, which likewise build an
    /// `InsertedPart::PositionRangeInv` only from a range. #1184.
    inv: bool,
    /// Whether the positions index the transcript (`c.`/`n.`/`r.`) rather than a
    /// genomic / absolute frame (`g./m./o.`). When set, a genomic-context
    /// compound accession (`NG_(NM_)`) must be reduced to its transcript before
    /// lookup. `n.` shares the `Direct` *fetch* path with `g./m./o.` (it needs
    /// no `cds_start` translation), and `r.` resolves to `Direct` for a
    /// non-coding transcript, so `kind` alone cannot distinguish either — this
    /// flag carries that bit.
    transcript_relative: bool,
}

/// Parse a cross-reference string (e.g.
/// `"NC_000022.10:g.42536337_42536382"`) into its components for
/// provider lookup. Returns a `CrossReferenceParse` (the inner accession,
/// the `InsCoordKind` used to fetch its bases, the one-based `start`/`end`,
/// the `transcript_relative` flag, and the `inv` orientation) on success.
///
/// Two payload forms are accepted:
///   - **Axis-qualified** (`<ACC>:<axis>.<range>`): supports g./m./o./c./n./r.
///     axes with simple positive-integer positions or ranges (no offsets,
///     no `*N`, no `?`, no `pter`/`qter` markers). A range may carry a trailing
///     `inv` (`A_Binv`), which inverts the fetched bases (#1184).
///   - **Axis-less native-frame** (`<ACC>:<range>`, #759): the bare positions
///     index the accession's *native frame* — the axis an explicit
///     description would default to for that accession
///     (`Accession::inferred_variant_type`, resolved via
///     `axisless_native_frame`): `c.` for a coding transcript (`NM_`/`ENST`),
///     `n.` for a non-coding one (`NR_`/`XR_`), `g.` for a genomic reference
///     (`NC_`/`NG_`). So `NM_003002.4:100_102` == `NM_003002.4:c.100_102`.
///     Protein / unknown native frames are out of scope.
///
/// Out-of-scope shapes return `None` so the caller can surface a focused
/// `FerroError::UnsupportedVariant`. #422.
///
/// **Axis-set rationale:** the parser's `parse_reference_location`
/// (`src/hgvs/parser/edit.rs`) accepts `g/c/m/n/r/p/o` axes for
/// `Reference` strings. This helper accepts the subset that has a
/// well-defined position-range fetch on the inner accession:
///   - `g.`/`m.`/`o.`: genomic position-range fetch (`Direct`).
///   - `c.`: CDS-relative; translated via the FOREIGN transcript's
///     `cds_start` inside `fetch_position_range_bases`.
///   - `n.`: non-coding transcript; positions are already transcript-
///     relative (`Direct`).
///
/// **`r.` (RNA):** supported — numbering follows the associated DNA frame
/// (CDS-relative for coding transcripts, transcript-relative for non-coding),
/// resolved coding-aware in `fetch_position_range_bases` via
/// `Transcript::is_coding()`. No `U→T` step is needed: the provider stores
/// transcript bases as DNA. #773.
///
/// **Deferred axis (`p.`):** protein is structurally invalid as a DNA-
/// insertion payload, so it stays `UnsupportedVariant`. If
/// `parse_reference_location` ever broadens to additional axes, keep the
/// match arm below in sync.
fn parse_cross_reference(reference: &str) -> Option<CrossReferenceParse> {
    let (acc_part, after_colon) = reference.split_once(':')?;
    if acc_part.is_empty() || after_colon.is_empty() {
        return None;
    }
    let bytes = after_colon.as_bytes();
    // An axis-qualified payload is `<axis>.<range>` (`g.`/`c.`/…); anything else
    // is treated as an axis-less native-frame payload (`<ACC>:<range>`, #759).
    let axis_qualified = bytes.len() >= 2 && bytes[1] == b'.';
    let (kind, transcript_relative, range_str) = if axis_qualified {
        // c. requires CDS translation via the foreign transcript's cds_start;
        // g./m./o./n. are direct position-range fetches on the inner accession.
        // r. is resolved coding-aware in fetch_position_range_bases; p. is
        // deferred — see the doc-comment above. #773.
        let kind = match bytes[0] {
            b'g' | b'm' | b'o' | b'n' => InsCoordKind::Direct,
            b'c' => InsCoordKind::Cds,
            b'r' => InsCoordKind::Rna,
            _ => return None,
        };
        // `c.`, `n.`, and `r.` are transcript-relative (their positions index the
        // transcript), so a genomic-context compound payload accession
        // (`NG_(NM_)`) must be reduced to its transcript before lookup.
        // `g./m./o.` are genomic / absolute and fetch on the named accession
        // unchanged. (`n.` shares the `Direct` *fetch* path with `g./m./o.` — it
        // needs no cds_start translation — so the kind alone can't distinguish
        // it; this flag carries that out. `r.` likewise resolves to `Direct` for
        // a non-coding transcript, so it relies on the flag too.)
        let transcript_relative = matches!(bytes[0], b'c' | b'n' | b'r');
        (kind, transcript_relative, &after_colon[2..])
    } else {
        // Axis-less payload `<ACC>:<range>` (#759): the bare positions index the
        // accession's *native frame*, i.e. the axis an explicit description would
        // default to for that accession — `c.` for a coding transcript (`NM_`,
        // `ENST`), `n.` for a non-coding one (`NR_`/`XR_`), `g.` for a genomic
        // reference (`NC_`/`NG_`). So `NM_003002.4:100_102` == `NM_003002.4:c.100_102`.
        // Protein / unknown native frames are out of scope (return None).
        if !bytes[0].is_ascii_digit() {
            return None;
        }
        let (kind, transcript_relative) = axisless_native_frame(acc_part)?;
        (kind, transcript_relative, after_colon)
    };
    // A trailing `inv` inverts the payload (`A_Binv`). Strip it up front so the
    // position parse below still sees a bare range, and record the orientation
    // for `resolve_cross_reference_bases` to apply. `strip_suffix` anchors to the
    // end, so an embedded `inv` (`1inv_3`) is left alone and still rejected by
    // the digits-only checks. #1184.
    let (range_str, inv) = match range_str.strip_suffix("inv") {
        Some(stripped) => (stripped, true),
        None => (range_str, false),
    };
    // A payload may name a range (`A_B`) or a single position (`A`, a one-base
    // copy with start == end).
    let (start_str, end_str) = match range_str.split_once('_') {
        Some((s, e)) => (s, e),
        // A single position has no orientation, so `Ainv` is not a meaningful
        // shape — and it is unrepresentable in the same-reference form, where
        // both parse paths build a `PositionRangeInv` only from a two-part
        // range. Keep it out of scope rather than complement one base. #1184.
        None if inv => return None,
        None => (range_str, range_str),
    };
    if start_str.is_empty() || end_str.is_empty() {
        return None;
    }
    // Reject any decoration: offsets (`+`/`-`), UTR markers (`*`),
    // unknown markers (`?`), pter/qter/cen. Only positive-integer
    // ranges in scope per the issue's acceptance criteria — the one
    // permitted suffix, a trailing `inv`, is already stripped above.
    if !start_str.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    if !end_str.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let start: u64 = start_str.parse().ok()?;
    let end: u64 = end_str.parse().ok()?;
    Some(CrossReferenceParse {
        accession: acc_part.to_string(),
        kind,
        start,
        end,
        transcript_relative,
        inv,
    })
}

/// Resolve a cross-reference string to its literal sequence via the
/// provider. Routes through `fetch_position_range_bases` for the
/// position-translation logic (`Cds` arm translates via the foreign
/// transcript's `cds_start`; `Direct` arm fetches directly). The
/// outer variant's accession is not used here — the inner reference
/// carries its own accession.
///
/// Cross-chromosome translocations (e.g.
/// `NC_000002.12:g.X_Y delins[NC_000011.10:g.A_B]`) are supported as
/// long as both accessions are present in the provider; otherwise the
/// inner lookup returns `FerroError::ReferenceNotFound` citing the
/// missing accession. #422.
/// The shape-based (provider-independent) `UnsupportedVariant` error for
/// a cross-reference whose form is out of scope (`parse_cross_reference`
/// returned `None`). Shared by `resolve_cross_reference_bases` and the
/// `detect_deferred_part` pre-scan so the same message is surfaced
/// whether the unsupported `ExternalRef` is reached during expansion or
/// caught up front (which is what makes the deferral order-independent).
fn cross_reference_shape_error(reference: &str) -> FerroError {
    FerroError::UnsupportedVariant {
        variant_type: format!(
            "ins[{}] cross-reference shape not yet supported. \
             Expansion currently covers g./m./o./c./n./r. axes with simple \
             positive-integer positions or ranges, optionally with a trailing \
             `inv` on a range. Out-of-scope: p. (protein — structurally invalid \
             as DNA-insertion payload), offsets (+/-), UTR markers (*N), \
             unknown markers (?), and pter/qter/cen decorations",
            reference,
        ),
    }
}

fn resolve_cross_reference_bases<P: ReferenceProvider>(
    reference: &str,
    provider: &P,
) -> Result<String, FerroError> {
    let CrossReferenceParse {
        accession,
        kind,
        start,
        end,
        transcript_relative,
        inv,
    } = parse_cross_reference(reference).ok_or_else(|| cross_reference_shape_error(reference))?;
    // A transcript-relative payload (`c.`, `n.`, or `r.`) indexes the transcript, so a
    // genomic-context compound accession (`NG_(NM_)`) must be reduced to its
    // transcript (`NM_`/`NR_`) before lookup — otherwise the fetch fails on the
    // compound string (`c.` at the cds_start lookup, `n.` at the sequence
    // lookup). Genomic / absolute payloads (`g./m./o.`) fetch on the named
    // accession unchanged.
    let accession = if transcript_relative {
        transcript_accession_of(&accession).unwrap_or(accession)
    } else {
        accession
    };
    let bases = fetch_position_range_bases(&accession, start, end, kind, provider)?;
    // A trailing `inv` on the payload range inserts the segment
    // reverse-complemented, exactly as the same-reference
    // `InsertedPart::PositionRangeInv` / `InsertedSequence::PositionRangeInv`
    // arms do. Applied after the fetch so the axis-specific position
    // translation is unaffected by the orientation. #1184.
    if inv {
        Ok(crate::sequence::reverse_complement(&bases))
    } else {
        Ok(bases)
    }
}

/// Resolve the `(InsCoordKind, transcript_relative)` an axis-less cross-reference
/// payload should use, from the accession's *native frame* (#759).
///
/// An axis-less `<ACC>:<range>` is interpreted in the axis an explicit
/// description would default to for that accession (`Accession::inferred_variant_type`):
/// - `c.` (coding transcript `NM_`/`ENST`) → `Cds` (CDS-translated), transcript-relative;
/// - `n.` (non-coding `NR_`/`XR_`) → `Direct` fetch, transcript-relative;
/// - `g.`/`m.`/`o.` (genomic `NC_`/`NG_`/…) → `Direct` fetch, not transcript-relative.
///
/// `transcript_relative` lets a compound `NG_(NM_)` accession reduce to its inner
/// transcript before lookup (the parsed accession of a compound is the inner one).
/// Returns `None` when the accession does not parse or its native frame is out of
/// scope (protein / unknown), so the cross-reference is rejected as unsupported.
fn axisless_native_frame(accession: &str) -> Option<(InsCoordKind, bool)> {
    let (rest, parsed) = crate::hgvs::parser::accession::parse_accession(accession).ok()?;
    if !rest.is_empty() {
        return None;
    }
    match parsed.inferred_variant_type() {
        Some("c") => Some((InsCoordKind::Cds, true)),
        Some("n") => Some((InsCoordKind::Direct, true)),
        Some("g") | Some("m") | Some("o") => Some((InsCoordKind::Direct, false)),
        _ => None,
    }
}

/// Reduce an accession string to its transcript accession (`NM_`/`NR_`/`ENST`),
/// dropping any `NG_(…)`/`LRG_(…)` genomic-context wrapper. Returns `None` when
/// the string does not parse as a single accession, so the caller can fall back
/// to the original string unchanged.
fn transcript_accession_of(accession: &str) -> Option<String> {
    let (rest, parsed) = crate::hgvs::parser::accession::parse_accession(accession).ok()?;
    if !rest.is_empty() {
        return None;
    }
    Some(parsed.transcript_accession())
}

/// Append bases of an `InsertedPart` to `out` if the part is a
/// flat-resolvable form. Returns `Err(FerroError::UnsupportedVariant)`
/// for the remaining deferred shape (intronic-offset CDS range).
/// Returns `Err(...)` for genuine provider lookup failures.
///
/// Returns `Ok(false)` if the part is not flat-resolvable in a way the
/// caller should treat as "leave the variant alone" (currently never —
/// any unhandled part is an error). Returns `Ok(true)` if the part was
/// resolved and appended.
fn append_part_bases<P: ReferenceProvider>(
    part: &InsertedPart,
    accession: &str,
    kind: InsCoordKind,
    provider: &P,
    out: &mut String,
) -> Result<bool, FerroError> {
    match part {
        InsertedPart::Literal(seq) => {
            // `Sequence` Display prints uppercase IUPAC bytes; this is
            // exactly what the flat `insXXX` form wants.
            out.push_str(&seq.to_string());
            Ok(true)
        }
        InsertedPart::PositionRange { start, end } => {
            let bases = fetch_position_range_bases(accession, *start, *end, kind, provider)?;
            out.push_str(&bases);
            Ok(true)
        }
        InsertedPart::PositionRangeInv { start, end } => {
            let bases = fetch_position_range_bases(accession, *start, *end, kind, provider)?;
            out.push_str(&crate::sequence::reverse_complement(&bases));
            Ok(true)
        }
        InsertedPart::CdsPositionRange(range) => {
            // `CdsPositionRange` carries any CDS-coordinate range whose
            // positions are not plain positive integers — intronic
            // offsets (`244-8_249`, `c.244+8_249`) or 3'UTR markers
            // (`*100_*200`). All such shapes need offset/UTR-aware
            // coordinate translation the same-reference lookup machinery
            // does not provide. The HGVS spec does not pin canonical
            // expansion for any of these shapes, so we defer with a
            // grep-friendly tag. The only real-world usage we've seen
            // pairs an intronic-offset range with `N[k]` padding (e.g.
            // `c.249_250ins[N[2800];244-8_249]`).
            Err(FerroError::UnsupportedVariant {
                variant_type: format!(
                    "ins[{}] CDS-offset range (intronic or UTR-marker) is spec-undefined and not yet supported by ferro; see follow-up",
                    range
                ),
            })
        }
        InsertedPart::Repeat { .. } => {
            // Bracketed `[A[10];T]` carries a Repeat part whose flat
            // form is well-defined (`A` ten times). We could expand
            // these too, but the issue scope is limited to bases /
            // ranges / inv — fall through to a leave-alone signal.
            Ok(false)
        }
        InsertedPart::ExternalRef(reference) => {
            // Cross-reference inside a Complex bracket (e.g.
            // `ins[A;NC_000022.11:g.100_200]`). Resolve via the
            // shared cross-reference resolver, which routes through
            // `fetch_position_range_bases` for c.-axis CDS translation
            // and direct position fetch for g./m./o./n. axes. The
            // resolved bases are appended to the surrounding flat
            // payload — `[A; cross-ref]` becomes `[A; <bases>]` and
            // ultimately a single `Literal("A" + <bases>)`. #422.
            let bases = resolve_cross_reference_bases(reference, provider)?;
            out.push_str(&bases);
            Ok(true)
        }
    }
}

/// Pre-scan parts for shape-based deferred shapes that must be surfaced
/// regardless of part order:
///
/// - intronic-offset CDS range (`CdsPositionRange`), and
/// - an out-of-scope cross-reference (`ExternalRef` whose form
///   `parse_cross_reference` cannot handle — p. axis, offsets, UTR
///   markers, `pter`/`qter`, etc.).
///
/// Both are *shape* deferrals (provider-independent), unlike a resolvable
/// `ExternalRef` or same-reference `PositionRange` whose lookup needs the
/// provider. Catching them up front means a leading non-flatten part
/// (e.g. `Repeat`) that short-circuits `expand_complex_parts` to
/// `Ok(None)` can no longer mask a later unsupported cross-reference —
/// the deferral is order-independent. #422.
fn detect_deferred_part(parts: &[InsertedPart]) -> Option<FerroError> {
    for part in parts {
        match part {
            InsertedPart::CdsPositionRange(range) => {
                return Some(FerroError::UnsupportedVariant {
                    variant_type: format!(
                        "ins[{}] CDS-offset range (intronic or UTR-marker) is spec-undefined and not yet supported by ferro; see follow-up",
                        range
                    ),
                });
            }
            InsertedPart::ExternalRef(reference) if parse_cross_reference(reference).is_none() => {
                return Some(cross_reference_shape_error(reference));
            }
            _ => {}
        }
    }
    None
}

/// Expand the bracketed parts of a `Complex` insertion to a single flat
/// literal sequence, or return `Ok(None)` if at least one part is not
/// expandable in scope (e.g. a `Repeat { base, count }` member).
///
/// Deferred shapes — currently only `CdsPositionRange` (intronic-offset
/// or UTR-marker range) — are reported as
/// `Err(FerroError::UnsupportedVariant)` even when an earlier part is
/// non-flatten, because preserving a `Complex` that contains an
/// unresolvable CDS-offset range would otherwise silently pass through.
/// `ExternalRef` was deferred prior to #422; the cross-reference
/// resolver now handles it via `append_part_bases`.
pub(crate) fn expand_complex_parts<P: ReferenceProvider>(
    parts: &[InsertedPart],
    accession: &str,
    kind: InsCoordKind,
    provider: &P,
) -> Result<Option<InsertedSequence>, FerroError> {
    // Defer surfacing first so any deferred shape errors even when an
    // earlier part is non-flatten (e.g. `Repeat`).
    if let Some(err) = detect_deferred_part(parts) {
        return Err(err);
    }

    let mut out = String::new();
    for part in parts {
        let resolved = append_part_bases(part, accession, kind, provider, &mut out)?;
        if !resolved {
            // A part we don't flatten (e.g. Repeat) keeps the whole
            // Complex shape intact — partial flattening would silently
            // drop semantic content. Shape-unsupported parts (CDS-offset
            // ranges, out-of-scope cross-references) are already surfaced
            // by `detect_deferred_part` above, so returning here cannot
            // mask them regardless of part order.
            return Ok(None);
        }
    }
    let seq = Sequence::from_str(&out).map_err(|_| FerroError::ConversionError {
        msg: format!(
            "expanded ins[...] payload contains non-IUPAC base in {:?}",
            out
        ),
    })?;
    Ok(Some(InsertedSequence::Literal(seq)))
}

/// Canonicalize an `InsertedSequence` to its flat literal form when
/// possible. Returns:
///
/// - `Ok(None)` — the payload is either already flat (`Literal`) or not
///   in scope of this canonicalization (e.g. `Count`, `Range`, `Repeat`,
///   `Named`, `Uncertain`, `Empty`, and bare `PositionRange[Inv]` /
///   `Complex` shapes that contain at least one non-flat part).
/// - `Ok(Some(new_seq))` — a flat `Literal` rewrite for callers to
///   substitute in place of the input.
/// - `Err(FerroError::UnsupportedVariant)` — a deferred shape that
///   the spec permits but ferro does not yet expand. Two categories:
///   - **Intronic-offset / UTR-marker CDS range** (`CdsPositionRange`):
///     spec-undefined; error carries the grep tag `"CDS-offset range"`
///     and `"follow-up"`.
///   - **Out-of-scope cross-reference** (p. axis, or `pter`/`qter`
///     decoration, etc.): error carries `"cross-reference"` and
///     describes the supported axis set. Simple positive-integer
///     ranges on g./m./o./c./n. (#422) and r. (#773) are now
///     expanded and no longer surface as errors, as is such a range
///     carrying a trailing `inv` (#1184).
/// - `Err(FerroError::...)` — genuine provider lookup failure
///   (`ReferenceNotFound`, `InvalidCoordinates`, etc.).
pub fn expand_inserted_sequence<P: ReferenceProvider>(
    seq: &InsertedSequence,
    accession: &str,
    kind: InsCoordKind,
    provider: &P,
) -> Result<Option<InsertedSequence>, FerroError> {
    match seq {
        // `Reference("ACC:coord_a_b")`: bare cross-reference (no
        // surrounding Complex bracket). Resolve via the shared
        // cross-reference resolver. #422.
        //
        // Same-accession optimization: when the inner accession
        // matches the outer variant's accession, the resolver still
        // routes through `fetch_position_range_bases` — the helper
        // already handles same-accession lookups efficiently, so no
        // extra branch is needed here.
        InsertedSequence::Reference(reference) => {
            let bases = resolve_cross_reference_bases(reference, provider)?;
            let seq = Sequence::from_str(&bases).map_err(|_| FerroError::ConversionError {
                msg: format!(
                    "expanded ins[{}] cross-reference payload contains non-IUPAC base in {:?}",
                    reference, bases
                ),
            })?;
            Ok(Some(InsertedSequence::Literal(seq)))
        }

        // Same-reference position range — expand directly.
        InsertedSequence::PositionRange { start, end } => {
            let bases = fetch_position_range_bases(accession, *start, *end, kind, provider)?;
            let bases_seq =
                Sequence::from_str(&bases).map_err(|_| FerroError::ConversionError {
                    msg: format!(
                        "expanded ins[{}..={}] payload contains non-IUPAC base in {:?}",
                        start, end, bases
                    ),
                })?;
            Ok(Some(InsertedSequence::Literal(bases_seq)))
        }
        InsertedSequence::PositionRangeInv { start, end } => {
            let bases = fetch_position_range_bases(accession, *start, *end, kind, provider)?;
            let rc = crate::sequence::reverse_complement(&bases);
            let rc_seq = Sequence::from_str(&rc).map_err(|_| FerroError::ConversionError {
                msg: format!(
                    "expanded ins[{}..={}inv] payload contains non-IUPAC base in {:?}",
                    start, end, rc
                ),
            })?;
            Ok(Some(InsertedSequence::Literal(rc_seq)))
        }

        // The common bracketed shape — emitted by the parser for
        // anything `ins[...]` whose payload isn't a bare reference,
        // bare count, or bare external accession.
        InsertedSequence::Complex(parts) => expand_complex_parts(parts, accession, kind, provider),

        // Already flat or out-of-scope for this canonicalization.
        // `UncertainRangeInv` has uncertain (parenthesized) breakpoints that
        // are not sequenced, so there are no exact bases to fetch — leave it
        // as-is (DNA/complex.md). `SpecialPositionRange` has pter/qter/cen
        // endpoints that have no concrete coordinate, so it likewise cannot be
        // expanded to literal bases.
        InsertedSequence::Literal(_)
        | InsertedSequence::Count(_)
        | InsertedSequence::Range(_, _)
        | InsertedSequence::Repeat { .. }
        | InsertedSequence::SequenceRepeat { .. }
        | InsertedSequence::Named(_)
        | InsertedSequence::UncertainRangeInv { .. }
        | InsertedSequence::SpecialPositionRange { .. }
        | InsertedSequence::Uncertain
        | InsertedSequence::Empty => Ok(None),
    }
}

/// Rewrite the inserted-sequence payload of an `Insertion` / `Delins` /
/// `DupIns` edit to a flat literal. Returns `Ok(None)` when the edit
/// is not one of those three kinds or the payload is already flat.
///
/// Mirrors [`canonicalize_conversion_to_delins`]'s `None`-means-no-op
/// contract so the call site can chain helpers.
pub fn canonicalize_insertion_expand<P: ReferenceProvider>(
    edit: &NaEdit,
    accession: &str,
    kind: InsCoordKind,
    provider: &P,
) -> Result<Option<NaEdit>, FerroError> {
    match edit {
        NaEdit::Insertion { sequence } => {
            match expand_inserted_sequence(sequence, accession, kind, provider)? {
                Some(new_seq) => Ok(Some(NaEdit::Insertion { sequence: new_seq })),
                None => Ok(None),
            }
        }
        NaEdit::Delins {
            sequence,
            deleted,
            deleted_length,
            ..
        } => match expand_inserted_sequence(sequence, accession, kind, provider)? {
            Some(new_seq) => Ok(Some(NaEdit::Delins {
                sequence: new_seq,
                deleted: deleted.clone(),
                deleted_length: *deleted_length,
                substitution_reference: None,
            })),
            None => Ok(None),
        },
        NaEdit::DupIns { sequence } => {
            match expand_inserted_sequence(sequence, accession, kind, provider)? {
                Some(new_seq) => Ok(Some(NaEdit::DupIns { sequence: new_seq })),
                None => Ok(None),
            }
        }
        _ => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalize::config::ShuffleDirection;

    #[test]
    fn three_prime_align_odd_di_tract_shifts_phase() {
        // ref "CTGTGTT": tract TGTG at [1,5) (count 2, unit TG). idx 5 is 'T' ==
        // rotated[0], so it slides one base 3': GTGT at [2,6), unit GT.
        let r = b"CTGTGTT";
        let (s, e, u) = three_prime_align_tract(r, 1, 5, b"TG");
        assert_eq!((s, e, u.as_slice()), (2, 6, b"GT".as_slice()));
    }

    #[test]
    fn three_prime_align_clean_tract_is_noop() {
        // ref "CGCGCA": GCGC at [1,5) bounded by 'A' (!= 'G') -> no shift.
        let r = b"CGCGCA";
        let (s, e, u) = three_prime_align_tract(r, 1, 5, b"GC");
        assert_eq!((s, e, u.as_slice()), (1, 5, b"GC".as_slice()));
    }

    #[test]
    fn three_prime_align_homopolymer_shifts_to_run_end() {
        // ref "GAAAAG": single 'A' window [1,2) inside AAAA slides to [4,5).
        let r = b"GAAAAG";
        let (s, e, u) = three_prime_align_tract(r, 1, 2, b"A");
        assert_eq!((s, e, u.as_slice()), (4, 5, b"A".as_slice()));
    }

    #[test]
    fn three_prime_align_tri_tract_clean_noop() {
        // ref "ATTGTTGA": TTGTTG at [1,7) bounded by 'A' -> no shift.
        let r = b"ATTGTTGA";
        let (s, e, u) = three_prime_align_tract(r, 1, 7, b"TTG");
        assert_eq!((s, e, u.as_slice()), (1, 7, b"TTG".as_slice()));
    }

    #[test]
    fn insertion_to_repeat_three_prime_aligns_odd_di_tract() {
        // ref CTGTGTT: tract TGTGT at 0-based [1,6); inserting TGTGTG (3 copies of
        // TG) makes 5 copies. The HGVS 3'-rule anchors GT (not TG) at the 3'-most
        // 2-copy reference window. Hand-verified: unit GT, count 5, 1-based 3_6.
        let r = b"CTGTGTT";
        let (base, count, start, end, unit) =
            insertion_to_repeat(r, 4, b"TGTGTG", false, None, ShuffleDirection::ThreePrime)
                .expect("repeat");
        assert_eq!(base, b'G');
        assert_eq!(count, 5);
        assert_eq!((start, end), (3, 6));
        assert_eq!(unit, b"GT".to_vec());
    }

    /// #1389: the codon gate must be answered about the window the repeat is
    /// **emitted** over, not the raw tract `find_tandem_extent` discovered.
    ///
    /// Same fixture as the test above, with the CDS pinned to exactly the emitted
    /// window. The raw `TG` tract is 1-based `2..=5`, which is not CDS-resident
    /// (`2 < 3`), so gating on it switches the codon rule off — but the emitted
    /// window is `3..=6`, wholly inside the CDS, with a two-base unit that
    /// `repeated.md`'s codon rule forbids. The literal `ins` must win.
    #[test]
    fn insertion_to_repeat_gates_on_the_emitted_window_not_the_raw_tract() {
        let r = b"CTGTGTT";
        assert_eq!(
            insertion_to_repeat(
                r,
                4,
                b"TGTGTG",
                true,
                Some((3, 6)),
                ShuffleDirection::ThreePrime
            ),
            None,
            "the emitted window 3..=6 is wholly CDS-resident with a 2-base unit, \
             so the codon gate must block repeat notation"
        );
    }

    /// The 5' direction of the same fixture, which must NOT change: it emits `TG`
    /// over `2..=5`, genuinely outside the CDS, so gating the codon rule off there
    /// is correct. Pinned so the #1389 fix cannot over-reach into a window whose
    /// residency answer was already right.
    #[test]
    fn insertion_to_repeat_five_prime_window_outside_the_cds_still_repeats() {
        let r = b"CTGTGTT";
        let (_, count, start, end, unit) = insertion_to_repeat(
            r,
            4,
            b"TGTGTG",
            true,
            Some((3, 6)),
            ShuffleDirection::FivePrime,
        )
        .expect("the 5' window lies outside the CDS, so the codon gate stays off");
        assert_eq!((start, end), (2, 5));
        assert_eq!(unit, b"TG".to_vec());
        assert_eq!(count, 5);
    }

    #[test]
    fn normalize_repeat_single_anchor_off_phase_di_repeat_expands() {
        // #852 Option 2: a reverse-complemented projected di-repeat lands at a
        // single anchor off the tract's phase boundary. ref "CTGTGTT": the literal
        // unit "TG" does not tile at anchor idx 2 (which reads "GT"), so without
        // the rotation retry `normalize_repeat` bails to `Unchanged`. It must
        // instead find the GT tract and emit the entire-range 3'-form.
        let r = b"CTGTGTT";
        let got = normalize_repeat(r, 2, 2, b"TG", 5, false, ShuffleDirection::ThreePrime);
        assert_eq!(
            got,
            RepeatNormResult::Repeat {
                start: 3,
                end: 6,
                sequence: b"GT".to_vec(),
                count: 5,
            }
        );
    }

    #[test]
    fn normalize_repeat_single_anchor_tract_end_di_repeat_expands() {
        // #852 Option 2: minus-strand projection lands the single anchor at the
        // tract's 3' END. ref "TGCGCA": GCGC tract at idx 1..5; single anchor at
        // the last base (idx 4) — a boundary-anchored count undercounts, so the
        // offset search must recover the full 2-copy tract (the `c.3GC[5]` shape).
        let r = b"TGCGCA";
        let got = normalize_repeat(r, 4, 4, b"GC", 5, false, ShuffleDirection::ThreePrime);
        assert_eq!(
            got,
            RepeatNormResult::Repeat {
                start: 2,
                end: 5,
                sequence: b"GC".to_vec(),
                count: 5,
            }
        );
    }

    #[test]
    fn normalize_repeat_single_anchor_in_phase_di_repeat_unaffected() {
        // The already-phase-aligned unit takes the direct path and yields the same
        // result (phase-independence / idempotency).
        let r = b"CTGTGTT";
        let got = normalize_repeat(r, 2, 2, b"GT", 5, false, ShuffleDirection::ThreePrime);
        assert_eq!(
            got,
            RepeatNormResult::Repeat {
                start: 3,
                end: 6,
                sequence: b"GT".to_vec(),
                count: 5,
            }
        );
    }

    #[test]
    fn insertion_to_repeat_agrees_with_normalize_repeat_phase() {
        // #866 cross-check. The two ins→repeat entry points reach the canonical
        // repeat window by DIFFERENT tract-finders (`find_tandem_extent` anchored
        // at an insertion point vs `count_tandem_repeats` + off-phase rotation
        // search anchored at an existing position) but MUST agree on the final
        // 3'-aligned window `(start, end, rotated_unit)` for the same reference
        // tract — both share `three_prime_align_tract` (#852). That agreement is
        // what prevents the #852 non-idempotency from recurring; the count
        // differs by construction (insertion adds copies) and is not compared.
        //
        // Covers both off-phase di-repeat shapes #852 fixed: a mid-tract phase
        // rotation (`CTGTGTT`, unit reported as `GT` not `TG`) and a tract-end
        // boundary anchor (`TGCGCA`, anchor at the last base).
        struct Case {
            ref_seq: &'static [u8],
            ins_pos: u64,
            ins_seq: &'static [u8],
            nr_pos: usize,
            nr_end: usize,
            nr_unit: &'static [u8],
            nr_count: u64,
        }
        let cases = [
            Case {
                ref_seq: b"CTGTGTT",
                ins_pos: 4,
                ins_seq: b"TGTGTG",
                nr_pos: 1,
                nr_end: 5,
                nr_unit: b"TG",
                nr_count: 5,
            },
            // Case 2: in-phase rotation `CGCG` (not the out-of-phase `GCGC`,
            // which #882's gate now rejects). Same aligned window (2,5,"GC").
            Case {
                ref_seq: b"TGCGCA",
                ins_pos: 1,
                ins_seq: b"CGCG",
                nr_pos: 4,
                nr_end: 4,
                nr_unit: b"GC",
                nr_count: 5,
            },
            // Case 3 (#866): tri-nucleotide `CAG` tract ref[2..11); ins abuts its 3'
            // end. Both paths must agree on window (3,11,"CAG").
            Case {
                ref_seq: b"GGCAGCAGCAGGG",
                ins_pos: 10,
                ins_seq: b"CAGCAG",
                nr_pos: 2,
                nr_end: 10,
                nr_unit: b"CAG",
                nr_count: 5,
            },
            // Case 4 (#866): multi-copy homopolymer `A` tract ref[3..8). Window (4,8,"A").
            Case {
                ref_seq: b"GGGAAAAAGGG",
                ins_pos: 7,
                ins_seq: b"AA",
                nr_pos: 3,
                nr_end: 7,
                nr_unit: b"A",
                nr_count: 7,
            },
            // Case 5 (#866): single-reference-copy homopolymer (`ref_count == 1`) — the
            // one structurally distinct delegated branch, where the re-expressed range
            // collapses to `end_pos == pos` with `unit_len == 1`. `count_tandem_repeats`
            // still succeeds, so the len>=2 rotation fallback is unreachable. Window
            // (2,2,"A").
            Case {
                ref_seq: b"GAG",
                ins_pos: 1,
                ins_seq: b"AA",
                nr_pos: 1,
                nr_end: 1,
                nr_unit: b"A",
                nr_count: 3,
            },
        ];
        for c in cases {
            let (_b, _c, ins_start, ins_end, ins_unit) = insertion_to_repeat(
                c.ref_seq,
                c.ins_pos,
                c.ins_seq,
                false,
                None,
                ShuffleDirection::ThreePrime,
            )
            .unwrap_or_else(|| panic!("insertion_to_repeat None for {:?}", c.ref_seq));
            match normalize_repeat(
                c.ref_seq,
                c.nr_pos,
                c.nr_end,
                c.nr_unit,
                c.nr_count,
                false,
                ShuffleDirection::ThreePrime,
            ) {
                RepeatNormResult::Repeat {
                    start,
                    end,
                    sequence,
                    ..
                } => {
                    assert_eq!(
                        (ins_start, ins_end, ins_unit),
                        (start, end, sequence),
                        "ins→repeat phase disagreement on {:?}",
                        std::str::from_utf8(c.ref_seq).unwrap()
                    );
                }
                other => panic!("expected Repeat for {:?}, got {other:?}", c.ref_seq),
            }
        }

        // Case 6 (#866): coding codon-gate agreement. A non-codon-aligned (`A`, len 1)
        // expansion in c. context is NOT repeat notation — the two paths must agree that
        // it degrades: `insertion_to_repeat` declines (→ literal `ins`), and
        // `normalize_repeat` routes the expansion to `Insertion`. Kept out of the shared
        // loop above (which unwraps a `Repeat`) since this case is intentionally None.
        assert!(
            insertion_to_repeat(
                b"CAAAAAC",
                5,
                b"AA",
                true,
                None,
                ShuffleDirection::ThreePrime
            )
            .is_none(),
            "coding non-codon-aligned expansion must not be repeat notation"
        );
        match normalize_repeat(
            b"CAAAAAC",
            1,
            5,
            b"A",
            7,
            true,
            ShuffleDirection::ThreePrime,
        ) {
            RepeatNormResult::Insertion {
                start,
                end,
                sequence,
            } => assert_eq!((start, end, sequence), (6, 7, b"AA".to_vec())),
            other => panic!("expected codon-gated Insertion, got {other:?}"),
        }
    }

    #[test]
    fn test_needs_normalization() {
        use crate::hgvs::edit::Base;

        assert!(needs_normalization(&NaEdit::Deletion {
            sequence: None,
            length: None,
        }));
        assert!(needs_normalization(&NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: None,
        }));
        // Inversions need normalization so the complementary-outer-bases
        // shortening rule fires (RULE 10 in normalize_tests.rs).
        assert!(needs_normalization(&NaEdit::Inversion {
            sequence: None,
            length: None,
        }));
        // Real substitutions are now routed for validation (never shuffling)
        // — see #1052.
        assert!(needs_normalization(&NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        }));
        // Degenerate `A>A` substitutions must enter normalization so they
        // can be rewritten to identity (`=`).
        assert!(needs_normalization(&NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::A,
        }));
    }

    #[test]
    fn test_is_real_substitution() {
        use crate::hgvs::edit::Base;

        assert!(is_real_substitution(&NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        }));
        assert!(!is_real_substitution(&NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::A,
        }));
        assert!(!is_real_substitution(&NaEdit::Deletion {
            sequence: None,
            length: None,
        }));
    }

    #[test]
    fn test_is_uncertain_real_substitution() {
        use crate::hgvs::edit::Base;
        use crate::hgvs::Mu;

        let real = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        };
        // Only the `Uncertain`-wrapped real substitution qualifies.
        assert!(!is_uncertain_real_substitution(
            &real,
            &Mu::Certain(real.clone())
        ));
        assert!(is_uncertain_real_substitution(
            &real,
            &Mu::Uncertain(real.clone())
        ));
        assert!(!is_uncertain_real_substitution(&real, &Mu::Unknown));

        // A degenerate (`ref == alt`) substitution is never a *real* sub, even
        // wrapped in `Uncertain` — it must still route to the identity rewrite.
        let identity = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::A,
        };
        assert!(!is_uncertain_real_substitution(
            &identity,
            &Mu::Uncertain(identity.clone())
        ));
    }

    #[test]
    fn test_insertion_is_duplication() {
        // Sequence: ATGATGATG
        // Inserting ATG at position 3 (after first ATG) is a dup
        let ref_seq = b"ATGATGATG";
        assert!(insertion_is_duplication(ref_seq, 3, b"ATG"));

        // Inserting TGA at position 3 is not a dup
        assert!(!insertion_is_duplication(ref_seq, 3, b"TGA"));

        // Inserting ATG at position 6 (after second ATG) is a dup
        assert!(insertion_is_duplication(ref_seq, 6, b"ATG"));
    }

    #[test]
    fn test_insertion_is_duplication_pos_beyond_ref() {
        // Regression: inverted-range insertions (start > end) like
        // NC_000011.10:g.5238138_5153222insTATTT produce a pos far
        // beyond the reference sequence length.  Must return false,
        // not panic.
        let ref_seq = b"ATGATGATG";
        assert!(!insertion_is_duplication(ref_seq, 95, b"TATTT"));
        assert!(!insertion_is_duplication(b"", 95, b"TATTT"));
        assert!(!insertion_is_duplication(b"", 0, b"A"));
    }

    #[test]
    fn test_deletion_stays_deletion() {
        // Deletions should ALWAYS stay as deletions, never become dups
        let ref_seq = b"ATGATGATG";
        let del_edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        // Even if the deleted sequence matches preceding sequence,
        // it should remain a deletion
        assert_eq!(
            get_canonical_form(&del_edit, ref_seq, 3, 6),
            CanonicalForm::Deletion
        );
    }

    #[test]
    fn test_canonicalize_delins() {
        use crate::hgvs::edit::Base;
        let ref_seq = b"ACGTACGT";

        // Identity: insert == ref
        assert!(matches!(
            canonicalize_delins(ref_seq, 3, 4, b"T"),
            DelinsCanonical::Identity
        ));
        assert!(matches!(
            canonicalize_delins(ref_seq, 1, 4, b"CGT"),
            DelinsCanonical::Identity
        ));
        assert!(matches!(
            canonicalize_delins(ref_seq, 0, 8, b"ACGTACGT"),
            DelinsCanonical::Identity
        ));

        // Substitution: 1->1, ref != alt (no trimming needed; trimmed position
        // == input position)
        assert!(matches!(
            canonicalize_delins(ref_seq, 3, 4, b"A"),
            DelinsCanonical::Substitution {
                position: 3,
                reference: Base::T,
                alternative: Base::A,
            }
        ));
        // 1-base complement (ref=A, alt=T) -> Substitution, NEVER Inversion (HGVS spec gate)
        assert!(matches!(
            canonicalize_delins(b"A", 0, 1, b"T"),
            DelinsCanonical::Substitution {
                position: 0,
                reference: Base::A,
                alternative: Base::T,
            }
        ));

        // Substitution after shared-suffix trim: ref CTAG (idx 0..4) -> TTAG.
        // Suffix TAG is shared, leaving C->T at position 0.
        assert!(matches!(
            canonicalize_delins(b"CTAG", 0, 4, b"TTAG"),
            DelinsCanonical::Substitution {
                position: 0,
                reference: Base::C,
                alternative: Base::T,
            }
        ));
        // Substitution after shared-prefix trim: ref CTAG -> CTAA. Prefix CTA
        // is shared, leaving G->A at position 3.
        assert!(matches!(
            canonicalize_delins(b"CTAG", 0, 4, b"CTAA"),
            DelinsCanonical::Substitution {
                position: 3,
                reference: Base::G,
                alternative: Base::A,
            }
        ));

        // Substitution after both-side trim leaving exactly 1 residual base.
        // ref window NACGN[1..4] = ACG, ins = ATG. Prefix `A` and suffix `G`
        // shared; trimmed range idx 2..3 (C) -> T. The both-ends control flow
        // is distinct from the prefix-only and suffix-only paths above.
        assert!(matches!(
            canonicalize_delins(b"NACGN", 1, 4, b"ATG"),
            DelinsCanonical::Substitution {
                position: 2,
                reference: Base::C,
                alternative: Base::T,
            }
        ));

        // Mutalyzer A9 canonical example. ref window NGTAN[1..4] = GTA, ins = GTT.
        // Greedy prefix `GT` matches; the residual A>T lands at the last
        // position of the input window (the 3'-most index when no suffix
        // shares). This is the byte-level analogue of the original report's
        // `c.5948_5950delinsGTT over GTA -> c.5950A>T`.
        assert!(matches!(
            canonicalize_delins(b"NGTAN", 1, 4, b"GTT"),
            DelinsCanonical::Substitution {
                position: 3,
                reference: Base::A,
                alternative: Base::T,
            }
        ));

        // Homopolymer residual: ref AAAA, ins AAAT. Multiple positions are
        // semantically equivalent (any A could be the one that mutated to T);
        // greedy prefix-trim consumes 3 As, locking the residual to the
        // 3'-most position. This pins the tie-breaking convention.
        assert!(matches!(
            canonicalize_delins(b"AAAA", 0, 4, b"AAAT"),
            DelinsCanonical::Substitution {
                position: 3,
                reference: Base::A,
                alternative: Base::T,
            }
        ));

        // Pure deletion after both-side trim: ref ACGT (idx 0..4) -> AT.
        // Prefix A and suffix T shared; deleted CG remains at idx 1..3.
        assert!(matches!(
            canonicalize_delins(b"ACGT", 0, 4, b"AT"),
            DelinsCanonical::Deletion { start: 1, end: 3 }
        ));

        // Pure insertion after both-side trim: ref ACT (idx 0..3) -> ACGT.
        // Prefix AC and suffix T shared; G inserted between idx 2 and idx 2.
        assert!(matches!(
            canonicalize_delins(b"ACT", 0, 3, b"ACGT"),
            DelinsCanonical::Insertion {
                after_index: 2,
                ref sequence,
            } if sequence == b"G"
        ));

        // Inversion: insert == revcomp(ref), no shortening
        // ref = CTA (idx 0..3), revcomp = TAG; A and C not complements -> stays at (0,3)
        assert!(matches!(
            canonicalize_delins(b"CTA", 0, 3, b"TAG"),
            DelinsCanonical::Inversion { start: 0, end: 3 }
        ));

        // Inversion with shortening: ref CTATG, revcomp CATAG
        // outer C/G are complement -> shortens to inner TAT (0-idx 1..4)
        assert!(matches!(
            canonicalize_delins(b"CTATG", 0, 5, b"CATAG"),
            DelinsCanonical::Inversion { start: 1, end: 4 }
        ));

        // Inversion revealed by shared-affix trim: ref ACGAGT (idx 0..6) -> ACTCGT.
        // Full-range revcomp(ACGAGT) = ACTCGT — so this is also a full-range
        // inversion, but more importantly the trim path classifies it via the
        // trimmed-range revcomp check on GA -> TC.
        assert!(matches!(
            canonicalize_delins(b"ACGAGT", 0, 6, b"ACTCGT"),
            DelinsCanonical::Inversion { start, end } if start < end
        ));

        // Palindrome: ref ATAT, revcomp ATAT == ref -> Identity (NOT Inversion)
        assert!(matches!(
            canonicalize_delins(b"ATAT", 0, 4, b"ATAT"),
            DelinsCanonical::Identity
        ));

        // Duplication: 1->2 doubling. Range fields equal the input.
        assert!(matches!(
            canonicalize_delins(b"GATG", 1, 2, b"AA"),
            DelinsCanonical::Duplication { start: 1, end: 2 }
        ));
        // Duplication: 3->6 doubling
        assert!(matches!(
            canonicalize_delins(b"NATGCN", 1, 4, b"ATGATG"),
            DelinsCanonical::Duplication { start: 1, end: 4 }
        ));

        // KeepAsDelins: only complement, not reverse (AAGC -> TTCG). No
        // shared affix, so trim is a no-op and (start, end, sequence) match
        // the input.
        assert!(matches!(
            canonicalize_delins(b"AAGC", 0, 4, b"TTCG"),
            DelinsCanonical::KeepAsDelins {
                start: 0,
                end: 4,
                ref sequence,
            } if sequence == b"TTCG"
        ));
        // KeepAsDelins: only reverse, not revcomp (AAGC -> CGAA)
        assert!(matches!(
            canonicalize_delins(b"AAGC", 0, 4, b"CGAA"),
            DelinsCanonical::KeepAsDelins {
                start: 0,
                end: 4,
                ref sequence,
            } if sequence == b"CGAA"
        ));
        // KeepAsDelins: length doesn't fit any rule
        assert!(matches!(
            canonicalize_delins(b"AAGC", 0, 4, b"GGG"),
            DelinsCanonical::KeepAsDelins {
                start: 0,
                end: 4,
                ref sequence,
            } if sequence == b"GGG"
        ));

        // KeepAsDelins after trim: ref AGGCT (idx 0..5) -> AAACT. Prefix A
        // and suffix CT shared; trimmed range idx 1..3 (GG) -> AA. revcomp(GG)
        // = CC != AA, so not an inversion. 2->2 stays as a (trimmed) delins.
        assert!(matches!(
            canonicalize_delins(b"AGGCT", 0, 5, b"AAACT"),
            DelinsCanonical::KeepAsDelins {
                start: 1,
                end: 3,
                ref sequence,
            } if sequence == b"AA"
        ));

        // Bounds / degenerate: empty insert. Returns the (untrimmed) input
        // since classification requires a non-empty insert.
        assert!(matches!(
            canonicalize_delins(b"AAGC", 0, 1, b""),
            DelinsCanonical::KeepAsDelins {
                start: 0,
                end: 1,
                ref sequence,
            } if sequence.is_empty()
        ));
        // start >= end
        assert!(matches!(
            canonicalize_delins(b"AAGC", 2, 2, b"X"),
            DelinsCanonical::KeepAsDelins {
                start: 2,
                end: 2,
                ref sequence,
            } if sequence == b"X"
        ));
        // OOB end
        assert!(matches!(
            canonicalize_delins(b"AAGC", 3, 5, b"X"),
            DelinsCanonical::KeepAsDelins {
                start: 3,
                end: 5,
                ref sequence,
            } if sequence == b"X"
        ));
    }

    #[test]
    fn test_insertion_becomes_dup() {
        use crate::hgvs::edit::{InsertedSequence, Sequence};
        use std::str::FromStr;

        let ref_seq = b"ATGATGATG";
        let ins_edit = NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::from_str("ATG").unwrap()),
        };
        // Inserting ATG after first ATG (pos 3) should become dup
        assert_eq!(
            get_canonical_form(&ins_edit, ref_seq, 3, 3),
            CanonicalForm::Duplication
        );

        // Inserting TGA should stay as insertion
        let ins_edit2 = NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::from_str("TGA").unwrap()),
        };
        assert_eq!(
            get_canonical_form(&ins_edit2, ref_seq, 3, 3),
            CanonicalForm::Insertion
        );
    }

    // =========================================================================
    // HOMOPOLYMER REPEAT TESTS
    // =========================================================================

    #[test]
    fn test_find_homopolymer_at() {
        // Sequence: GGGAAAAAGGG (4 A's at positions 3-6, 0-indexed)
        let ref_seq = b"GGGAAAAAGGG";

        // Find homopolymer at position 4 (middle of A's)
        let result = find_homopolymer_at(ref_seq, 4);
        assert!(result.is_some());
        let analysis = result.unwrap();
        assert!(analysis.is_homopolymer);
        assert_eq!(analysis.base, Some(b'A'));
        assert_eq!(analysis.ref_start, 3);
        assert_eq!(analysis.ref_end, 8); // exclusive
        assert_eq!(analysis.ref_count, 5);

        // Position in G's at start
        let result = find_homopolymer_at(ref_seq, 1);
        assert!(result.is_some());
        let analysis = result.unwrap();
        assert_eq!(analysis.base, Some(b'G'));
        assert_eq!(analysis.ref_count, 3);

        // Single base (no repeat) - should return None
        let single_seq = b"ATGC";
        assert!(find_homopolymer_at(single_seq, 0).is_none());
    }

    #[test]
    fn test_insertion_to_repeat() {
        // Sequence: GGGAAAAAGGG (5 A's at positions 3-7, 0-indexed inclusive)
        let ref_seq = b"GGGAAAAAGGG";

        // `pos` is the 0-based index of the base immediately 5' of the
        // insertion point. pos=7 means insert between index 7 and 8 — i.e.,
        // immediately after the last A in the homopolymer. Inserting AA here
        // should become A[7] (5 ref + 2 inserted).
        let result =
            insertion_to_repeat(ref_seq, 7, b"AA", false, None, ShuffleDirection::ThreePrime);
        assert!(result.is_some());
        let (base, count, start, end, _unit) = result.unwrap();
        assert_eq!(base, b'A');
        assert_eq!(count, 7); // 5 original + 2 inserted
        assert_eq!(start, 4); // 1-indexed start of A region in reference
        assert_eq!(end, 8); // 1-indexed end of A region in reference (per HGVS, positions refer to reference tract)

        // Single-copy inserts (added_copies < 2) are duplications, not repeats.
        let result =
            insertion_to_repeat(ref_seq, 7, b"A", false, None, ShuffleDirection::ThreePrime);
        assert!(result.is_none());

        // Inserting T (non-matching) should return None
        let result =
            insertion_to_repeat(ref_seq, 7, b"T", false, None, ShuffleDirection::ThreePrime);
        assert!(result.is_none());

        // Inserting mixed bases should return None
        let result =
            insertion_to_repeat(ref_seq, 7, b"AT", false, None, ShuffleDirection::ThreePrime);
        assert!(result.is_none());

        // Multi-base tandem unit: insert ACAC into ACAC tract → AC[4].
        // ref_seq has ACAC at indices 0-3. Insert ACAC just before index 4
        // (pos=3, between A at index 3 and base at index 4). Expected:
        // 2 ref AC units + 2 inserted AC units = AC[4] at ref indices 0..3.
        let ref_ac = b"ACACGGG";
        let result = insertion_to_repeat(
            ref_ac,
            3,
            b"ACAC",
            false,
            None,
            ShuffleDirection::ThreePrime,
        );
        assert!(result.is_some());
        let (base, count, start, end, unit) = result.unwrap();
        assert_eq!(base, b'A');
        assert_eq!(count, 4);
        assert_eq!(start, 1); // 1-indexed
        assert_eq!(end, 4); // 1-indexed
        assert_eq!(unit, b"AC");
    }

    #[test]
    fn test_duplication_to_repeat() {
        // Sequence: GGGAAAAAGGG (5 A's at positions 3-7, 0-indexed)
        let ref_seq = b"GGGAAAAAGGG";

        // Duplicating 2 A's (positions 3-5, 0-indexed) should become A[7]
        let result = duplication_to_repeat(ref_seq, 3, 5, false);
        assert!(result.is_some());
        match result.unwrap() {
            DupToRepeatResult::Homopolymer {
                base, count, start, ..
            } => {
                assert_eq!(base, b'A');
                assert_eq!(count, 7); // 5 original + 2 duplicated
                assert_eq!(start, 4); // 1-indexed start
            }
            _ => panic!("Expected Homopolymer result"),
        }

        // Duplicating non-homopolymer region should return None (if not a tandem repeat)
        // ATGCXYZ has no repeats, so duplicating ATG should return None
        let non_repeat_seq = b"ATGCXYZ";
        let result = duplication_to_repeat(non_repeat_seq, 0, 3, false);
        assert!(result.is_none());
    }

    #[test]
    fn test_duplication_to_tandem_repeat() {
        // Sequence with GCA repeats
        // String: AAAAAGCAGCAGCAGCAGCAGCAGCAGCAAAAA (33 chars)
        // 5 A's + 8 GCAs (24 chars) + 4 A's = 33 chars
        // 8 GCA repeats at 0-indexed positions 5-28 (exclusive 29)
        let ref_seq = b"AAAAAGCAGCAGCAGCAGCAGCAGCAGCAAAAA";

        // Duplicating single GCA (one copy) should NOT become repeat notation
        // It should stay as a simple dup per HGVS rules
        let result = duplication_to_repeat(ref_seq, 5, 8, false);
        assert!(result.is_none(), "Single-copy dup should not become repeat");

        // Duplicating GCAGCA (2 copies) SHOULD become repeat notation
        let result = duplication_to_repeat(ref_seq, 5, 11, false);
        assert!(result.is_some());
        match result.unwrap() {
            DupToRepeatResult::TandemRepeat {
                unit,
                count,
                start,
                end,
            } => {
                assert_eq!(unit, b"GCA");
                assert_eq!(count, 10); // 8 original + 2 duplicated
                assert_eq!(start, 6); // 1-indexed
                assert_eq!(end, 29); // 1-indexed end of 8 GCAs (0-indexed 28 + 1)
            }
            _ => panic!("Expected TandemRepeat result"),
        }
    }

    // =========================================================================
    // TANDEM REPEAT COUNTING TESTS
    // =========================================================================

    #[test]
    fn test_count_tandem_repeats_basic() {
        // Sequence: GGGCATCATCATGGG (3 CAT repeats at positions 3-11)
        let ref_seq = b"GGGCATCATCATGGG";

        let result = count_tandem_repeats(ref_seq, 3, b"CAT");
        assert!(result.is_some());
        let (count, start, end) = result.unwrap();
        assert_eq!(count, 3);
        assert_eq!(start, 3);
        assert_eq!(end, 12); // exclusive

        // Try from middle of the repeat
        let result = count_tandem_repeats(ref_seq, 6, b"CAT");
        assert!(result.is_some());
        let (count, start, end) = result.unwrap();
        assert_eq!(count, 3);
        assert_eq!(start, 3);
        assert_eq!(end, 12);
    }

    #[test]
    fn test_count_tandem_repeats_single_base() {
        // Sequence: GGGAAAAAAGGG (6 A's)
        let ref_seq = b"GGGAAAAAAGGG";

        let result = count_tandem_repeats(ref_seq, 5, b"A");
        assert!(result.is_some());
        let (count, start, end) = result.unwrap();
        assert_eq!(count, 6);
        assert_eq!(start, 3);
        assert_eq!(end, 9);
    }

    #[test]
    fn test_count_tandem_repeats_no_match() {
        let ref_seq = b"GGGAAAAAAGGG";
        // XYZ doesn't appear in the sequence
        let result = count_tandem_repeats(ref_seq, 5, b"XYZ");
        assert!(result.is_none());
    }

    // =========================================================================
    // REPEAT NORMALIZATION TESTS
    // =========================================================================

    #[test]
    fn test_normalize_repeat_to_deletion() {
        // Sequence: GGGCATCATCATCATGGG (4 CAT repeats at positions 3-14, 0-indexed)
        // Specifying CAT[1]: ref_count=4, specified=1, k=3 >= 2, post=1 >= 1 → B2 → Repeat
        let ref_seq = b"GGGCATCATCATCATGGG";

        let result = normalize_repeat(
            ref_seq,
            3,
            3,
            b"CAT",
            1,
            false,
            ShuffleDirection::ThreePrime,
        );
        match result {
            RepeatNormResult::Repeat {
                sequence, count, ..
            } => {
                assert_eq!(sequence, b"CAT");
                assert_eq!(count, 1, "Should reflect specified count of 1");
            }
            _ => panic!("Expected Repeat (B2), got {:?}", result),
        }
    }

    #[test]
    fn test_normalize_repeat_to_duplication() {
        // Sequence: GGGCATCATGGG (2 CAT repeats at positions 3-8)
        // Specifying CAT[3] (ref is 2, so 2+1=3) should become duplication
        let ref_seq = b"GGGCATCATGGG";

        let result = normalize_repeat(
            ref_seq,
            3,
            3,
            b"CAT",
            3,
            false,
            ShuffleDirection::ThreePrime,
        );
        match result {
            RepeatNormResult::Duplication {
                start,
                end,
                sequence,
            } => {
                assert_eq!(sequence, b"CAT");
                assert_eq!(end - start + 1, 3, "Should duplicate 3 bases (1 CAT)");
            }
            _ => panic!("Expected Duplication, got {:?}", result),
        }
    }

    #[test]
    fn test_normalize_repeat_stays_repeat() {
        // Sequence: GGGCATCATGGG (2 CAT repeats)
        // Specifying CAT[5] (ref is 2, 5 > 2+1) should stay as repeat
        let ref_seq = b"GGGCATCATGGG";

        let result = normalize_repeat(
            ref_seq,
            3,
            3,
            b"CAT",
            5,
            false,
            ShuffleDirection::ThreePrime,
        );
        match result {
            RepeatNormResult::Repeat {
                count, sequence, ..
            } => {
                assert_eq!(sequence, b"CAT");
                assert_eq!(count, 5);
            }
            _ => panic!("Expected Repeat, got {:?}", result),
        }
    }

    #[test]
    fn test_normalize_repeat_unchanged() {
        // Sequence: GGGCATCATGGG (2 CAT repeats)
        // Specifying CAT[2] (same as ref) should be unchanged
        let ref_seq = b"GGGCATCATGGG";

        let result = normalize_repeat(
            ref_seq,
            3,
            3,
            b"CAT",
            2,
            false,
            ShuffleDirection::ThreePrime,
        );
        assert!(matches!(result, RepeatNormResult::Unchanged));
    }

    #[test]
    fn test_normalize_repeat_empty_unit_is_unchanged() {
        // Pre-refactor, `normalize_repeat` delegated emptiness handling to
        // `count_tandem_repeats`, which returns `None` on an empty unit, yielding
        // `Unchanged`. The B2 canonicalization step now calls
        // `smallest_repeat_unit` first, which returns the empty slice unchanged
        // for an empty input — leaving a divide-by-zero on `repeat_unit.len() /
        // canonical_unit.len()` unless we guard up front. This test pins the
        // pre-refactor contract.
        let ref_seq = b"GGGCATCATGGG";
        let result = normalize_repeat(ref_seq, 3, 3, b"", 1, false, ShuffleDirection::ThreePrime);
        assert!(matches!(result, RepeatNormResult::Unchanged));
    }

    #[test]
    fn test_normalize_repeat_canonicalizes_non_minimal_unit() {
        // Caller-spelled `ATAT[1]` over an `AT[4]` tract is a contraction:
        // canonical unit AT, ref_count=4, specified_canonical=2, k=2 → B2
        // emits `AT[2]` at the canonical tract span. Without canonicalization
        // this would fall to a 1-unit (ATAT) reduction → deletion (k<2).
        let ref_seq = b"GGGATATATATGGG"; // AT-tract at indices 3..11 (4 AT)

        let result = normalize_repeat(
            ref_seq,
            3,
            3,
            b"ATAT",
            1,
            false,
            ShuffleDirection::ThreePrime,
        );
        match result {
            RepeatNormResult::Repeat {
                start,
                end,
                sequence,
                count,
            } => {
                assert_eq!(sequence, b"AT", "Should emit canonical (smallest) unit");
                assert_eq!(count, 2, "Specified ATAT[1] = 2 canonical AT copies");
                assert_eq!(start, 4, "Canonical tract starts at HGVS pos 4");
                assert_eq!(end, 11, "Canonical tract ends at HGVS pos 11");
            }
            _ => panic!("Expected canonical AT[2] Repeat, got {:?}", result),
        }
    }

    // =========================================================================
    // INVERSION SHORTENING TESTS
    // =========================================================================

    #[test]
    fn test_complement() {
        assert_eq!(complement_base(b'A'), Some(b'T'));
        assert_eq!(complement_base(b'T'), Some(b'A'));
        assert_eq!(complement_base(b'G'), Some(b'C'));
        assert_eq!(complement_base(b'C'), Some(b'G'));
        assert_eq!(complement_base(b'N'), Some(b'N')); // N stays N
                                                       // The contract that replaced the old `_ => base` fallback (#1318): an
                                                       // unmodelled byte has no complement to report, so callers cannot mistake
                                                       // "no answer" for "self-complementary".
        assert_eq!(complement_base(b'X'), None);
    }

    #[test]
    fn test_shorten_inversion_basic() {
        // Sequence: ATGCAT (positions 0-5)
        // A(0) is complement of T(5) - cancel
        // T(1) is complement of A(4) - cancel
        // G(2) is complement of C(3) - cancel
        // All bases cancel -> becomes identity
        let seq = b"ATGCAT";
        let result = shorten_inversion(seq, 0, 6);
        assert!(
            result.is_none(),
            "Fully complementary inversion should become identity"
        );
    }

    #[test]
    fn test_shorten_inversion_partial() {
        // Sequence: ATGGAT (positions 0-5)
        // A(0) is complement of T(5) - cancel
        // T(1) is complement of A(4) - cancel
        // G(2) is NOT complement of G(3) - stop
        // Result: positions 2-4
        let seq = b"ATGGAT";
        let result = shorten_inversion(seq, 0, 6);
        assert!(result.is_some());
        let (s, e) = result.unwrap();
        assert_eq!(s, 2);
        assert_eq!(e, 4);
    }

    #[test]
    fn test_shorten_inversion_no_change() {
        // Sequence: GGCC (positions 0-3)
        // G(0) is complement of C(3) - cancel
        // G(1) is complement of C(2) - cancel
        // All cancel -> identity
        let seq = b"GGCC";
        let result = shorten_inversion(seq, 0, 4);
        assert!(result.is_none());

        // Sequence: GATT (positions 0-3)
        // G(0) is NOT complement of T(3) - no cancellation
        let seq2 = b"GATT";
        let result2 = shorten_inversion(seq2, 0, 4);
        assert!(result2.is_some());
        let (s, e) = result2.unwrap();
        assert_eq!(s, 0);
        assert_eq!(e, 4);
    }

    #[test]
    fn complement_models_the_iupac_ambiguity_codes() {
        // Ambiguity codes complement by complementing their member set.
        for (base, expected) in [
            (b'R', b'Y'), // A/G   -> C/T
            (b'Y', b'R'), // C/T   -> A/G
            (b'K', b'M'), // G/T   -> A/C
            (b'M', b'K'), // A/C   -> G/T
            (b'B', b'V'), // C/G/T -> A/C/G
            (b'V', b'B'), // A/C/G -> C/G/T
            (b'D', b'H'), // A/G/T -> A/C/T
            (b'H', b'D'), // A/C/T -> A/G/T
        ] {
            assert_eq!(
                complement_base(base),
                Some(expected),
                "complement({}) should be {}",
                base as char,
                expected as char
            );
            // Lowercase complements to the same uppercase symbol as A/C/G/T/U do.
            assert_eq!(complement_base(base.to_ascii_lowercase()), Some(expected));
        }

        // Only W, S and N are genuinely self-complementary.
        for base in *b"WSN" {
            assert_eq!(
                complement_base(base),
                Some(base),
                "{} is self-complementary",
                base as char
            );
        }
    }

    #[test]
    fn shorten_inversion_keeps_an_ambiguous_one_base_residue() {
        // "ART": A(0) cancels against T(2), leaving R at the centre. R is not
        // self-complementary — its reverse complement is Y — so the residue must
        // survive as a one-base range for the caller to describe. Treating it as
        // identity would silently discard the variant, which is #1249 all over
        // again for the IUPAC codes.
        let result = shorten_inversion(b"ART", 0, 3);
        assert_eq!(result, Some((1, 2)));
        assert_eq!(
            complementary_substitution(b'R'),
            Some((crate::hgvs::edit::Base::R, crate::hgvs::edit::Base::Y)),
            "a lone R inverts to Y"
        );
    }

    #[test]
    fn shorten_inversion_never_returns_a_backwards_range() {
        // A self-complementary centre base used to cancel against *itself* —
        // `ref_seq[s]` and `ref_seq[e - 1]` are the same byte once the span is
        // one base wide — driving `e` below `s` and yielding `Some((2, 1))`.
        // These spans are identity and must report `None`.
        for seq in [b"ANT", b"AWT", b"AST"] {
            assert_eq!(
                shorten_inversion(seq, 0, 3),
                None,
                "{} reduces to a self-complementary centre, i.e. identity",
                std::str::from_utf8(seq).unwrap()
            );
        }

        // The invariant, stated directly: whatever the sequence, a returned
        // range is non-empty and correctly ordered.
        for seq in [
            b"ART".as_slice(),
            b"ANT".as_slice(),
            b"RR".as_slice(),
            b"RY".as_slice(),
            b"NNN".as_slice(),
            b"ATGCAT".as_slice(),
        ] {
            if let Some((s, e)) = shorten_inversion(seq, 0, seq.len()) {
                assert!(
                    s < e,
                    "{} produced a backwards range ({s}, {e})",
                    std::str::from_utf8(seq).unwrap()
                );
            }
        }
    }

    #[test]
    fn shorten_inversion_cancels_ambiguity_codes_only_when_they_pair() {
        // revcomp("RY") == "RY", so the pair genuinely cancels to identity.
        assert_eq!(shorten_inversion(b"RY", 0, 2), None);

        // revcomp("RR") == "YY", so nothing cancels and the inversion stands.
        // Before the IUPAC codes were modelled, `complement(R) == R` made this
        // look like a complementary pair and collapsed it to identity.
        assert_eq!(shorten_inversion(b"RR", 0, 2), Some((0, 2)));
    }

    #[test]
    fn shorten_inversion_treats_soft_masked_bases_like_uppercase() {
        // `complement` always answers uppercase, so every `complement(x) == y`
        // test used to fail on a soft-masked reference: `complement(b'a')` is
        // `b'T'`, never `b't'`. That made lowercase spans behave differently
        // from the identical uppercase span at both comparison sites — the
        // outer-pair loop and the lone-centre check.
        for (lower, upper) in [
            (b"at".as_slice(), b"AT".as_slice()), // outer pair cancels to nothing
            (b"awt".as_slice(), b"AWT".as_slice()), // cancels to a self-complementary centre
            (b"ast".as_slice(), b"AST".as_slice()),
            (b"ant".as_slice(), b"ANT".as_slice()),
            (b"art".as_slice(), b"ART".as_slice()), // centre survives: revcomp(r) is y
            (b"rr".as_slice(), b"RR".as_slice()),   // nothing cancels
            (b"atgcat".as_slice(), b"ATGCAT".as_slice()),
        ] {
            assert_eq!(
                shorten_inversion(lower, 0, lower.len()),
                shorten_inversion(upper, 0, upper.len()),
                "{} and {} must shorten identically",
                std::str::from_utf8(lower).unwrap(),
                std::str::from_utf8(upper).unwrap()
            );
        }

        // The two headline cases, spelled out rather than left to the parity
        // check above: a soft-masked self-complementary centre is identity, and
        // a soft-masked non-self-complementary centre still survives.
        assert_eq!(shorten_inversion(b"awt", 0, 3), None);
        assert_eq!(shorten_inversion(b"art", 0, 3), Some((1, 2)));
    }

    /// The delins path must answer identically for a soft-masked reference and
    /// the same sequence uppercased (#1318).
    ///
    /// #1250 fixed the *inversion-shortening* comparisons and was deliberately
    /// scope-limited; the delins→inv typing path kept its raw byte tests. On a
    /// soft-masked reference `deleted` arrives lowercase while `inserted` comes
    /// from the author's uppercase description, so `is_revcomp` and the
    /// shared-affix/identity comparisons all saw a mismatch that is not one.
    ///
    /// Parity between the two cases is the assertion, rather than a pinned
    /// output shape: whatever the canonical form of a given delins is, the
    /// reference's masking must not change it.
    #[test]
    fn canonicalize_delins_treats_soft_masked_bases_like_uppercase() {
        // (reference, start, end, inserted) — each exercises a different one of
        // the case-sensitive comparison sites.
        let cases: &[(&[u8], usize, usize, &[u8])] = &[
            // is_revcomp: `atgcat` reverse-complements to `ATGCAT`, so this is
            // an inversion — but only if the comparison folds case.
            (b"atgcat", 0, 6, b"ATGCAT"),
            // A longer revcomp that survives outer-pair shortening.
            (b"aggcct", 0, 6, b"AGGCCT"),
            // The Duplication branch (`inserted_seq.len() == 2 * deleted.len()`),
            // which is checked before trimming and has its own pair of
            // comparisons. Nothing above reaches it: every other case here has
            // an insert the same length as the deleted range, so the branch's
            // length gate is false and its case-folding is never exercised.
            (b"acgt", 0, 4, b"ACGTACGT"),
            // shared_affix_lengths: a shared prefix/suffix that differs only in
            // case must still be trimmed.
            (b"acgtacgt", 0, 8, b"ACGTTTGT"),
            // The identity-run detector inside the decomposition loop.
            (b"acgtacgt", 0, 8, b"ACGAACGT"),
            // Whole-range identity: insert equals the deleted reference, modulo
            // case.
            (b"acgtacgt", 0, 8, b"ACGTACGT"),
            // A plain substitution reached through the trimming path.
            (b"acgt", 0, 4, b"ACTT"),
        ];
        for (lower, start, end, inserted) in cases {
            let upper: Vec<u8> = lower.to_ascii_uppercase();
            assert_eq!(
                canonicalize_delins(lower, *start, *end, inserted),
                canonicalize_delins(&upper, *start, *end, inserted),
                "`{}` and `{}` with insert `{}` must canonicalize identically",
                std::str::from_utf8(lower).unwrap(),
                std::str::from_utf8(&upper).unwrap(),
                std::str::from_utf8(inserted).unwrap(),
            );
        }
    }

    /// `decompose_delins` must hold the same soft-mask parity, and did not.
    ///
    /// The identity test at the top of its loop was case-folded while the
    /// mismatch-run scan below it was not — the same rule written twice, fixed
    /// in one copy. On a soft-masked reference the run then swallowed positions
    /// that the identity test calls identities, so `has_identity` never became
    /// true and the function's final gate returned `None` for a span whose
    /// uppercase twin decomposes. Measured before the fix:
    /// `decompose_delins(b"acgt", 0, 4, b"TCGA")` was `None` against
    /// `Some(4)` for `b"ACGT"`.
    #[test]
    fn decompose_delins_treats_soft_masked_bases_like_uppercase() {
        let cases: &[(&[u8], &[u8])] = &[
            // An identity sitting *inside* what the run scan saw as a mismatch
            // run — the shape the unfixed copy swallowed.
            (b"acgt", b"TCGA"),
            (b"acgtac", b"TCGTAC"),
            // A palindrome, so the inversion branch is reached on both sides.
            (b"atgcat", b"TTGCAA"),
        ];
        for (lower, inserted) in cases {
            let upper: Vec<u8> = lower.to_ascii_uppercase();
            assert_eq!(
                decompose_delins(lower, 0, lower.len(), inserted),
                decompose_delins(&upper, 0, upper.len(), inserted),
                "`{}` and `{}` with insert `{}` must decompose identically",
                std::str::from_utf8(lower).unwrap(),
                std::str::from_utf8(&upper).unwrap(),
                std::str::from_utf8(inserted).unwrap(),
            );
        }
    }

    /// The same parity for the insertion→duplication test, which compares a
    /// reference slice against the author's inserted sequence and so carried
    /// the identical defect: `insertion_is_duplication(b"acgtacgt", 4, b"ACGT")`
    /// answered `false` where the uppercase reference answered `true`, silently
    /// leaving a duplication spelled as an insertion.
    #[test]
    fn insertion_is_duplication_treats_soft_masked_bases_like_uppercase() {
        for (lower, pos, inserted) in [
            (&b"acgtacgt"[..], 4u64, &b"ACGT"[..]), // 3' duplication
            (&b"aaggtt"[..], 2u64, &b"AA"[..]),     // 5' duplication
        ] {
            let upper: Vec<u8> = lower.to_ascii_uppercase();
            assert_eq!(
                insertion_is_duplication(lower, pos, inserted),
                insertion_is_duplication(&upper, pos, inserted),
                "`{}` at {pos} with insert `{}` must classify identically",
                std::str::from_utf8(lower).unwrap(),
                std::str::from_utf8(inserted).unwrap(),
            );
        }
    }

    /// Unifying on `complement_base` changed what happens to a byte outside the
    /// modelled alphabet, and that change is deliberate — pinned here so it is a
    /// decision rather than a side effect.
    ///
    /// `complement` returned such a byte unchanged, so it read as
    /// self-complementary and an inversion over it collapsed to identity — the
    /// #1249 class of silent loss. `complement_base` answers `None`, so the span
    /// survives as a residue and the caller leaves the inversion as authored
    /// rather than claiming an identity it cannot justify.
    #[test]
    fn an_unmodelled_byte_no_longer_reads_as_self_complementary() {
        assert!(!is_self_complementary(b'X'));
        // Was `None` (identity) before #1318; the span now survives.
        assert_eq!(shorten_inversion(b"X", 0, 1), Some((0, 1)));
        assert_eq!(shorten_inversion(b"XX", 0, 2), Some((0, 2)));
        // No typed substitution exists for it, so the caller must refuse.
        assert_eq!(complementary_substitution(b'X'), None);
        // The genuinely self-complementary symbols are unaffected.
        for base in *b"WSN" {
            assert!(is_self_complementary(base), "{} ", base as char);
        }
        // And a real complementary pair still collapses.
        assert_eq!(shorten_inversion(b"AT", 0, 2), None);
    }

    /// The headline case, spelled out: a soft-masked delins that genuinely *is*
    /// a reverse complement must be typed as an inversion.
    ///
    /// `atgc` deliberately, not a palindrome like `aggcct`: a palindrome is its
    /// own reverse complement, so the insert equals the deleted bases and the
    /// *identity* branch is the correct answer — such a case cannot show that
    /// inversion typing works.
    #[test]
    fn a_soft_masked_delins_that_is_a_revcomp_is_typed_as_an_inversion() {
        assert!(
            matches!(
                canonicalize_delins(b"atgc", 0, 4, b"GCAT"),
                DelinsCanonical::Inversion { .. }
            ),
            "a lowercase reference whose delins is a true reverse complement \
             must still be recognised as an inversion",
        );
        // The upper-case reference must of course agree.
        assert_eq!(
            canonicalize_delins(b"atgc", 0, 4, b"GCAT"),
            canonicalize_delins(b"ATGC", 0, 4, b"GCAT"),
        );
    }

    /// `is_revcomp` itself, at the byte level.
    #[test]
    fn is_revcomp_folds_case_on_both_sides() {
        assert!(is_revcomp(b"atgcat", b"ATGCAT"));
        assert!(is_revcomp(b"ATGCAT", b"atgcat"));
        assert!(is_revcomp(b"atgcat", b"atgcat"));
        // Still discriminating — this is not a fold that makes everything true.
        assert!(!is_revcomp(b"atgcat", b"ATGCAA"));
        assert!(!is_revcomp(b"acgt", b"ACG"));
    }

    #[test]
    fn complementary_substitution_never_describes_a_base_as_itself() {
        // A soft-masked self-complementary base used to escape the identity
        // branch — `complement(b'w')` is `b'W'`, which is not `b'w'` — and then
        // `Base::from_char` uppercased both sides, yielding the degenerate
        // `Some((W, W))`. Normalization returned that straight out of the
        // `Inversion` arm as `<pos>W>W`, bypassing the `ref == alt` → identity
        // rewrite that only runs on the original input edit.
        for base in *b"WSNwsn" {
            assert_eq!(
                complementary_substitution(base),
                None,
                "{} is self-complementary, so inverting it substitutes nothing",
                base as char
            );
        }

        // Bases that do have a distinct complement are unaffected, in either case.
        for (base, reference, alternative) in [
            (b'a', crate::hgvs::edit::Base::A, crate::hgvs::edit::Base::T),
            (b'A', crate::hgvs::edit::Base::A, crate::hgvs::edit::Base::T),
            (b'r', crate::hgvs::edit::Base::R, crate::hgvs::edit::Base::Y),
            (b'R', crate::hgvs::edit::Base::R, crate::hgvs::edit::Base::Y),
        ] {
            assert_eq!(
                complementary_substitution(base),
                Some((reference, alternative)),
                "{} inverts to its complement",
                base as char
            );
        }

        // A byte outside the modelled alphabet has no complement to substitute
        // in, upper or lower case.
        for base in *b"Xx" {
            assert_eq!(complementary_substitution(base), None);
        }
    }

    // =========================================================================
    // EXTEND_TANDEM_TRACT TESTS
    // =========================================================================

    #[test]
    fn test_extend_tandem_tract_homopolymer() {
        // ref: "TTAAAATT", anchor [2,4) (the first AA), unit "A"
        let ref_seq = b"TTAAAATT";
        let tract = extend_tandem_tract(ref_seq, 2..4, b"A").unwrap();
        assert_eq!(tract.start, 2);
        assert_eq!(tract.end, 6);
        assert_eq!(tract.ref_count, 4);
    }

    #[test]
    fn test_extend_tandem_tract_multi_base_unit() {
        // ref: "TTGCAGCAGCATT", anchor [5,8) (middle GCA), unit "GCA"
        let ref_seq = b"TTGCAGCAGCATT";
        let tract = extend_tandem_tract(ref_seq, 5..8, b"GCA").unwrap();
        assert_eq!(tract.start, 2);
        assert_eq!(tract.end, 11);
        assert_eq!(tract.ref_count, 3);
    }

    #[test]
    fn test_extend_tandem_tract_anchor_at_5prime_edge() {
        // ref: "AAAATT", anchor [0,2), unit "A"
        let ref_seq = b"AAAATT";
        let tract = extend_tandem_tract(ref_seq, 0..2, b"A").unwrap();
        assert_eq!(tract.start, 0);
        assert_eq!(tract.end, 4);
        assert_eq!(tract.ref_count, 4);
    }

    #[test]
    fn test_extend_tandem_tract_anchor_at_3prime_edge() {
        // ref: "TTAAAA", anchor [4,6), unit "A"
        let ref_seq = b"TTAAAA";
        let tract = extend_tandem_tract(ref_seq, 4..6, b"A").unwrap();
        assert_eq!(tract.start, 2);
        assert_eq!(tract.end, 6);
        assert_eq!(tract.ref_count, 4);
    }

    #[test]
    fn test_extend_tandem_tract_zero_width_anchor() {
        // Zero-width anchor (insertion-point semantics): the anchor itself is
        // trivially unit-periodic; the helper just walks left/right.
        // ref: "TTAAAATT", anchor [4,4), unit "A"
        let ref_seq = b"TTAAAATT";
        let tract = extend_tandem_tract(ref_seq, 4..4, b"A").unwrap();
        assert_eq!(tract.start, 2);
        assert_eq!(tract.end, 6);
        assert_eq!(tract.ref_count, 4);
    }

    #[test]
    fn test_extend_tandem_tract_anchor_not_unit_periodic() {
        // ref: "TTAGAATT", anchor [2,4) is "AG" — not "A"-periodic.
        let ref_seq = b"TTAGAATT";
        assert!(extend_tandem_tract(ref_seq, 2..4, b"A").is_none());
    }

    #[test]
    fn test_extend_tandem_tract_anchor_length_not_multiple_of_unit() {
        // ref: "AAA", anchor [0,3), unit "AA" — len 3 % 2 == 1.
        let ref_seq = b"AAA";
        assert!(extend_tandem_tract(ref_seq, 0..3, b"AA").is_none());
    }

    // =========================================================================
    // DELETION_TO_REPEAT TESTS
    // =========================================================================

    #[test]
    fn test_deletion_to_repeat_homopolymer_two_removed() {
        // ref "TTAAAAATT" (5 A's at indices 2..7). Delete 2 A's at [4..6).
        // After 3' shift, del lands at [5..7). post-shift slice "AA", unit "A", k=2.
        // Tract [2..7), ref_count=5, post_count=3 → A[3] at HGVS [3..7] (1-based).
        let ref_seq = b"TTAAAAATT";
        let r = deletion_to_repeat(ref_seq, 5, 7, false, None).expect("should fire");
        assert_eq!(r.unit, b"A");
        assert_eq!(r.count, 3);
        assert_eq!(r.start, 3); // 1-based HGVS
        assert_eq!(r.end, 7);
    }

    #[test]
    fn test_deletion_to_repeat_multi_base_tandem_two_removed() {
        // ref "TTGCAGCAGCATT" (3 GCAs at [2..11)). Delete 2 GCAs at [5..11).
        // post-shift slice "GCAGCA", unit "GCA", k=2. Tract [2..11), ref_count=3,
        // post_count=1 → GCA[1] at HGVS [3..11] (1-based).
        let ref_seq = b"TTGCAGCAGCATT";
        let r = deletion_to_repeat(ref_seq, 5, 11, false, None).expect("should fire");
        assert_eq!(r.unit, b"GCA");
        assert_eq!(r.count, 1);
        assert_eq!(r.start, 3);
        assert_eq!(r.end, 11);
    }

    #[test]
    fn test_deletion_to_repeat_one_unit_removed_returns_none() {
        // k=1 → stays as del.
        // ref "TTAAAAATT", delete 1 A at [6..7).
        let ref_seq = b"TTAAAAATT";
        assert!(deletion_to_repeat(ref_seq, 6, 7, false, None).is_none());
    }

    #[test]
    fn test_deletion_to_repeat_full_tract_removal_returns_none() {
        // post_count == 0 → stays as del.
        // ref "TTAATT", delete both A's at [2..4). ref_count=2, k=2, post_count=0.
        let ref_seq = b"TTAATT";
        assert!(deletion_to_repeat(ref_seq, 2, 4, false, None).is_none());
    }

    #[test]
    fn test_deletion_to_repeat_non_tandem_returns_none() {
        // Single isolated unit, ref_count < 2.
        // ref "TTGCATT", delete "GCA" at [2..5). smallest_repeat_unit("GCA")="GCA",
        // ref_count=1, returns None.
        let ref_seq = b"TTGCATT";
        assert!(deletion_to_repeat(ref_seq, 2, 5, false, None).is_none());
    }

    #[test]
    fn test_deletion_to_repeat_finer_periodicity() {
        // ref "TTATATATATATT" (5 ATs at [2..12)). Delete "ATAT" at [8..12).
        // smallest_repeat_unit("ATAT")="AT", k=2, ref_count=5, post_count=3 → AT[3].
        let ref_seq = b"TTATATATATATT";
        let r = deletion_to_repeat(ref_seq, 8, 12, false, None).expect("should fire");
        assert_eq!(r.unit, b"AT");
        assert_eq!(r.count, 3);
        assert_eq!(r.start, 3); // 1-based HGVS [3..12]
        assert_eq!(r.end, 12);
    }

    // =========================================================================
    // CANONICALIZATION TESTS
    // =========================================================================

    #[test]
    fn test_canonicalize_deletion_with_length() {
        use crate::hgvs::edit::Sequence;
        use std::str::FromStr;

        // del12 -> del (remove explicit length)
        let edit = NaEdit::Deletion {
            sequence: None,
            length: Some(12),
        };
        let canonical = canonicalize_edit(&edit);
        assert!(matches!(
            canonical,
            NaEdit::Deletion {
                sequence: None,
                length: None
            }
        ));

        // delATG -> del (remove explicit sequence)
        let edit = NaEdit::Deletion {
            sequence: Some(Sequence::from_str("ATG").unwrap()),
            length: None,
        };
        let canonical = canonicalize_edit(&edit);
        assert!(matches!(
            canonical,
            NaEdit::Deletion {
                sequence: None,
                length: None
            }
        ));
    }

    #[test]
    fn test_canonicalize_duplication_with_length() {
        use crate::hgvs::edit::Sequence;
        use std::str::FromStr;

        // dup12 -> dup
        let edit = NaEdit::Duplication {
            sequence: None,
            length: Some(12),
            uncertain_extent: None,
        };
        let canonical = canonicalize_edit(&edit);
        assert!(matches!(
            canonical,
            NaEdit::Duplication {
                sequence: None,
                length: None,
                ..
            }
        ));

        // dupATG -> dup
        let edit = NaEdit::Duplication {
            sequence: Some(Sequence::from_str("ATG").unwrap()),
            length: None,
            uncertain_extent: None,
        };
        let canonical = canonicalize_edit(&edit);
        assert!(matches!(
            canonical,
            NaEdit::Duplication {
                sequence: None,
                length: None,
                ..
            }
        ));
    }

    #[test]
    fn test_should_canonicalize() {
        // Deletion with length should be canonicalized
        assert!(should_canonicalize(&NaEdit::Deletion {
            sequence: None,
            length: Some(12)
        }));

        // Deletion without length shouldn't be canonicalized
        assert!(!should_canonicalize(&NaEdit::Deletion {
            sequence: None,
            length: None
        }));

        // Real substitutions stay on the no-op fast path
        use crate::hgvs::edit::Base;
        assert!(!should_canonicalize(&NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G
        }));

        // Degenerate `A>A` substitutions enter the no-reference canonicalize
        // path so they can be rewritten to identity (`=`) even when the
        // provider has no transcript / genomic sequence loaded.
        assert!(should_canonicalize(&NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::A
        }));

        // Delins with explicit deleted sequence enters the no-reference
        // canonicalize path so the strip-deleted-bases rule fires
        // (closes #338). `Delins { sequence, deleted: None, deleted_length: None }`
        // is already canonical and must NOT re-enter the pass.
        use crate::hgvs::edit::{InsertedSequence, Sequence};
        use std::str::FromStr;
        let trimmed_seq = InsertedSequence::Literal(Sequence::from_str("AC").unwrap());
        assert!(should_canonicalize(&NaEdit::Delins {
            sequence: trimmed_seq.clone(),
            deleted: Some(Sequence::from_str("CA").unwrap()),
            deleted_length: None,
            substitution_reference: None,
        }));
        assert!(should_canonicalize(&NaEdit::Delins {
            sequence: trimmed_seq.clone(),
            deleted: None,
            deleted_length: Some(2),
            substitution_reference: None,
        }));
        assert!(!should_canonicalize(&NaEdit::Delins {
            sequence: trimmed_seq,
            deleted: None,
            deleted_length: None,
            substitution_reference: None,
        }));
    }

    #[test]
    fn test_canonicalize_edit_degenerate_substitution_to_identity() {
        // A4: c.100A>A is "not allowed" per HGVS v21 spec; canonicalize to `=`.
        use crate::hgvs::edit::Base;

        let degenerate = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::A,
        };
        assert_eq!(canonicalize_edit(&degenerate), NaEdit::position_identity());

        // Real SNVs pass through unchanged.
        let real_sub = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        };
        assert_eq!(canonicalize_edit(&real_sub), real_sub);
    }

    // =========================================================================
    // Codon-frame gate tests (#81 B1)
    //
    // HGVS spec (docs/recommendations/DNA/repeated.md, §Notes):
    // > using a coding DNA reference sequence ("c." description) a Repeated
    // > sequence variant description can be used only for repeat units with
    // > a length which is a multiple of 3, i.e. which can not affect the
    // > reading frame.
    // =========================================================================

    #[test]
    fn test_insertion_to_repeat_codon_frame_gate_blocks_a_in_coding() {
        // Reference: 5-A homopolymer flanked by Cs. Insert AA → would normally
        // emit A[7], but unit_len=1 is not a multiple of 3 in c. context,
        // so the gate returns None.
        let ref_seq = b"CAAAAAC";
        let result =
            insertion_to_repeat(ref_seq, 5, b"AA", true, None, ShuffleDirection::ThreePrime);
        assert!(
            result.is_none(),
            "is_coding=true + unit_len=1 must return None"
        );
    }

    #[test]
    fn test_insertion_to_repeat_codon_frame_gate_blocks_at_in_coding() {
        // Reference: AT[3] tandem flanked by Cs. Insert ATAT → unit_len=2,
        // gate blocks in coding context.
        let ref_seq = b"CATATATC";
        let result = insertion_to_repeat(
            ref_seq,
            6,
            b"ATAT",
            true,
            None,
            ShuffleDirection::ThreePrime,
        );
        assert!(
            result.is_none(),
            "is_coding=true + unit_len=2 must return None"
        );
    }

    #[test]
    fn test_insertion_to_repeat_codon_frame_gate_passes_cag_in_coding() {
        // Reference: CAG[3] tandem. Insert CAGCAG → unit_len=3, codon-aligned,
        // gate passes; result is Some(...) carrying CAG[5].
        let ref_seq = b"CCAGCAGCAGT";
        let result = insertion_to_repeat(
            ref_seq,
            9,
            b"CAGCAG",
            true,
            None,
            ShuffleDirection::ThreePrime,
        );
        assert!(
            result.is_some(),
            "is_coding=true + unit_len=3 must allow rewrite"
        );
        let (_first, count, _start, _end, unit) = result.unwrap();
        assert_eq!(count, 5, "expected CAG[5]");
        assert_eq!(unit, b"CAG");
    }

    #[test]
    fn test_insertion_to_repeat_gate_no_op_in_genomic() {
        // Same A-homopolymer case as the blocking test, but is_coding=false
        // → gate is a no-op, repeat rewrite proceeds.
        let ref_seq = b"CAAAAAC";
        let result =
            insertion_to_repeat(ref_seq, 5, b"AA", false, None, ShuffleDirection::ThreePrime);
        assert!(result.is_some(), "is_coding=false must not gate");
    }

    #[test]
    fn test_deletion_to_repeat_codon_frame_gate_blocks_a_in_coding() {
        // 5-A tract, delete 2 A's. Span 2..4 covers two of the As so that
        // without the codon-frame gate the function would return Some(A[3]);
        // with the gate (coding) it must return None. Span 2..3 would also
        // return None via the unrelated `k < 2` early exit, so it would not
        // discriminate the gate.
        let ref_seq = b"CAAAAAC";
        let result = deletion_to_repeat(ref_seq, 2, 4, true, None);
        assert!(
            result.is_none(),
            "is_coding=true + unit_len=1 must return None"
        );
    }

    #[test]
    fn test_deletion_to_repeat_codon_frame_gate_passes_cag_in_coding() {
        // ref "CCAGCAGCAGT": 3-CAG tract at indices 1..10. Delete 2 CAGs at
        // [1..7) (6 bases CAGCAG). With codon-aligned unit, gate passes.
        let ref_seq = b"CCAGCAGCAGT";
        let result = deletion_to_repeat(ref_seq, 1, 7, true, None);
        assert!(
            result.is_some(),
            "is_coding=true + unit_len=3 must allow rewrite"
        );
    }

    #[test]
    fn test_duplication_to_repeat_codon_frame_gate_routes_a_to_gated_insertion() {
        // 4-A tract at indices 1..5 (0-based). Duplicate 2 A's at positions 1..3.
        // Under the gate, structural conditions for repeat notation are met
        // but unit_len=1 in c. is forbidden, so the result routes to a
        // GatedInsertion that the caller renders as `ins<dup_seq>`.
        let ref_seq = b"CAAAAC";
        let result = duplication_to_repeat(ref_seq, 1, 3, true);
        match result {
            Some(DupToRepeatResult::GatedInsertion { sequence, .. }) => {
                assert_eq!(sequence, b"AA", "sequence is the duplicated literal");
            }
            other => panic!("expected GatedInsertion, got {:?}", other),
        }
    }

    #[test]
    fn test_duplication_to_repeat_codon_frame_gate_passes_cag_in_coding() {
        // 3-CAG tract. Duplicate 2 CAGs.
        let ref_seq = b"CCAGCAGCAGT";
        let result = duplication_to_repeat(ref_seq, 1, 7, true);
        assert!(
            result.is_some(),
            "is_coding=true + unit_len=3 must allow rewrite"
        );
    }

    #[test]
    fn test_normalize_repeat_codon_frame_gate_routes_contraction_to_deletion() {
        // 5-A tract, specified A[3] in coding → must NOT emit Repeat; emits Deletion.
        let ref_seq = b"CAAAAAC";
        let result = normalize_repeat(ref_seq, 1, 1, b"A", 3, true, ShuffleDirection::ThreePrime);
        match result {
            RepeatNormResult::Deletion { .. } => {}
            other => panic!("expected Deletion under gate, got {:?}", other),
        }
    }

    #[test]
    fn test_normalize_repeat_codon_frame_gate_routes_expansion_to_insertion() {
        // 5-A tract, specified A[8] in coding → must NOT emit Repeat; emits Insertion.
        let ref_seq = b"CAAAAAC";
        let result = normalize_repeat(ref_seq, 1, 1, b"A", 8, true, ShuffleDirection::ThreePrime);
        match result {
            RepeatNormResult::Insertion { sequence, .. } => {
                assert_eq!(sequence, b"AAA", "3 extra A's");
            }
            other => panic!("expected Insertion under gate, got {:?}", other),
        }
    }

    #[test]
    fn test_normalize_repeat_codon_frame_gate_passes_through_dup_branch() {
        // 5-A tract, specified A[6] in coding → +1 copy = dup, gate doesn't change this.
        let ref_seq = b"CAAAAAC";
        let result = normalize_repeat(ref_seq, 1, 1, b"A", 6, true, ShuffleDirection::ThreePrime);
        match result {
            RepeatNormResult::Duplication { .. } => {}
            other => panic!("expected Duplication, got {:?}", other),
        }
    }

    // =========================================================================
    // ISSUE #160 / #165: decompose_delins tests
    // =========================================================================

    fn sub_at(position: usize, r: char, a: char) -> DelinsSubedit {
        DelinsSubedit::Substitution {
            position,
            reference: crate::hgvs::edit::Base::from_char(r).unwrap(),
            alternative: crate::hgvs::edit::Base::from_char(a).unwrap(),
        }
    }
    fn inv_at(start: usize, end: usize) -> DelinsSubedit {
        DelinsSubedit::Inversion { start, end }
    }
    fn ident_at(position: usize, b: char) -> DelinsSubedit {
        DelinsSubedit::IdentityAt {
            position,
            base: crate::hgvs::edit::Base::from_char(b).unwrap(),
        }
    }

    #[test]
    fn decompose_leading_revcomp_subrun_stays_delins() {
        // ref=TCC, alt=GAG: the 2-nt prefix TC->GA satisfies revcomp(TC)=GA,
        // but the whole run TCC->GAG is a single contiguous change and
        // revcomp(TCC)=GGA != GAG, so it is NOT an inversion. A reverse-
        // complement sub-run may not be carved out (issue #1034): emit three
        // substitutions, which the caller re-merges into one delins → None.
        // This is the exact issue #160 row-2 example, now spec-corrected.
        let result = decompose_delins(b"TCC", 0, 3, b"GAG");
        assert_eq!(result, None);
    }

    #[test]
    fn decompose_trailing_revcomp_subrun_stays_delins() {
        // ref=AAG, alt=GCT: the 2-nt suffix AG->CT satisfies revcomp(AG)=CT,
        // but the whole contiguous run revcomp(AAG)=CTT != GCT is not an
        // inversion. No sub-run carve-out (issue #1034) → None.
        let result = decompose_delins(b"AAG", 0, 3, b"GCT");
        assert_eq!(result, None);
    }

    #[test]
    fn decompose_full_span_inv_returns_none() {
        // ref=GCT, alt=AGC: revcomp(GCT)=AGC, full span is inv. Single-
        // element result with one Inversion → trigger doesn't fire (already
        // handled by canonicalize_delins's full-span check).
        let result = decompose_delins(b"GCT", 0, 3, b"AGC");
        assert_eq!(result, None);
    }

    #[test]
    fn decompose_no_inv_returns_none() {
        // ref=AT, alt=GC: no inv possible (revcomp(AT)=AT != GC).
        let result = decompose_delins(b"AT", 0, 2, b"GC");
        assert_eq!(result, None);
    }

    #[test]
    fn decompose_contiguous_run_with_multiple_revcomp_subruns_stays_delins() {
        // ref=AGACC, alt=CTTGG: the sub-runs AG->CT (revcomp(AG)=CT) and
        // CC->GG (revcomp(CC)=GG) each look like inversions in isolation, but
        // they belong to ONE contiguous mismatch run (all five positions
        // differ). revcomp(AGACC)=GGTCT != CTTGG, so the whole run is not an
        // inversion and no sub-run may be carved out (issue #1034). Emits
        // five substitutions → None.
        let result = decompose_delins(b"AGACC", 0, 5, b"CTTGG");
        assert_eq!(result, None);
    }

    #[test]
    fn decompose_codon_frame_merge_emits_sub_identity_sub() {
        // ref=TAG, alt=AAC: T>A at pos 0, identity at pos 1 (the
        // synthesized middle from issue #79), G>C at pos 2. The
        // decomposition itself yields `[Sub; Identity; Sub]`; preserving
        // the spec's `general.md:35-38` codon-frame exception when the
        // surrounding variant lives in a CDS triplet is the caller's
        // responsibility (`build_split_variants` reads this pattern back
        // out as a 3-base delins for in-codon `c.` variants). Other
        // coord systems decompose to two separate subs.
        let result = decompose_delins(b"TAG", 0, 3, b"AAC");
        assert_eq!(
            result,
            Some(vec![
                sub_at(0, 'T', 'A'),
                ident_at(1, 'A'),
                sub_at(2, 'G', 'C'),
            ])
        );
    }

    #[test]
    fn decompose_complement_only_returns_none() {
        // complement(AC)=TG matches at each position but isn't reversed.
        // revcomp(AC)=GT != TG.
        let result = decompose_delins(b"AC", 0, 2, b"TG");
        assert_eq!(result, None);
    }

    #[test]
    fn decompose_reverse_only_returns_none() {
        // reverse(AC)=CA matches but isn't complemented.
        let result = decompose_delins(b"AC", 0, 2, b"CA");
        assert_eq!(result, None);
    }

    #[test]
    fn decompose_unequal_length_returns_none() {
        let result = decompose_delins(b"AC", 0, 2, b"GTT");
        assert_eq!(result, None);
    }

    #[test]
    fn decompose_offset_start_propagates_position() {
        // A full-run inv (AG->CT) separated by an identity from a trailing
        // non-inv run (CC->GT), placed at offset 100. Positions in the result
        // are 0-indexed offsets into the input ref_seq slice, so the inv and
        // substitution positions must all carry the +100 offset.
        let mut seq = vec![b'A'; 200];
        seq[100..105].copy_from_slice(b"AGACC");
        let result = decompose_delins(&seq, 100, 105, b"CTAGT");
        assert_eq!(
            result,
            Some(vec![
                inv_at(100, 102),
                ident_at(102, 'A'),
                sub_at(103, 'C', 'G'),
                sub_at(104, 'C', 'T'),
            ])
        );
    }

    #[test]
    fn decompose_palindromic_full_span_returns_none() {
        // ref=GCTA, alt=TAGC: full-span inv (revcomp(GCTA)=TAGC). Single-
        // element result — trigger doesn't fire.
        let result = decompose_delins(b"GCTA", 0, 4, b"TAGC");
        assert_eq!(result, None);
    }

    #[test]
    fn decompose_palindromic_inv_subspan_skipped_in_isolation() {
        // ATATC -> ATATG: the ATAT prefix is its own revcomp (palindrome) but
        // `shorten_inversion` collapses it to identity, so it must NOT be
        // emitted as an inv. The scan falls through to per-position
        // `IdentityAt` for each palindrome base and a single trailing
        // `Sub` for C>G. Under issue #165's loosened gate the emission
        // commits because at least one `IdentityAt` is present (item A10
        // sub-only branch).
        //
        // In the full pipeline this input is trimmed away by
        // `canonicalize_delins` before reaching `decompose_delins`; the
        // unit-test case here exercises the function in isolation.
        let result = decompose_delins(b"ATATC", 0, 5, b"ATATG");
        assert_eq!(
            result,
            Some(vec![
                ident_at(0, 'A'),
                ident_at(1, 'T'),
                ident_at(2, 'A'),
                ident_at(3, 'T'),
                sub_at(4, 'C', 'G'),
            ])
        );
    }

    #[test]
    fn decompose_inv_run_bounded_by_identities() {
        // CTATGC -> CATAGG:
        //   pos 0 C==C : identity (bounds the run on the left)
        //   pos 1-3 TAT -> ATA : revcomp(TAT)=ATA ✓ full-run inv
        //   pos 4 G==G : identity (bounds the run on the right)
        //   pos 5 C -> G : substitution
        // The inv run [1..4) is maximal — it is bounded by identities on both
        // sides, so no complementary-outer-pair shortening is needed (a
        // peelable outer pair would imply an identity, which is already a run
        // boundary here). Regression note: the old greedy scan folded the
        // flanking identities into a [0..5) revcomp window and shortened it
        // back to [1..4); the run-based scan reaches [1..4) directly and
        // reports the flanking identities explicitly.
        let result = decompose_delins(b"CTATGC", 0, 6, b"CATAGG");
        assert_eq!(
            result,
            Some(vec![
                ident_at(0, 'C'),
                inv_at(1, 4),
                ident_at(4, 'G'),
                sub_at(5, 'C', 'G'),
            ])
        );
    }

    #[test]
    fn decompose_inv_run_with_identity_in_middle() {
        // Construct a 5-base pattern with inv at left, identity in middle,
        // sub at right: ref=AGXCC, alt=CTXGT for some X where alt[2]==X.
        // Pick X='A': ref=AGACC, alt=CTAGT
        //   inv(0,2): revcomp(AG)=CT ✓
        //   identity(2): A == A ✓
        //   No inv (3,5): revcomp(CC)=GG, alt[3..5]=GT ✗
        //   sub(3): C>G; sub(4): C>T.
        // Verify the algorithm threads through: [Inv, Identity, Sub, Sub].
        let result = decompose_delins(b"AGACC", 0, 5, b"CTAGT");
        assert_eq!(
            result,
            Some(vec![
                inv_at(0, 2),
                ident_at(2, 'A'),
                sub_at(3, 'C', 'G'),
                sub_at(4, 'C', 'T'),
            ])
        );
    }

    #[test]
    fn decompose_sub_only_with_two_identity_gap_emits_split() {
        // ref=ACGT, alt=TCGG: pos 0 sub (A>T), pos 1-2 identity (C, G),
        // pos 3 sub (T>G). The two-identity gap exercises the issue
        // #165 sub-only branch with the gate firing on `has_identity`.
        let result = decompose_delins(b"ACGT", 0, 4, b"TCGG");
        assert_eq!(
            result,
            Some(vec![
                sub_at(0, 'A', 'T'),
                ident_at(1, 'C'),
                ident_at(2, 'G'),
                sub_at(3, 'T', 'G'),
            ])
        );
    }

    #[test]
    fn decompose_length_one_returns_none() {
        // n=1 is below the precondition (`n < 2`). A 1-base "delins" is
        // a substitution range; never a delins shape.
        let result = decompose_delins(b"A", 0, 1, b"T");
        assert_eq!(result, None);
    }

    #[test]
    fn decompose_all_identity_emits_identities_in_isolation() {
        // Pure identity emissions (ref == alt) have no edits to split
        // out. The `has_identity` gate still fires and the function
        // returns `Some([Identity; Identity; Identity])`; downstream
        // `build_split_variants` would emit zero variants. Such inputs
        // never reach `decompose_delins` in the real pipeline —
        // `canonicalize_delins` collapses them to `Identity` first —
        // so the test pins the isolated-call behavior only.
        let result = decompose_delins(b"AAA", 0, 3, b"AAA");
        assert_eq!(
            result,
            Some(vec![ident_at(0, 'A'), ident_at(1, 'A'), ident_at(2, 'A'),])
        );
    }

    #[test]
    fn decompose_out_of_bounds_returns_none() {
        // `start >= end` and `end > ref_seq.len()` both fail the bounds
        // precondition.
        assert_eq!(decompose_delins(b"AAA", 2, 2, b""), None);
        assert_eq!(decompose_delins(b"AAA", 0, 5, b"TTTTT"), None);
    }

    #[test]
    fn decompose_non_iupac_in_identity_position_returns_none() {
        // Sub-position non-IUPAC bytes already abandon the
        // decomposition. After issue #165's review, identity-position
        // non-IUPAC bytes do the same: the whole decomposition is
        // discarded so the caller keeps the original delins and the
        // next round-trip lands on the same string. Without this guard
        // an identity at a non-IUPAC byte would coerce to `Base::N` in
        // the emitted codon-frame triplet's alt sequence, diverging
        // from the input.
        let result = decompose_delins(b"AXG", 0, 3, b"TXC");
        assert_eq!(result, None);
    }

    #[test]
    fn canonicalize_conversion_same_reference_emits_position_range() {
        use crate::hgvs::edit::InsertedSequence;
        let edit = NaEdit::Conversion {
            source: "42536337_42536382".to_string(),
        };
        let got = canonicalize_conversion_to_delins(&edit).expect("expected Some");
        match got {
            NaEdit::Delins {
                sequence: InsertedSequence::PositionRange { start, end },
                ..
            } => {
                assert_eq!(start, 42536337);
                assert_eq!(end, 42536382);
            }
            other => panic!("expected Delins{{PositionRange}}, got {:?}", other),
        }
    }

    #[test]
    fn canonicalize_conversion_cross_reference_emits_bracketed_reference() {
        use crate::hgvs::edit::InsertedSequence;
        let edit = NaEdit::Conversion {
            source: "NM_000089.1:c.789_1011".to_string(),
        };
        let got = canonicalize_conversion_to_delins(&edit).expect("expected Some");
        match &got {
            NaEdit::Delins {
                sequence: InsertedSequence::Reference(s),
                ..
            } => {
                assert_eq!(s, "NM_000089.1:c.789_1011");
            }
            other => panic!("expected Delins{{Reference}}, got {:?}", other),
        }
        // Display of Reference adds the SVD-WG009 brackets.
        assert_eq!(format!("{}", got), "delins[NM_000089.1:c.789_1011]");
    }

    #[test]
    fn canonicalize_conversion_returns_none_for_non_conversion() {
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        assert!(canonicalize_conversion_to_delins(&edit).is_none());
    }

    #[test]
    fn strip_outer_brackets_handles_edge_cases() {
        use super::strip_outer_brackets;
        // Simple matched pair: strip.
        assert_eq!(
            strip_outer_brackets("[NC_000001.11:g.1_2]"),
            "NC_000001.11:g.1_2"
        );
        // No leading bracket: unchanged.
        assert_eq!(
            strip_outer_brackets("NC_000001.11:g.1_2"),
            "NC_000001.11:g.1_2"
        );
        // Two adjacent groups: outer `[` closes before the end. Stripping
        // would yield `A][B`; leave the source unchanged so the caller
        // doesn't fabricate an invalid reference payload.
        assert_eq!(strip_outer_brackets("[A][B]"), "[A][B]");
        // Trailing bracket only: leave unchanged.
        assert_eq!(strip_outer_brackets("[A]B"), "[A]B");
        // Nested matched brackets: strip the outermost pair.
        assert_eq!(strip_outer_brackets("[A[10];T]"), "A[10];T");
        // Degenerate empty pair: strip to empty (caller decides what to
        // do with the empty string).
        assert_eq!(strip_outer_brackets("[]"), "");
        // Single bracket only: unchanged.
        assert_eq!(strip_outer_brackets("["), "[");
        // Reversed brackets: unchanged (would be unbalanced).
        assert_eq!(strip_outer_brackets("]A["), "]A[");
    }

    #[test]
    fn canonicalize_conversion_strips_outer_brackets_from_cross_reference() {
        // Regression: `con[NC_...:g.X_Y]` parses with the brackets included
        // in `source` (the conversion parser captures the whole rest of the
        // token). Without stripping, `Display` of `Reference("[NC_...]")`
        // emits `delins[[NC_...]]` -- doubled brackets that then leak into
        // downstream error messages (`ins[[NC_...]] cross-reference ...`).
        // Canonicalization must drop a single matching outer `[...]` pair
        // before storing the source in `Reference`.
        use crate::hgvs::edit::InsertedSequence;
        let edit = NaEdit::Conversion {
            source: "[NC_000022.11:g.17178616_17178886]".to_string(),
        };
        let got = canonicalize_conversion_to_delins(&edit).expect("expected Some");
        match &got {
            NaEdit::Delins {
                sequence: InsertedSequence::Reference(s),
                ..
            } => {
                assert_eq!(s, "NC_000022.11:g.17178616_17178886");
            }
            other => panic!(
                "expected Delins{{Reference}} without doubled brackets, got {:?}",
                other
            ),
        }
        // Single pair of brackets on round-trip display.
        assert_eq!(
            format!("{}", got),
            "delins[NC_000022.11:g.17178616_17178886]"
        );
    }

    #[test]
    fn resolve_cross_reference_extracts_transcript_from_genomic_context_cds_payload() {
        // A cross-reference `c.` payload whose accession is a genomic-context
        // compound `NG_(NM_)` must resolve the *transcript* (`NM_`) for the CDS
        // lookup, not the whole compound string — otherwise the CDS translation
        // fails with "transcript NG_…(NM_…) has no CDS start". `NM_000088.3` in
        // the mock has cds_start=1 and sequence "ATGCCCAAG…", so c.1_3 == "ATG".
        let provider = crate::reference::mock::MockProvider::with_test_data();
        let bases = resolve_cross_reference_bases("NG_999999.9(NM_000088.3):c.1_3", &provider)
            .expect("genomic-context compound c. cross-reference should resolve");
        assert_eq!(bases, "ATG");
    }

    #[test]
    fn resolve_cross_reference_accepts_a_single_position_payload() {
        // A cross-reference payload may name a single position (no `_` range),
        // e.g. `NM_…:c.2` — a one-base copy. `NM_000088.3` is "ATG…" with
        // cds_start=1, so c.2 == "T". Without single-position support the
        // payload is rejected as an out-of-scope cross-reference shape.
        let provider = crate::reference::mock::MockProvider::with_test_data();
        let bases = resolve_cross_reference_bases("NM_000088.3:c.2", &provider)
            .expect("single-position c. cross-reference should resolve");
        assert_eq!(bases, "T");
    }

    #[test]
    fn resolve_cross_reference_reduces_transcript_accession_for_noncoding_payload() {
        // `n.` is transcript-relative like `c.` (it is classified `Direct` only
        // because it needs no cds_start translation). A genomic-context compound
        // payload addressed in `n.` must therefore also reduce to the transcript
        // accession — `NM_000088.3` is "ATG…", so n.2 == "T". Without the
        // reduction the lookup fails on the compound `NG_…(NM_…)` string.
        let provider = crate::reference::mock::MockProvider::with_test_data();
        let bases = resolve_cross_reference_bases("NG_999999.9(NM_000088.3):n.2", &provider)
            .expect("genomic-context compound n. cross-reference should resolve");
        assert_eq!(bases, "T");
    }

    #[test]
    fn resolve_cross_reference_accepts_a_single_position_direct_axis_payload() {
        // A `Direct`-axis (here `n.`) single-position payload — exercises the
        // single-position path on the Direct branch (no cds_start translation),
        // complementing the `c.` single-position test above. n.2 == "T".
        let provider = crate::reference::mock::MockProvider::with_test_data();
        let bases = resolve_cross_reference_bases("NM_000088.3:n.2", &provider)
            .expect("single-position n. cross-reference should resolve");
        assert_eq!(bases, "T");
    }

    #[test]
    fn resolve_cross_reference_resolves_axisless_coding_payload_in_cds_frame() {
        // An axis-less payload `<ACC>:<range>` (no `c.`/`n.`/`g.` prefix) is
        // interpreted in the accession's native frame (#759). For a *coding*
        // transcript that is `c.` (CDS), not `n.` — `NM_…:100` == `NM_…:c.100`.
        // `NM_001234.1` is "AAAAATGCCC…" with cds_start=5, so the frames differ:
        //   axis-less 1_3 == c.1_3 == "ATG"   (CDS-translated)
        //   n.1_3 == "AAA"                     (sequence-native, NOT used here)
        let provider = crate::reference::mock::MockProvider::with_test_data();
        let bases = resolve_cross_reference_bases("NM_001234.1:1_3", &provider)
            .expect("axis-less coding cross-reference should resolve");
        assert_eq!(bases, "ATG");
        // Same as the explicit `c.` form, and distinct from the `n.` form.
        assert_eq!(
            resolve_cross_reference_bases("NM_001234.1:c.1_3", &provider).unwrap(),
            "ATG"
        );
        assert_eq!(
            resolve_cross_reference_bases("NM_001234.1:n.1_3", &provider).unwrap(),
            "AAA"
        );
    }

    #[test]
    fn resolve_cross_reference_resolves_axisless_noncoding_payload_in_transcript_frame() {
        // A non-coding accession's native frame is `n.` (sequence-native). The
        // mock has no NR_; emulate by checking that the coding/non-coding split
        // is driven by inferred_variant_type: NR_ would be Direct/transcript.
        // Here we assert the helper classifies a coding NM_ as CDS-framed.
        assert_eq!(
            axisless_native_frame("NM_001234.1"),
            Some((InsCoordKind::Cds, true))
        );
        assert_eq!(
            axisless_native_frame("NR_000001.1"),
            Some((InsCoordKind::Direct, true))
        );
        assert_eq!(
            axisless_native_frame("NC_000001.11"),
            Some((InsCoordKind::Direct, false))
        );
        // Compound NG_(NM_) parses to the inner transcript → CDS-framed.
        assert_eq!(
            axisless_native_frame("NG_012337.1(NM_003002.4)"),
            Some((InsCoordKind::Cds, true))
        );
        // Protein / unparseable → out of scope.
        assert_eq!(axisless_native_frame("NP_000079.2"), None);
    }

    #[test]
    fn resolve_cross_reference_accepts_an_axisless_single_position() {
        // A single axis-less position (no `_` range) — a one-base copy in the
        // coding (c.) frame. `NM_001234.1` c.2 == transcript pos 6 == "T".
        let provider = crate::reference::mock::MockProvider::with_test_data();
        let bases = resolve_cross_reference_bases("NM_001234.1:2", &provider)
            .expect("axis-less single-position cross-reference should resolve");
        assert_eq!(bases, "T");
    }

    #[test]
    fn resolve_cross_reference_rna_coding_is_cds_relative() {
        // r. numbering follows c. (CDS-relative) for a coding transcript, NOT
        // the transcript/n. frame. NM_001234.1 has cds_start=5 over
        // "AAAAATGCCCAAG…", so c.2 == transcript pos 6 == "T", while n.2 ==
        // transcript pos 2 == "A". r.2 must equal c.2 ("T"), not n.2 ("A"). #773.
        let provider = crate::reference::mock::MockProvider::with_test_data();
        let r = resolve_cross_reference_bases("NM_001234.1:r.2", &provider)
            .expect("coding r. cross-reference should resolve");
        let c = resolve_cross_reference_bases("NM_001234.1:c.2", &provider)
            .expect("c.2 should resolve");
        let n = resolve_cross_reference_bases("NM_001234.1:n.2", &provider)
            .expect("n.2 should resolve");
        assert_eq!(r, "T");
        assert_eq!(r, c, "coding r.N must match c.N");
        assert_ne!(r, n, "coding r.N must differ from n.N (off by 5'UTR)");
    }

    #[test]
    fn resolve_cross_reference_rna_coding_range_is_cds_relative() {
        // r.1_3 == c.1_3 == transcript pos 5..7 == "ATG" (cds_start=5), not
        // n.1_3 == "AAA".
        let provider = crate::reference::mock::MockProvider::with_test_data();
        let r = resolve_cross_reference_bases("NM_001234.1:r.1_3", &provider)
            .expect("coding r. range should resolve");
        assert_eq!(r, "ATG");
    }

    #[test]
    fn resolve_cross_reference_rna_noncoding_is_transcript_relative() {
        // A non-coding transcript has no CDS, so r. is transcript-relative
        // (== n.). NR_000123.1 is "ACGTACGTACGT", so r.2 == n.2 == "C" and
        // r.1_4 == "ACGT". #773.
        let provider = crate::reference::mock::MockProvider::with_test_data();
        let single = resolve_cross_reference_bases("NR_000123.1:r.2", &provider)
            .expect("non-coding r. cross-reference should resolve");
        assert_eq!(single, "C");
        let range = resolve_cross_reference_bases("NR_000123.1:r.1_4", &provider)
            .expect("non-coding r. range should resolve");
        assert_eq!(range, "ACGT");
    }

    #[test]
    fn resolve_cross_reference_rna_reduces_genomic_context_compound() {
        // A genomic-context compound payload addressed in r. must reduce to its
        // transcript before lookup (r. is transcript/CDS-relative). The inner
        // NM_001234.1 has cds_start=5, so r.2 == c.2 == "T". #773.
        let provider = crate::reference::mock::MockProvider::with_test_data();
        let bases = resolve_cross_reference_bases("NG_999999.9(NM_001234.1):r.2", &provider)
            .expect("genomic-context compound r. cross-reference should resolve");
        assert_eq!(bases, "T");
    }

    #[test]
    fn resolve_cross_reference_inv_suffix_reverse_complements_the_payload() {
        // #1184: a trailing `inv` on a cross-reference range inverts the fetched
        // bases, mirroring the same-reference `PositionRangeInv` arm. `NM_000088.3`
        // has cds_start=1 over "ATGCCCAAG…", so c.1_3 == "ATG" and c.1_3inv is its
        // reverse complement "CAT".
        let provider = crate::reference::mock::MockProvider::with_test_data();
        let forward = resolve_cross_reference_bases("NM_000088.3:c.1_3", &provider)
            .expect("forward c. cross-reference should resolve");
        assert_eq!(forward, "ATG");
        let inverted = resolve_cross_reference_bases("NM_000088.3:c.1_3inv", &provider)
            .expect("inv-suffixed c. cross-reference should resolve");
        assert_eq!(inverted, "CAT");
        assert_eq!(
            inverted,
            crate::sequence::reverse_complement(&forward),
            "inv payload must be the reverse complement of the forward payload"
        );
    }

    #[test]
    fn resolve_cross_reference_inv_suffix_applies_on_every_axis() {
        // The suffix is stripped before axis-specific position translation, so it
        // composes with each axis's frame rather than replacing it. NR_000123.1 is
        // "ACGTACGTACGT" (non-coding: n./r. are transcript-relative), so 1_5 is
        // "ACGTA" — deliberately NOT a palindrome, so this would fail if the
        // suffix were stripped without inverting.
        let provider = crate::reference::mock::MockProvider::with_test_data();
        for reference in [
            "NR_000123.1:n.1_5inv",
            "NR_000123.1:r.1_5inv",
            // Axis-less native frame (#759) — non-coding accession resolves to n.
            "NR_000123.1:1_5inv",
        ] {
            let bases = resolve_cross_reference_bases(reference, &provider)
                .unwrap_or_else(|e| panic!("{reference} should resolve, got {e:?}"));
            assert_eq!(bases, "TACGT", "{reference} must reverse-complement ACGTA");
        }
    }

    #[test]
    fn resolve_cross_reference_inv_requires_a_two_part_range() {
        // Both same-reference parse paths (`parse_external_ref_part`,
        // `parse_simple_count`) only ever build a `PositionRangeInv` from an
        // `A_Binv` range — a bare `Ainv` is a parse error there. The
        // cross-reference form mirrors that: a single position with an `inv`
        // suffix stays an out-of-scope shape rather than silently complementing
        // one base.
        let provider = crate::reference::mock::MockProvider::with_test_data();
        assert!(parse_cross_reference("NM_000088.3:c.2inv").is_none());
        let err = resolve_cross_reference_bases("NM_000088.3:c.2inv", &provider)
            .expect_err("single-position inv must stay unsupported");
        assert!(
            matches!(err, FerroError::UnsupportedVariant { .. }),
            "expected the shape error, got {err:?}"
        );
    }

    #[test]
    fn resolve_cross_reference_inv_does_not_loosen_decoration_rejection() {
        // Stripping `inv` must not become a general-purpose suffix strip: the
        // decorations that were out of scope before are still out of scope with
        // an `inv` on the end, and a lone/embedded `inv` is not a range.
        for reference in [
            // pter/qter keep their own rejection (these appear in the spec corpus).
            "NC_000014.9:g.29745314_qterinv",
            "NC_000011.10:g.pter_15825266inv",
            // Offsets and UTR markers stay rejected.
            "NM_000088.3:c.244-8_249inv",
            "NM_000088.3:c.100_*200inv",
            // Not a range / not a position at all.
            "NM_000088.3:c.inv",
            "NM_000088.3:c.1_inv",
            // `inv` must be a *trailing* suffix, not matched mid-string.
            "NM_000088.3:c.1inv_3",
        ] {
            assert!(
                parse_cross_reference(reference).is_none(),
                "{reference} must remain an out-of-scope cross-reference shape"
            );
        }
    }

    #[test]
    fn canonicalize_conversion_overflow_falls_back_to_reference() {
        use crate::hgvs::edit::InsertedSequence;
        // 21-digit number overflows u64. Falls back to Reference.
        let edit = NaEdit::Conversion {
            source: "123456789012345678901_2".to_string(),
        };
        let got = canonicalize_conversion_to_delins(&edit).expect("expected Some");
        assert!(matches!(
            got,
            NaEdit::Delins {
                sequence: InsertedSequence::Reference(_),
                ..
            }
        ));
    }

    #[test]
    fn canonicalize_conversion_zero_position_falls_back_to_reference() {
        use crate::hgvs::edit::InsertedSequence;
        // HGVS positions are 1-based; `0_0` is not a valid range, so we must
        // not emit a `PositionRange`. Preserve the original numeric source as
        // a `Reference` so the parser sees the same string on round-trip.
        let edit = NaEdit::Conversion {
            source: "0_0".to_string(),
        };
        let got = canonicalize_conversion_to_delins(&edit).expect("expected Some");
        match got {
            NaEdit::Delins {
                sequence: InsertedSequence::Reference(s),
                ..
            } => assert_eq!(s, "0_0"),
            other => panic!("expected Delins{{Reference}} for 0_0, got {:?}", other),
        }
    }

    #[test]
    fn canonicalize_conversion_reversed_range_falls_back_to_reference() {
        use crate::hgvs::edit::InsertedSequence;
        // `10_2` violates `start <= end`. Don't promote it to PositionRange.
        let edit = NaEdit::Conversion {
            source: "10_2".to_string(),
        };
        let got = canonicalize_conversion_to_delins(&edit).expect("expected Some");
        match got {
            NaEdit::Delins {
                sequence: InsertedSequence::Reference(s),
                ..
            } => assert_eq!(s, "10_2"),
            other => panic!("expected Delins{{Reference}} for 10_2, got {:?}", other),
        }
    }

    #[test]
    fn canonicalize_conversion_zero_start_falls_back_to_reference() {
        use crate::hgvs::edit::InsertedSequence;
        // `0_5`: start < 1 violates HGVS 1-based position rule.
        let edit = NaEdit::Conversion {
            source: "0_5".to_string(),
        };
        let got = canonicalize_conversion_to_delins(&edit).expect("expected Some");
        assert!(
            matches!(
                got,
                NaEdit::Delins {
                    sequence: InsertedSequence::Reference(_),
                    ..
                }
            ),
            "expected Delins{{Reference}} fallback for 0_5"
        );
    }

    // =========================================================================
    // INSERTION_TO_DUPLICATION TESTS (issue #132)
    // =========================================================================

    #[test]
    fn test_insertion_to_duplication_homopolymer_matched_3prime() {
        // ref "TTAAATT": insert A at pos 3 (between idx 3 and idx 4 = inside
        // the A-tract). unit "A", only one rotation. tract [2..5), ref_count=3.
        // Most-3' A dup → dup_start_idx = tract_end-1 = 4 → HGVS [5..5].
        let ref_seq = b"TTAAATT";
        let r = insertion_to_duplication(ref_seq, 3, b"A", ShuffleDirection::ThreePrime)
            .expect("should fire");
        assert_eq!(r.unit, b"A");
        assert_eq!(r.start, 5);
        assert_eq!(r.end, 5);
    }

    #[test]
    fn test_insertion_to_duplication_homopolymer_matched_5prime() {
        // Same ref/tract as the 3prime case, but direction=5prime picks
        // the most-5' dup position. tract at idx [2..5), ref_count=3:
        // dup_start_idx = ref_start = 2 → HGVS [3..3]. Closes #340 (i).
        let ref_seq = b"TTAAATT";
        let r = insertion_to_duplication(ref_seq, 3, b"A", ShuffleDirection::FivePrime)
            .expect("should fire");
        assert_eq!(r.unit, b"A");
        assert_eq!(r.start, 3);
        assert_eq!(r.end, 3);
    }

    #[test]
    fn test_insertion_to_duplication_dinucleotide_5prime() {
        // Same ref / alt as the cyclic_rotation_two_base test but under
        // FivePrime: the GT tract spans idx [2..8), ref_count=3. The
        // 5'-most unit-aligned slot starts at ref_start=2, end=3 →
        // HGVS [3..4]. Locks in unit alignment of `ref_start` under
        // rotation — a regression that picked mid-unit would land at
        // an odd HGVS pair like [3..3] or [4..5].
        let ref_seq = b"ACGTGTGTAC";
        let r = insertion_to_duplication(ref_seq, 2, b"TG", ShuffleDirection::FivePrime)
            .expect("should fire");
        assert_eq!(r.unit, b"GT");
        assert_eq!(r.start, 3);
        assert_eq!(r.end, 4);
    }

    #[test]
    fn test_insertion_to_duplication_single_copy_tract_is_direction_invariant() {
        // ref_count == 1 (a single unit, no tandem run). Under either
        // direction the 5'-most and 3'-most slots coincide
        // (`ref_start == tract_end_idx - unit.len()`), so the helper
        // must return the same answer regardless of `ShuffleDirection`.
        // Locks in the intentional no-op of the new arm for length-1
        // tracts (Sonnet review feedback).
        let ref_seq = b"TTGCATT";
        let three = insertion_to_duplication(ref_seq, 2, b"G", ShuffleDirection::ThreePrime)
            .expect("should fire");
        let five = insertion_to_duplication(ref_seq, 2, b"G", ShuffleDirection::FivePrime)
            .expect("should fire");
        assert_eq!(three.ref_count, 1);
        assert_eq!(five.ref_count, 1);
        assert_eq!(three.start, five.start);
        assert_eq!(three.end, five.end);
    }

    #[test]
    fn test_insertion_to_duplication_cyclic_rotation_two_base() {
        // ref "ACGTGTGTAC": GT tract at indices [2..8), ref_count=3.
        // Insert TG at pos=2 (between idx 2=G and idx 3=T — at 5' edge of
        // the tract). alt "TG" is the r=1 rotation of canonical "GT".
        // Rotation iteration over `smallest_repeat_unit("TG") = "TG"`:
        //   r=0 "TG": anchor=3 ref[3..5]="TG" matches; tract [3..7),
        //     ref_count=2. (anchor=1 "CG" no; anchor=2 "GT" no.)
        //   r=1 "GT": anchor=2 ref[2..4]="GT" matches; tract [2..8),
        //     ref_count=3.
        // "GT" wins with the larger run. Most-3' GT dup: tract end=8,
        // dup_start_idx=6, dup_end_idx=7, HGVS [7..8].
        //
        // The non-tract flanking ('A' on both sides) ensures the TG-phase
        // run is shorter than the GT-phase run; this mirrors the actual
        // padded-sequence behavior in `MockProvider`-driven integration
        // tests (the issue #132 reproducer at `g.258_259insTG`).
        let ref_seq = b"ACGTGTGTAC";
        let r = insertion_to_duplication(ref_seq, 2, b"TG", ShuffleDirection::ThreePrime)
            .expect("should fire");
        assert_eq!(r.unit, b"GT");
        assert_eq!(r.start, 7);
        assert_eq!(r.end, 8);
    }

    #[test]
    fn test_insertion_to_duplication_offphase_slides_through_partial_flank() {
        // #864: an ODD-length di-repeat run ("TGTGTGTGT", 9 bases) ends one base
        // past the last full 2-copy window. The dup must still reach the run's
        // true 3' end (3'rule), which lands on the GT phase — NOT stop at the
        // abutting TG-phase boundary. ref "ACTGTGTGTGTAC": insert TG at pos=1
        // (between C@1 and the run). `three_prime_align_tract` slides through the
        // trailing base → dup unit "GT" at 1-based [10..11].
        //
        // Mirrors the real, mutalyzer-confirmed locus
        // NC_000001.11:g.5010037_5010038insTG → g.5010045_5010046dup (same
        // ACTGTGTGTGTAC context). Before #864 this stopped one base 5' (the
        // under-shifted, mutalyzer-divergent form).
        let ref_seq = b"ACTGTGTGTGTAC";
        let r = insertion_to_duplication(ref_seq, 1, b"TG", ShuffleDirection::ThreePrime)
            .expect("should fire");
        assert_eq!(r.unit, b"GT");
        assert_eq!(r.start, 10);
        assert_eq!(r.end, 11);
    }

    /// #1355: an out-of-phase rotation with a longer tract must not discard a
    /// valid in-phase repeat.
    ///
    /// The same defect #1349 fixed in `insertion_to_duplication`, one function
    /// over. `insertion_to_repeat` scored every cyclic rotation, kept the single
    /// highest-`ref_count` winner, and only then applied the #882 phase gate — so
    /// a rotation that abuts a longer but out-of-phase tract took the whole result
    /// down with it.
    ///
    /// On `GTACGTCCCCTCCTCCTCCTCCCGGGG` with a 2-copy `TCCTCC` inserted at 7:
    ///   - rotation `TCC` abuts a one-copy in-phase tract at `[5, 8)` — valid;
    ///   - rotation `CCT` abuts a four-copy tract at index 8 — out of phase, and
    ///     it wins on `ref_count`.
    #[test]
    fn an_out_of_phase_rotation_does_not_discard_a_valid_repeat() {
        let r = b"GTACGTCCCCTCCTCCTCCTCCCGGGG";
        // The candidate that was being thrown away is genuinely valid.
        assert!(
            tandem_edit_preserves_insertion(r, 7, b"TCCTCC", 5, 8, b"TCC", 3),
            "the discarded `TCC` candidate must pass the #882 gate"
        );
        // Both directions must find it. `find_tandem_extent` scores the rotations
        // identically either way, so the loser-takes-all bug was direction-blind.
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            let emitted = insertion_to_repeat(r, 7, b"TCCTCC", false, None, direction)
                .unwrap_or_else(|| panic!("{direction:?}: a valid in-phase repeat must be found"));
            let (_, count, start, end, unit) = emitted;
            assert_eq!(
                (unit.as_slice(), count, start, end),
                (b"TCC".as_slice(), 3, 6, 8),
                "{direction:?}: expected the in-phase `TCC[3]` over the 1-based tract 6_8"
            );
        }
    }

    /// The emitted repeat must decode back to the input insertion, whichever
    /// rotation won. This is the property the phase gate exists to guarantee, and
    /// asserting it separately means a future selection change cannot satisfy the
    /// pinned values above while breaking the contract underneath them.
    #[test]
    fn every_emitted_repeat_decodes_to_its_input_insertion() {
        let r = b"GTACGTCCCCTCCTCCTCCTCCCGGGG";
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            let (_, count, start, end, unit) =
                insertion_to_repeat(r, 7, b"TCCTCC", false, None, direction).expect("fires");
            let start_idx = hgvs_pos_to_index(start);
            let end_excl = hgvs_pos_to_index(end) + 1;
            assert!(
                tandem_edit_preserves_insertion(r, 7, b"TCCTCC", start_idx, end_excl, &unit, count),
                "{direction:?}: emitted `{}[{count}]` over [{start_idx}, {end_excl}) does not \
                 decode to the input insertion",
                String::from_utf8_lossy(&unit)
            );
        }
    }

    #[test]
    fn tandem_edit_preserves_insertion_in_phase_dup_true() {
        // GACACACATAGGT: insCA at the A|C cut (ins_pos=1) equals dup of the
        // 3'-most CA (region [6,8), one unit -> total 2).
        let r = b"GACACACATAGGT";
        assert!(tandem_edit_preserves_insertion(r, 1, b"CA", 6, 8, b"CA", 2));
    }

    #[test]
    fn tandem_edit_preserves_insertion_out_of_phase_dup_false() {
        // Same tract, insAC at the same cut: out of phase -> reject.
        let r = b"GACACACATAGGT";
        assert!(!tandem_edit_preserves_insertion(
            r, 1, b"AC", 6, 8, b"CA", 2
        ));
    }

    #[test]
    fn tandem_edit_preserves_insertion_in_phase_repeat_true() {
        // TGCGCA: in-phase insCGCG at ins_pos=1 equals GC[4] over region [1,5).
        let r = b"TGCGCA";
        assert!(tandem_edit_preserves_insertion(
            r, 1, b"CGCG", 1, 5, b"GC", 4
        ));
    }

    #[test]
    fn tandem_edit_preserves_insertion_out_of_phase_repeat_false() {
        // TGCGCA: out-of-phase insGCGC at ins_pos=1 vs GC[4] over [1,5) -> reject.
        let r = b"TGCGCA";
        assert!(!tandem_edit_preserves_insertion(
            r, 1, b"GCGC", 1, 5, b"GC", 4
        ));
    }

    #[test]
    fn tandem_edit_preserves_insertion_bounds_guards_no_panic() {
        let r = b"ACGT";
        // ins_pos == len (would panic on ref[..=ins_pos])
        assert!(!tandem_edit_preserves_insertion(r, 4, b"A", 0, 1, b"A", 2));
        // region_end_excl > len
        assert!(!tandem_edit_preserves_insertion(r, 1, b"A", 0, 5, b"A", 2));
        // empty unit
        assert!(!tandem_edit_preserves_insertion(r, 1, b"A", 0, 1, b"", 2));
        // region_start > region_end_excl
        assert!(!tandem_edit_preserves_insertion(r, 1, b"A", 2, 1, b"A", 2));
    }

    #[test]
    fn insertion_to_duplication_rejects_out_of_phase_882() {
        // GACACACATAGGT, insAC at the A|C cut (pos=1) is out of phase with the
        // AC-periodic tract; today it wrongly becomes a dup. Must be None now.
        let r = b"GACACACATAGGT";
        assert!(insertion_to_duplication(r, 1, b"AC", ShuffleDirection::ThreePrime).is_none());
    }

    #[test]
    fn insertion_to_duplication_accepts_in_phase_882() {
        // insCA at the same cut IS a legitimate tandem dup (3'-most CA).
        let r = b"GACACACATAGGT";
        let got = insertion_to_duplication(r, 1, b"CA", ShuffleDirection::ThreePrime)
            .expect("in-phase insCA should still convert to dup");
        assert_eq!(got.unit, b"CA");
        assert_eq!(got.start, 7);
        assert_eq!(got.end, 8);
    }

    /// Regression for the `FivePrime` branch of `insertion_to_duplication`: the
    /// odd-length off-phase partial flank that `five_prime_align_tract` exists to
    /// cross.
    ///
    /// In `ACTGTGTGTGTAC` the alternating run spans indices 2..=10
    /// (`T G T G T G T G T`). The `GT` tract abutting the cut starts at index 3,
    /// but the run continues ONE base further 5' in the *opposite* phase — the
    /// stray `T` at index 2. Stopping at `ref_start` would emit HGVS 4_5, which
    /// is not a fixed point (the dup re-shuffles one base further 5' on the next
    /// pass). The aligner must rotate `GT` → `TG` and land on 3_4.
    #[test]
    fn insertion_to_duplication_five_prime_crosses_off_phase_partial_flank() {
        let r = b"ACTGTGTGTGTAC";
        let got = insertion_to_duplication(r, 10, b"GT", ShuffleDirection::FivePrime)
            .expect("5' insGT inside the GT run must convert to a dup");
        assert_eq!(
            got.unit, b"TG",
            "unit rotates right as the window slides one base 5'"
        );
        assert_eq!(
            (got.start, got.end),
            (3, 4),
            "the 5'-most slot lies across the off-phase T at index 2"
        );

        // The emitted dup must describe the SAME haplotype as the input insertion.
        let inserted = {
            let mut s = r.to_vec();
            s.splice(11..11, b"GT".iter().copied());
            s
        };
        let duplicated = {
            let (s, e) = (got.start as usize - 1, got.end as usize); // 1-based incl → 0-based excl
            let mut v = r.to_vec();
            let copy = r[s..e].to_vec();
            v.splice(e..e, copy);
            v
        };
        assert_eq!(
            String::from_utf8_lossy(&inserted),
            String::from_utf8_lossy(&duplicated),
            "dup and insertion spellings must reconstruct byte-identical sequences"
        );
    }

    /// Regression for #1349: a higher-scoring rotation that fails the #882
    /// equivalence check must not shadow a lower-scoring rotation that passes it.
    ///
    /// In `GTACGTCCCCTCCTC` the insertion point is between indices 7 and 8
    /// (`pos = 7`) and the alt `TCC` equals the immediately preceding tract at
    /// indices 5..8 — a textbook single-copy duplication. But a *different*
    /// rotation of the same unit, `CCT`, has a two-copy tract at indices 8..14,
    /// so it wins the `ref_count` comparison outright. That tract is out of phase
    /// with the insertion, so it fails `tandem_edit_preserves_insertion` — and
    /// before this fix the whole function returned `None`, discarding the correct
    /// `TCC` candidate that was never reconsidered.
    ///
    /// Downstream that `None` cost the 5' path its `ins → dup` promotion, and the
    /// CDS-start clamp turned the surviving insertion into a boundary `delins`
    /// that was not a fixed point (`c.2_3insTCC` → `c.1delinsCCTC` → `c.-1_2dup`).
    #[test]
    fn insertion_to_duplication_falls_back_past_a_non_preserving_rotation() {
        let r = b"GTACGTCCCCTCCTC";
        let got = insertion_to_duplication(r, 7, b"TCC", ShuffleDirection::FivePrime)
            .expect("the in-phase TCC rotation is a valid dup and must not be discarded");
        assert_eq!(got.unit, b"TCC");
        assert_eq!(
            (got.start, got.end),
            (6, 8),
            "the dup covers the preceding TCC tract at 0-based 5..8"
        );

        // The emitted dup must describe the SAME haplotype as the input insertion.
        let inserted = {
            let mut s = r.to_vec();
            s.splice(8..8, b"TCC".iter().copied());
            s
        };
        let duplicated = {
            let (s, e) = (got.start as usize - 1, got.end as usize); // 1-based incl → 0-based excl
            let mut v = r.to_vec();
            let copy = r[s..e].to_vec();
            v.splice(e..e, copy);
            v
        };
        assert_eq!(
            String::from_utf8_lossy(&inserted),
            String::from_utf8_lossy(&duplicated),
            "dup and insertion spellings must reconstruct byte-identical sequences"
        );
    }

    /// The 3' direction reaches the same fallback: the winning rotation is chosen
    /// by `ref_count` before any equivalence check, so the out-of-phase `CCT`
    /// tract shadows the correct `TCC` one in both directions.
    #[test]
    fn insertion_to_duplication_falls_back_past_a_non_preserving_rotation_three_prime() {
        let r = b"GTACGTCCCCTCCTC";
        let got = insertion_to_duplication(r, 7, b"TCC", ShuffleDirection::ThreePrime)
            .expect("the in-phase TCC rotation is a valid dup and must not be discarded");
        assert_eq!(
            (got.unit.as_slice(), got.start, got.end),
            (&b"TCC"[..], 6, 8)
        );
    }

    /// #1846: a payload whose smallest repeat unit is longer than the reference
    /// window cannot convert — no rotation can match a tract that does not fit —
    /// so both conversions decline. Pins the *output* the linear early exit
    /// preserves, on inputs small enough to read at a glance.
    #[test]
    fn insertion_conversions_decline_when_the_unit_outgrows_the_window_1846() {
        // Window shorter than the (aperiodic, single-copy) payload: dup path.
        let ref_seq = b"ACGT";
        assert!(
            insertion_to_duplication(ref_seq, 1, b"ACGTAC", ShuffleDirection::ThreePrime).is_none()
        );
        assert!(
            insertion_to_duplication(ref_seq, 1, b"ACGTAC", ShuffleDirection::FivePrime).is_none()
        );
        // Two copies of a 6-base unit that is itself longer than the window:
        // repeat path (`added_copies == 2`).
        assert!(insertion_to_repeat(
            ref_seq,
            1,
            b"ACGTAAACGTAA",
            false,
            None,
            ShuffleDirection::ThreePrime
        )
        .is_none());
    }

    /// #1846 boundary guard: a unit whose length exactly EQUALS the reference
    /// window must still convert — `find_tandem_extent` accepts anchor 0 there
    /// (`0 + unit.len() == ref_seq.len()`, not `>`). This pins that the early
    /// exit is strict `>` and not `>=`; a `>=` would silently decline this real
    /// single-copy duplication while passing every other test.
    #[test]
    fn insertion_to_duplication_converts_when_the_unit_exactly_fills_the_window_1846() {
        // `ACGT` inserted at the 3' end duplicates the whole 4-base window.
        let ref_seq = b"ACGT";
        let got = insertion_to_duplication(ref_seq, 3, b"ACGT", ShuffleDirection::ThreePrime)
            .expect("a single-copy dup of the whole window must still convert");
        assert_eq!(got.unit, b"ACGT");
        assert_eq!((got.start, got.end), (1, 4));
    }

    /// #1846: `insertion_to_duplication` must not scan every rotation of the
    /// expanded insert payload. When the smallest repeat unit is longer than the
    /// reference window, `find_tandem_extent` skips every anchor (it needs
    /// `anchor + unit.len() <= ref_seq.len()`), so no rotation can match and the
    /// result is *determined* to be `None` — but the pre-#1846 loop still
    /// allocated and filled a payload-sized `Vec` on each of `payload.len()`
    /// iterations, i.e. Θ(L²) bytes copied to reach a `None` already fixed. A
    /// spec-legal `ins[ACC:g.A_B]` expands to that whole named range, so this
    /// never terminated on a megabase cross-reference insert.
    ///
    /// Measured on the unfixed code, this doubled from ~1.0 s at L=200 kb to
    /// ~5.7 s at L=400 kb (≈4-6x per doubling); the early exit makes it O(L).
    /// The wall-clock bound is a loud backstop against a quadratic regression,
    /// generous by ~three orders of magnitude over the linear cost (a few ms) so
    /// it cannot flake on a loaded machine.
    #[test]
    fn insertion_to_duplication_does_not_rescan_a_megabase_payload_1846() {
        let ref_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        // Aperiodic, single-copy: `smallest_repeat_unit` returns the whole
        // payload, so `added_copies == 1` and the rotation loop is the hot path.
        let mut payload = vec![b'A'; 2_000_000];
        payload[0] = b'C';
        let start = std::time::Instant::now();
        let got = insertion_to_duplication(ref_seq, 30, &payload, ShuffleDirection::ThreePrime);
        let elapsed = start.elapsed();
        assert!(
            got.is_none(),
            "a payload longer than the reference window cannot be a duplication"
        );
        assert!(
            elapsed < std::time::Duration::from_secs(5),
            "insertion_to_duplication took {elapsed:?} on a 2 Mb payload — the O(payload^2) \
             rotation scan is back"
        );
    }

    /// #1846 sibling: `insertion_to_repeat` carries the identical rotation loop,
    /// reached when the payload is two or more copies of its unit. Two copies of
    /// a 1 Mb aperiodic unit exercises it, and the same early exit applies — the
    /// unit outgrows the window, so every rotation's `find_tandem_extent`
    /// returns `None` and the result is a determined `None` in O(unit) rather
    /// than Θ(unit²).
    #[test]
    fn insertion_to_repeat_does_not_rescan_a_megabase_payload_1846() {
        let ref_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        // Two copies of a 1 Mb aperiodic unit: `smallest_repeat_unit` returns the
        // 1 Mb unit and `added_copies == 2`, so the rotation loop runs.
        let unit_len = 1_000_000;
        let mut unit = vec![b'A'; unit_len];
        unit[0] = b'C';
        let mut payload = unit.clone();
        payload.extend_from_slice(&unit);
        let start = std::time::Instant::now();
        let got = insertion_to_repeat(
            ref_seq,
            30,
            &payload,
            false,
            None,
            ShuffleDirection::ThreePrime,
        );
        let elapsed = start.elapsed();
        assert!(
            got.is_none(),
            "a repeat unit longer than the reference window cannot match a tract"
        );
        assert!(
            elapsed < std::time::Duration::from_secs(5),
            "insertion_to_repeat took {elapsed:?} on a 2 Mb payload — the O(payload^2) \
             rotation scan is back"
        );
    }

    #[test]
    fn insertion_to_repeat_rejects_out_of_phase_882() {
        // TGCGCA: insGCGC at the G|C cut (pos=1) is out of phase with the GC
        // tract; today it wrongly becomes GC[4]. Must be None now.
        let r = b"TGCGCA";
        assert!(
            insertion_to_repeat(r, 1, b"GCGC", false, None, ShuffleDirection::ThreePrime).is_none()
        );
    }

    #[test]
    fn insertion_to_repeat_accepts_in_phase_882() {
        // In-phase rotation insCGCG reaches the same GC[4] window and is preserved.
        let r = b"TGCGCA";
        let (_b, count, start, end, unit) =
            insertion_to_repeat(r, 1, b"CGCG", false, None, ShuffleDirection::ThreePrime)
                .expect("in-phase should convert");
        assert_eq!(count, 4);
        assert_eq!((start, end), (2, 5));
        assert_eq!(unit, b"GC");
    }

    #[test]
    fn test_insertion_to_duplication_no_adjacent_tract() {
        // ref "ACGTACGT": no tandem of "X"; insert "X" at any pos returns None.
        let ref_seq = b"ACGTACGT";
        assert!(
            insertion_to_duplication(ref_seq, 3, b"X", ShuffleDirection::ThreePrime,).is_none()
        );
    }

    #[test]
    fn test_insertion_to_duplication_empty_or_oob() {
        let ref_seq = b"TTAAATT";
        assert!(insertion_to_duplication(ref_seq, 3, b"", ShuffleDirection::ThreePrime,).is_none());
        assert!(insertion_to_duplication(b"", 0, b"A", ShuffleDirection::ThreePrime,).is_none());
    }

    #[test]
    fn test_insertion_to_duplication_rejects_multi_copy() {
        // alt is 2 copies of unit A → not a 1-copy ins; helper returns None
        // (caller routes 2+ copies to insertion_to_repeat instead).
        let ref_seq = b"TTAAATT";
        assert!(
            insertion_to_duplication(ref_seq, 3, b"AA", ShuffleDirection::ThreePrime,).is_none()
        );
    }

    #[test]
    fn test_insertion_to_duplication_rejects_out_of_phase_alt_at_interior_cut() {
        // #882: same reference layout as
        // `test_insertion_to_duplication_cyclic_rotation_two_base`, but the alt
        // is the OUT-OF-PHASE rotation "GT" (not the in-phase "TG"). At the
        // interior G|T cut (pos=2, between idx2=G and idx3=T of the GTGTGT
        // tract), inserting "TG" extends the tandem (-> ACGTGTGTGTAC, a real
        // dup), but inserting "GT" yields ACGGTTGTGTAC — which is NOT a tandem
        // expansion. Pre-#882 this was wrongly emitted as g.7_8dup (decoding to
        // the different ACGTGTGTGTAC). The gate now rejects it; it stays `ins`,
        // matching mutalyzer. (This test previously asserted the corrupt dup
        // under the misnomer `phase_matched_first_base`.)
        let ref_seq = b"ACGTGTGTAC";
        assert!(
            insertion_to_duplication(ref_seq, 2, b"GT", ShuffleDirection::ThreePrime).is_none()
        );
    }

    #[test]
    fn test_canonicalize_edit_delins_strips_explicit_deleted_seq() {
        use crate::hgvs::edit::{InsertedSequence, Sequence};
        use std::str::FromStr;
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::from_str("TTCC").unwrap()),
            deleted: Some(Sequence::from_str("ATG").unwrap()),
            deleted_length: None,
            substitution_reference: None,
        };
        let canonical = canonicalize_edit(&edit);
        match canonical {
            NaEdit::Delins {
                sequence: _,
                deleted,
                deleted_length,
                ..
            } => {
                assert_eq!(deleted, None);
                assert_eq!(deleted_length, None);
            }
            other => panic!("expected NaEdit::Delins, got {other:?}"),
        }
    }

    #[test]
    fn test_canonicalize_edit_delins_strips_explicit_deleted_count() {
        use crate::hgvs::edit::{InsertedSequence, Sequence};
        use std::str::FromStr;
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::from_str("TA").unwrap()),
            deleted: None,
            deleted_length: Some(3),
            substitution_reference: None,
        };
        let canonical = canonicalize_edit(&edit);
        match canonical {
            NaEdit::Delins {
                sequence: _,
                deleted,
                deleted_length,
                ..
            } => {
                assert_eq!(deleted, None);
                assert_eq!(deleted_length, None);
            }
            other => panic!("expected NaEdit::Delins, got {other:?}"),
        }
    }

    #[test]
    fn rna_uracil_to_thymine_rewrites_literal_u_in_insertion() {
        use crate::hgvs::edit::InsertedSequence;
        let edit = NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::from_str("UU").unwrap()),
        };
        let mapped = rna_uracil_to_thymine(&edit).expect("expected Some for a u-bearing edit");
        assert_eq!(
            mapped,
            NaEdit::Insertion {
                sequence: InsertedSequence::Literal(Sequence::from_str("TT").unwrap()),
            }
        );
    }

    #[test]
    fn rna_uracil_to_thymine_rewrites_delins_deleted_and_inserted() {
        use crate::hgvs::edit::InsertedSequence;
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::from_str("AU").unwrap()),
            deleted: Some(Sequence::from_str("UG").unwrap()),
            deleted_length: None,
            substitution_reference: None,
        };
        let mapped = rna_uracil_to_thymine(&edit).expect("expected Some");
        assert_eq!(
            mapped,
            NaEdit::Delins {
                sequence: InsertedSequence::Literal(Sequence::from_str("AT").unwrap()),
                deleted: Some(Sequence::from_str("TG").unwrap()),
                deleted_length: None,
                substitution_reference: None,
            }
        );
    }

    #[test]
    fn rna_uracil_to_thymine_is_noop_without_u() {
        use crate::hgvs::edit::{Base, InsertedSequence};
        // A DNA-only insertion (no U anywhere) must return None.
        assert!(rna_uracil_to_thymine(&NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::from_str("ACGT").unwrap()),
        })
        .is_none());
        // A non-literal-base shape returns None too.
        assert!(rna_uracil_to_thymine(&NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        })
        .is_none());
    }

    // -----------------------------------------------------------------
    // `inverted_adjacent_copy_span` and its coincidence floor.
    // -----------------------------------------------------------------

    mod inverted_adjacent_copy {
        use super::*;

        fn revcomp(s: &str) -> String {
            crate::sequence::reverse_complement(s)
        }

        /// A contig with a sharp AT→GC composition boundary, the junction
        /// 30 bases into the GC half. Both halves are non-repetitive enough
        /// that the 10-base span ending at the junction is neither
        /// self-complementary nor half of a tandem repeat.
        fn composition_boundary_contig() -> (String, usize) {
            let at = "ATTAATAAT".repeat(50);
            let gc = "GGCCGCGGC".repeat(50);
            let gap_idx = at.len() + 30;
            (format!("{at}{gc}"), gap_idx)
        }

        /// The verdict must not depend on how much reference the caller
        /// happened to hand in.
        ///
        /// `ref_seq` is the **normalization window**, whose extent is
        /// `NormalizeConfig::window_size` (a public field, default 100) grown
        /// by `normalize_in_grown_window` until the 3'-shift settles. None of
        /// that is a property of the variant or the reference, so keying the
        /// coincidence floor on the composition of the whole slice made the
        /// emitted description a function of a caller-set knob. Measured on
        /// this contig before the fix, one variant, one reference, one input:
        ///
        /// ```text
        /// len=10 payload=GCCGCCGCGG  w10..w80 = lit    w120..w300 = INV
        /// ```
        ///
        /// [`COMPOSITION_HALF_WIDTH`] makes the floor key on a fixed window
        /// centred on the junction, so every supplied extent that *contains*
        /// that window agrees. Extents narrower than it still clamp — there is
        /// no more sequence to read — which is why the assertion starts at the
        /// half-width and why the half-width is set well inside the default
        /// `window_size`.
        #[test]
        fn the_verdict_is_independent_of_how_much_reference_is_supplied() {
            let (contig, gap_idx) = composition_boundary_contig();
            let span = &contig[gap_idx - 10..gap_idx];
            let payload = revcomp(span);
            assert_ne!(payload, span, "case must not be self-complementary");

            // Slice the same contig to a series of extents around the
            // junction, exactly as a different `window_size` would.
            let mut verdicts = Vec::new();
            for half in [COMPOSITION_HALF_WIDTH, 80, 100, 120, 200, 300] {
                let start = gap_idx.saturating_sub(half);
                let end = (gap_idx + half).min(contig.len());
                let sliced = &contig.as_bytes()[start..end];
                // `gap` is expressed in the frame the slice is indexed in.
                let local_gap = (gap_idx - start) as u64;
                verdicts.push((
                    half,
                    inverted_adjacent_copy_span(sliced, local_gap, payload.as_bytes()).is_some(),
                ));
            }

            let first = verdicts[0].1;
            assert!(
                verdicts.iter().all(|(_, v)| *v == first),
                "the emitted form depends on the supplied window extent: {verdicts:?}"
            );
        }

        /// Both candidate spans match exactly when the `2 * len` bases
        /// straddling the junction are a **tandem repeat of period `len`** —
        /// not when the neighbourhood is palindromic.
        ///
        /// Candidate 0 is `A = ref[gap-len..gap]`, candidate 1 is
        /// `B = ref[gap..gap+len]`; both reverse-complement-match one payload
        /// iff `revcomp(A) == revcomp(B)` iff `A == B`. Nothing about that
        /// requires self-complementarity, and tandem repeats of this length
        /// are common in real sequence while palindromes of it are not — so
        /// this tie-break is exercised far more often than a palindrome-only
        /// reading would suggest. The 5'-abutting span wins, deterministically.
        #[test]
        fn a_tandem_repeat_makes_both_candidates_match_and_the_five_prime_span_wins() {
            // A 12-mer that is NOT self-complementary, laid down twice.
            let unit = "GGCTAACGTTCA";
            assert_ne!(revcomp(unit), unit, "the unit must not be palindromic");
            // Pad either side so the coincidence floor sees a normal window.
            let filler = "ACGTTGCAACGT".repeat(20);
            let contig = format!("{filler}{unit}{unit}{filler}");
            let gap_idx = filler.len() + unit.len();
            let payload = revcomp(unit);

            let span =
                inverted_adjacent_copy_span(contig.as_bytes(), gap_idx as u64, payload.as_bytes())
                    .expect("a 12-base inverted copy beside its source must be derived");

            // 1-based inclusive, the span ENDING at the junction.
            assert_eq!(
                span,
                ((gap_idx - unit.len() + 1) as u64, gap_idx as u64),
                "with both candidates matching, the 5'-abutting span must win"
            );
        }

        /// A single base is outside the concept of an inversion
        /// (`DNA/inversion.md` requires "more than one nucleotide").
        #[test]
        fn a_single_base_payload_is_declined() {
            let contig = "ACGTTGCAACGT".repeat(30);
            let gap_idx = 120usize;
            let span = &contig[gap_idx - 1..gap_idx];
            let payload = revcomp(span);
            assert!(inverted_adjacent_copy_span(
                contig.as_bytes(),
                gap_idx as u64,
                payload.as_bytes()
            )
            .is_none());
        }

        /// A self-complementary payload is a duplication, and
        /// `general.md:55` ranks duplication above insertion — so this
        /// declines and leaves the case to `insertion_is_duplication`.
        /// The guard `continue`s rather than returning, so the 3' candidate
        /// still gets its turn.
        #[test]
        fn a_self_complementary_payload_is_left_to_the_duplication_rule() {
            let contig = "ACGTTGCAACGT".repeat(30);
            let gap_idx = 120usize;
            let span = &contig[gap_idx - 2..gap_idx];
            // `GC`/`AT`-style two-mers are their own reverse complement.
            if revcomp(span) == span {
                assert!(inverted_adjacent_copy_span(
                    contig.as_bytes(),
                    gap_idx as u64,
                    span.as_bytes()
                )
                .is_none());
            }
        }

        /// Fail-safes: a window with no ACGT, and a payload carrying a
        /// non-ACGT base, both yield probability `1.0` and therefore decline.
        #[test]
        fn non_acgt_input_fails_safe_to_declining() {
            assert_eq!(chance_match_probability(b"NNNNNNNN", b"ACGT"), 1.0);
            assert_eq!(chance_match_probability(b"ACGTACGT", b"ACNT"), 1.0);
        }
    }
}
