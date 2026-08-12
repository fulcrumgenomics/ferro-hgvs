//! Apply a variant to its reference, and derive an encoding-invariant key from
//! the result (#1159).
//!
//! # Why this exists
//!
//! [`crate::spdi::hgvs_to_spdi`] converts one description into one SPDI triple.
//! That is a *transliteration*: it preserves however the caller chose to
//! partition the change, so two descriptions of the same edit convert to
//! different triples. `TEMPLATE:g.8_14delinsGATTA` and
//! `TEMPLATE:g.[8A>G;9G>A;11C>T;13_14del]` denote the same resulting sequence and
//! transliterate to one triple and four.
//!
//! Callers who want a key that is equal **iff two descriptions denote the same
//! edit** therefore cannot use it, and #1159 records what they do instead:
//! normalize both and compare the strings.
//!
//! That was failing when #1159 was filed, for the reasons #1157 and #1158
//! recorded. **Both are now closed and fixed on `main`** — #1157 by #1160/#1161
//! and #1158 by #1237/#1341 — and the normalizer is measurably confluent on
//! #1157's own shape today. So the argument for this module is no longer "the
//! normalize-and-compare route is broken"; it is that the route couples a key to
//! the normalizer, and every future confluence fix then churns stored keys.
//! `canonical_spdi` has no normalizer dependency at all, which is the property
//! worth having.
//!
//! [`canonical_spdi`] answers the question directly instead: apply every member
//! to the reference, then read the difference back out of the *resulting bases*.
//! Partitioning cannot survive that, because it is not represented in the
//! sequence.
//!
//! # What the key does and does not guarantee
//!
//! **Guaranteed:** two variants on one accession with the same resulting
//! sequence produce the same `SpdiVariant`, whatever their spelling, member count
//! or member order. That is the property #1159 asks for and what the tests pin.
//!
//! **Also guaranteed, as of the 3' padding:** the change is rolled to its
//! 3'-most equivalent position. A blunt trim over a window wide enough to
//! contain the roll *is* the maximal 3' shift, so `g.3_5del` and `g.4_6del` —
//! the same deletion spelled one base apart in a tract — now key identically.
//! Before the window was padded they did not, and that was the defect: the
//! window was exactly the members' own span, leaving the trim no room to move.
//!
//! **Also guaranteed, as of the alphabet fold:** the key does not depend on the
//! provider's spelling conventions. Its bases are emitted uppercase, and `U` is
//! folded to `T` on the `r.` axis, so a soft-masked FASTA and an uppercase one
//! key the same locus identically. Stated here because it is a guarantee a
//! *caller* has, while the mechanism is a private helper (`key_alphabet`) that
//! rustdoc does not publish — read that helper for the one boundary the fold
//! does not reach, a provider serving a uracil-spelled transcript.
//!
//! **Still not claimed:** byte-compatibility with the SPDI specification's own
//! canonical form. The form targeted here is **3'-maximal**, per `general.md:34`;
//! NCBI's is extended in *both* directions over the repeat, so the two differ on
//! a rolled indel — `9:G:` here against `8:GG:G` there. Equality across spellings
//! is the contract, not agreement with another implementation.
//!
//! **Not claimed either:** that a no-op carries a canonical position. `g.3A>A`
//! and `g.5T>T` both change nothing but key as `3::` and `5::`. Giving that one
//! answer means choosing one, and nothing here needs it.

use crate::error::FerroError;
use crate::hgvs::variant::HgvsVariant;
use crate::reference::ReferenceProvider;
use crate::spdi::convert::{apply_alphabet, AlphabetMode};
use crate::spdi::SpdiVariant;

/// A variant applied to its reference.
///
/// The window is the union of every member's footprint, so `reference` and
/// `resulting` are directly comparable strings of the same locus before and
/// after. `resulting` is what #1159 calls the ground truth for equivalence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AppliedVariant {
    /// The accession every member acts on.
    pub accession: String,
    /// 0-based start of the window both strings cover.
    pub start: u64,
    /// The reference bases over the window.
    pub reference: String,
    /// The same window after every member has been applied.
    pub resulting: String,
}

/// Largest window, in bases, that [`apply_to_reference`] will fetch.
///
/// A variant spanning more than this is declined rather than served: the point
/// of the primitive is a comparable local window, and a caller asking about a
/// megabase deletion wants a different tool. Mirrors the equivalence checker's
/// own bound, for the same reason.
pub const MAX_APPLY_WINDOW: u64 = 100_000;

/// Longest repeat tract [`canonical_spdi`] will roll a change 3' through.
///
/// Named separately from [`MAX_APPLY_WINDOW`] because it bounds a different
/// thing, and saying so is the point. `MAX_APPLY_WINDOW` answers "how much
/// reference will I read for this variant", which the caller can predict from
/// the description alone. This one answers "how long a tract will I tolerate
/// *around* it", which the caller cannot predict — it is a property of the
/// variant's neighbourhood in the reference.
///
/// That makes a previously total function partial on a new axis, which is the
/// strongest argument against extending the window at all. The mitigation is
/// this: its own name, a generous size, and its own error variant, so the
/// failure reads as "this variant sits in a repeat tract longer than N" rather
/// than as a generic decline the caller must guess the cause of.
///
/// 32 KB is far above any tract a key is meaningful for — the longest known
/// pathogenic expansions run to a few tens of kilobases, and a change that rolls
/// that far is a structural event rather than the local edit this key describes.
pub const MAX_SHIFT_TRACT: u64 = 32_768;

/// Apply every member of `variant` to `provider`'s reference and return the
/// before/after window.
///
/// # Errors
///
/// Declines — rather than guessing — when a single resulting sequence is not
/// well defined: a non-cis allele (trans / mosaic / chimeric / unknown phase, all
/// of which describe more than one molecule), a null or unknown allele, members
/// on different accessions, members whose edits SPDI cannot represent, members
/// that overlap (applying them depends on order, so there is no one answer), a
/// stated deletion that disagrees with the reference, or a window wider than
/// [`MAX_APPLY_WINDOW`].
pub fn apply_to_reference<P: ReferenceProvider + ?Sized>(
    variant: &HgvsVariant,
    provider: &P,
) -> Result<AppliedVariant, FerroError> {
    apply_to_reference_padded(variant, provider, 0).map(|(applied, _)| applied)
}

/// As [`apply_to_reference`], but fetching `pad_3prime` extra reference bases
/// past the members' own span, and reporting whether the full pad was served.
///
/// The pad is what gives [`canonical_spdi`]'s blunt trim room to roll a change
/// 3'. **Only the 3' side is padded, and that is not an oversight.** Prepending
/// `n` bases of 5' flank adds exactly `n` to the common prefix, so
/// `position = start + prefix` is unchanged — left context cannot move the key,
/// and fetching it would be cost with no effect.
///
/// The returned flag says the window is **final**: it ends where the pad asked,
/// or where the sequence itself ends. False means the contig ran out first, so
/// a change still touching the 3' edge has genuinely nowhere left to roll.
///
/// A provider that can serve neither the wider span nor a length to bound it by
/// is a third case, and it declines rather than reporting either — the caller
/// could not otherwise tell "the change stops here" from "I could not read far
/// enough to find out", and would key off a window that may have cut the roll
/// short.
pub(crate) fn apply_to_reference_padded<P: ReferenceProvider + ?Sized>(
    variant: &HgvsVariant,
    provider: &P,
    pad_3prime: u64,
) -> Result<(AppliedVariant, bool), FerroError> {
    // `UnsupportedVariant` carries only a `variant_type` string, so the variant
    // and the reason both go in it — a caller that cannot apply a description
    // needs to know which one and why, and a bare type name would say neither.
    let decline = |reason: &str| FerroError::UnsupportedVariant {
        variant_type: format!("{variant}: cannot apply to reference — {reason}"),
    };

    let (accession, triples) = variant_edit_triples(variant, provider).ok_or_else(|| {
        decline(
            "no single resulting sequence is defined for it — it is a multi-molecule \
             or null allele, spans more than one accession, or carries an edit SPDI \
             cannot represent",
        )
    })?;

    let (mut start, mut end) = (u64::MAX, u64::MIN);
    for triple in &triples {
        start = start.min(triple.position);
        // Checked: `position` is provider-derived and `deletion` caller-derived,
        // so nothing upstream bounds their sum. A wrap here would silently
        // shrink the window and key the variant off the wrong bases.
        let triple_end = triple
            .position
            .checked_add(triple.deletion.len() as u64)
            .ok_or_else(|| decline("its span overflows the coordinate space"))?;
        end = end.max(triple_end);
    }
    if start > end {
        return Err(decline("its members name no reference span"));
    }
    if end - start > MAX_APPLY_WINDOW {
        return Err(decline(&format!(
            "it spans {} bases, more than the {MAX_APPLY_WINDOW}-base limit",
            end - start
        )));
    }

    // Clamp the pad to the contig when its length is knowable. A failed lookup
    // is not fatal: `get_sequence_length` carries a trait default that always
    // errors (`provider.rs:424`), so an out-of-tree provider that serves bases
    // perfectly well may not implement it, and refusing there would regress a
    // caller that works today. The padded fetch is simply attempted and, if the
    // provider cannot serve it, the unpadded window is used instead.
    let requested_end = end.saturating_add(pad_3prime);
    let known_length = provider.get_sequence_length(&accession).ok();
    let wanted_end = known_length.map_or(requested_end, |length| requested_end.min(length));

    let reference = fetch_window(provider, &accession, start, wanted_end).ok_or_else(|| {
        // Falling back to the unpadded window here would be the wrong kind of
        // safe. The caller cannot then tell "the change genuinely stops here"
        // from "I could not read far enough to find out", and would key the
        // variant off a window that may have cut the roll short — two spellings
        // either side of the cut keying differently, which is the exact
        // non-convergence this padding removes. Declining says so instead.
        if wanted_end > end {
            decline(
                "its reference window could not be widened far enough to settle where the \
                 change ends — the provider served neither the wider span nor a length to \
                 bound it by",
            )
        } else {
            decline("its reference window could not be read")
        }
    })?;

    // The roll is complete when the window ends where we asked, or where the
    // sequence itself ends. Only a *known* length makes the second case
    // trustworthy: without one, a short window is indistinguishable from a
    // provider that simply stopped, which is why the fetch above declines
    // rather than falling back.
    let window_is_final =
        wanted_end == requested_end || known_length.is_some_and(|length| wanted_end == length);

    let resulting = apply_triples(&reference, start, &triples).ok_or_else(|| {
        decline(
            "its members overlap, or a stated reference base disagrees with the \
             reference",
        )
    })?;

    Ok((
        AppliedVariant {
            accession,
            start,
            reference,
            resulting,
        },
        window_is_final,
    ))
}

/// An encoding-invariant SPDI key for `variant`, derived from the bases it
/// results in rather than from how it was written (#1159).
///
/// Two descriptions on one accession denoting the same resulting sequence give
/// the same triple, whatever their spelling, member count or member order. See
/// the module docs for what this does and does not claim.
///
/// # Errors
///
/// As [`apply_to_reference`].
pub fn canonical_spdi<P: ReferenceProvider + ?Sized>(
    variant: &HgvsVariant,
    provider: &P,
) -> Result<SpdiVariant, FerroError> {
    /// First pad tried once a change is found to reach the window's 3' edge.
    /// Doubling from here makes a homopolymer cost `O(log tract)` fetches, and
    /// covers every non-repeat case in one.
    const FIRST_PAD: u64 = 64;

    let alphabet = key_alphabet(variant);
    let mut pad = 0u64;
    loop {
        let (applied, window_is_final) = apply_to_reference_padded(variant, provider, pad)?;
        let (offset, deletion, insertion) =
            trim_common_flanks(applied.reference.as_bytes(), applied.resulting.as_bytes());

        // A description that changes nothing has no block to roll, and the
        // 3'-edge test below is vacuously true for it — every base matches, so
        // the trimmed block sits at the window's end however wide the window is.
        // Without this the loop would widen to the accession's end chasing an
        // edge it can never leave, then decline a variant that is merely inert.
        //
        // It does **not** make two no-ops agree, and measuring rather than
        // assuming is what showed that: `g.3A>A` keys as `3::` and `g.5T>T` as
        // `5::`, each at its own position. That is a real residual — a key
        // meaning "nothing changed" arguably should not carry a position at all
        // — but it is pre-existing, orthogonal to the 3'-shift this change is
        // about, and giving it an answer means choosing one, so it is left
        // stated rather than silently picked.
        let is_no_op = deletion.is_empty() && insertion.is_empty();

        // The blunt trim *is* the maximal 3' shift, but only over a window wide
        // enough to contain the roll. While the trimmed block still runs to the
        // window's 3' edge, the shift may have been cut off by the window rather
        // than by the sequence, so widen and ask again.
        let reaches_edge = offset + deletion.len() == applied.reference.len();

        if is_no_op || !reaches_edge || !window_is_final {
            // Folded through `apply_alphabet`, not emitted raw: see
            // `key_alphabet` for why a raw window makes the key depend on the
            // provider's case convention rather than on the bases.
            return Ok(SpdiVariant {
                sequence: applied.accession,
                position: applied.start + offset as u64,
                deletion: apply_alphabet(&String::from_utf8_lossy(deletion), alphabet),
                insertion: apply_alphabet(&String::from_utf8_lossy(insertion), alphabet),
            });
        }

        // Decline rather than answer over a truncated window. Returning the
        // clamped key would not be a smaller answer, it would be a *different*
        // one: two spellings either side of the cap key differently, which is
        // exactly the non-convergence this function exists to remove.
        if pad >= MAX_SHIFT_TRACT {
            return Err(FerroError::UnsupportedVariant {
                variant_type: format!(
                    "{variant}: cannot derive a stable key — it sits in a repeat tract \
                     still running past {MAX_SHIFT_TRACT} bases 3' of it, so how far the \
                     change shifts depends on how much reference is read"
                ),
            });
        }
        pad = if pad == 0 {
            FIRST_PAD
        } else {
            pad.saturating_mul(2)
        };
    }
}

/// The alphabet convention [`canonical_spdi`] renders a key's bases in.
///
/// **A key must not depend on the provider's spelling conventions**, and without
/// this fold it did. The two halves of a key come from two places: the inserted
/// payload is carried in by the triples, which [`crate::spdi::hgvs_to_spdi`] has
/// already run through [`apply_alphabet`] with the member's own axis, while the
/// deleted bases — and any inserted bases a `dup` or `inv` reads *out of* the
/// reference — come from the fetched window verbatim. So one key could hold two
/// conventions at once, and two providers serving the same sequence could key it
/// differently.
///
/// Measured before this fold, on a soft-masked copy of the module's own fixture
/// contig: `g.3_7delinsGGCTA` keyed as `2:ATTAC:GGCTA` against the uppercase
/// contig and as `2:attac:GGCTA` against the lowercase one. Real genomic FASTAs
/// lowercase repeat tracts, so that is the common case rather than a corner:
/// `SpdiKey`'s whole contract is that equal keys mean equal bases, and case
/// carries no biological meaning. `trim_common_flanks` was already
/// case-insensitive, so only the *emitted* strings were affected — the position
/// and the block boundaries were right all along, which is what made this
/// invisible to every test using an uppercase fixture.
/// `a_soft_masked_reference_keys_the_same_as_an_uppercase_one` pins it.
///
/// [`AlphabetMode::Rna`] additionally rewrites `U` to `T`, the same convention
/// the `r.` axis's triples already carry. **That reconciles a uracil-spelled
/// provider only where the U-bearing reference bases reach the emitted block
/// through [`trim_common_flanks`], and the boundary is worth stating rather than
/// generalising from** — the fold runs *after* both comparisons that decide the
/// block, and each is `eq_ignore_ascii_case` only, which does not relate `U` to
/// `T`. Measured on one 40-base transcript served both ways:
///
/// ```text
///                        T-spelled provider   U-spelled provider
///   r.[3a>g;7c>a]        2:ATTAC:GTTAA        2:ATTAC:GTTAA   agrees — the fold
///   r.3a>g               2:A:G                2:A:G           agrees
///   r.14dup              14::T                14::T           agrees
///   r.13_14insu          14::T                13::T           position shifts
///   r.3_5del             3:TTA:               declined        apply_triples
///   r.3_7delinsggcua     2:ATTAC:GGCTA        declined        apply_triples
/// ```
///
/// The two failing shapes are bounded by named code, not by luck. A **stated**
/// deletion is validated against the reference in [`apply_triples`] before any
/// fold, so `T` against `U` reads as a ref-mismatch and the variant is declined
/// outright. An insertion or `dup` rolling 3' has its common-prefix scan stopped
/// at the first `U`/`T` disagreement, so the roll is cut short and the *position*
/// moves. Neither is reachable from a RefSeq-spelled provider — `refseq.md` and
/// [`apply_alphabet`] both record that transcript sequences are stored as DNA —
/// which is why this is left stated rather than fixed here: making the two
/// comparisons alphabet-aware is a change to the `n.`/`c.` paths as well.
/// `the_alphabet_fold_reconciles_case_and_bounds_the_uracil_case` pins every row
/// above.
///
/// Note this deliberately does **not** fold [`AppliedVariant`], whose whole
/// point is to hand back the reference window as it was served; a caller
/// inspecting soft-masking there is asking a legitimate question. Nothing
/// compares a folded key against an unfolded window: [`EquivalenceChecker`]
/// compares two [`apply_triples`] results with each other, and
/// `examples/dump_normalized_corpus.rs` compares `reference` against the frame
/// it generated — neither reads a [`SpdiVariant`] field.
///
/// [`EquivalenceChecker`]: crate::equivalence::EquivalenceChecker
fn key_alphabet(variant: &HgvsVariant) -> AlphabetMode {
    fn is_rna(variant: &HgvsVariant) -> bool {
        matches!(variant, HgvsVariant::Rna(_))
    }
    match variant {
        // `any` rather than `all`, and the distinction is **reachable** — not, as
        // this comment once claimed, ruled out by the members sharing one
        // accession. `[NR_X.1:r.3a>g;NR_X.1:n.7C>A]` parses, is one cis allele on
        // one accession, and keys; measured on a U-spelled provider it keys
        // `2:ATTAC:GTTAA` under `any` and `2:AUUAC:GUUAA` under `all`. So `any`
        // is what stops one key from carrying two alphabets, which is the whole
        // point of the fold. Pinned by
        // `a_mixed_axis_cis_allele_folds_on_any_rna_member`.
        HgvsVariant::Allele(allele) if allele.variants.iter().any(is_rna) => AlphabetMode::Rna,
        single if is_rna(single) => AlphabetMode::Rna,
        _ => AlphabetMode::Dna,
    }
}

/// Strip the bases the two windows share at each end, returning the offset of the
/// remaining block and the two differing slices.
///
/// This is what makes the key window-independent: a wider window differs from a
/// narrower one only by flanking bases that are identical on both sides, and
/// those are exactly what is removed. Case-insensitive, because reference FASTAs
/// are often soft-masked and case carries no biological meaning — the same
/// reasoning `apply_triples` already applies to stated deletions.
fn trim_common_flanks<'a>(reference: &'a [u8], resulting: &'a [u8]) -> (usize, &'a [u8], &'a [u8]) {
    let max_prefix = reference.len().min(resulting.len());
    let mut prefix = 0;
    while prefix < max_prefix && reference[prefix].eq_ignore_ascii_case(&resulting[prefix]) {
        prefix += 1;
    }
    let mut suffix = 0;
    while suffix < max_prefix - prefix
        && reference[reference.len() - 1 - suffix]
            .eq_ignore_ascii_case(&resulting[resulting.len() - 1 - suffix])
    {
        suffix += 1;
    }
    (
        prefix,
        &reference[prefix..reference.len() - suffix],
        &resulting[prefix..resulting.len() - suffix],
    )
}

/// Apply SPDI triples to `reference` — the bases spanning the interbase
/// interval that begins at `win_start` — and return the edited sequence.
///
/// Triples are applied from the 3' end (descending position) so that an earlier
/// splice never shifts the coordinates of a later one. Each triple's stated
/// deleted bases are validated against the actual reference bases at that span;
/// if they disagree (a ref-mismatched input — e.g. `c.5A>G` where the reference
/// base is not `A`), we cannot faithfully reconstruct the edit, so we decline
/// (`None`) rather than assert a sequence equivalence we cannot trust. Also
/// returns `None` if any triple falls outside the window (a defensive guard;
/// callers build the window to cover every triple), or if two triples overlap
/// (see [`triples_are_disjoint`]).
pub(crate) fn apply_triples(
    reference: &str,
    win_start: u64,
    triples: &[SpdiVariant],
) -> Option<String> {
    let ref_bytes = reference.as_bytes();
    let mut bytes = ref_bytes.to_vec();
    let mut ordered: Vec<&SpdiVariant> = triples.iter().collect();
    ordered.sort_by_key(|t| std::cmp::Reverse(t.position));
    if !triples_are_disjoint(&ordered) {
        return None;
    }
    for t in ordered {
        let rel = t.position.checked_sub(win_start)? as usize;
        let end = rel.checked_add(t.deletion.len())?;
        if end > ref_bytes.len() {
            return None;
        }
        // Validate the stated deletion against the original reference span.
        // (Checked against `ref_bytes`, not the mutated `bytes`: descending
        // order means every already-applied splice sits strictly 3' of `rel`,
        // so this span is untouched either way.)
        if !ref_bytes[rel..end].eq_ignore_ascii_case(t.deletion.as_bytes()) {
            return None;
        }
        // The splice targets the mutated buffer, whose length no longer matches
        // the reference once a length-changing edit has been applied. The
        // disjointness guard above already makes `end <= bytes.len()` hold;
        // bound it explicitly anyway so a future change cannot turn a logic
        // slip back into an out-of-bounds panic (#1244).
        if end > bytes.len() {
            return None;
        }
        bytes.splice(rel..end, t.insertion.bytes());
    }
    String::from_utf8(bytes).ok()
}

/// Whether no two of `ordered` claim the same reference base.
///
/// `ordered` must be sorted by descending position, as [`apply_triples`] leaves
/// it. That descending application order is what lets each stated deletion be
/// validated against the pristine reference: every already-applied splice sits
/// strictly 3' of the next one, so the span about to be read is untouched. The
/// argument holds only while the triples are disjoint — overlapping ones both
/// invalidate that validation and can index past the end of the shrinking
/// buffer, which is the out-of-bounds panic of #1244.
///
/// Declining is also the honest answer semantically: an allele whose members
/// claim the same base has no single well-defined resulting sequence, so there
/// is nothing to compare. The caller uses the comparison only to *upgrade* a
/// `NotEquivalent` verdict, so a decline never invents an equivalence.
///
/// Two triples that merely abut are disjoint, and so is any number of pure
/// insertions at one interbase position — an insertion deletes nothing and
/// therefore claims no base.
fn triples_are_disjoint(ordered: &[&SpdiVariant]) -> bool {
    // Walk 5' -> 3' (the reverse of `ordered`) carrying the furthest 3' reach
    // of every triple seen so far. Comparing against the running maximum rather
    // than the immediate predecessor is what makes this complete: one long
    // triple can span several shorter ones that do not touch each other.
    let mut reach: Option<u64> = None;
    for t in ordered.iter().rev() {
        let Some(end) = t.position.checked_add(t.deletion.len() as u64) else {
            return false;
        };
        if reach.is_some_and(|r| r > t.position) {
            return false;
        }
        reach = Some(reach.map_or(end, |r| r.max(end)));
    }
    true
}

/// The SPDI triples that make up `variant`'s resulting sequence, and the single
/// accession they act on.
///
/// `None` when a single resulting sequence is undefined or cannot be derived —
/// see [`apply_to_reference`]'s error list, which this decides. Use
/// [`variant_edit_triples_reason`] when the *cause* matters.
pub(crate) fn variant_edit_triples<P: ReferenceProvider + ?Sized>(
    variant: &HgvsVariant,
    provider: &P,
) -> Option<(String, Vec<SpdiVariant>)> {
    variant_edit_triples_reason(variant, provider).ok()
}

/// Why a description yields no SPDI triple set.
///
/// The split exists because the two are **not equally severe**, and a caller
/// that cannot tell them apart draws the wrong conclusion from a decline. See
/// [`compare_denoted_sequences`], where reading every decline as
/// [`Self::SelfContradictory`] produced 328 false alarms across the test suite
/// in a single run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NoTriples {
    /// The applier could not transliterate it: a shape with no single resulting
    /// sequence by construction (a trans allele, a null allele, members on
    /// different accessions), an edit SPDI cannot carry, a position it cannot
    /// resolve, or reference data the provider does not hold.
    ///
    /// **A limit of the applier or of the shape, never a verdict on the
    /// description.** The last cause is the common one and the one that misleads:
    /// `g.1000delC` states its deleted base and so converts with no provider at
    /// all, while the `g.1000del` a normalizer emits for it must read the
    /// reference — so on a provider holding no bases the input converts and the
    /// output does not, through no fault of the normalization.
    Untransliterable,
    /// The description states an edit set with no defined result **whatever the
    /// reference says** — members claiming the same base, or two insertions at
    /// one interbase with no stated order.
    ///
    /// A real fault in a description that claims to denote one sequence.
    SelfContradictory,
}

/// [`variant_edit_triples`], reporting why it declined.
pub(crate) fn variant_edit_triples_reason<P: ReferenceProvider + ?Sized>(
    variant: &HgvsVariant,
    provider: &P,
) -> Result<(String, Vec<SpdiVariant>), NoTriples> {
    use crate::hgvs::variant::AllelePhase;

    let members: Vec<&HgvsVariant> = match variant {
        HgvsVariant::Allele(allele) => {
            // A single resulting sequence is only well-defined for a cis allele —
            // every edit applied to the same molecule.
            if allele.phase != AllelePhase::Cis {
                return Err(NoTriples::Untransliterable);
            }
            allele.variants.iter().collect()
        }
        HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => {
            return Err(NoTriples::Untransliterable)
        }
        single => vec![single],
    };
    if members.is_empty() {
        return Err(NoTriples::Untransliterable);
    }

    let mut accession: Option<String> = None;
    let mut triples = Vec::with_capacity(members.len());
    for member in members {
        let spdi =
            crate::spdi::hgvs_to_spdi(member, provider).map_err(|_| NoTriples::Untransliterable)?;
        match &accession {
            None => accession = Some(spdi.sequence.clone()),
            Some(acc) if *acc != spdi.sequence => return Err(NoTriples::Untransliterable),
            Some(_) => {}
        }
        triples.push(spdi);
    }

    // Two insertions at one interbase have no order between them, and this path
    // *publishes* a key — so it must decline rather than let the application
    // order pick a winner. Measured before this guard, on the fixture contig:
    //
    //     g.[5_6insA;5_6insC]  ->  5::CA
    //     g.[5_6insC;5_6insA]  ->  8::CA
    //
    // One variant, two keys, decided by nothing the description states. HGVS
    // spells "insert both, in this order" as a single ordered compound payload
    // (`ins[A;C]`, general.md:79), so a caller who means an order has a way to
    // say it and a caller who wrote two members has not said one.
    //
    // **Deliberately here and not in `triples_are_disjoint`.** That predicate is
    // shared with `EquivalenceChecker`, where the permissive reading is correct:
    // a decline there only forgoes upgrading a `NotEquivalent` verdict, so it can
    // never invent an equivalence, and two pinned tests depend on it. One
    // predicate cannot serve both a checker that may guess and a key that may
    // not.
    let mut zero_width: Vec<u64> = triples
        .iter()
        .filter(|t| t.deletion.is_empty())
        .map(|t| t.position)
        .collect();
    zero_width.sort_unstable();
    if zero_width.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(NoTriples::SelfContradictory);
    }

    accession
        .map(|acc| (acc, triples))
        .ok_or(NoTriples::Untransliterable)
}

/// Why two descriptions' denoted sequences could not be compared (#1615).
///
/// Every variant here is a limit of the *comparison*, not a verdict on either
/// description. They are enumerated rather than collapsed to one "skip" because
/// a skip that reads as a pass is the failure mode the denoted-sequence oracle
/// exists to remove: a caller that cannot say which of these it hit cannot tell
/// a clean run from a run that compared nothing.
///
/// Marked `#[non_exhaustive]` so a new decline reason is additive rather than a
/// breaking change, matching `SpdiParseError` and `ConversionError` in this
/// module. The set is demonstrably still growing:
/// [`Self::UnresolvableSpecialPosition`] exists only until `hgvs_to_spdi` stops
/// resolving `pter` silently — which as of #1643 is true of the **transcript**
/// axis only, the genomic half now being refused outright — and #1618/#1619 are
/// two more disagreements in flight.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum NotComparable {
    /// The **input** denotes no single sequence, so there is no baseline to
    /// compare the output against.
    ///
    /// A multi-molecule or null allele, members on different accessions, an edit
    /// SPDI cannot represent, a self-contradictory member set, or a stated
    /// reference base that disagrees with the reference (the `REFSEQ_MISMATCH`
    /// inputs normalization exists to correct). Blaming normalization for any of
    /// those would make a fire mean two different things — the discipline
    /// `assert_reparseable` and `assert_in_bounds` already apply.
    InputDenotesNoSequence,
    /// The **output** cannot be transliterated, though it is not
    /// self-contradictory. A limit of the applier or of the provider, not of the
    /// description — see [`NoTriples::Untransliterable`], which is exactly the
    /// asymmetry that makes this its own verdict rather than a fire.
    OutputUntransliterable,
    /// One side is an allele wrapped in the predicted marker `[(…)]`, whose
    /// members are uncertain by construction.
    UncertainAllele,
    /// One side names `pter`/`qter`/`cen`.
    ///
    /// Those carry no numeric coordinate — `GenomePos::pter()` and
    /// `CdsPos::pter()` both set `base: 0` — and comparing against a position
    /// that is not the one meant would report every correct `pter`
    /// normalization as a corruption. So `names_a_special_position` declines
    /// the pair on **both** axes, before either side is transliterated.
    ///
    /// **Only the transcript half is still resolved silently.** #1643 closed
    /// the genomic one: `hgvs_to_spdi` now refuses a `g.`/`m.`/`o.` special
    /// position with `ConversionError::InvalidPosition` (see
    /// `convert::reject_unresolvable_genomic_position`), so on that axis this
    /// decline has become belt-and-braces rather than the only thing standing
    /// between the oracle and a false fault. On the transcript axis `CdsPos`
    /// still flattens onto `base: 0` and `hgvs_to_spdi` still reads it as a
    /// coordinate — measured: `NM_003002.2:c.pterdel` transliterates to a
    /// deletion of the sequence's LAST base — which is the half this variant
    /// now exists for, and it retires when that half is fixed.
    UnresolvableSpecialPosition,
    /// The two descriptions name different accessions, so their sequences are
    /// not comparable at all. Normalization can do this legitimately (the #785
    /// transcript-version substitution).
    AccessionChanged,
    /// The union of the two descriptions' spans is wider than
    /// [`MAX_APPLY_WINDOW`].
    WindowTooWide,
    /// The provider could not serve the union window.
    ReferenceUnreadable,
}

impl std::fmt::Display for NotComparable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let reason = match self {
            Self::InputDenotesNoSequence => "the input denotes no single sequence",
            Self::OutputUntransliterable => "the output cannot be transliterated to SPDI",
            Self::UncertainAllele => "one side is a predicted `[(…)]` allele",
            Self::UnresolvableSpecialPosition => "one side names pter/qter/cen",
            Self::AccessionChanged => "the two descriptions name different accessions",
            Self::WindowTooWide => "their union spans more than MAX_APPLY_WINDOW bases",
            Self::ReferenceUnreadable => "the provider could not serve the union window",
        };
        f.write_str(reason)
    }
}

/// The outcome of comparing the sequences two descriptions denote (#1615).
///
/// Deliberately **not** `#[non_exhaustive]`, unlike its
/// [`NotComparable`] payload: `issue_1615_denoted_sequence_oracle::
/// the_oracle_fires_on_every_recorded_defect` matches this enum exhaustively
/// from outside the crate, which is what makes a new verdict fail to compile
/// rather than be swallowed by a wildcard arm. A new *decline reason* is
/// additive and belongs in `NotComparable`; a new *verdict* is a change to what
/// the oracle can say, and every caller should have to answer for it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DenotedSequenceComparison {
    /// Both descriptions apply, and to the same bases.
    Agree,
    /// Both apply, and to **different** bases.
    Differ {
        /// The accession both descriptions act on.
        accession: String,
        /// 0-based start of the window both strings cover.
        start: u64,
        /// The reference bases over that window.
        reference: String,
        /// The window after applying the first description.
        from_input: String,
        /// The window after applying the second.
        from_output: String,
    },
    /// The **output** is self-contradictory although the input denotes a
    /// sequence: its members claim the same base, or two of its insertions share
    /// one interbase with no stated order.
    ///
    /// Deliberately *not* a [`NotComparable`] variant. #1281's `g.[1del;1del]`
    /// denotes nothing at all, which is strictly worse than denoting the wrong
    /// thing, so folding it in with the skips would hide the more severe defect
    /// behind the milder one. Equally deliberately, it is **narrower** than "the
    /// output produced no triples" — see [`NotComparable::OutputUntransliterable`].
    OutputContradictsItself,
    /// Neither verdict is available, for a stated reason.
    NotComparable(NotComparable),
}

/// Whether `variant` is an allele wrapped in the predicted marker `[(…)]`.
fn is_uncertain_allele(variant: &HgvsVariant) -> bool {
    matches!(variant, HgvsVariant::Allele(allele) if allele.uncertain)
}

/// Whether `variant` names a chromosome-arm or centromere special position.
///
/// Only [`crate::hgvs::location::GenomePos`] and
/// [`crate::hgvs::location::CdsPos`] carry the marker, so only the axes built on
/// those are walked; the rest cannot express one. See
/// [`NotComparable::UnresolvableSpecialPosition`] for why the answer matters.
fn names_a_special_position(variant: &HgvsVariant) -> bool {
    use crate::hgvs::interval::{Interval, UncertainBoundary};
    use crate::hgvs::location::{CdsPos, GenomePos};

    fn boundary_is_special<T>(boundary: &UncertainBoundary<T>, special: fn(&T) -> bool) -> bool {
        match boundary {
            UncertainBoundary::Single(mu) => mu.inner().is_some_and(special),
            UncertainBoundary::Range { start, end } => {
                start.inner().is_some_and(special) || end.inner().is_some_and(special)
            }
        }
    }
    fn interval_is_special<T>(interval: &Interval<T>, special: fn(&T) -> bool) -> bool {
        boundary_is_special(&interval.start, special) || boundary_is_special(&interval.end, special)
    }
    fn genomic(interval: &Interval<GenomePos>) -> bool {
        interval_is_special(interval, |p| p.special.is_some())
    }

    match variant {
        HgvsVariant::Allele(allele) => allele.variants.iter().any(names_a_special_position),
        HgvsVariant::Genome(v) => genomic(&v.loc_edit.location),
        HgvsVariant::Mt(v) => genomic(&v.loc_edit.location),
        HgvsVariant::Circular(v) => genomic(&v.loc_edit.location),
        HgvsVariant::Cds(v) => {
            interval_is_special(&v.loc_edit.location, |p: &CdsPos| p.special.is_some())
        }
        _ => false,
    }
}

/// Compare the sequences `input` and `output` denote against one reference
/// window (#1615).
///
/// This is the primitive behind the denoted-sequence seam oracle, and it is
/// **independent of the normalizer**: it reaches the bases through
/// [`variant_edit_triples_reason`] and [`apply_triples`], the same SPDI-splicing
/// walk [`apply_to_reference`] uses, so nothing here can agree with
/// normalization merely because normalization produced it. (`EquivalenceChecker`
/// cannot serve this role — it normalizes both sides, which is circular.)
///
/// Both descriptions are applied over the **union** of their two spans, one
/// fetch, so a 3'-shift is compared where it belongs: `g.3_4del` and `g.7_8del`
/// in a tract denote the same bases over any window containing both, and a
/// per-description window would give each its own frame and make them look
/// different. The comparison is ASCII-case-insensitive for the reason
/// [`trim_common_flanks`] gives — reference FASTAs are often soft-masked and
/// case carries no biological meaning, while one side's payload may be
/// reference-derived (a `dup`) and the other's a literal.
///
/// # Declining is the default; only two outcomes are faults
///
/// [`DenotedSequenceComparison::Differ`] and
/// [`DenotedSequenceComparison::OutputContradictsItself`] are the faults.
/// Everything else the applier cannot do is a [`NotComparable`] with a stated
/// reason, and that asymmetry is load-bearing rather than cautious: an earlier
/// revision reported *any* untranslatable output as a fault and raised **328**
/// false alarms over the test suite, essentially all of them one shape — an
/// input that states its own deleted bases (`g.1000delC`, `g.[1000G>A;1001A>C]`)
/// against an output that does not (`g.1000del`, `g.1000_1001delinsAC`), on a
/// provider holding no reference at that locus. The input converted with no
/// provider and the output could not, and nothing was wrong with either.
pub fn compare_denoted_sequences<P: ReferenceProvider + ?Sized>(
    input: &HgvsVariant,
    output: &HgvsVariant,
    provider: &P,
) -> DenotedSequenceComparison {
    use DenotedSequenceComparison as Outcome;

    if names_a_special_position(input) || names_a_special_position(output) {
        return Outcome::NotComparable(NotComparable::UnresolvableSpecialPosition);
    }
    // An allele wrapped in the predicted marker `[(…)]` states that its members
    // are *uncertain*, so "which bases does it denote" is not a question it
    // answers. `normalize` agrees and leaves such members where the caller put
    // them — `sort_cis_members_by_genomic_order` and the sibling clamps all gate
    // on cis-and-not-uncertain — so its output may legitimately still overlap.
    //
    // Checked here rather than in `variant_edit_triples_reason` on purpose: that
    // function also backs the public `apply_to_reference` and `canonical_spdi`,
    // and narrowing what *they* accept is a separate decision from what this
    // comparison is willing to adjudicate.
    if is_uncertain_allele(input) || is_uncertain_allele(output) {
        return Outcome::NotComparable(NotComparable::UncertainAllele);
    }

    let Ok((accession, input_triples)) = variant_edit_triples_reason(input, provider) else {
        return Outcome::NotComparable(NotComparable::InputDenotesNoSequence);
    };
    let (output_accession, output_triples) = match variant_edit_triples_reason(output, provider) {
        Ok(pair) => pair,
        Err(NoTriples::SelfContradictory) => return Outcome::OutputContradictsItself,
        Err(NoTriples::Untransliterable) => {
            return Outcome::NotComparable(NotComparable::OutputUntransliterable)
        }
    };
    if accession != output_accession {
        return Outcome::NotComparable(NotComparable::AccessionChanged);
    }

    let (mut start, mut end) = (u64::MAX, u64::MIN);
    for triple in input_triples.iter().chain(&output_triples) {
        start = start.min(triple.position);
        // Checked for the same reason `apply_to_reference_padded` checks it: a
        // wrap would silently shrink the window and compare the wrong bases.
        let Some(triple_end) = triple.position.checked_add(triple.deletion.len() as u64) else {
            return Outcome::NotComparable(NotComparable::WindowTooWide);
        };
        end = end.max(triple_end);
    }
    if start > end || end - start > MAX_APPLY_WINDOW {
        return Outcome::NotComparable(NotComparable::WindowTooWide);
    }

    let Some(reference) = fetch_window(provider, &accession, start, end) else {
        return Outcome::NotComparable(NotComparable::ReferenceUnreadable);
    };
    let from_input = match splice_denoted_sequence(&reference, start, &input_triples) {
        Ok(bases) => bases,
        Err(_) => return Outcome::NotComparable(NotComparable::InputDenotesNoSequence),
    };
    let from_output = match splice_denoted_sequence(&reference, start, &output_triples) {
        Ok(bases) => bases,
        // Only one of the two refusals is the description's own fault.
        // Overlapping members claim the same base whatever the reference holds;
        // a stated base that disagrees with the reference is a claim this
        // comparison cannot adjudicate — and is exactly what an input carrying a
        // `REFSEQ_MISMATCH` looks like.
        Err(SpliceFailure::Overlapping) => return Outcome::OutputContradictsItself,
        Err(SpliceFailure::StatedBasesMismatch) => {
            return Outcome::NotComparable(NotComparable::OutputUntransliterable)
        }
    };

    if same_bases(&from_input, &from_output) {
        Outcome::Agree
    } else {
        Outcome::Differ {
            accession,
            start,
            reference,
            from_input,
            from_output,
        }
    }
}

/// Why [`splice_denoted_sequence`] refused a triple set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SpliceFailure {
    /// Two members claim the same reference base, so the result depends on an
    /// order the description does not state. A fault in the description itself.
    Overlapping,
    /// A member's stated deleted bases disagree with the reference.
    StatedBasesMismatch,
}

/// Splice `triples` into `reference` — the bases beginning at interbase
/// `win_start` — and return the resulting sequence.
///
/// Deliberately **not** [`apply_triples`], which is shared with
/// `EquivalenceChecker` and whose [`triples_are_disjoint`] guard reports a pure
/// insertion flush against the 5' end of a deletion as an overlap. That
/// combination is well defined — the payload lands at the junction, then the
/// span is removed — and it is a shape cis alleles reach constantly: measured
/// over one suite run, reading it as an overlap produced **233** false alarms,
/// nearly all of them `g.[…;261_262insA;262_264A[5]]`-shaped. (Its verdict there
/// is also input-order dependent, since the sort breaks position ties by input
/// order.) Changing that predicate is out of scope — two pinned tests depend on
/// its permissiveness in the other direction — so this walk is the one used for
/// the comparison.
///
/// The walk is `tests/it/common/cis_apply_oracle.rs`'s, the applier every
/// sibling-crossing test already rests on: 3'→5' so an applied splice never
/// moves a later one's coordinates, with **longer deletions first** among
/// members at one position so a span-claiming member is always applied before
/// the zero-length one flush against it. Coincident insertions are refused
/// upstream, by [`variant_edit_triples_reason`].
fn splice_denoted_sequence(
    reference: &str,
    win_start: u64,
    triples: &[SpdiVariant],
) -> Result<String, SpliceFailure> {
    let reference = reference.as_bytes();
    let mut ordered: Vec<&SpdiVariant> = triples.iter().collect();
    ordered.sort_by_key(|t| {
        (
            std::cmp::Reverse(t.position),
            std::cmp::Reverse(t.deletion.len()),
        )
    });

    let mut edited = reference.to_vec();
    let mut claimed_from = reference.len();
    for triple in ordered {
        // The caller builds the window from these very triples, so an
        // out-of-window position is unreachable; treat it defensively as a
        // stated-bases problem rather than panicking on the slice.
        let Some(start) = triple
            .position
            .checked_sub(win_start)
            .and_then(|offset| usize::try_from(offset).ok())
        else {
            return Err(SpliceFailure::StatedBasesMismatch);
        };
        let Some(end) = start.checked_add(triple.deletion.len()) else {
            return Err(SpliceFailure::StatedBasesMismatch);
        };
        if end > reference.len() {
            return Err(SpliceFailure::StatedBasesMismatch);
        }
        if end > claimed_from {
            return Err(SpliceFailure::Overlapping);
        }
        if !same_bases_bytes(&reference[start..end], triple.deletion.as_bytes()) {
            return Err(SpliceFailure::StatedBasesMismatch);
        }
        // The splice targets the mutated buffer, whose length no longer matches
        // the reference once a length-changing edit has been applied. The
        // descending walk and the overlap check above already make
        // `end <= edited.len()` hold — every applied splice sits at or past
        // `claimed_from`, so the prefix this indexes into is untouched — but bound
        // it explicitly anyway, so a future change cannot turn a logic slip back
        // into the out-of-bounds panic of #1244.
        if end > edited.len() {
            return Err(SpliceFailure::StatedBasesMismatch);
        }
        edited.splice(start..end, triple.insertion.bytes());
        // Unconditionally, including for a pure insertion — an insertion
        // *interior* to a member 5' of it is an overlap even though it claims no
        // base of its own, and that is precisely the overlap-conflicting allele
        // (`g.[2_3del;2_3insAA]`) this repository declines to canonicalize. The
        // flush case the tie-break above exists for is already safe: the
        // span-claiming member is applied first, so the zero-length one sits
        // exactly *at* `claimed_from` rather than past it.
        claimed_from = start;
    }
    String::from_utf8(edited).map_err(|_| SpliceFailure::StatedBasesMismatch)
}

/// The base `b` is compared as, folding case and the RNA alphabet.
///
/// Case, for the reason [`trim_common_flanks`] gives: reference FASTAs are often
/// soft-masked and case carries no biological meaning. `U` to `T` because the
/// two sides of a comparison need not agree on the alphabet — an `r.` payload is
/// transliterated to RNA by `hgvs_to_spdi` while the reference window it is
/// spliced into is served as DNA, so `r.6_8dupugc` and `r.6_8dup` denote one
/// sequence and would otherwise read as `UGC` against `TGC`.
fn canonical_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'U' => b'T',
        other => other,
    }
}

/// Whether two base strings denote the same sequence under [`canonical_base`].
fn same_bases(left: &str, right: &str) -> bool {
    same_bases_bytes(left.as_bytes(), right.as_bytes())
}

fn same_bases_bytes(left: &[u8], right: &[u8]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(l, r)| canonical_base(*l) == canonical_base(*r))
}

/// Read `[start, end)` from `accession`, or `None` if the provider cannot serve
/// exactly that many bases.
///
/// Tries the genomic accessor first and falls back to the generic one, so a
/// contig and a transcript are both reachable.
pub(crate) fn fetch_window<P: ReferenceProvider + ?Sized>(
    provider: &P,
    accession: &str,
    start: u64,
    end: u64,
) -> Option<String> {
    let bases = provider
        .get_genomic_sequence(accession, start, end)
        .or_else(|_| provider.get_sequence(accession, start, end))
        .ok()?;
    if bases.len() as u64 != end - start {
        return None;
    }
    Some(bases)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;
    use crate::reference::MockProvider;

    /// A 40-base contig whose bases are known, so the expected triples can be
    /// written by hand rather than read back out of the code under test.
    ///
    ///  1-based: 1234567890...
    ///           GGATTACAGGCATTAGCCTGAGGATTACAGGCATTAGCCT
    fn provider() -> MockProvider {
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_KEY.1", "GGATTACAGGCATTAGCCTGAGGATTACAGGCATTAGCCT");
        // A second contig, so "members on different accessions" can be tested as
        // the real thing rather than as an absent-reference failure.
        provider.add_genomic_sequence("NC_OTHER.1", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        provider
    }

    fn key(descriptor: &str) -> SpdiVariant {
        let variant = parse_hgvs(descriptor).expect("fixture must parse");
        canonical_spdi(&variant, &provider())
            .unwrap_or_else(|e| panic!("`{descriptor}` must canonicalize: {e}"))
    }

    /// #1159's whole point: a spanning `delins` and its decomposed cis allele are
    /// the same edit, and must produce the same key.
    ///
    /// `hgvs_to_spdi` cannot do this — it transliterates the partitioning, so it
    /// answers with one triple for the first and four for the second. The key is
    /// derived from the resulting bases, where partitioning is not represented.
    #[test]
    fn a_spanning_delins_and_its_decomposition_share_one_key() {
        // Reference 3..=7 is `ATTAC`; replace it with `GGCTA` piecewise. Every
        // one of the five bases changes, so the decomposition is five members —
        // A>G, T>G, T>C, A>T, C>A.
        let spanning = key("NC_KEY.1:g.3_7delinsGGCTA");
        let decomposed = key("NC_KEY.1:g.[3A>G;4T>G;5T>C;6A>T;7C>A]");
        assert_eq!(
            spanning, decomposed,
            "the same edit written two ways must give one key"
        );
        // And the key really is the minimal changed block, not the input's span.
        assert_eq!(spanning.sequence, "NC_KEY.1");
        assert_eq!(
            (spanning.deletion.as_str(), spanning.insertion.as_str()),
            ("ATTAC", "GGCTA")
        );
    }

    /// Member order must not matter, since a cis allele is a set of edits on one
    /// molecule.
    ///
    /// Note what this test alone cannot show: it uses two substitutions at 3 and
    /// 7, which are disjoint whichever order they are applied in, so
    /// `apply_triples` sorts them to one sequence and the assertion holds even
    /// for an implementation that keyed off application order. The case that
    /// discriminates is two insertions at one interbase, and it is
    /// `coincident_insertions_are_declined_rather_than_ordered` below.
    #[test]
    fn member_order_does_not_change_the_key() {
        assert_eq!(key("NC_KEY.1:g.[3A>G;7C>A]"), key("NC_KEY.1:g.[7C>A;3A>G]"));
    }

    /// The pad doubles until it contains the roll, and a tract that outruns the
    /// cap is declined rather than answered over a truncated window.
    ///
    /// Both halves need a tract longer than `FIRST_PAD`, so this builds its own
    /// contig rather than using the 40-base fixture: 300 `A`s, which forces at
    /// least one doubling (64 -> 128 -> 256 -> 512), and a tract past
    /// `MAX_SHIFT_TRACT`, which must decline.
    ///
    /// Without the doubling the first pad would silently cap the roll and the
    /// two spellings would key apart — the same defect one order of magnitude
    /// further out, which is exactly the failure a single fixed pad would hide.
    #[test]
    fn the_pad_grows_until_it_contains_the_roll_then_declines() {
        let homopolymer = |run: usize| {
            let mut provider = MockProvider::new();
            provider.add_genomic_sequence("NC_RUN.1", format!("C{}C", "A".repeat(run)));
            provider
        };

        // A 300-base `A` run at 2..=301. Deleting any one `A` leaves the same
        // sequence, so every spelling must key identically — which needs a
        // window wider than one `FIRST_PAD`.
        let provider = homopolymer(300);
        let key_of = |descriptor: &str| {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            canonical_spdi(&variant, &provider).expect("must canonicalize")
        };
        let first = key_of("NC_RUN.1:g.2del");
        assert_eq!(first, key_of("NC_RUN.1:g.150del"));
        assert_eq!(first, key_of("NC_RUN.1:g.301del"));
        assert_eq!(
            (first.position, first.deletion.as_str()),
            // 0-based, so the 3'-most `A` of a run occupying 1-based 2..=301.
            (300, "A"),
            "the deletion must roll to the 3' end of the run"
        );

        // Past the cap it declines, and says why rather than returning a key
        // derived from however much reference it happened to read.
        let huge = homopolymer(MAX_SHIFT_TRACT as usize + 1_000);
        let variant = parse_hgvs("NC_RUN.1:g.2del").expect("fixture must parse");
        let error = canonical_spdi(&variant, &huge)
            .expect_err("a tract past the cap has no window-independent key");
        let message = error.to_string();
        assert!(
            message.contains("repeat tract"),
            "the decline must name the tract, not read as a generic failure; got: {message}"
        );
    }

    /// Two insertions at one interbase have no order between them, so the key
    /// declines rather than inventing one.
    ///
    /// Before this guard the two spellings produced **different keys** —
    /// `5::CA` and `8::CA` — which is worse than declining: a caller comparing
    /// keys would read one variant as two. HGVS spells "insert both, in this
    /// order" as a single ordered payload (`ins[A;C]`, `general.md:79`), so a
    /// caller who means an order has a way to say it.
    ///
    /// The decline is asserted on both spellings, because a guard that fired on
    /// only one would leave exactly the asymmetry it exists to remove.
    #[test]
    fn coincident_insertions_are_declined_rather_than_ordered() {
        for descriptor in [
            "NC_KEY.1:g.[5_6insA;5_6insC]",
            "NC_KEY.1:g.[5_6insC;5_6insA]",
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            assert!(
                canonical_spdi(&variant, &provider()).is_err(),
                "`{descriptor}` has no order between its two insertions, so no single \
                 key describes it — declining is the answer, not picking one"
            );
        }
        // The single-member spelling states its own order, so it still keys.
        //
        // Note the payload comes back **rotated**, `CA` at 8 rather than `AC` at
        // 5: the reference reads `AC` at 6-7, so inserting `AC` at 5|6 is the
        // same edit as inserting `CA` at 7|8, and the 3' roll takes it there.
        // Pinned rather than left as "some key" because a rotation is exactly
        // the kind of thing that looks like a bug to the next reader, and
        // because it is the property that makes the roll worth having.
        let ordered = key("NC_KEY.1:g.5_6insAC");
        assert_eq!((ordered.position, ordered.insertion.as_str()), (8, "CA"));
    }

    /// The same change spelled at two positions in a tract keys identically.
    ///
    /// This is what the 3' padding buys, and it is the property #1159 needs:
    /// without a window wider than the members' own span the blunt trim has no
    /// room to roll, so each spelling keys where it was written. Measured before
    /// the padding: `g.3_5del` gave `2:ATT:` and `g.4_6del` gave `3:TTA:`; `g.1del`
    /// gave `0:G:` and `g.2del` gave `1:G:`.
    ///
    /// The pairs are asserted equal *and* pinned to the 3'-most form, because
    /// equality alone is satisfied by two spellings agreeing on a wrong answer.
    #[test]
    fn one_change_spelled_two_ways_in_a_tract_keys_once() {
        // `GG` at 1-2: deleting either `G` leaves the same sequence.
        let first = key("NC_KEY.1:g.1del");
        assert_eq!(first, key("NC_KEY.1:g.2del"));
        assert_eq!((first.position, first.deletion.as_str()), (1, "G"));

        // `GG` at 9-10, away from the contig edge.
        let inner = key("NC_KEY.1:g.9del");
        assert_eq!(inner, key("NC_KEY.1:g.10del"));
        assert_eq!((inner.position, inner.deletion.as_str()), (9, "G"));

        // `ATT` at 3-5 and `TTA` at 4-6 remove the same three bases, because
        // `ref[2] == ref[5] == 'A'`.
        let rolled = key("NC_KEY.1:g.3_5del");
        assert_eq!(rolled, key("NC_KEY.1:g.4_6del"));
        assert_eq!((rolled.position, rolled.deletion.as_str()), (3, "TTA"));

        // An insertion and the `dup` that denotes it are one variant.
        let inserted = key("NC_KEY.1:g.13_14insT");
        assert_eq!(inserted, key("NC_KEY.1:g.14dup"));
        assert_eq!((inserted.position, inserted.insertion.as_str()), (14, "T"));
    }

    /// A net deletion and a net insertion round-trip through the same path, and
    /// the key trims to the block that actually changed.
    #[test]
    fn the_key_is_the_minimal_changed_block() {
        // `g.3_5del` removes three bases. The block is reported one position 3'
        // of where it was authored, as `TTA` rather than `ATT`, because
        // `ref[2] == ref[5] == 'A'` and the change therefore rolls: deleting
        // `ATT` at 2 and deleting `TTA` at 3 leave the same sequence.
        //
        // That roll is the whole point — before the window was padded 3' this
        // returned `2:ATT:` here while `g.4_6del`, the same deletion spelled one
        // base over, returned `3:TTA:`. Two keys, one variant. Which of the two
        // survives is the 3'-maximal one, per `general.md:34`, and it is the
        // form the sibling spelling already keyed to, so converging here costs
        // no shipped key that was not already ambiguous.
        let deletion = key("NC_KEY.1:g.3_5del");
        assert_eq!(deletion.deletion, "TTA");
        assert_eq!(deletion.insertion, "");
        // A pure insertion deletes nothing.
        let insertion = key("NC_KEY.1:g.5_6insCCC");
        assert_eq!(insertion.deletion, "");
        assert_eq!(insertion.insertion, "CCC");
    }

    /// Two genuinely different edits must not collide — a key that made
    /// everything equal would pass every test above.
    #[test]
    fn different_edits_get_different_keys() {
        assert_ne!(key("NC_KEY.1:g.3A>G"), key("NC_KEY.1:g.3A>C"));
        assert_ne!(key("NC_KEY.1:g.3A>G"), key("NC_KEY.1:g.4T>G"));
        assert_ne!(key("NC_KEY.1:g.3_5del"), key("NC_KEY.1:g.3_6del"));
    }

    /// `apply_to_reference` returns the window before and after, and the two
    /// differ exactly by the edit.
    #[test]
    fn apply_to_reference_returns_both_windows() {
        let variant = parse_hgvs("NC_KEY.1:g.3_7delinsGGCTA").unwrap();
        let applied = apply_to_reference(&variant, &provider()).expect("applies");
        assert_eq!(applied.accession, "NC_KEY.1");
        assert_eq!(applied.start, 2, "0-based interbase start of g.3");
        assert_eq!(applied.reference, "ATTAC");
        assert_eq!(applied.resulting, "GGCTA");
    }

    /// The shapes with no single resulting sequence must decline, not guess.
    ///
    /// Each is a distinct reason, and a caller needs the decline rather than an
    /// answer derived from one arbitrary reading of an ambiguous input. Every
    /// fixture is `expect`ed to parse, so one that stopped parsing fails loudly
    /// instead of being skipped — a skipped fixture asserts nothing.
    #[test]
    fn shapes_without_one_resulting_sequence_decline() {
        for (descriptor, why) in [
            (
                "NC_KEY.1:g.[3_5del;4T>G]",
                "overlapping members — applying them depends on order",
            ),
            (
                "NC_KEY.1:g.[3A>G(;)7C>A]",
                "trans phase — two molecules, not one sequence",
            ),
            (
                "[NC_KEY.1:g.3A>G;NC_OTHER.1:g.7T>G]",
                "members on different accessions",
            ),
        ] {
            let variant = parse_hgvs(descriptor)
                .unwrap_or_else(|e| panic!("fixture `{descriptor}` must parse: {e}"));
            assert!(
                apply_to_reference(&variant, &provider()).is_err(),
                "`{descriptor}` must decline ({why})"
            );
            assert!(
                canonical_spdi(&variant, &provider()).is_err(),
                "`{descriptor}` must decline for the key too ({why})"
            );
        }
    }

    /// The control for the test above: a single-member variant on a served
    /// accession applies. Without this, every decline assertion would also be
    /// satisfied by an `apply_to_reference` that declined *everything*.
    #[test]
    fn a_plain_single_member_variant_applies() {
        let variant = parse_hgvs("NC_KEY.1:g.3A>G").expect("must parse");
        assert!(apply_to_reference(&variant, &provider()).is_ok());
        assert!(canonical_spdi(&variant, &provider()).is_ok());
    }

    /// An unknown accession declines rather than panicking.
    #[test]
    fn an_unreadable_reference_declines() {
        let variant = parse_hgvs("NC_ABSENT.1:g.3A>G").unwrap();
        assert!(apply_to_reference(&variant, &provider()).is_err());
    }

    /// The 40-base fixture served as a transcript, so the `r.` axis is reachable;
    /// `spell` chooses the DNA or the uracil spelling of the same molecule.
    fn transcript_provider(spell_uracil: bool) -> MockProvider {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};

        const DNA: &str = "GGATTACAGGCATTAGCCTGAGGATTACAGGCATTAGCCT";
        let sequence = if spell_uracil {
            DNA.replace('T', "U")
        } else {
            DNA.to_string()
        };
        let length = sequence.len() as u64;
        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript::new(
            "NR_KEY.1".to_string(),
            Some("SYNTH".to_string()),
            Strand::Plus,
            sequence.clone(),
            None,
            None,
            vec![Exon::with_genomic(1, 1, length, 1, length)],
            Some("chr_key".to_string()),
            Some(1),
            Some(length),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        ));
        provider.add_genomic_sequence("chr_key", sequence);
        provider
    }

    /// The alphabet fold, in both directions: what it reconciles, and the
    /// boundary it does not reach.
    ///
    /// **Case is reconciled unconditionally.** Real genomic FASTAs lowercase
    /// repeat tracts, so the same accession can be served soft-masked by one
    /// provider and uppercase by another; case carries no biological meaning, so
    /// the key must not move. Every shape whose key can carry reference-derived
    /// bases is covered — a stated deletion, an unspelled one, a `dup` and an
    /// `inv` (both of which read their inserted bases *out of* the reference),
    /// and a `delins` mixing both directions.
    ///
    /// **`U`/`T` is reconciled only part of the way, and that is the point of
    /// asserting it.** [`key_alphabet`] folds *after* the two comparisons that
    /// decide the emitted block, and both are `eq_ignore_ascii_case`, which does
    /// not relate `U` to `T`. So a uracil-spelled provider agrees with a
    /// DNA-spelled one exactly when the U-bearing bases reach the block through
    /// [`trim_common_flanks`]; where a `U` sits in a **stated** deletion
    /// [`apply_triples`] reads a ref-mismatch and declines, and where it stops an
    /// insertion's common-prefix scan the roll is cut short and the position
    /// moves. Both rows are pinned as the current answer rather than left to be
    /// rediscovered as a surprise — a uracil-spelled transcript provider is not a
    /// RefSeq spelling, so neither is a defect on any supported path, but a doc
    /// claiming blanket reconciliation would be wrong and this is what keeps it
    /// honest.
    #[test]
    fn the_alphabet_fold_reconciles_case_and_bounds_the_uracil_case() {
        for descriptor in [
            "NC_KEY.1:g.3_7delinsGGCTA",
            "NC_KEY.1:g.3_5del",
            "NC_KEY.1:g.3_5dup",
            "NC_KEY.1:g.3_7inv",
            "NC_KEY.1:g.3A>G",
            "NC_KEY.1:g.5_6insAC",
        ] {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            let mut masked = MockProvider::new();
            masked.add_genomic_sequence(
                "NC_KEY.1",
                "GGATTACAGGCATTAGCCTGAGGATTACAGGCATTAGCCT".to_ascii_lowercase(),
            );
            let upper =
                canonical_spdi(&variant, &provider()).expect("the uppercase reference keys");
            let lower = canonical_spdi(&variant, &masked).expect("the soft-masked reference keys");
            assert_eq!(
                upper, lower,
                "`{descriptor}` must key identically on a soft-masked reference"
            );
            assert!(
                upper.deletion.chars().all(|c| !c.is_ascii_lowercase())
                    && upper.insertion.chars().all(|c| !c.is_ascii_lowercase()),
                "`{descriptor}` emitted lowercase bases: {upper}"
            );
        }

        let dna = transcript_provider(false);
        let uracil = transcript_provider(true);
        let keyed = |descriptor: &str, provider: &MockProvider| {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            canonical_spdi(&variant, provider)
                .ok()
                .map(|spdi| spdi.to_string())
        };

        // Reconciled: the U-bearing bases reach the block through the trim.
        for descriptor in [
            "NR_KEY.1:r.[3a>g;7c>a]",
            "NR_KEY.1:r.3a>g",
            "NR_KEY.1:r.14dup",
        ] {
            assert_eq!(
                keyed(descriptor, &dna),
                keyed(descriptor, &uracil),
                "`{descriptor}` is a shape the fold does reconcile"
            );
        }
        assert_eq!(
            keyed("NR_KEY.1:r.[3a>g;7c>a]", &uracil).as_deref(),
            Some("NR_KEY.1:2:ATTAC:GTTAA"),
            "and it reconciles by folding, not by refusing both sides"
        );

        // Not reconciled: a stated deletion is validated before the fold.
        for descriptor in [
            "NR_KEY.1:r.3_5del",
            "NR_KEY.1:r.3_7delinsggcua",
            "NR_KEY.1:r.3_7inv",
        ] {
            assert!(
                keyed(descriptor, &dna).is_some(),
                "`{descriptor}` keys on the RefSeq spelling"
            );
            assert_eq!(
                keyed(descriptor, &uracil),
                None,
                "`{descriptor}` is declined on a uracil provider — `apply_triples` \
                 validates the stated deletion case-insensitively, not alphabet-insensitively"
            );
        }

        // Not reconciled: the roll's common-prefix scan stops at the first U/T
        // disagreement, so the position moves rather than the bases.
        assert_eq!(
            (
                keyed("NR_KEY.1:r.13_14insu", &dna).as_deref(),
                keyed("NR_KEY.1:r.13_14insu", &uracil).as_deref()
            ),
            (Some("NR_KEY.1:14::T"), Some("NR_KEY.1:13::T")),
            "the 3' roll is cut short on a uracil provider"
        );
    }

    /// A cis allele mixing an `r.` member with an `n.` one on a single accession
    /// is reachable, so [`key_alphabet`]'s `any`-versus-`all` choice is a real
    /// decision and not a formality.
    ///
    /// The comment on that arm used to say the members "must share one accession
    /// to key at all, so in practice they share one axis". They need not: the
    /// description below parses, is one cis allele on one accession, and keys.
    /// With `any` the whole key is folded and agrees with the DNA-spelled
    /// provider; with `all` it would emit `2:AUUAC:GUUAA`, one key carrying two
    /// alphabets, which is exactly what the fold exists to prevent. The `n.`-only
    /// spelling is asserted alongside as the control that the `r.` member is what
    /// selects the fold.
    #[test]
    fn a_mixed_axis_cis_allele_folds_on_any_rna_member() {
        let uracil = transcript_provider(true);
        let keyed = |descriptor: &str| {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            canonical_spdi(&variant, &uracil)
                .unwrap_or_else(|e| panic!("`{descriptor}` must key: {e}"))
                .to_string()
        };
        assert_eq!(
            keyed("[NR_KEY.1:r.3a>g;NR_KEY.1:n.7C>A]"),
            "NR_KEY.1:2:ATTAC:GTTAA",
            "one `r.` member folds the whole key"
        );
        assert_eq!(
            keyed("NR_KEY.1:n.[3A>G;7C>A]"),
            "NR_KEY.1:2:AUUAC:GUUAA",
            "the same allele with no `r.` member is not folded — the control that \
             the `r.` member is what selects `AlphabetMode::Rna`"
        );
    }
}
