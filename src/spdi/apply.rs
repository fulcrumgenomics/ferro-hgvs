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
//! normalize both and compare the strings — which fails, because normalization
//! is not encoding-invariant for complex indels (#1157) and
//! `EquivalenceChecker` was affected too (#1158).
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
//! **Not claimed:** that the triple is "maximally shifted" in the sense the SPDI
//! specification uses. The result is blunt-trimmed — common leading and trailing
//! bases removed — which is deterministic and window-independent (a wider window
//! only contributes flanks that trim away), but it does not roll an indel through
//! a repeat to a canonical extreme. Equality is the contract here, not
//! byte-compatibility with another implementation's canonical form.

use crate::error::FerroError;
use crate::hgvs::variant::HgvsVariant;
use crate::reference::ReferenceProvider;
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
        end = end.max(triple.position + triple.deletion.len() as u64);
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

    let reference = fetch_window(provider, &accession, start, end)
        .ok_or_else(|| decline("its reference window could not be read"))?;
    let resulting = apply_triples(&reference, start, &triples).ok_or_else(|| {
        decline(
            "its members overlap, or a stated reference base disagrees with the \
             reference",
        )
    })?;

    Ok(AppliedVariant {
        accession,
        start,
        reference,
        resulting,
    })
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
    let applied = apply_to_reference(variant, provider)?;
    let (offset, deletion, insertion) =
        trim_common_flanks(applied.reference.as_bytes(), applied.resulting.as_bytes());
    Ok(SpdiVariant {
        sequence: applied.accession,
        position: applied.start + offset as u64,
        deletion: String::from_utf8_lossy(deletion).into_owned(),
        insertion: String::from_utf8_lossy(insertion).into_owned(),
    })
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
/// see [`apply_to_reference`]'s error list, which this decides.
pub(crate) fn variant_edit_triples<P: ReferenceProvider + ?Sized>(
    variant: &HgvsVariant,
    provider: &P,
) -> Option<(String, Vec<SpdiVariant>)> {
    use crate::hgvs::variant::AllelePhase;

    let members: Vec<&HgvsVariant> = match variant {
        HgvsVariant::Allele(allele) => {
            // A single resulting sequence is only well-defined for a cis allele —
            // every edit applied to the same molecule.
            if allele.phase != AllelePhase::Cis {
                return None;
            }
            allele.variants.iter().collect()
        }
        HgvsVariant::NullAllele | HgvsVariant::UnknownAllele => return None,
        single => vec![single],
    };
    if members.is_empty() {
        return None;
    }

    let mut accession: Option<String> = None;
    let mut triples = Vec::with_capacity(members.len());
    for member in members {
        let spdi = crate::spdi::hgvs_to_spdi(member, provider).ok()?;
        match &accession {
            None => accession = Some(spdi.sequence.clone()),
            Some(acc) if *acc != spdi.sequence => return None,
            Some(_) => {}
        }
        triples.push(spdi);
    }
    accession.map(|acc| (acc, triples))
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
    #[test]
    fn member_order_does_not_change_the_key() {
        assert_eq!(key("NC_KEY.1:g.[3A>G;7C>A]"), key("NC_KEY.1:g.[7C>A;3A>G]"));
    }

    /// A net deletion and a net insertion round-trip through the same path, and
    /// the key trims to the block that actually changed.
    #[test]
    fn the_key_is_the_minimal_changed_block() {
        // `g.3_5del` removes `ATT`. Trimming leaves exactly those three bases
        // deleted and nothing inserted.
        let deletion = key("NC_KEY.1:g.3_5del");
        assert_eq!(deletion.deletion, "ATT");
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
}
