//! Reference-backed verification that a rewrite preserved the sequence a cis
//! allele denotes.
//!
//! Structural reasoning about when a rewrite is safe keeps failing on this
//! class of problem. A 3'-shift is a rotation, so pulling a member back toward
//! where it started *looks* unconditionally sequence-preserving — and it is,
//! for that member alone. It stops being true the moment a sibling edits a base
//! inside the tract the member rotated through, and no property of the edit
//! kinds or span lengths distinguishes the two cases. Only the resulting
//! sequence can.
//!
//! So this module answers exactly one question, by construction rather than by
//! inference: applied to their shared reference, do these two member lists
//! produce the same bases?
//!
//! # Why not reuse `equivalence::checker`
//!
//! [`crate::equivalence`] has a sibling of [`apply_triples`] serving its
//! `SequenceMatch` rung. It is deliberately not reused here:
//!
//! * it compares two whole [`HgvsVariant`]s, while the clamp holds member
//!   *lists* that are not yet an allele;
//! * its applier assumes disjoint triples and indexes out of bounds when that
//!   fails (#1244) — and overlapping members are precisely the input the clamp
//!   starts from, so a caller here must be able to ask about them safely.
//!
//! Folding the two into one applier is worth doing once #1244 lands; keeping
//! them separate now avoids coupling this fix to that one.

use crate::hgvs::variant::HgvsVariant;
use crate::reference::ReferenceProvider;
use crate::spdi::convert::hgvs_to_spdi;
use crate::spdi::SpdiVariant;

/// Widest reference window this module will fetch and compare, in bases.
///
/// A cis allele whose members span more than this is not something the clamp
/// can usefully verify, and fetching it would cost more than the rewrite saves.
/// Mirrors the equivalence checker's own cap.
const MAX_COMPARE_WINDOW: u64 = 100_000;

/// Whether `before` and `after` denote the same sequence against `provider`.
///
/// Returns `None` — "cannot tell" — rather than guessing, whenever the question
/// is not well posed or the answer cannot be derived: a member that does not
/// project to SPDI, members spread across more than one accession, a window
/// wider than [`MAX_COMPARE_WINDOW`], a reference the provider will not serve,
/// or a member list whose own edits overlap (which has no single resulting
/// sequence at all).
///
/// Callers must treat `None` as failure to verify, never as success.
pub(crate) fn same_resulting_sequence<P: ReferenceProvider>(
    provider: &P,
    before: &[HgvsVariant],
    after: &[HgvsVariant],
) -> Option<bool> {
    let (before_accession, before_triples) = project(provider, before)?;
    let (after_accession, after_triples) = project(provider, after)?;
    if before_accession != after_accession {
        return None;
    }

    // One window covering both lists, so the two edited sequences are built
    // over identical bases and can be compared directly.
    let (start, end) = window(before_triples.iter().chain(after_triples.iter()))?;
    let reference = fetch(provider, &before_accession, start, end)?;

    let before_sequence = apply_triples(&reference, start, &before_triples)?;
    let after_sequence = apply_triples(&reference, start, &after_triples)?;
    // Case-insensitive: reference FASTAs are often soft-masked, and case
    // carries no biological meaning here.
    Some(before_sequence.eq_ignore_ascii_case(&after_sequence))
}

/// Project every member to an SPDI triple on one shared accession.
fn project<P: ReferenceProvider>(
    provider: &P,
    members: &[HgvsVariant],
) -> Option<(String, Vec<SpdiVariant>)> {
    if members.is_empty() {
        return None;
    }
    let mut accession: Option<String> = None;
    let mut triples = Vec::with_capacity(members.len());
    for member in members {
        let triple = hgvs_to_spdi(member, provider).ok()?;
        match &accession {
            None => accession = Some(triple.sequence.clone()),
            // Members on different accessions do not share a sequence to
            // compare, so the question is not well posed.
            Some(existing) if *existing != triple.sequence => return None,
            Some(_) => {}
        }
        triples.push(triple);
    }
    accession.map(|accession| (accession, triples))
}

/// The 0-based interbase interval covering every triple, or `None` if it is
/// wider than [`MAX_COMPARE_WINDOW`].
fn window<'a>(triples: impl Iterator<Item = &'a SpdiVariant>) -> Option<(u64, u64)> {
    let (mut start, mut end) = (u64::MAX, u64::MIN);
    for triple in triples {
        start = start.min(triple.position);
        end = end.max(triple.position.checked_add(triple.deletion.len() as u64)?);
    }
    if start > end || end - start > MAX_COMPARE_WINDOW {
        return None;
    }
    Some((start, end))
}

/// Reference bases for the 0-based half-open interval `[start, end)`.
fn fetch<P: ReferenceProvider>(
    provider: &P,
    accession: &str,
    start: u64,
    end: u64,
) -> Option<String> {
    let bases = provider
        .get_genomic_sequence(accession, start, end)
        .or_else(|_| provider.get_sequence(accession, start, end))
        .ok()?;
    // A short read would silently shift every downstream splice.
    (bases.len() as u64 == end - start).then_some(bases)
}

/// Splice `triples` into `reference` and return the edited bases.
///
/// Applies from the 3' end so an earlier splice never shifts a later one's
/// coordinates, which also lets each stated deletion be validated against the
/// pristine reference span.
///
/// Returns `None` when the triples overlap (no single resulting sequence
/// exists, and splicing them would index past the shrinking buffer), when a
/// triple falls outside the window, or when a stated deletion disagrees with
/// the reference — a ref-mismatched input cannot be reconstructed faithfully.
fn apply_triples(reference: &str, window_start: u64, triples: &[SpdiVariant]) -> Option<String> {
    let reference_bytes = reference.as_bytes();
    let mut edited = reference_bytes.to_vec();

    let mut ordered: Vec<&SpdiVariant> = triples.iter().collect();
    ordered.sort_by_key(|triple| std::cmp::Reverse(triple.position));
    if !are_disjoint(&ordered) {
        return None;
    }

    for triple in ordered {
        let start = triple.position.checked_sub(window_start)? as usize;
        let end = start.checked_add(triple.deletion.len())?;
        if end > reference_bytes.len() {
            return None;
        }
        if !reference_bytes[start..end].eq_ignore_ascii_case(triple.deletion.as_bytes()) {
            return None;
        }
        edited.splice(start..end, triple.insertion.bytes());
    }
    String::from_utf8(edited).ok()
}

/// Whether no two of `ordered` claim the same reference base.
///
/// `ordered` must be sorted by descending position. Compares each triple
/// against the furthest 3' reach of everything 5' of it rather than against its
/// immediate neighbour, so one long triple spanning several shorter ones that
/// do not touch each other is still caught.
///
/// Abutting triples are disjoint, and so is any number of pure insertions at
/// one interbase position: an insertion deletes nothing and claims no base.
fn are_disjoint(ordered: &[&SpdiVariant]) -> bool {
    let mut reach: Option<u64> = None;
    for triple in ordered.iter().rev() {
        let Some(end) = triple.position.checked_add(triple.deletion.len() as u64) else {
            return false;
        };
        if reach.is_some_and(|furthest| furthest > triple.position) {
            return false;
        }
        reach = Some(reach.map_or(end, |furthest| furthest.max(end)));
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;
    use crate::reference::MockProvider;

    /// `TAATATATATAATATATATT` on contig `T`, the reference from #1254.
    fn provider() -> MockProvider {
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("T", "TAATATATATAATATATATT");
        provider
    }

    fn members(descriptions: &[&str]) -> Vec<HgvsVariant> {
        descriptions
            .iter()
            .map(|d| parse_hgvs(d).expect("parse"))
            .collect()
    }

    #[test]
    fn identical_member_lists_match() {
        let list = members(&["T:g.3_4del", "T:g.9del"]);
        assert_eq!(
            same_resulting_sequence(&provider(), &list, &list),
            Some(true)
        );
    }

    #[test]
    fn a_rotation_within_a_tract_matches() {
        // Deleting `AT` anywhere inside the (AT) tract at 3..=11 denotes the
        // same sequence, which is the whole premise of the 3'-shift.
        assert_eq!(
            same_resulting_sequence(
                &provider(),
                &members(&["T:g.3_4del"]),
                &members(&["T:g.9_10del"]),
            ),
            Some(true)
        );
    }

    #[test]
    fn a_move_past_the_tract_boundary_does_not_match() {
        // The tract ends at 11 (`11=A, 12=A` breaks the repeat), so a deletion
        // placed past it is a different edit. This is the discrimination the
        // clamp relies on, and the case no positional rule catches.
        assert_eq!(
            same_resulting_sequence(
                &provider(),
                &members(&["T:g.3_4del", "T:g.9del"]),
                &members(&["T:g.12_14del"]),
            ),
            Some(false)
        );
    }

    #[test]
    fn overlapping_members_are_undecidable() {
        // Two members claiming position 13 have no single resulting sequence,
        // so the answer is "cannot tell" rather than a verdict.
        let overlapping = members(&["T:g.12_14del", "T:g.13T>A"]);
        assert_eq!(
            same_resulting_sequence(&provider(), &members(&["T:g.3_4del"]), &overlapping),
            None
        );
        assert_eq!(
            same_resulting_sequence(&provider(), &overlapping, &members(&["T:g.3_4del"])),
            None
        );
    }

    #[test]
    fn abutting_members_are_decidable() {
        // Abutting is not overlapping: an edit ending where the next begins
        // still has a well-defined resulting sequence.
        let abutting = members(&["T:g.3_4del", "T:g.5A>G"]);
        assert_eq!(
            same_resulting_sequence(&provider(), &abutting, &abutting),
            Some(true)
        );
    }

    #[test]
    fn members_on_different_accessions_are_undecidable() {
        let mut provider = provider();
        provider.add_genomic_sequence("U", "TAATATATATAATATATATT");
        assert_eq!(
            same_resulting_sequence(
                &provider,
                &members(&["T:g.3_4del"]),
                &members(&["U:g.3_4del"]),
            ),
            None
        );
    }

    #[test]
    fn an_unservable_reference_is_undecidable() {
        // No bases for this accession, so nothing can be reconstructed.
        assert_eq!(
            same_resulting_sequence(
                &MockProvider::new(),
                &members(&["T:g.3_4del"]),
                &members(&["T:g.9_10del"]),
            ),
            None
        );
    }

    #[test]
    fn an_empty_member_list_is_undecidable() {
        assert_eq!(
            same_resulting_sequence(&provider(), &[], &members(&["T:g.3_4del"])),
            None
        );
    }
}
