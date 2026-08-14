//! Resolving an `n.` interval's start onto the flat transcript coordinate the
//! [`CdotMapper`](crate::data::cdot::CdotMapper) API takes.
//!
//! `CdotMapper::tx_to_cds` and `tx_to_genome` accept a bare `u64`, so every
//! caller has to flatten a [`TxPos`] down to one number itself. Two handlers
//! did that inline, with the same expression, and the same defect: they read
//! `base` and dropped `downstream`, which is the `n.*N` marker.

use crate::hgvs::interval::TxInterval;

/// Resolve the start of an `n.` interval to a flat transcript coordinate.
///
/// # Refuses `n.*N`
///
/// [`TxPos::downstream`](crate::hgvs::location::TxPos::downstream) marks the
/// `n.*N` notation, which names a nucleotide **past the transcript's last
/// base**. It is not a transcript coordinate and has no `u64` to be flattened
/// to, so this refuses rather than answering:
///
/// - `background/numbering.md:52` gives the whole of `n.` numbering as "`n.1`,
///   `n.2`, `n.3`, ..., etc., from the first to the last nucleotide of the
///   reference sequence" — no `*` zone exists on this axis;
/// - `:54` states that "it is **not** allowed to describe variants in
///   nucleotides beyond the boundaries of a transcript reference sequence,
///   using that transcript reference sequence".
///
/// This is the service-side twin of the refusal
/// [`CoordinateMapper::tx_to_cds`](crate::convert::CoordinateMapper::tx_to_cds)
/// makes, and it exists because that refusal cannot reach here: these handlers
/// go through `CdotMapper`, whose entry points take the flattened `u64` and so
/// never see the flag. Dropping it made `n.5` and `n.*5` — two different
/// nucleotides — resolve to the same `5`, so the handlers reported a `c.` or
/// `g.` position, or a whole VCF record, for the wrong base with no diagnostic.
///
/// # The spelling is unreachable today; the contract still belongs here
///
/// **Do not read this as a live service defect.** `n.*N` is refused at parse in
/// *every* mode by #1748 (`E1003 InvalidPosition`; see
/// [`crate::hgvs::noncoding_zones`]), and both handlers obtain their variant
/// from `parse_hgvs_lenient` on the request string — so no HTTP request can
/// deliver a `TxPos` carrying this flag. Measured on this head: `parse_hgvs`
/// and `parse_hgvs_lenient` both refuse `NM_TEST.1:n.*5del`.
///
/// It is guarded anyway for the reason the PR's own `tx_to_cds` refusal is:
/// [`TxPos::downstream`](crate::hgvs::location::TxPos::downstream) remains
/// public API, a Rust caller can construct one, and a contract that holds only
/// because a *different* layer currently refuses the spelling is one grammar
/// change away from not holding. Claiming the route is exploitable over HTTP
/// would be overstating it; claiming it therefore needs no guard would be the
/// mirror error.
///
/// # Two other in-band markers this deliberately does NOT read
///
/// Both are pre-existing and neither is widened or narrowed here; they are
/// named so the omissions are visible rather than looking like oversights.
///
/// - **`TxPos::offset`.** `n.5+10` is intronic, and `tx_to_genome(5)` answers
///   for exonic base 5. That is the same class of silent substitution as the
///   `downstream` drop, on a far more common shape, so correcting it is a
///   behaviour change for the service that owes its own measurement and
///   disclosure rather than riding along here.
/// - **An absent inner position.** `inner()` returns `None` for an unknown
///   (`?`) endpoint and for a range boundary (`(1_5)`), and the historical
///   behaviour is to substitute `1`. Preserved verbatim.
pub(crate) fn resolve_tx_start(location: &TxInterval) -> Result<u64, String> {
    match location.start.inner() {
        Some(pos) if pos.downstream => Err(format!(
            "n.*{} lies beyond the transcript's last base and has no transcript \
             coordinate: the n. axis numbers only n.1..n.<length> \
             (background/numbering.md:52) and a variant outside a transcript's \
             boundaries may not be described against that transcript \
             (background/numbering.md:54)",
            pos.base
        )),
        Some(pos) => Ok(pos.base as u64),
        None => Ok(1),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::interval::UncertainBoundary;
    use crate::hgvs::location::TxPos;
    use crate::hgvs::uncertainty::Mu;

    /// The defect in its bare form: `n.5` and `n.*5` are different nucleotides,
    /// so the resolution must not answer the same thing for both.
    ///
    /// Before the fix both sides read `5` — a *collapse*. The assertion is
    /// written as a divergence rather than as a pinned error string so it
    /// cannot be satisfied by a message change.
    #[test]
    fn a_downstream_start_does_not_resolve_to_its_in_transcript_twin() {
        let plain = resolve_tx_start(&TxInterval::point(TxPos::new(5)));
        let downstream = resolve_tx_start(&TxInterval::point(TxPos::downstream(5)));

        assert_eq!(plain, Ok(5), "control: n.5 is transcript coordinate 5");
        assert!(
            downstream.is_err(),
            "n.*5 names a nucleotide past the transcript's last base \
             (background/numbering.md:52, :54) and has no transcript coordinate; \
             resolution answered {downstream:?} instead of refusing"
        );
        assert_ne!(
            plain, downstream,
            "n.5 and n.*5 are different nucleotides and must not resolve alike"
        );
    }

    /// The refusal is on the flag, not on the base, and it names the notation
    /// so a caller can tell this decline from a bounds one.
    #[test]
    fn every_downstream_start_is_refused_by_name() {
        for base in [1, 5, 31, 4000] {
            let err = resolve_tx_start(&TxInterval::point(TxPos::downstream(base)))
                .expect_err("n.*{base} must be refused");
            assert!(
                err.contains("n.*"),
                "the decline must name the notation it refused, got: {err}"
            );
        }
    }

    /// An offset on a downstream position is refused for the same reason — the
    /// anchor it is measured from is itself unnumberable.
    #[test]
    fn a_downstream_start_carrying_an_offset_is_refused() {
        assert!(
            resolve_tx_start(&TxInterval::point(TxPos::downstream_with_offset(5, 10))).is_err(),
            "n.*5+10 must be refused like n.*5"
        );
    }

    /// The two shapes the doc comment names as deliberately unchanged. Pinned
    /// so that "completing" the fix is a visible decision rather than a silent
    /// one, and so the `1` substitution is not mistaken for a fix side-effect.
    #[test]
    fn the_two_unread_markers_keep_their_historical_answers() {
        assert_eq!(
            resolve_tx_start(&TxInterval::point(TxPos::with_offset(5, 10))),
            Ok(5),
            "an intronic offset is still dropped; changing that is its own change"
        );
        assert_eq!(
            resolve_tx_start(&TxInterval {
                start: UncertainBoundary::Single(Mu::Unknown),
                end: UncertainBoundary::Single(Mu::Unknown),
            }),
            Ok(1),
            "an absent inner position still substitutes 1"
        );
    }
}
