//! Sequence-first canonicalization on the transcript axes (`c.`, `n.`, `r.`).
//!
//! The pass was originally gated to `g.`/`m.` because the transcript axes layer
//! extra rules on a raw re-derivation — the CDS/UTR and transcript-end
//! insertion clamps, the exon-junction exception to the 3' rule, the coding
//! one-amino-acid exception to the separation rule, and the `r.` `U`/`T`
//! alphabet. Those are now covered by the pass's own refusals rather than by an
//! axis gate, so these tests pin that the canonicalization genuinely *fires*
//! here and is not merely declining everywhere.

use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};
use ferro_hgvs::reference::transcript::Strand;

/// Core `CAATT`-bearing sequence; the interesting positions are the two changed
/// bases separated by an unchanged one, the shape of #1231/#1233.
const CORE: &str = "GGGCAATTGGGCCCAAATTTGGG";

/// The description after the `r.` prefix. The accession itself contains an
/// upper-case `T` (`NM_TEST.1`), so a thymine check has to look only at the
/// edit, which the `r.` formatter renders in lower case.
fn edit_of(description: &str) -> &str {
    description
        .split_once("r.")
        .map(|(_, edit)| edit)
        .unwrap_or(description)
}

/// `c.` — two changes separated by one unchanged base must be described
/// individually, and the spanning delins spelling must reach the same string.
#[test]
fn cds_axis_converges_on_separated_changes() {
    let provider = || SyntheticBuilder::cds(CORE, 1, 21, Strand::Plus).build();
    // Core positions 4..8 are CAATT; change 4 and 6 (C>T, A>G), leaving 5
    // unchanged between them.
    // c.4, c.5, c.6 are one codon, so the coding one-amino-acid exception
    // (`general.md:35`) applies and the canonical form is the delins — the
    // opposite of the `g.` answer for the same shape. Both spellings must still
    // reach it.
    let split = normalize_to_string(provider(), "NM_TEST.1:c.[4C>T;6A>G]");
    let spanning = normalize_to_string(provider(), "NM_TEST.1:c.4_6delinsTAG");
    assert_eq!(
        split, spanning,
        "the two spellings of one c. variant must converge"
    );
    assert_eq!(split, "NM_TEST.1:c.4_6delinsTAG");
}

/// The canonicalization really does fire on the coding axis: across a codon
/// boundary the exception does not apply, so the separation rule splits the
/// spanning delins — and the split spelling is already canonical.
#[test]
fn cds_axis_splits_across_a_codon_boundary() {
    let provider = || SyntheticBuilder::cds(CORE, 1, 21, Strand::Plus).build();
    // c.6 is in codon 2, c.8 in codon 3, so the two changes do not share a codon.
    let spanning = normalize_to_string(provider(), "NM_TEST.1:c.6_8delinsGTA");
    let split = normalize_to_string(provider(), "NM_TEST.1:c.[6A>G;8T>A]");
    assert_eq!(
        spanning, split,
        "the two spellings of one c. variant must converge"
    );
    assert_eq!(split, "NM_TEST.1:c.[6A>G;8T>A]");
}

/// `n.` — same property on a non-coding transcript, where there is no CDS and
/// so no codon-frame exception to interact with.
#[test]
fn noncoding_axis_converges_on_separated_changes() {
    let provider = || SyntheticBuilder::noncoding(CORE, Strand::Plus).build();
    let split = normalize_to_string(provider(), "NR_TEST.1:n.[4C>T;6A>G]");
    let spanning = normalize_to_string(provider(), "NR_TEST.1:n.4_6delinsTAG");
    assert_eq!(
        split, spanning,
        "the two spellings of one n. variant must converge"
    );
    assert_eq!(split, "NR_TEST.1:n.[4C>T;6A>G]");
}

/// `r.` — the RNA axis renders alt bases as `U`, never `T`. A re-derivation that
/// read bases from the transcript would emit `T` and silently change the
/// alphabet; this pins that it does not.
#[test]
fn rna_axis_converges_and_keeps_the_u_alphabet() {
    // A CDS-bearing transcript, so `r.` is CDS-relative and the codon exception
    // applies — the non-coding `NR_TEST.1` fixture would be transcript-relative
    // with no reading frame at all.
    let provider = || SyntheticBuilder::cds(CORE, 1, 21, Strand::Plus).build();
    // `r.` is CDS-relative, so the codon exception applies here too.
    let split = normalize_to_string(provider(), "NM_TEST.1:r.[4c>u;6a>g]");
    let spanning = normalize_to_string(provider(), "NM_TEST.1:r.4_6delinsuag");
    assert_eq!(
        split, spanning,
        "the two spellings of one r. variant must converge"
    );
    assert_eq!(split, "NM_TEST.1:r.4_6delinsuag");
    assert!(
        !edit_of(&split).contains('t'),
        "the r. axis must not emit thymine, got `{split}`"
    );
}

/// **#1241 regression.** A *non-coding* `r.` description has no reading frame,
/// so the coding one-amino-acid exception must not apply to it and the two
/// spellings must converge on the separated two-member form.
///
/// The exception used to be keyed off the axis alone (`CisKind::Cds | Rna`),
/// which says nothing about whether the transcript is coding — so it stamped a
/// codon frame onto `NR_TEST.1` and held `r.4_6delinsuag` and `r.[4c>u;6a>g]`
/// apart as two stable fixed points. The gate is now the transcript's `CDS`,
/// not the axis letter.
///
/// The oracle is the `n.` axis over the *same* transcript, the *same*
/// positions, and the *same* change — a different internal code path, since
/// `n.` never consults `cds_start` at all. Whatever `r.` decides here it must
/// agree with `noncoding_axis_converges_on_separated_changes`, modulo the
/// `U`/`T` alphabet.
#[test]
fn noncoding_rna_axis_converges_without_a_reading_frame() {
    let provider = || SyntheticBuilder::rna(CORE, Strand::Plus).build();
    let spanning = normalize_to_string(provider(), "NR_TEST.1:r.4_6delinsuag");
    let split = normalize_to_string(provider(), "NR_TEST.1:r.[4c>u;6a>g]");
    assert_eq!(
        spanning, split,
        "the two spellings of one non-coding r. variant must converge (#1241)"
    );
    assert_eq!(split, "NR_TEST.1:r.[4c>u;6a>g]");
    // Independent path: the `n.` spelling of this very change over the same
    // fixture, which has no CDS lookup to get wrong.
    let noncoding = || SyntheticBuilder::noncoding(CORE, Strand::Plus).build();
    assert_eq!(
        normalize_to_string(noncoding(), "NR_TEST.1:n.4_6delinsTAG"),
        "NR_TEST.1:n.[4C>T;6A>G]",
        "the n. oracle must reach the same partition"
    );
    assert!(
        !edit_of(&split).contains('t'),
        "the r. axis must not emit thymine, got `{split}`"
    );
    // Still individually idempotent once converged.
    for form in [&spanning, &split] {
        assert_eq!(&normalize_to_string(provider(), form), form);
    }
}

/// A *coding* transcript keeps its reading frame on the `r.` axis: the same
/// shape that separates on `NR_TEST.1` above stays a delins on `NM_TEST.1`.
/// This is the discriminating pair for #1241's fix — a gate that simply dropped
/// `r.` from the exception would make both converge on the split form and this
/// test would catch it.
#[test]
fn coding_rna_axis_keeps_its_reading_frame() {
    let provider = || SyntheticBuilder::cds(CORE, 1, 21, Strand::Plus).build();
    assert_eq!(
        normalize_to_string(provider(), "NM_TEST.1:r.[4c>u;6a>g]"),
        "NM_TEST.1:r.4_6delinsuag"
    );
}

/// The pass genuinely *fires* on `c.` — a mixed-type allele reduces to the
/// substitutions it actually describes, which the per-member pipeline cannot do.
///
/// This is #1233's shape on the coding axis. It is the discriminating test for
/// this PR: every other case in this file is one the per-member pipeline already
/// answered, so they pass with the axis gate in place and prove nothing about it.
/// Here the gate is decisive — with it, the two spellings of this one variant
/// settle on *different* stable strings:
///
/// ```text
///                        gated (before)        this PR
///   c.[4_5insT;7del]  ->  c.[4_5insT;8del]      c.[5A>T;7T>A]
///   c.[5A>T;7T>A]     ->  c.[5A>T;7T>A]         c.[5A>T;7T>A]
/// ```
///
/// c.5 and c.7 sit in different codons (codon 2 is c.4-6), so the one-amino-acid
/// exception does not apply and the separation rule stands.
#[test]
fn cds_axis_reduces_a_mixed_type_allele_to_substitutions() {
    let provider = || SyntheticBuilder::cds(CORE, 1, 21, Strand::Plus).build();
    let mixed = normalize_to_string(provider(), "NM_TEST.1:c.[4_5insT;7del]");
    let subs = normalize_to_string(provider(), "NM_TEST.1:c.[5A>T;7T>A]");
    assert_eq!(
        mixed, subs,
        "the ins/del and substitution spellings of one c. variant must converge"
    );
    assert_eq!(mixed, "NM_TEST.1:c.[5A>T;7T>A]");
}

/// The same discriminator on `r.`, which additionally pins the alphabet.
///
/// The re-derivation reads a DNA reference and writes `Base::to_u8()`, so a
/// canonicalization that reached the output without the `U`/`T` folding would
/// either mis-compare every uracil or emit `T` on an `r.` description. The
/// expected string is fully lower-case with `u`, so both failures are visible.
#[test]
fn rna_axis_reduces_a_mixed_type_allele_and_keeps_the_u_alphabet() {
    let provider = || SyntheticBuilder::cds(CORE, 1, 21, Strand::Plus).build();
    let mixed = normalize_to_string(provider(), "NM_TEST.1:r.[4_5insu;7del]");
    let subs = normalize_to_string(provider(), "NM_TEST.1:r.[5a>u;7u>a]");
    assert_eq!(
        mixed, subs,
        "the ins/del and substitution spellings of one r. variant must converge"
    );
    assert_eq!(mixed, "NM_TEST.1:r.[5a>u;7u>a]");
    let edit = edit_of(&mixed);
    assert!(
        !edit.contains('T') && !edit.contains('A'),
        "an r. description must carry no DNA bases; got `{mixed}`"
    );
}

/// An allele whose members overlap must not be canonicalized on any axis.
///
/// Since #395/#486/#1004 ferro declines to *resolve* an overlap conflict: it
/// preserves the input's own resolution and reports `W5002`
/// (`OverlapConflictingEdits`), which strict mode promotes to an error. The
/// sequence-first pass is defined as "re-derive the description from the
/// resulting sequence", and an overlap conflict is precisely the case where
/// there is no well-defined resulting sequence to derive from — so it must
/// decline rather than launder the conflict into a tidy, non-overlapping
/// string.
///
/// `c.[4_10inv;5_6insAA]` is an insertion strictly interior to an inversion
/// span (the `#486` shape, on the transcript axis this PR opens). Before the
/// gate it came back as a single `c.5_9delinsCAAAATT` — two visibly colliding
/// members silently fused into one well-formed edit.
///
/// The refusal keys off the edit geometry, not off the `OverlapConflict`
/// warning: a gate reading a detector that only recognised `NaEdit::Insertion`
/// would decline on the first pass and accept on the second — costing
/// idempotency. `detect_insertion_overlaps` registers a `dup` as a junction
/// occupant as well, so the diagnostic survives whichever spelling the pipeline
/// settles on; this test pins both halves.
///
/// **Settled form updated by #1284.** The pipeline used to leave this as
/// `c.[5_9inv;5_6dup]`, and this test pinned that. It no longer does:
/// `respell_colliding_duplications` repairs a duplication whose span collides
/// with a sibling's bases into the equivalent insertion, and #1284 lifted that
/// repair's genomic-only gate onto the CDS-relative axes. So the `c.` spelling
/// now settles where the `g.` spelling of the same shape already did —
/// `g.[261_265inv;262_263insAA]`, measured — and the axis divergence #1284
/// describes is gone. Everything this test actually guards is unchanged: the
/// allele stays multi-member, strict still rejects it as `W5002`, the authored
/// member order is still preserved, and the result is still idempotent.
#[test]
fn overlap_conflicting_allele_is_not_canonicalized() {
    use ferro_hgvs::error::FerroError;
    use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

    const CONFLICT: &str = "NM_TEST.1:c.[4_10inv;5_6insAA]";
    let provider = || SyntheticBuilder::cds(CORE, 1, 21, Strand::Plus).build();

    // The conflict is real: strict mode still rejects it as W5002.
    let variant = parse_hgvs(CONFLICT).expect("parse");
    let err = Normalizer::with_config(provider(), NormalizeConfig::strict())
        .normalize(&variant)
        .expect_err("strict mode must reject an insertion interior to an inversion");
    match err {
        FerroError::InvalidCoordinates { msg } => assert!(
            msg.contains("W5002") || msg.contains("OverlapConflictingEdits"),
            "expected W5002 / OverlapConflictingEdits; got: {msg}"
        ),
        other => panic!("unexpected error variant: {other:?}"),
    }

    // Lenient mode must leave the conflict as the per-member pipeline spells
    // it: still two visibly-colliding members, not one re-derived edit.
    let normalized = normalize_to_string(provider(), CONFLICT);
    assert!(
        normalized.starts_with("NM_TEST.1:c.[") && normalized.contains(';'),
        "an overlap-conflicting allele must stay a multi-member allele, got `{normalized}`"
    );
    // The per-member pipeline shifts the inversion 3' and leaves the interior
    // insertion at its own junction; neither member is re-derived away.
    assert_eq!(
        normalized, "NM_TEST.1:c.[5_9inv;6_7insAA]",
        "the per-member pipeline must settle both members here — the inversion \
         3'-shifted, the interior insertion respelled at its own junction — \
         rather than the sequence-first pass re-deriving the allele as one edit"
    );

    // Idempotent on the nose: normalizing the output must return it byte for
    // byte. This is the property that broke while the overlap detector was
    // spelling-sensitive — pass one emitted the conflict with the members
    // deliberately left in authored order (the `#395` verbatim contract), then
    // pass two failed to see a conflict in the settled spelling, un-gated the
    // genomic-order sort, and reordered them.
    let twice = normalize_to_string(provider(), &normalized);
    assert_eq!(
        twice, normalized,
        "normalizing an overlap-conflicting allele must be idempotent"
    );

    // The `g.`-axis parity the doc comment above rests on, asserted rather than
    // recited. This is the independent code path — the repair was genomic-only
    // before #1284, so the `g.` answer was reached without any of the
    // CDS-relative conversion under test — and it is what makes the `c.` pin
    // above more than a self-consistency check. Left unasserted, a move in the
    // `g.` repair would leave this test green and the stated rationale silently
    // false.
    //
    // `SyntheticBuilder::cds` serves its genomic sequence under `chr_synth`
    // (`TX_CONTIG`) at `PAD_OFFSET`, so `c.N` is `g.N + 256` here: the same
    // shape authored as `g.[260_266inv;261_262insAA]`.
    let genomic = normalize_to_string(provider(), "chr_synth:g.[260_266inv;261_262insAA]");
    assert_eq!(
        genomic, "chr_synth:g.[261_265inv;262_263insAA]",
        "the g. spelling of the same shape must settle where the doc comment says"
    );
    // And the two axes now agree member-for-member, which is the divergence
    // #1284 set out to remove: strip the 256-base padding offset off the
    // genomic answer and it is the `c.` answer.
    assert_eq!(
        genomic
            .replace("chr_synth:g.", "")
            .replace("261_265", "5_9")
            .replace("262_263", "6_7"),
        normalized.replace("NM_TEST.1:c.", ""),
        "the c. and g. spellings of one shape must no longer diverge"
    );

    // The `dup` spelling of the same conflict converges onto the settled `ins`
    // spelling rather than being a second fixed point of its own (#1284): the
    // collision repair rewrites it, in whichever member order it was authored.
    // Convergence is the point — one variant, one canonical string — and it is
    // reached in a single pass, so idempotency is preserved either way.
    for (dup_spelling, converges_to) in [
        (
            "NM_TEST.1:c.[5_9inv;5_6dup]",
            "NM_TEST.1:c.[5_9inv;6_7insAA]",
        ),
        (
            "NM_TEST.1:c.[5_6dup;5_9inv]",
            "NM_TEST.1:c.[6_7insAA;5_9inv]",
        ),
    ] {
        assert_eq!(
            normalize_to_string(provider(), dup_spelling),
            converges_to,
            "a `dup` interior to an `inv` must settle on the `ins` spelling of the same conflict"
        );
        assert_eq!(
            normalize_to_string(provider(), converges_to),
            converges_to,
            "and that settled spelling must be a fixed point"
        );
    }

    // Both spellings are recognised as the same conflict by strict mode, and
    // the authored member order survives in each.
    for settled in [
        "NM_TEST.1:c.[5_9inv;5_6dup]",
        "NM_TEST.1:c.[5_6dup;5_9inv]",
        "NM_TEST.1:c.[5_9inv;6_7insAA]",
        "NM_TEST.1:c.[6_7insAA;5_9inv]",
    ] {
        // Strict mode must reject it too — and for the *same* reason as the
        // `ins` spelling, not merely reject it for some other one.
        let parsed = parse_hgvs(settled).expect("parse");
        let err = Normalizer::with_config(provider(), NormalizeConfig::strict())
            .normalize(&parsed)
            .expect_err("strict mode must reject a duplication interior to an inversion");
        match err {
            FerroError::InvalidCoordinates { msg } => assert!(
                msg.contains("W5002") || msg.contains("OverlapConflictingEdits"),
                "expected W5002 / OverlapConflictingEdits for `{settled}`; got: {msg}"
            ),
            other => panic!("unexpected error variant for `{settled}`: {other:?}"),
        }
    }
}

// ---------------------------------------------------------------------------
// The derived terminal-insertion clamp, on each transcript axis
// ---------------------------------------------------------------------------
//
// #1205 and #1217 pin this on `g.` and `m.`. It is not axis-specific — the
// derivation asks the provider for the sequence's own length — so it is live on
// `c.`, `n.` and `r.` too, and was unpinned there.
//
// The probe tail is #1205's contig, `C`x4 `A`x7 `G`x4: `insGGA` shifts twice
// through the `G` tract and then mismatches, which is what leaves it an
// insertion resting on the bound rather than a repeat or a `dup`. Each allele
// carries a second member because the transcript axes are only reachable here as
// a cis allele — `is_splittable_single_member` still gates *lone* members to
// `g.`/`m.`.

/// `n.` — the transcript body is the whole transcript, so its last base is the
/// sequence's last base and the clamp lands on `n.<len>`.
///
/// This is the one of the three the derivation actually owns: with the clamp
/// removed from the piece builder this comes back as `n.15_16insAGG`, naming a
/// position past the end of a 15-base transcript. `normalize_tx`'s own clamp
/// cannot repair that, for the same reason `normalize_genome`'s cannot on `g.` —
/// it keys on where a *shuffle* comes to rest, and the payload has nowhere left
/// to shift.
#[test]
fn noncoding_axis_clamps_a_derived_terminal_insertion() {
    let provider = || SyntheticBuilder::noncoding("CCCCAAAAAAAGGGG", Strand::Plus).build();
    assert_eq!(
        normalize_to_string(provider(), "NR_TEST.1:n.[3C>T;13_14insGGA]"),
        "NR_TEST.1:n.[3C>T;15delinsGAGG]",
    );
}

/// `r.` — same shape, lower-case RNA alphabet, and the clamp must not leak a
/// `T` into an axis that spells `u`.
#[test]
fn rna_axis_clamps_a_derived_terminal_insertion() {
    let provider = || SyntheticBuilder::rna("ccccaaaaaaagggg", Strand::Plus).build();
    let out = normalize_to_string(provider(), "NR_TEST.1:r.[3c>u;13_14insgga]");
    assert_eq!(out, "NR_TEST.1:r.[3c>u;15delinsgagg]");
    assert!(
        !edit_of(&out).contains('t') && !edit_of(&out).contains('T'),
        "the r. axis must not render a thymine: {out}",
    );
}

/// `c.` — and the answer to "what happens when the transcript's last base is in
/// the 3'UTR", which is where a naive `w_lo + offset` would be suspect.
///
/// The transcript here has a real 5'UTR *and* a real 3'UTR (`cds_start = 3`,
/// `cds_end = 17`, length 19), so `delta` is non-zero and the derived insertion
/// comes to rest on `c.*2` — two bases past the CDS.
///
/// **It renders correctly, and the reason is worth recording.** The piece
/// builder places every piece on the body region at `w_lo + offset`, which is
/// the `AxisFrame` identity `sequence = axis + delta`. That identity has a
/// discontinuity at the **5'** end of the CDS — `c.-1` is `cds_start - 1`, not
/// `cds_start - 2`, because the `c.` axis has no zero — which is why
/// `region_sequence_delta` exists. Going *up* from the CDS there is no
/// discontinuity: `c.<n>` for `n` past the CDS length denotes exactly the base
/// `c.*(n - cds_len)` does, so the overflow is arithmetically sound and the
/// formatter renders the canonical `*` form. Nothing here needs a per-region
/// offset.
///
/// Unlike the `n.` case above, the derived clamp is not load-bearing for this
/// one — `normalize_cds` re-clamps it, because a `c.` position outside the CDS
/// is representable and its shuffle can still be recomputed. It is pinned
/// anyway: it is the only place the 3'UTR rendering is asserted at all.
#[test]
fn cds_axis_clamps_a_derived_insertion_at_the_transcript_end() {
    let provider = || SyntheticBuilder::cds("CTCCCCAAAAAAAGGGGGG", 3, 17, Strand::Plus).build();
    assert_eq!(
        normalize_to_string(provider(), "NM_TEST.1:c.[2C>T;14_15insGGGA]"),
        "NM_TEST.1:c.[2C>T;*2delinsGAGGG]",
    );
}
