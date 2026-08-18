//! #1536 — a lone `delins` whose span crosses a CDS boundary was never
//! re-derived from the sequence it denotes, so two spellings of one variant
//! reached two different fixed points.
//!
//! The identical block placed entirely inside the CDS converged correctly. The
//! discriminator was nothing but where the start or stop codon happens to fall.
//!
//! # Mechanism — CORRECTED, and the correction is the point
//!
//! This module's first form, and #1536 itself, named `merge::join_pos`
//! (`src/normalize/merge.rs`) as the cause: it refuses any interval whose two
//! ends are in different regions, so `collect_canonical_edits` returns `None`
//! and the sequence-first pass is skipped. **That refusal is real and it is not
//! what this issue was.** Measured with `normalize_with_diagnostics`:
//!
//! ```text
//! c.-4_4delinsCCTGATCG   -> unchanged    [CrossAxisVariantNotShuffled{5utr,cds}]
//! c.17_*4delinsTAAGCTAC  -> unchanged    [CrossAxisVariantNotShuffled{cds,3utr}]
//! c.*1_*8delinsACCGTAAG  -> c.*1_*8inv   []
//! c.-10_-3delinsCGGACAGC -> c.-10_-3inv  []
//! ```
//!
//! The last two rows are the disproof. A lone `delins` sitting **entirely** in
//! either UTR is refused by `collect_canonical_edits` just as hard — its region
//! is not `body_region(Cds)` — and it comes back typed as `inv` anyway. The
//! `delins -> inv` re-typing for a lone member belongs to the *per-member*
//! pipeline (`rules`' single-span typing); `merge::anchor_for_piece` states
//! outright that it does not type an inversion, and a whole-block reverse
//! complement partitions to a single piece identical to the input, so
//! `canonicalize_from_sequence` returns `None` on `rebuilt == variants` whatever
//! `join_pos` decides.
//!
//! What actually produced the two fixed points is the **#350 cross-axis bail**
//! in `normalize_cds` (`src/normalize/mod.rs`): when the two endpoints fall in
//! different axis regions it returns `canonicalize_cds_variant(variant)` — no
//! shuffle, and no canonicalization either — plus
//! `CrossAxisVariantNotShuffled`. Its stated ground is that *the shuffle* has no
//! defined semantics across an axis boundary, which is a claim about moving a
//! member, not about re-typing one against the bases under its own span.
//!
//! # The fix
//!
//! A third carve-out from that bail, beside #402's zero-width insertion and
//! #918's whole-CDS del/dup: a cross-axis `Delins`/`Inversion` runs the ordinary
//! per-member pipeline with the shuffle bound pinned to its own input span, so
//! it is re-typed and affix-trimmed but cannot be carried off the footprint it
//! arrived on. If the trim leaves it inside a single region, `normalize_cds`
//! recurses once and the ordinary 3' rule resumes — without that step the pass
//! emits `c.*1_*4del`, which is not its own fixed point.
//!
//! # What is measured here
//!
//! Two fixtures. [`PLACEMENTS`] is the original 20-base one, which holds the
//! block, the bases and the edit constant and moves only the CDS annotation.
//! [`walk`] is the wider one this fix was developed against: a 40-mer with
//! **ten real 5'UTR bases and ten real 3'UTR bases**, so "crosses a CDS
//! boundary" and "reaches the transcript end" are no longer the same statement —
//! the confound the 20-base fixture and `cis_confluence_adjudication`'s
//! `CDS_START = 1` / `CDS_END = 63` fixture both carry.
//!
//! Walking one 8-base whole-block reverse complement across all 33 placements on
//! that transcript partitions perfectly on the boundary and on nothing else:
//! before the fix, all 14 crossing placements diverged and all 19 non-crossing
//! ones converged — **including the three that sit wholly in a UTR and the one
//! that runs to the last base of the transcript**. Re-annotating the same
//! transcript with `CDS 11..=40`, so that `cds_end` *is* the transcript end,
//! makes the tx 27..34 placement converge where `CDS 11..=30` made it diverge.
//! Same bases, same payload, only the stop codon moved.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// 20 bases, with the 8-base block `AATGCACA` at transcript positions 5..=12.
/// `TGTGCATT` is its exact whole-block reverse complement, so `delins` of that
/// payload denotes precisely `inv` of that span — nothing about the pair
/// depends on where the CDS is.
const CORE: &str = "GCTGAATGCACATCCGTAGC";

/// 40 bases, the wide fixture. Long enough to carry ten 5'UTR bases and ten
/// 3'UTR bases either side of a 20-base CDS, which is what separates the two
/// competing explanations — see the module docs.
const WIDE_CORE: &str = "GCTGTCCGATCAGGTAATGCACATCCGTAGCTTACGGTAC";

fn provider(cds_start: u64, cds_end: u64) -> MockProvider {
    transcript_provider(CORE, cds_start, cds_end)
}

fn transcript_provider(core: &str, cds_start: u64, cds_end: u64) -> MockProvider {
    let mut provider = MockProvider::new();
    let length = core.len() as u64;
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        core.to_string(),
        Some(cds_start),
        Some(cds_end),
        vec![Exon::new(1, 1, length)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

fn normalize(provider: &MockProvider, input: &str, direction: ShuffleDirection) -> String {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
    Normalizer::with_config(
        provider.clone(),
        NormalizeConfig::default().with_direction(direction),
    )
    .normalize(&variant)
    .unwrap_or_else(|e| panic!("{input} must normalize: {e}"))
    .to_string()
}

/// The three placements of one block: `(cds_start, cds_end, delins, inv)`.
///
/// The `inv` spelling is the same variant written the other way. Confluence
/// means the two reach one output; the defect was that across a boundary they
/// did not.
const PLACEMENTS: [(u64, u64, &str, &str); 3] = [
    // Entirely inside the CDS — the control, which converged even before the fix.
    (
        1,
        20,
        "NM_TEST.1:c.5_12delinsTGTGCATT",
        "NM_TEST.1:c.5_12inv",
    ),
    // Across the stop codon.
    (
        1,
        8,
        "NM_TEST.1:c.5_*4delinsTGTGCATT",
        "NM_TEST.1:c.5_*4inv",
    ),
    // Across the start codon.
    (
        8,
        20,
        "NM_TEST.1:c.-3_5delinsTGTGCATT",
        "NM_TEST.1:c.-3_5inv",
    ),
];

/// The half that always worked: inside the CDS, the two spellings converge.
///
/// It is the control the whole issue rests on — "the identical block placed
/// entirely inside the CDS converges correctly" — so if it ever stops holding,
/// the defect below has been mis-stated and the record needs rewriting rather
/// than the fix landing.
#[test]
fn a_delins_inside_the_cds_converges_with_its_inv_spelling() {
    let (cds_start, cds_end, delins, inv) = PLACEMENTS[0];
    let provider = provider(cds_start, cds_end);
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        let from_delins = normalize(&provider, delins, direction);
        assert_eq!(
            from_delins,
            normalize(&provider, inv, direction),
            "{delins} and {inv} denote one variant inside the CDS and must converge"
        );
        // Named rather than merely equal: the issue's claim is not only that
        // they agree but that the `delins` is *re-derived from the sequence*
        // and typed `inv`. Two spellings that both froze would satisfy an
        // equality check and prove nothing.
        assert_eq!(
            from_delins, inv,
            "{delins} must be re-derived and typed as an inversion"
        );
    }
}

/// #1536: across a CDS boundary they now converge too.
///
/// Measured on `320e98dc`, before the fix, with the geometry in the module docs
/// — the same block, the same payload, only the CDS annotation moved:
///
/// ```text
/// CDS 1..=20   c.5_12delinsTGTGCATT   -> c.5_12inv     (converged)
/// CDS 1..=8    c.5_*4delinsTGTGCATT   -> unchanged     (did not)
/// CDS 8..=20   c.-3_5delinsTGTGCATT   -> unchanged     (did not)
/// ```
///
/// So `c.5_*4inv` and `c.5_*4delinsTGTGCATT` were both fixed points: two
/// answers for one variant, discriminated by nothing but where the stop codon
/// falls. This was `#[ignore]`d and RED against #1536 until the #350 carve-out
/// landed; it is un-ignored in the same change that closes it.
#[test]
fn a_delins_crossing_a_cds_boundary_converges_with_its_inv_spelling() {
    let mut divergent = Vec::new();
    for (cds_start, cds_end, delins, inv) in PLACEMENTS.iter().skip(1) {
        let provider = provider(*cds_start, *cds_end);
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            let from_delins = normalize(&provider, delins, direction);
            let from_inv = normalize(&provider, inv, direction);
            if from_delins != from_inv {
                divergent.push(format!(
                    "  CDS {cds_start}..={cds_end}: {delins} -> {from_delins}, but {inv} -> {from_inv}"
                ));
            }
        }
    }
    assert!(
        divergent.is_empty(),
        "#1536: {} boundary-crossing spelling pair(s) reached two fixed points for one \
         variant:\n{}",
        divergent.len(),
        divergent.join("\n")
    );
}

/// The premise, asserted separately so a fix cannot be credited to a change of
/// geometry.
///
/// If the block ever stops being an exact whole-block reverse complement, the
/// guard above becomes vacuous — `delins` and `inv` would no longer denote one
/// variant and their divergence would be correct. Checked from the bases
/// themselves rather than trusted from the constant.
#[test]
fn the_payload_is_the_exact_reverse_complement_of_the_block() {
    let block = &CORE[4..12];
    assert_eq!(block, "AATGCACA");
    assert_eq!(reverse_complement(block), "TGTGCATT");
}

// ---------------------------------------------------------------------------
// The wide fixture: a real UTR at both ends, so the boundary and the sequence
// end are separable.
// ---------------------------------------------------------------------------

fn reverse_complement(block: &str) -> String {
    block
        .chars()
        .rev()
        .map(|b| match b {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => panic!("unexpected base {other}"),
        })
        .collect()
}

/// The `c.` spelling of a 1-based transcript position under the given CDS.
fn cds_spelling(tx: i64, cds_start: i64, cds_end: i64) -> String {
    if tx < cds_start {
        format!("{}", tx - cds_start)
    } else if tx <= cds_end {
        format!("{}", tx - cds_start + 1)
    } else {
        format!("*{}", tx - cds_end)
    }
}

/// Normalize both spellings of the 8-base whole-block reverse complement at
/// transcript `start..=start + 7`, returning `(delins_output, inv_output)`.
fn walk(cds_start: i64, cds_end: i64, start: i64) -> (String, String) {
    let end = start + 7;
    let block = &WIDE_CORE[(start - 1) as usize..end as usize];
    let payload = reverse_complement(block);
    let s = cds_spelling(start, cds_start, cds_end);
    let e = cds_spelling(end, cds_start, cds_end);
    let provider = transcript_provider(WIDE_CORE, cds_start as u64, cds_end as u64);
    (
        normalize(
            &provider,
            &format!("NM_TEST.1:c.{s}_{e}delins{payload}"),
            ShuffleDirection::ThreePrime,
        ),
        normalize(
            &provider,
            &format!("NM_TEST.1:c.{s}_{e}inv"),
            ShuffleDirection::ThreePrime,
        ),
    )
}

/// Every placement of one block on a transcript with a real UTR at both ends
/// converges, crossing or not.
///
/// 33 placements — the 8-base block at every transcript start from 1 to 33 on a
/// 40-mer with `CDS 11..=30`. Before the fix this was 14 divergent and 19
/// converged, partitioning exactly on whether the block crossed `cds_start` or
/// `cds_end`.
///
/// The three placements wholly inside the 5'UTR and the three wholly inside the
/// 3'UTR are the load-bearing rows: they converged *before* the fix, which is
/// what rules out "outside the CDS" as the explanation and forces the answer to
/// be the boundary itself.
#[test]
fn every_placement_of_one_block_converges_on_a_transcript_with_two_real_utrs() {
    let (cds_start, cds_end) = (11i64, 30i64);
    let mut divergent = Vec::new();
    for start in 1..=33i64 {
        let (from_delins, from_inv) = walk(cds_start, cds_end, start);
        if from_delins != from_inv {
            divergent.push(format!(
                "  tx {start}..={}: delins -> {from_delins}, but inv -> {from_inv}",
                start + 7
            ));
        }
    }
    assert!(
        divergent.is_empty(),
        "#1536: {} of 33 placements reached two fixed points for one variant:\n{}",
        divergent.len(),
        divergent.join("\n")
    );
}

/// The confound, resolved: the discriminator is the CDS boundary, not the end
/// of the transcript.
///
/// The 20-base fixture above cannot tell the two apart — with `CDS 1..=8` on a
/// 20-mer the block's far end is in the 3'UTR but nowhere near the sequence end,
/// which is suggestive and not conclusive. `cis_confluence_adjudication`'s
/// fixture is worse: `CDS_END = 63` on a 64-base core makes `c.64` simultaneously
/// the last base of the transcript and the first base of the 3'UTR, so the two
/// hypotheses are literally the same row.
///
/// Here they are separated in both directions on one transcript:
///
/// * tx 33..=40 ends on the transcript's **last base** and converged even before
///   the fix, because it lies wholly inside the 3'UTR — so reaching the sequence
///   end is not sufficient to break it;
/// * tx 27..=34 crosses `cds_end` under `CDS 11..=30` and diverged, and sits
///   wholly inside the CDS under `CDS 11..=40` — where `cds_end` *is* the
///   transcript end — and converged. Same bases, same payload, only the
///   annotation moved.
#[test]
fn the_discriminator_is_the_cds_boundary_and_not_the_transcript_end() {
    // Reaching the transcript's last base, wholly inside the 3'UTR.
    let (from_delins, from_inv) = walk(11, 30, 33);
    assert_eq!(
        from_delins, from_inv,
        "a block running to the transcript end but not crossing a boundary must converge"
    );
    assert_eq!(
        from_delins, "NM_TEST.1:c.*3_*10inv",
        "and must be re-derived as an inversion"
    );

    // The same block, once crossing `cds_end` and once not, distinguished only
    // by where `cds_end` is annotated.
    let (crossing_delins, crossing_inv) = walk(11, 30, 27);
    assert_eq!(crossing_delins, crossing_inv);
    assert_eq!(crossing_delins, "NM_TEST.1:c.17_*4inv");

    let (interior_delins, interior_inv) = walk(11, 40, 27);
    assert_eq!(interior_delins, interior_inv);
    assert_eq!(
        interior_delins, "NM_TEST.1:c.17_24inv",
        "with `cds_end` at the transcript end the same block is CDS-interior"
    );
}

/// The re-typed member is shuffled once it stops straddling, and the output is
/// its own fixed point.
///
/// The carve-out pins the shuffle to the input's own span, because a straddling
/// member has no defined most-3' position. The shared-affix trim can dissolve
/// the straddle, and at that point the pin has no ground left: `normalize_cds`
/// recurses once and the ordinary 3' rule resumes.
///
/// Without the recursion this shape emits `c.*1_*4del`, which sits wholly in
/// the 3'UTR and therefore shuffles again on the next pass — an output that is
/// not its own fixed point, which `spec_conformance_axis` counts as a rank-1
/// conformance regression rather than a representation choice. It caught exactly
/// two rows (`s00-c3m-cds-end-del-del-p2-sep1`, `s00-c3p-cds-end-del-del-p2-sep1`).
///
/// # The input matters, and the first one chosen here did not reach the recursion
///
/// This test asserted `once == twice` on `c.20_*7delinsCC`, and **passed with the
/// recursion block disabled** — so it was pinning nothing, for two review rounds.
/// Measured at the time: that input's `ref` (`GCTTACGG`) and payload (`CC`) shared
/// neither a prefix nor a suffix, so the trim removed nothing, the member came back
/// on the footprint it arrived on, and the recursion's `new_variant != *variant`
/// gate never opened. The mechanism above was right; the demonstration was of a
/// different shape.
///
/// **That input is no longer available as the demonstration it was, and the reason
/// is this change.** It still shares no affix — but sharing no affix no longer
/// implies standing still, because `merge::ExtendedBody` folds `c.*N` onto
/// `cds_len + N`, so a footprint that *crosses* `cds_end` now reaches
/// `canonicalize_from_sequence` and is partitioned from the bases instead of handed
/// back: `c.20_*7delinsCC` -> `c.[20del;*2_*4del;*6_*7del]`.
/// [`the_input_that_taught_us_this_test_can_be_vacuous`] keeps the finding and pins
/// the new value rather than loosening to fit it.
///
/// **The trap it recorded has not gone away; it has moved to the side of the CDS
/// this change deliberately does not reach.** `ExtendedBody`'s fold is one-sided —
/// 3'UTR only, `c.-N` left to #1816 — so on the CDS-start side the shared-affix trim
/// is still the only thing that can move a straddling member, and "untrimmable
/// therefore declines" still holds there.
/// [`a_cds_start_straddle_is_the_surviving_untrimmable_decline`] is where that is
/// now anchored, and it is the test to read before choosing an input for this one.
///
/// `c.20_*4delinsG` does reach it. `ref` is `GCTTA` and the payload `G`, so the
/// trim eats the one CDS-side base and leaves `c.*1_*4del` — wholly in the 3'UTR,
/// and one base short of 3'-most, because transcript 35 repeats transcript 31:
///
/// ```text
/// with the recursion:     c.20_*4delinsG -> c.*2_*5del    (a fixed point)
/// without the recursion:  c.20_*4delinsG -> c.*1_*4del    (NOT a fixed point)
/// ```
///
/// Both halves are asserted, so neither a lost shuffle nor a lost fixed point can
/// pass, and the un-recursed form is asserted *not* to be a fixed point so the
/// pair cannot go vacuous the way it did before. Verified by forcing the
/// recursion's gate closed and re-running: this test fails, and 3'
/// `non_idempotent_outputs` in `spec_conformance_axis` reads 6 rather than 4.
#[test]
fn a_retyped_member_that_stops_straddling_is_shuffled_and_is_a_fixed_point() {
    // `CDS 11..=30` on the 40-mer, so `c.*1` is transcript 31.
    let provider = transcript_provider(WIDE_CORE, 11, 30);
    // The form the trim produces, before anything shuffles it.
    let untrimmed = "NM_TEST.1:c.*1_*4del";

    let once = normalize(
        &provider,
        "NM_TEST.1:c.20_*4delinsG",
        ShuffleDirection::ThreePrime,
    );
    assert_ne!(
        once, untrimmed,
        "the trim dissolved the straddle and nothing shuffled the result — the \
         recursion in `normalize_cds` did not run"
    );
    assert_eq!(
        once, "NM_TEST.1:c.*2_*5del",
        "the re-typed member must be shuffled to its 3'-most position once it stops \
         straddling"
    );
    let twice = normalize(&provider, &once, ShuffleDirection::ThreePrime);
    assert_eq!(
        once, twice,
        "the re-typed output must be its own fixed point, or the next pass moves it again"
    );

    // The control that makes the `assert_ne!` above measure the recursion rather
    // than a coincidence: `c.*1_*4del` really is the un-shuffled form of the same
    // answer, and really is not a fixed point — so emitting it would be the
    // rank-1 regression, not a representation choice.
    let from_untrimmed = normalize(&provider, untrimmed, ShuffleDirection::ThreePrime);
    assert_eq!(
        from_untrimmed, once,
        "{untrimmed} must reach the same answer, or it is not the un-recursed form of it"
    );
    assert_ne!(
        from_untrimmed, untrimmed,
        "{untrimmed} must NOT be a fixed point — were it one, the recursion would be \
         unobservable on this input and the assertions above would be vacuous again"
    );
}

/// The recursion returns what the intermediate's own normalization returns —
/// warnings included — so it cannot leak a warning that the direct route does not.
///
/// **Written as an equivalence, not as an absence, and that is the whole point.**
/// A bare `assert!(no AXIS_CLAMP_APPLIED)` here would be satisfied by a build that
/// never emits the warning at all, which is the vacuity the test directly above
/// exists to remember. Comparing the two routes is falsifiable instead: the
/// straddling input reaches its answer *through* the recursion in `normalize_cds`,
/// the intermediate reaches it without recursing, and if the recursion kept the
/// `AxisClampApplied` it accumulates before re-entering, these two warning lists
/// would differ and this test would fail.
///
/// The warning is genuinely live on the recursing path — an assertion at the
/// re-entry line fires on `a_retyped_member_that_stops_straddling_is_shuffled_and_is_a_fixed_point`
/// carrying `AxisClampApplied { direction: "3prime", clamp_kind: "cds_end" }`.
/// Dropping it is correct rather than lossy: its contract is to explain why a
/// position did *not* shift further, and here it did (`c.*1_*4del` ->
/// `c.*2_*5del`). See the comment at that line in `normalize_cds`.
#[test]
fn the_recursion_emits_no_stale_clamp_warning() {
    let provider = transcript_provider(WIDE_CORE, 11, 30);
    let normalize_with_warnings = |input: &str| {
        let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
        let result = Normalizer::with_config(
            provider.clone(),
            NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
        )
        .normalize_with_diagnostics(&variant)
        .unwrap_or_else(|e| panic!("{input} must normalize: {e}"));
        let mut codes: Vec<&'static str> = result.warnings.iter().map(|w| w.code()).collect();
        codes.sort_unstable();
        (result.result.to_string(), codes)
    };

    // Straddles the CDS/3'UTR boundary, so the carve-out re-types it and the
    // re-entry runs.
    let (recursed, recursed_codes) = normalize_with_warnings("NM_TEST.1:c.20_*4delinsG");
    // The form the trim produces; single-axis, so it never reaches the carve-out.
    let (direct, direct_codes) = normalize_with_warnings("NM_TEST.1:c.*1_*4del");

    assert_eq!(
        recursed, direct,
        "the recursion must reach the same answer as normalizing the intermediate, \
         or the two routes are not comparable and the warning check below is moot"
    );
    assert_ne!(
        recursed, "NM_TEST.1:c.*1_*4del",
        "the answer must have shifted past the intermediate — that is what makes a \
         retained `AxisClampApplied` a false statement rather than a true one"
    );
    assert_eq!(
        recursed_codes, direct_codes,
        "the recursing route must not carry warnings the direct route does not: \
         recursed={recursed_codes:?} direct={direct_codes:?}"
    );
}

/// The input the test above used to run on — kept for the finding, re-pinned for
/// the behaviour.
///
/// **The finding is the valuable part and it is unchanged.** This input is the
/// reason that test read as coverage for two review rounds while guarding nothing:
/// `c.20_*7delinsCC` shares no affix between `ref` (`GCTTACGG`) and payload (`CC`),
/// so the #1651 carve-out re-typed nothing, the member came back exactly as
/// authored, and any `once == twice` assertion on it was satisfied by the decline
/// rather than by the recursion.
///
/// **What this change falsifies is the premise, not the finding.** The assertion
/// below used to read "an untrimmable straddling delins must be preserved", and that
/// is the precise statement `merge::ExtendedBody` makes false: folding `c.*N` onto
/// `cds_len + N` admits a footprint that *crosses* `cds_end` to
/// `canonicalize_from_sequence`, which partitions it from the bases. Untrimmable no
/// longer implies preserved on this side of the CDS — sharing no affix stopped being
/// sufficient the moment something other than the trim could move the member.
///
/// The test's own standing instruction was that if this ever started moving, "the
/// note on the test above needs rewriting rather than this pin loosening". Both
/// notes are rewritten and the pin is kept, re-pointed at the derived split, so
/// whatever moves it next is still visible here.
///
/// Verified by mutation rather than asserted: with the extension forced off
/// (`ExtendedBody::OFF` in `canonicalize_from_sequence_with_rule`) this input goes
/// straight back to being preserved and this test fails, while
/// [`a_cds_start_straddle_is_the_surviving_untrimmable_decline`] stays green — which
/// is what says the two rows are on opposite sides of the one-sided fold.
///
/// The vacuity trap this test used to hold survives on the CDS-start side, and is
/// re-anchored by [`a_cds_start_straddle_is_the_surviving_untrimmable_decline`].
#[test]
fn the_input_that_taught_us_this_test_can_be_vacuous() {
    let provider = transcript_provider(WIDE_CORE, 11, 30);
    let taught_us = "NM_TEST.1:c.20_*7delinsCC";
    // Untrimmable read off the bases, not trusted from the string: `c.20`..`c.*7`
    // is transcript 30..=37 under `CDS 11..=30`.
    let reference = &WIDE_CORE[29..37];
    assert_eq!(reference, "GCTTACGG");
    assert_ne!(reference.as_bytes()[0], b'C', "no shared prefix with `CC`");
    assert_ne!(reference.as_bytes()[7], b'C', "no shared suffix with `CC`");
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize(&provider, taught_us, direction),
            "NM_TEST.1:c.[20del;*2_*4del;*6_*7del]",
            "the 3'UTR extension partitions this straddling delins from the bases, so \
             sharing no affix no longer keeps it still and it is no longer the decline \
             that made the test above vacuous"
        );
    }
}

/// The vacuity trap, re-anchored on the side of the CDS this change does not reach.
///
/// [`the_input_that_taught_us_this_test_can_be_vacuous`] used to hold this role and
/// cannot any more: its input now derives a three-member split. The trap itself is
/// unchanged — a straddling `delins` that the sequence-first derivation cannot reach
/// and whose affix trim removes nothing comes back exactly as authored, so an
/// `once == twice` assertion on it is satisfied by the decline rather than by the
/// recursion — so it needs an input the derivation still cannot reach.
///
/// **`merge::ExtendedBody`'s fold is deliberately one-sided**, 3'UTR only: `c.-N`
/// would need the window origin itself to move, which changes behaviour for every
/// CDS-start-adjacent variant, and #1650 leaves it to #1816 with its own
/// measurement. So across `cds_start` the shared-affix trim in the #1651 carve-out
/// is still the *only* mechanism that can move a straddling member, and the property
/// the marker above recorded holds there exactly: preserved if and only if
/// untrimmable.
///
/// That biconditional is what is asserted, over every straddling placement this
/// fixture admits, rather than a count taken elsewhere and quoted here. Both classes
/// are asserted non-empty, so a generator change that stopped producing either one
/// cannot pass as a result.
///
/// **It has exactly one exception and the exception is excluded by name, not by
/// shrinking the corpus around it.** A whole-span reverse complement is re-typed
/// `delins` -> `inv` by the per-member single-span typing, which is region-blind and
/// is the mechanism the module docs above spend their first section separating from
/// #1536; on this fixture it is `c.-1_1delinsGA` over `TC`. Nothing else in the
/// enumeration disagrees.
///
/// **Both named rows are asserted, in both directions, so this cannot go vacuous the
/// way its predecessor did.** `c.-1_2delinsGG` shares neither affix with its
/// reference `TCA` and must stand still; `c.-1_2delinsTG` shares its first base, so
/// the trim eats it, the member stops straddling and it must move. A build in which
/// the carve-out had stopped running altogether would satisfy the first assertion and
/// fail the second.
///
/// When #1816 extends the fold to the 5'UTR this test is what will go red, and that
/// is its point: the 3'UTR half of the same statement moved silently enough to sit in
/// a committed pin for two review rounds.
#[test]
fn a_cds_start_straddle_is_the_surviving_untrimmable_decline() {
    let provider = transcript_provider(WIDE_CORE, 11, 30);

    // `c.-1`..`c.2` is transcript 10..=12 under `CDS 11..=30`. Read off the bases
    // rather than trusted from the constant, for the same reason
    // `the_payload_is_the_exact_reverse_complement_of_the_block` is its own test.
    assert_eq!(&WIDE_CORE[9..12], "TCA");
    let declined = "NM_TEST.1:c.-1_2delinsGG";
    let trimmed = "NM_TEST.1:c.-1_2delinsTG";
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize(&provider, declined, direction),
            declined,
            "`TCA` and `GG` share neither affix, so the trim removes nothing, the \
             member comes back on the footprint it arrived on, and the recursion never \
             runs — which is what makes an idempotency assertion on it vacuous"
        );
        assert_eq!(
            normalize(&provider, trimmed, direction),
            "NM_TEST.1:c.1_2delinsG",
            "`TCA` and `TG` share a prefix, so the trim eats it and the member stops \
             straddling — without this the assertion above would be satisfied by a \
             carve-out that had stopped running at all"
        );
    }

    // The property those two rows are instances of, over every placement across
    // `cds_start` this fixture admits: preserved iff untrimmable. The contrast with
    // `cds_end` is stated by the test above, which pins one row there that breaks it;
    // it is not re-derived here.
    let mut counterexamples = Vec::new();
    let (mut untrimmable_rows, mut trimmable_rows) = (0usize, 0usize);
    for start in 1..=10i64 {
        for end in 11..=30i64 {
            let reference = &WIDE_CORE[(start - 1) as usize..end as usize];
            let s = cds_spelling(start, 11, 30);
            let e = cds_spelling(end, 11, 30);
            for payload in one_and_two_base_payloads() {
                // The one mechanism beside the trim that moves a straddling member,
                // excluded by name rather than by narrowing the corpus until it
                // vanished. A whole-span reverse complement is re-typed `delins` ->
                // `inv` by the per-member single-span typing, which is region-blind
                // — the module docs above disentangle exactly that from #1536, with
                // `c.*1_*8delinsACCGTAAG -> c.*1_*8inv` as the disproof. Here it is
                // `c.-1_1delinsGA` over `TC`, the only placement on this fixture the
                // biconditional below does not hold for.
                if payload == reverse_complement(reference) {
                    continue;
                }
                let untrimmable = reference.as_bytes()[0] != payload.as_bytes()[0]
                    && reference.as_bytes()[reference.len() - 1]
                        != payload.as_bytes()[payload.len() - 1];
                if untrimmable {
                    untrimmable_rows += 1;
                } else {
                    trimmable_rows += 1;
                }
                let input = format!("NM_TEST.1:c.{s}_{e}delins{payload}");
                let preserved = normalize(&provider, &input, ShuffleDirection::ThreePrime) == input;
                if preserved != untrimmable {
                    counterexamples.push(format!(
                        "  {input} (ref {reference}): untrimmable={untrimmable} \
                         preserved={preserved}"
                    ));
                }
            }
        }
    }
    assert!(
        untrimmable_rows > 0 && trimmable_rows > 0,
        "both classes must be populated, or the biconditional below is `0 of 0`: \
         untrimmable={untrimmable_rows} trimmable={trimmable_rows}"
    );
    assert!(
        counterexamples.is_empty(),
        "across `cds_start` the shared-affix trim is the only mechanism that can move \
         a straddling member, so preserved and untrimmable must agree on every \
         placement. {} disagree:\n{}",
        counterexamples.len(),
        counterexamples.join("\n")
    );
}

/// Every one- and two-base `delins` payload, for the enumeration above.
fn one_and_two_base_payloads() -> Vec<String> {
    let mut payloads = Vec::with_capacity(20);
    for first in ['A', 'C', 'G', 'T'] {
        payloads.push(first.to_string());
        for second in ['A', 'C', 'G', 'T'] {
            payloads.push(format!("{first}{second}"));
        }
    }
    payloads
}

/// The multi-member spelling of one variant does **not** yet converge with its
/// `inv` spelling once the members straddle `cds_end` — deferred to #1816.
///
/// # RED BY DESIGN AGAINST #1816 SINCE #1703 — the canonical form flipped, and
/// the multi-member spelling cannot yet reach it across the boundary
///
/// #1650 (`merge::ExtendedBody`, folding `c.*N` onto `cds_len + N`) made this
/// converge, on the **split** — `general.md:56` ranked the two substitution
/// competitors above the inversion, so `c.13_*4inv` reduced to
/// `c.[13_15delinsCTA;*2T>G;*4G>T]` and the multi-member spelling was already
/// there.
///
/// `rulings[whole-span-reverse-complement-types-as-inv]` (`DNA/inversion.md:5`,
/// 2026-08-13, #1703) **overturns #1230's competitor-type ranking**: an exact
/// whole-span reverse complement is typed `inv` uniformly, so the canonical form
/// of this variant is now `c.13_*4inv`. Two spellings reach it —
/// `c.13_*4inv` stays `inv`, and the single `c.13_*4delinsCTACGGAT` is typed
/// `inv` by the per-member pipeline — and
/// `the_boundary_class_converges_on_the_inv_from_a_single_block` pins that.
///
/// The MULTI-MEMBER spelling does not converge, and the cause is a pre-existing
/// limitation the inversion ruling explicitly leaves untouched ("where a
/// whole-span reverse complement's boundaries fall … untouched"). To merge the
/// three members into one `inv`, `canonicalize_from_sequence` derives the whole
/// block correctly — `coalesce_inversion_runs` produces the single
/// `13_*4inv` piece — but `rebuild_members` → `render_on_its_own_region` then
/// **refuses a straddling anchor**: `Anchor` carries one `Region` for both
/// endpoints, and a member running from `c.<cds_len-1>` to `c.*2` has none. That
/// refusal, and its fix (widen `Anchor` to a region per endpoint), are the
/// deferred **#1816**. So the group falls back to the per-member pipeline, which
/// normalizes each member in isolation and cannot merge them — leaving the split.
///
/// It is `#[ignore]`d, red by design against #1816, exactly like this module's
/// other capability-gap guards (the #1536 and #1650 doc comments above). The
/// CDS-interior twin — where no anchor straddles — DOES converge on the `inv`
/// today and is pinned live by `the_multi_member_spelling_converges_inside_the_cds`,
/// so the mechanism is proven and only the straddling-anchor rendering is owed.
#[test]
#[ignore = "#1816: rendering a merged whole-span inv member that straddles the \
            CDS/3'UTR boundary needs a per-endpoint-region Anchor; until then the \
            multi-member spelling cannot converge on the c.13_*4inv the ruling makes canonical"]
fn a_multi_member_spelling_converges_with_its_inv_spelling_across_a_boundary() {
    let provider = transcript_provider(WIDE_CORE, 11, 26);
    let allele = "NM_TEST.1:c.[13_15delinsCTA;*2T>G;*4G>T]";
    let inv = "NM_TEST.1:c.13_*4inv";
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize(&provider, allele, direction),
            normalize(&provider, inv, direction),
            "{allele} and {inv} denote one variant and must converge"
        );
    }
}

/// The form the single-block spellings converge **on**, measured rather than
/// inferred.
///
/// # RE-PINNED BY #1703 — the canonical is the `inv`, not the split
///
/// This test used to assert all three spellings derive to the SPLIT
/// `c.[13_15delinsCTA;*2T>G;*4G>T]`, on the reasoning that "two of the three
/// competing members are substitutions and substitution outranks inversion
/// (`general.md:56`)". That is the #1230 competitor-type ranking, and
/// `rulings[whole-span-reverse-complement-types-as-inv]` (`DNA/inversion.md:5`,
/// 2026-08-13, #1703) **overturns it**: an exact whole-span reverse complement is
/// typed `inv` uniformly, whatever the competing partition is made of. So the
/// canonical form of this variant is `c.13_*4inv`, and this pins that the two
/// SINGLE-BLOCK spellings reach it — `c.13_*4inv` stays `inv`, and the single
/// `c.13_*4delinsCTACGGAT` is typed `inv` by the per-member pipeline.
///
/// The MULTI-MEMBER spelling `c.[13_15delinsCTA;*2T>G;*4G>T]` is deliberately not
/// in this set: it cannot yet reach the `inv` across the CDS/3'UTR boundary — a
/// straddling-`Anchor` limitation deferred to #1816 — and that gap is pinned,
/// red-by-design, by `a_multi_member_spelling_converges_with_its_inv_spelling_across_a_boundary`.
#[test]
fn the_boundary_class_converges_on_the_inv_from_a_single_block() {
    let provider = transcript_provider(WIDE_CORE, 11, 26);
    let expected = "NM_TEST.1:c.13_*4inv";
    for spelling in ["NM_TEST.1:c.13_*4inv", "NM_TEST.1:c.13_*4delinsCTACGGAT"] {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            assert_eq!(
                normalize(&provider, spelling, direction),
                expected,
                "{spelling} is an exact whole-span reverse complement and must type as one `inv`"
            );
        }
    }
}

/// A pair whose triplet reaches past `cds_end` does not merge — **recorded as a
/// measurement, NOT as a demonstration of `merge::apply_coding_codon_exception`'s
/// `within_cds` guard.**
///
/// The distinction is the whole value of this test, and it cost two attempts to
/// get right. `general.md:35`'s exception is a conjunction — "two variants
/// separated by one nucleotide, **together affecting one amino acid**" — and the
/// second conjunct is unstatable in the 3'UTR, so `within_cds` was added with the
/// exception's own reasoning behind it. But a guard is only coverage if some
/// input reaches it, and **neither input tried here does**:
///
/// | input | why it does not reach the guard |
/// |---|---|
/// | `c.[*1C>G;*3T>A]` under `CDS 11..=25` | sits wholly in the 3'UTR, so the straddle gate in `canonicalize_from_sequence` declines the group before any partition exists |
/// | `c.[16C>G;*2T>A]` under `CDS 11..=26` | crosses `cds_end` and IS admitted, yet still does not merge with the guard forced open |
///
/// **Both were checked against the sabotage — `within_cds` forced to `true` — and
/// both stayed green.** So `within_cds` is defensive: correct by the clause's own
/// reading, and not currently demonstrated to change any answer. It is left in
/// place rather than removed because `same_codon` is defined for every positive
/// integer and would silently answer for a codon that does not exist the moment
/// some other change makes the line reachable; but it must not be described
/// anywhere as guarded, and this test must not be read as guarding it.
///
/// What the assertion below *is* worth: it pins that the boundary-crossing pair
/// stays two members, which is `general.md:34`'s plain rule and the answer the
/// extension must not disturb.
///
/// Whoever makes the line reachable owes this test a real sabotage check.
#[test]
fn a_pair_whose_triplet_reaches_past_the_stop_codon_stays_two_members() {
    let provider = transcript_provider(WIDE_CORE, 11, 26);
    // `c.16` is the last CDS base (tx 26) and `c.*2` is tx 28, with tx 27
    // unchanged between them — so the group crosses `cds_end`.
    let split = "NM_TEST.1:c.[16C>G;*2T>A]";
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize(&provider, split, direction),
            split,
            "a pair whose triplet reaches past the stop codon affects no single amino \
             acid, so `general.md:34` governs and the members stay individual"
        );
    }
}

/// The in-CDS positive control: the identical bases at the identical transcript
/// positions, with only the stop codon moved so all three land inside the CDS,
/// **do** merge.
///
/// This one is not vacuous — it fails if the codon exception stops firing at all
/// — which is why it is kept while its sibling above is downgraded to a
/// measurement.
#[test]
fn the_same_pair_one_codon_earlier_still_merges() {
    let provider = transcript_provider(WIDE_CORE, 11, 30);
    // tx 26..28 is `c.16`..`c.18` here, and 16..18 is one codon.
    assert_eq!(
        normalize(
            &provider,
            "NM_TEST.1:c.[16C>G;18T>A]",
            ShuffleDirection::ThreePrime
        ),
        "NM_TEST.1:c.16_18delinsGGA",
        "inside the CDS the same pair is `general.md:35`'s exception and merges"
    );
}

/// The live control for the `#[ignore]`d
/// `a_multi_member_spelling_converges_with_its_inv_spelling_across_a_boundary`, so
/// that guard's ignore cannot quietly become vacuous: the identical design
/// converges when the members sit **inside** the CDS, where no anchor straddles
/// the CDS/3'UTR boundary. Post-#1703 both spellings converge on the `inv`
/// (`c.13_20inv`), which is what proves the merge mechanism works and isolates the
/// boundary guard's failure to the straddling-`Anchor` limitation (#1816) alone.
///
/// If this ever goes red, the merge itself has regressed and the boundary guard's
/// record needs rewriting rather than un-ignoring.
#[test]
fn the_multi_member_spelling_converges_inside_the_cds() {
    let provider = transcript_provider(WIDE_CORE, 11, 30);
    let allele = "NM_TEST.1:c.[13_15delinsCTA;18T>G;20G>T]";
    let inv = "NM_TEST.1:c.13_20inv";
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize(&provider, allele, direction),
            normalize(&provider, inv, direction),
            "{allele} and {inv} denote one variant inside the CDS and must converge"
        );
    }
}

/// A reference mismatch must survive the recursion, or strict mode stops refusing.
///
/// The recursion in `normalize_cds` re-enters on the re-typed variant and returns its
/// result. It used to `return` that result outright, discarding **every** warning the
/// outer pass had accumulated — and `warnings` there comes from `normalize_na_edit`,
/// which validates the reference and pushes `RefSeqMismatch`. Since
/// `NormalizeResult::has_ref_mismatch` is what strict mode refuses on, a straddling
/// delins whose stated deleted bases are wrong was silently accepted.
///
/// **The measurement that preceded this test is why it exists.** An assertion at the
/// re-entry line fired on exactly one of 9 956 tests, carrying only
/// `AxisClampApplied` — and that was read as licence to drop both vectors. It was a
/// claim about the corpus: no row builds a straddling retype with mismatched explicit
/// bases, so this class never reached the line and its loss was invisible. Filter on
/// what a warning *asserts*, not on what the corpus happens to produce.
///
/// Reference at the straddling span (tx 30..=34 under `CDS 11..=30`) is `GCTTA`, so
/// `delAAAAA` is a genuine mismatch rather than a coincidence of this fixture.
#[test]
fn a_reference_mismatch_survives_the_recursion_and_strict_mode_still_refuses() {
    let provider = transcript_provider(WIDE_CORE, 11, 30);
    let input = "NM_TEST.1:c.20_*4delAAAAAinsG";
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));

    let lenient = Normalizer::with_config(
        provider.clone(),
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    let result = lenient
        .normalize_with_diagnostics(&variant)
        .unwrap_or_else(|e| panic!("{input} must normalize in the default mode: {e}"));
    assert!(
        result.has_ref_mismatch(),
        "the reference mismatch must survive the re-typing recursion — without it \
         strict mode has nothing to refuse on. warnings={:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>()
    );

    let strict = Normalizer::with_config(
        provider,
        NormalizeConfig::default()
            .with_direction(ShuffleDirection::ThreePrime)
            .with_error_mode(ErrorMode::Strict),
    );
    assert!(
        strict.normalize(&variant).is_err(),
        "strict mode must refuse a straddling delins whose stated deleted bases do \
         not match the reference"
    );
}
