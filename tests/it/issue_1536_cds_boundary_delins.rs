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
/// **The trap it recorded was on the CDS-start side, until #1816 closed that too.**
/// `ExtendedBody`'s fold was one-sided — 3'UTR only, `c.-N` left to #1816 — so on
/// the CDS-start side the shared-affix trim was once the only thing that could move
/// a straddling member. #1816's 5'UTR arm removed that asymmetry, so both sides now
/// derive from the sequence;
/// [`a_cds_start_straddle_now_derives_from_the_sequence`] pins the CDS-start side,
/// and it is the test to read alongside this one.
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
/// straight back to being preserved and this test fails.
///
/// #1816 has since extended the same fold to the 5'UTR, so the CDS-start side is no
/// longer the surviving decline this note once pointed at:
/// [`a_cds_start_straddle_now_derives_from_the_sequence`] is its post-fold form —
/// the mirror of this test rather than its complement.
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

/// The 5'UTR mirror of [`the_input_that_taught_us_this_test_can_be_vacuous`]:
/// once #1816 extends the sequence-first fold to the 5'UTR, a `delins` straddling
/// `cds_start` is partitioned from the bases across the boundary instead of being
/// handed back.
///
/// This test used to assert the opposite — that such a straddle was the
/// *surviving untrimmable decline*, because `merge::ExtendedBody`'s fold was
/// 3'UTR only and across `cds_start` the #1651 shared-affix trim was the sole
/// mechanism that could move a straddling member, so "untrimmable therefore
/// declines" held. #1816's 5'UTR fold ([`unfold_extended_body`]'s `FivePrimeUtr`
/// arm, `cross-zone-c-positions-order-by-transcript-coordinate`) removes that
/// asymmetry: the derivation reaches a footprint crossing `cds_start` exactly as
/// it reaches one crossing `cds_end`.
///
/// `c.-4_4delinsCC` is the mirror of the 3'UTR `c.20_*7delinsCC`: it derives to
/// `c.[-3_-1del;2_4del]`, two deletion members on OPPOSITE sides of `cds_start`
/// with `c.1` unchanged between them. That interior gap is the proof the
/// sequence-first partition ran — the affix trim can only remove shared ends, it
/// can never split one member into two around an unchanged interior base.
///
/// The whole class is then checked for the one property that must hold over every
/// placement: normalization is meaning-preserving. Each non-`inv` straddle's
/// output is applied to the reference with the independent applier
/// (`cis_apply_oracle`, not the normalizer) and must denote exactly the bases its
/// input does — a wrong coordinate in the fold would change the denoted sequence
/// and be caught here. `changed > 0` keeps it from passing vacuously on a fold
/// that declined everything, the way its decline-form predecessor could.
#[test]
fn a_cds_start_straddle_now_derives_from_the_sequence() {
    use crate::common::cis_apply_oracle::apply_with;
    let provider = transcript_provider(WIDE_CORE, 11, 30);

    // The mirror of the 3'UTR `c.20_*7delinsCC -> c.[20del;*2_*4del;*6_*7del]`: a
    // straddle across `cds_start` derives to a boundary-crossing split, `c.1`
    // unchanged between the two members — a shape the affix trim alone cannot make.
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize(&provider, "NM_TEST.1:c.-4_4delinsCC", direction),
            "NM_TEST.1:c.[-3_-1del;2_4del]",
            "the 5'UTR fold partitions this straddling delins across `cds_start`"
        );
    }

    // The two rows the decline form pinned, at their post-fold values. The first
    // is now the sequence-canonical single straddling `delins` (a clean 3->2 net
    // deletion with no interior coincidence to split), emitted as the derived
    // answer rather than handed back as a decline; the second still trims its
    // shared `T` prefix and shuffles wholly into the CDS. `c.-1`..`c.2` is
    // transcript 10..=12 under `CDS 11..=30`, read off the bases below.
    assert_eq!(&WIDE_CORE[9..12], "TCA");
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize(&provider, "NM_TEST.1:c.-1_2delinsGG", direction),
            "NM_TEST.1:c.-1_2delinsGG",
            "`TCA` -> `GG` has no interior unchanged base, so its canonical form is \
             the single straddling `delins` — now the derived answer, not a decline"
        );
        assert_eq!(
            normalize(&provider, "NM_TEST.1:c.-1_2delinsTG", direction),
            "NM_TEST.1:c.1_2delinsG",
            "`TCA` -> `TG` shares a `T` prefix; the trim eats it, the member stops \
             straddling, and it shuffles inside the CDS"
        );
    }

    // The property that must hold over EVERY non-`inv` straddle placement across
    // `cds_start`: the output denotes exactly the bases the input does. A whole-
    // span reverse complement is excluded by name — it re-types `delins` -> `inv`
    // through the region-blind per-member typing, a separate mechanism
    // (`c.-1_1delinsGA -> c.-1_1inv` on this fixture).
    let mut checked = 0usize;
    let mut changed = 0usize;
    let mut violations = Vec::new();
    for start in 1..=10i64 {
        for end in 11..=30i64 {
            let reference = &WIDE_CORE[(start - 1) as usize..end as usize];
            let s = cds_spelling(start, 11, 30);
            let e = cds_spelling(end, 11, 30);
            for payload in one_and_two_base_payloads() {
                if payload == reverse_complement(reference) {
                    continue;
                }
                let input = format!("NM_TEST.1:c.{s}_{e}delins{payload}");
                let output = normalize(&provider, &input, ShuffleDirection::ThreePrime);
                checked += 1;
                if output != input {
                    changed += 1;
                }
                let denoted_in = apply_with(&provider, WIDE_CORE, &input);
                let denoted_out = apply_with(&provider, WIDE_CORE, &output);
                if denoted_in.is_none() || denoted_in != denoted_out {
                    violations.push(format!(
                        "  {input} -> {output} (in {denoted_in:?}, out {denoted_out:?})"
                    ));
                }
            }
        }
    }
    assert!(
        changed > 0,
        "the 5'UTR fold must move at least some straddle, or this is a decline test \
         wearing a convergence test's clothes: checked={checked} changed={changed}"
    );
    assert!(
        violations.is_empty(),
        "every straddle output must denote exactly its input's bases; {} of {checked} \
         do not:\n{}",
        violations.len(),
        violations.join("\n")
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

/// The multi-member spelling of one variant converges with its `inv` spelling
/// once the members straddle `cds_end` — **closed by #1816**.
///
/// # CLOSED BY #1816 — the straddling anchor is now built, one region per endpoint
///
/// #1650 (`merge::ExtendedBody`, folding `c.*N` onto `cds_len + N`) first made
/// this converge on the **split**, because `general.md:56` then ranked the two
/// substitution competitors above the inversion. `rulings[whole-span-reverse-
/// complement-types-as-inv]` (`DNA/inversion.md:5`, 2026-08-13, #1703)
/// **overturns #1230's competitor-type ranking**: an exact whole-span reverse
/// complement is typed `inv` uniformly, so the canonical form of this variant is
/// `c.13_*4inv`, and `the_boundary_class_converges_on_the_inv_from_a_single_block`
/// pins the two single-block spellings reaching it.
///
/// The MULTI-MEMBER spelling used to be left behind: `canonicalize_from_sequence`
/// derives the whole block correctly (`coalesce_inversion_runs` produces the
/// single `13_*4inv` piece), but `rebuild_members` → `render_on_its_own_region`
/// **refused a straddling anchor** — `Anchor` carried one `Region` for both
/// endpoints, and a member running from `c.<cds_len-1>` to `c.*2` had none. #1816
/// gives `Anchor` a region per endpoint (`region` for the start, `end_region` for
/// the end), so `render_on_its_own_region` now *builds* the straddling member and
/// `build_naedit` names each endpoint in its own zone — `c.13_*4`, `*` on the end
/// alone. The group merges, the per-member pipeline types the whole-span
/// reverse complement as `inv`, and the two spellings converge. Its CDS-interior
/// twin `the_multi_member_spelling_converges_inside_the_cds` remains the live
/// control proving the merge mechanism when nothing straddles.
#[test]
fn a_multi_member_spelling_converges_with_its_inv_spelling_across_a_boundary() {
    let provider = transcript_provider(WIDE_CORE, 11, 26);
    let allele = "NM_TEST.1:c.[13_15delinsCTA;*2T>G;*4G>T]";
    let inv = "NM_TEST.1:c.13_*4inv";
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        let from_allele = normalize(&provider, allele, direction);
        assert_eq!(
            from_allele,
            normalize(&provider, inv, direction),
            "{allele} and {inv} denote one variant and must converge"
        );
        // Pin the exact canonical form here, not only transitively through the
        // sibling `the_boundary_class_...`: a mutation making both spellings equal
        // but wrong would satisfy the convergence check above on its own.
        assert_eq!(
            from_allele, inv,
            "the multi-member spelling must reach the c.13_*4inv the ruling makes canonical"
        );
    }
}

/// The 5'UTR mirror of
/// [`a_multi_member_spelling_converges_with_its_inv_spelling_across_a_boundary`]:
/// a multi-member spelling straddling `cds_start` converges with its `inv`
/// spelling once #1816 extends the sequence-first fold to the 5'UTR.
///
/// `c.[-3_-1delinsCTG;1_3delinsATC]` and `c.-3_3inv` denote one variant — the
/// whole-span reverse complement of `GATCAG` at transcript 8..=13 under
/// `CDS 11..=26` (both apply to `...CCTGATCG...`, checked with the independent
/// applier). The SINGLE-BLOCK spellings already converge on `c.-3_3inv` via the
/// region-blind per-member `inv` typing, exactly as on the 3'UTR side. The
/// MULTI-MEMBER spelling was left split for the mirror-image reason #2186's 3'UTR
/// half addressed: [`unfold_extended_body`] produced no `FivePrimeUtr` endpoint,
/// so `render_on_its_own_region` could not build an anchor straddling `cds_start`.
/// #1816's 5'UTR fold makes the merged block's endpoints expressible across the
/// boundary; the group merges and the whole-span reverse complement types as one
/// `inv`. Its 3'UTR twin
/// [`a_multi_member_spelling_converges_with_its_inv_spelling_across_a_boundary`]
/// is the control on the far side of the CDS.
#[test]
fn a_five_prime_multi_member_spelling_converges_with_its_inv_spelling_across_a_boundary() {
    let provider = transcript_provider(WIDE_CORE, 11, 26);
    let allele = "NM_TEST.1:c.[-3_-1delinsCTG;1_3delinsATC]";
    let inv = "NM_TEST.1:c.-3_3inv";
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        let from_allele = normalize(&provider, allele, direction);
        assert_eq!(
            from_allele,
            normalize(&provider, inv, direction),
            "{allele} and {inv} denote one variant and must converge"
        );
        // Pin the exact canonical, not only the convergence: a mutation making
        // both spellings equal but wrong would satisfy the check above on its own.
        assert_eq!(
            from_allele, inv,
            "the multi-member spelling must reach the c.-3_3inv the ruling makes canonical"
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
/// The MULTI-MEMBER spelling `c.[13_15delinsCTA;*2T>G;*4G>T]` is not in this
/// set because it has its own guard: #1816 gave `Anchor` a region per endpoint,
/// so it now reaches the same `inv` across the CDS/3'UTR boundary, and
/// `a_multi_member_spelling_converges_with_its_inv_spelling_across_a_boundary`
/// pins that convergence directly.
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

/// The 5'UTR mirror of `a_pair_whose_triplet_reaches_past_the_stop_codon_stays_two_members`,
/// and — unlike that sibling — a **real** sabotage check, not a defensive one.
///
/// #1816's 5'UTR fold makes `merge::apply_coding_codon_exception` reachable at the
/// CDS-*start* seam. A substitution pair separated by one nucleotide whose triplet
/// straddles `cds_start` is now admitted (the group derives from the sequence
/// rather than being refused), so the codon-eligibility check runs on it. Measured
/// with the guard instrumented: the triplet arrives at internal `triplet_start == 0`
/// with `within_cds` true, and `same_codon(0, 2)` returns `false` on its `< 1`
/// guard — so the pair correctly stays two members.
///
/// **That `same_codon` lower-bound guard is load-bearing here, and this test proves
/// it.** With the `< 1` guard removed, `same_codon` computes `(0 - 1) / 3 == (2 - 1)
/// / 3` → `0 == 0` → `true`, and this exact pair wrongly merges to
/// `c.-2_1delinsGTA` — applying `general.md:35`'s "together affecting one amino
/// acid" exception to a codon that does not exist, since `c.-2` is 5'UTR and carries
/// no reading frame. This is the sabotage check the 3' sibling's doc-comment says
/// the line is owed "the moment some other change makes the line reachable"; #1816's
/// 5' fold is that change. (`within_cds`, the 3' upper guard, remains defensive —
/// this test does not reach or demonstrate it.)
#[test]
fn a_pair_whose_triplet_straddles_the_cds_start_stays_two_members() {
    let provider = transcript_provider(WIDE_CORE, 11, 26);
    // c.-2 = tx9, c.-1 = tx10 (the unchanged middle base), c.1 = tx11: the triplet
    // crosses `cds_start`, so `general.md:34` governs and the members stay individual.
    let split = "NM_TEST.1:c.[-2A>G;1C>A]";
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize(&provider, split, direction),
            split,
            "a substitution pair whose triplet straddles `cds_start` affects no single \
             amino acid (the 5'UTR carries no reading frame), so `general.md:34` governs \
             and the members stay individual"
        );
    }
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
