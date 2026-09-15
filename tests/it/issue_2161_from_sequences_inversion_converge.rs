//! #2161 Path 1 increment 1 — a geometry-varying convergence guard for
//! flanked/interior inversion typing on the derivation surface.
//!
//! `derive_block_members` now runs its placement-and-coalesce chain through
//! `place_direction_symmetrically` — both shuffle directions, member-minimal
//! kept — with `coalesce_inversion_runs` as the final coalesce pass and a
//! separation-zero re-merge after it, mirroring `canonicalize_from_sequence`. So
//! `from_sequences` types a flanked inversion the way `normalize` does. Before,
//! the derivation surface shifted a single direction and reached only the weak
//! member-level `retype_inversions`, which reads a reverse-complement *member*
//! but cannot see an inversion the partitioner split into `[del;dup]`:
//! `NC_TEST.1:g.[14del;21_23inv]` derived as `g.[14del;21del;24dup]` while
//! `normalize` retyped it to `g.[14del;21_23inv]`.
//!
//! The point of this file is the *corpus*, not the reproducer. The abandoned
//! skip-renormalize PR (#2161) passed 8k real ClinVar cis alleles because they
//! lack this geometry — "a corpus zero is a claim about the corpus". So this
//! generator deliberately *varies* the axis that corpus could not: it places a
//! genuine inversion span (revcomp ≠ identity) beside sibling `del`/`sub`
//! members at separations 0–8, on both the 5' and 3' side, in alleles of 2–4
//! members, over four structured contigs (a general contig, an inverted-repeat
//! contig, a tandem tract, and a GC palindrome), in both shuffle directions.
//!
//! # The properties this file pins
//!
//! `rederive(recommended_form = false)` (the derivation alone) and
//! `rederive(recommended_form = true)` (derivation followed by `normalize`)
//! should agree once the derivation produces the canonical form — that agreement
//! is the redundancy #2161 is about, checked here without the second partition
//! ever being removed. This file pins the increments that close the
//! inversion-bearing classes of that gap, not full convergence (repeat notation
//! remains):
//!
//! - [`from_sequences_types_a_flanked_inversion_like_normalize`] — the flanked
//!   `[del;dup]`-fragmentation port (increment 1).
//! - [`from_sequences_decomposes_a_merged_sub_flanked_inversion`] — the over-merge
//!   follow-up: a merged equal-length `sub`+`inv` `delins` split back to its
//!   canonical members (`split_canonical_delins`).
//! - `the_geometry_corpus_slice_{00..23}_converges_and_the_gate_holds` — the
//!   corpus is built in `SLICE_N` parallel slices, because a whole-corpus build
//!   is a single armed test nextest's `--partition hash:` cannot split across
//!   shards, and it was the slowest test in `test-oracle`. Each slice checks, on
//!   its 1/N, that the corpus is fully built (exact per-slice count), that
//!   convergence and inversion-typing hold their floors, and that the gate moved
//!   nothing. `slice_floors_preserve_the_global_guarantees` asserts the per-slice
//!   pins still sum to the former whole-corpus floors, so slicing did not weaken
//!   coverage.
//!
//! The guard is proven able to fail: on the pre-fix tree every curated row above
//! diverges (the `NC_TEST.1:g.[14del;21_23inv]` reproducer is one of them); the
//! shipped fix converges at least 60,731 of the 62,824 corpus rows (the sum of the
//! per-slice `CONVERGED_FLOOR` pins; 60,974 on the current tree).
//!
//! # The Drop and its gate (the #2161 Path 1 payoff)
//!
//! Once the derivation surface converges, the second block partition inside the
//! final `normalize()` is redundant, so `rederive(recommended_form = true)` skips
//! it unless a cheap AST-only gate (`repartition_gate`) says it could change the
//! output. This file guards that the SHIPPED gated path is byte-identical to the
//! full pre-drop path:
//!
//! - The `the_geometry_corpus_slice_*` tests each check this over their share of
//!   the corpus: the raw unconditional skip diverges on 742 shapes corpus-wide
//!   (`SKIP_FLOOR` sums to a `>= 72` non-vacuity floor), and the SHIPPED gated
//!   path closes every one (`gated_diverged` empty, strictly per row).
//! - [`the_gate_is_correct_on_dense_member_alleles`] +
//!   [`report_gate_fallback_on_dense_member_alleles`] — correctness and the gate's
//!   fallback rate over a densely-spaced multi-member cis corpus (the distribution
//!   the geometry corpus does not exercise), built from the public issue
//!   reproducers' member patterns.
//! - [`a_lone_delins_is_a_fixed_point_of_the_repartition`] +
//!   [`a_lone_inv_is_a_fixed_point_of_the_repartition`] — the measured fixed-point
//!   probes that justify the gate skipping every lone member.

use ferro_hgvs::{
    parse_hgvs, FromSequencesOptions, JsonProvider, NormalizeConfig, Normalizer, ShuffleDirection,
};
use std::io::Write;

/// A genome-capable provider over one contig named `NC_TEST.1` carrying `seq`.
/// Genomic on purpose: `from_sequences` emits `g.` and refuses a
/// transcript/protein accession.
fn provider(seq: &str) -> JsonProvider {
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { "NC_TEST.1": seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

/// `rederive` over a parsed description, or `None` if the description does not
/// parse or the surface declines it (overlapping members, an unbuildable
/// window). A decline is *counted* by the caller, never silently treated as a
/// pass.
fn rederive(
    nz: &Normalizer<JsonProvider>,
    desc: &str,
    direction: ShuffleDirection,
    recommended_form: bool,
) -> Option<String> {
    let variant = parse_hgvs(desc).ok()?;
    let options = FromSequencesOptions::default().with_direction(direction);
    Some(
        nz.rederive(&variant, &options, recommended_form)
            .ok()?
            .to_string(),
    )
}

/// The contigs. Each carries a structure that makes an interior inversion
/// non-trivial to type, and enough flank (≥ 12 bases either side of the region
/// the generator uses) that the derivation window is never starved.
fn contigs() -> Vec<(&'static str, &'static str)> {
    vec![
        // General mixed contig — the reproducer's own.
        (
            "general",
            "GCTAGCATGCATGCGTACAGTCGATCGATCAAAAAATGCAGTCAGTGGATCCGATTACGATCAGCT",
        ),
        // Inverted repeat: `TCGATCGA` near the start and its reverse complement
        // downstream, so an inversion's revcomp can coincide with real bases
        // elsewhere — the shape that makes a fragmented `[del;dup]` look
        // plausible.
        (
            "inverted_repeat",
            "AACCGGTTAATCGATCGATTGCACGTACGTGCAATCGATCGATTAACCGGTTAACCGGTTAACCGG",
        ),
        // Tandem tract: an `ATAT…` run in the middle where an interior inversion
        // sits inside a repeat and could shift.
        (
            "tandem",
            "GGCCAATTGGCCATATATATATATATGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATT",
        ),
        // GC palindrome region `GGGCGCGCCC` embedded in AT flanks.
        (
            "gc_palindrome",
            "ATATATATATGGGCGCGCCCATATATATATGGGCGCGCCCATATATATATGGGCGCGCCCATATAT",
        ),
    ]
}

/// A single generated member, as an HGVS fragment plus the 1-based reference
/// span it occupies (`[start, end]`, inclusive) so the caller can reject
/// overlaps before building a string that would not parse or would be declined.
struct Member {
    text: String,
    start: u64,
    end: u64,
}

/// A `del` of one base at 1-based `pos`.
fn del(pos: u64) -> Member {
    Member {
        text: format!("{pos}del"),
        start: pos,
        end: pos,
    }
}

/// A substitution at 1-based `pos`, choosing an alt base different from the
/// reference. `None` if the reference base is not one of A/C/G/T.
fn sub(seq: &[u8], pos: u64) -> Option<Member> {
    let refb = *seq.get((pos - 1) as usize)?;
    let alt = match refb {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        b'T' => b'A',
        _ => return None,
    };
    Some(Member {
        text: format!("{pos}{}>{}", refb as char, alt as char),
        start: pos,
        end: pos,
    })
}

/// An inversion of the 1-based span `[a, b]`, but only if its reverse complement
/// actually differs from the span (a palindromic span is unchanged and is not an
/// inversion at all).
fn inv(seq: &[u8], a: u64, b: u64) -> Option<Member> {
    let span = &seq[(a - 1) as usize..b as usize];
    let rc: Vec<u8> = span
        .iter()
        .rev()
        .map(|&x| match x {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => other,
        })
        .collect();
    if rc == span {
        return None;
    }
    Some(Member {
        text: format!("{a}_{b}inv"),
        start: a,
        end: b,
    })
}

/// Assemble members into a cis description if they are strictly ordered and
/// non-overlapping; `None` otherwise. Members must arrive sorted by start.
fn assemble(members: &[Member]) -> Option<String> {
    for pair in members.windows(2) {
        if pair[0].end >= pair[1].start {
            return None;
        }
    }
    let body: Vec<&str> = members.iter().map(|m| m.text.as_str()).collect();
    Some(if body.len() == 1 {
        format!("NC_TEST.1:g.{}", body[0])
    } else {
        format!("NC_TEST.1:g.[{}]", body.join(";"))
    })
}

/// One outcome of putting a generated variant through both surfaces.
struct Outcome {
    input: String,
    direction: ShuffleDirection,
    /// `Some((derive, recommend))` when both surfaces produced a string.
    compared: Option<(String, String)>,
    /// The #2161 Path 1 measurement: `Some((full, skip, gated, gate_fired))` for
    /// the derived form put through three normalization paths, plus the actual
    /// gate decision —
    /// - `full`: FULL normalization, second block partition included (the pre-drop
    ///   `rederive(true)`);
    /// - `skip`: second partition unconditionally SKIPPED
    ///   (`normalize_skipping_repartition`) — shows where the raw Drop diverges;
    /// - `gated`: the SHIPPED path (`normalize_gated_repartition`) — skip unless
    ///   the cheap `repartition_gate` fires, in which case run the full
    ///   re-partition;
    /// - `gate_fired`: whether `repartition_gate` ran the re-partition (the true
    ///   fallback signal, read from `repartition_gate_fires`).
    ///
    /// `full != skip` is a raw-Drop divergence (the residual the gate exists for);
    /// the gate is correct iff `full == gated` for **every** outcome. `None` when
    /// the derivation declined or its output did not re-parse — *counted*, never
    /// silently passed.
    drop_compared: Option<(String, String, String, bool)>,
}

/// Generate the corpus and run every variant through both surfaces — the whole
/// corpus, or one modulo-slice of it.
///
/// Every generated variant is a flanked inversion — a genuine `inv` span beside
/// at least one `del`/`sub` sibling — so the count of outcomes is the geometry
/// this file guards, floored below.
///
/// `slice = Some((k, n))` does the expensive `normalize` work only for the
/// outcomes whose global generation index `i` satisfies `i % n == k`, and skips
/// it for the rest — so `n` slices partition the corpus's cost evenly and their
/// row-sets tile it exactly, once each. `slice = None` builds the whole corpus
/// (the diagnostic and the pin meta-check use this). The generation index is
/// assigned in the deterministic iteration order, **independently of `slice`**,
/// so a row's slice membership is stable no matter which slice is asked for.
///
/// Splitting the build this way is what lets nextest's `--partition hash:` (by
/// test name) distribute the armed pass across shards: a single test can never be
/// split, but N named slice-tests can. Coverage is unchanged — every row is still
/// built armed in whichever slice owns it, and both properties are still checked
/// on every row.
fn run_corpus_slice(slice: Option<(usize, usize)>) -> Vec<Outcome> {
    let mut outcomes = Vec::new();
    // Global generation index, incremented once per (input, direction) that
    // reaches the normalize step, in deterministic order and independent of
    // `slice`. This is the number the modulo partition is taken over.
    let mut idx = 0usize;

    for (_name, seq) in contigs() {
        let bytes = seq.as_bytes();
        let len = seq.len() as u64;
        // Two normalizers, one per shuffle direction, so the gate comparison in
        // the direction loop below normalizes in the SAME direction the
        // derivation used — otherwise every `nz.normalize*` call runs 3' and the
        // 5' half of the gate is never exercised.
        let nz_3p = Normalizer::with_config(
            provider(seq),
            NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
        );
        let nz_5p = Normalizer::with_config(
            provider(seq),
            NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
        );

        // Interior inversion spans, kept ≥ 12 bases from either end.
        for a in 13..=(len - 12) {
            for l in 2..=6u64 {
                let b = a + l - 1;
                if b > len - 12 {
                    continue;
                }
                let Some(inv_member) = inv(bytes, a, b) else {
                    continue;
                };

                // Sibling positions on the 5' and 3' side, separations 0..=8.
                for sep in 0..=8u64 {
                    // 5' sibling immediately before `a - sep` (separation = the
                    // count of unchanged bases between the sibling and the inv).
                    let five = a.checked_sub(sep + 1).filter(|&p| p >= 8);
                    // 3' sibling after the inv.
                    let three = (b + sep < len - 8).then_some(b + sep + 1);

                    for &(use5, use3, use_sub) in &[
                        (true, false, false),
                        (false, true, false),
                        (true, true, false),
                        (true, false, true),
                        (false, true, true),
                        (true, true, true),
                    ] {
                        // A sibling at `pos`, on whichever side — a `sub` when
                        // `use_sub` (so substitution siblings are exercised on
                        // BOTH the 5' and 3' side, not only the 5'), else a `del`.
                        // `None` only when a `sub`'s reference base is not A/C/G/T,
                        // which never happens on these contigs.
                        let sibling = |pos: u64| {
                            if use_sub {
                                sub(bytes, pos)
                            } else {
                                Some(del(pos))
                            }
                        };

                        let mut members = Vec::new();
                        if use5 {
                            // A requested 5' sibling with no valid position must
                            // skip the WHOLE combination — building it without the
                            // sibling silently relabels it as a different shape and
                            // double-counts it in the landscape census.
                            let Some(p) = five else { continue };
                            match sibling(p) {
                                Some(m) => members.push(m),
                                None => continue,
                            }
                        }
                        // The inversion always participates.
                        members.push(Member {
                            text: inv_member.text.clone(),
                            start: inv_member.start,
                            end: inv_member.end,
                        });
                        if use3 {
                            let Some(p) = three else { continue };
                            match sibling(p) {
                                Some(m) => members.push(m),
                                None => continue,
                            }
                        }
                        members.sort_by_key(|m| m.start);
                        if members.len() < 2 {
                            continue;
                        }
                        let Some(input) = assemble(&members) else {
                            continue;
                        };

                        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime]
                        {
                            // Assign this outcome's stable global index, then skip
                            // the expensive work unless it belongs to the asked-for
                            // slice. The index advances regardless, so the same row
                            // always lands in the same slice.
                            let this_idx = idx;
                            idx += 1;
                            if let Some((k, n)) = slice {
                                if this_idx % n != k {
                                    continue;
                                }
                            }
                            let nz = if matches!(direction, ShuffleDirection::ThreePrime) {
                                &nz_3p
                            } else {
                                &nz_5p
                            };
                            let derive = rederive(nz, &input, direction, false);
                            let recommend = rederive(nz, &input, direction, true);
                            // The drop-relevant comparison: put the SETTLED
                            // derivation (`derive`) through the full normalization
                            // and through the re-partition-skipping one, and hold
                            // the two against each other. This is the property the
                            // drop actually needs — that the second block partition
                            // is a no-op on `normalize_core`'s output — as opposed
                            // to `compared`, which asks the stricter, not-required
                            // question of whether the derivation already equals the
                            // normalized form.
                            let drop_compared = derive.as_ref().and_then(|d| {
                                let variant = parse_hgvs(d).ok()?;
                                let full = nz.normalize(&variant).ok()?.to_string();
                                // `normalize_core`'s output — what the gate inspects.
                                let skip_variant =
                                    nz.normalize_skipping_repartition(&variant).ok()?;
                                let gate_fired = nz.repartition_gate_fires(&skip_variant);
                                // The SHIPPED path: skip unless the gate fires.
                                let gated =
                                    nz.normalize_gated_repartition(&variant).ok()?.to_string();
                                Some((full, skip_variant.to_string(), gated, gate_fired))
                            });
                            outcomes.push(Outcome {
                                input: input.clone(),
                                direction,
                                compared: derive.zip(recommend),
                                drop_compared,
                            });
                        }
                    }
                }
            }
        }
    }

    outcomes
}

/// The whole corpus. Kept for the ignored diagnostic and any full-corpus reader;
/// the process-cached [`geometry_corpus`] wraps it.
fn run_corpus() -> Vec<Outcome> {
    run_corpus_slice(None)
}

/// The geometry corpus is hermetic and deterministic, so its size is exact — the
/// per-slice `EXPECTED_COMPARED` pins tile it, and the slice guards assert against
/// this single source of truth.
///
/// (There is no longer a process-cached `geometry_corpus()` wrapper: under
/// nextest's process-per-test isolation an `OnceLock` cache never crossed tests,
/// so each consumer rebuilt the whole corpus — the EXPENSIVE armed
/// normalize/rederive work included. The slices split that expensive work 1/N
/// via [`run_corpus_slice`] instead — one full pass of the costly work total,
/// distributed and parallel. The cheap enumerating loop itself still runs in
/// every slice process; only the per-row normalization it guards is slice-gated.)
const GEOMETRY_CORPUS_SIZE: usize = 62_824;

// ===========================================================================
// Sliced geometry-corpus guards
// ===========================================================================
//
// The corpus's expensive armed build (`run_corpus_slice`) is split into
// `SLICE_N` named slice-tests so nextest's `--partition hash:` (by test name)
// distributes it across shards — a single test can never be split. Each slice
// builds its 1/N armed and checks, on its own rows, exactly what the two former
// whole-corpus tests checked on all of them:
//
//   - `gated_diverged` empty              — the SHIPPED gate's safety property,
//                                            strictly per-row, the real guard;
//   - `drop_compared == EXPECTED_COMPARED[k]` — non-vacuity for that guard: every
//                                            row produced a drop comparison;
//   - `compared == EXPECTED_COMPARED[k]`  — the corpus size/geometry, exact;
//   - `converged  >= CONVERGED_FLOOR[k]`  — convergence regression floor;
//   - `recommend_has_inv >= INV_FLOOR[k]` — non-vacuity: real inversions exist;
//   - `skip_diverged     >= SKIP_FLOOR[k]`— non-vacuity: the gate is exercised.
//
// The four non-per-row floors were whole-corpus in the former tests; here they
// are per-slice, and `slice_floors_preserve_the_global_guarantees` (a cheap
// test, no corpus build) asserts the pins still sum to at least the former
// global floors — so coverage is a superset: per-slice floors additionally
// catch a regression concentrated in one slice. Re-tune with the ignored
// `report_slice_pins` diagnostic, and heed the census-pin cautions in CONTRIBUTING.md
// ("Assert the property. Measure the count. Never let a count BE the property").

/// Number of slices the geometry corpus's build is partitioned into.
const SLICE_N: usize = 24;

/// Exact per-slice outcome count; sums to [`GEOMETRY_CORPUS_SIZE`]. Guards that
/// the generator still emits this slice's geometry unchanged.
const EXPECTED_COMPARED: [usize; SLICE_N] = [
    2618, 2618, 2618, 2618, 2618, 2618, 2618, 2618, 2618, 2618, 2618, 2618, 2618, 2618, 2618, 2618,
    2617, 2617, 2617, 2617, 2617, 2617, 2617, 2617,
];

/// Per-slice convergence floor (measured minus a ~20-row margin). Sum ≥ 60_731,
/// the former whole-corpus floor.
const CONVERGED_FLOOR: [usize; SLICE_N] = [
    2571, 2495, 2595, 2533, 2491, 2446, 2596, 2519, 2502, 2477, 2607, 2548, 2570, 2494, 2593, 2531,
    2495, 2452, 2599, 2524, 2507, 2478, 2607, 2552,
];

/// Per-slice "the normalize output carries an inv" floor (measured minus ~30).
/// Sum ≥ 39_760, the former whole-corpus floor.
const INV_FLOOR: [usize; SLICE_N] = [
    1642, 1642, 1764, 1764, 1411, 1411, 1729, 1729, 1520, 1520, 1874, 1874, 1656, 1656, 1768, 1768,
    1431, 1431, 1737, 1737, 1516, 1516, 1869, 1869,
];

/// Per-slice raw-skip divergence floor — a loose non-vacuity floor (the gate is
/// exercised), zero where a slice's geometries produce no skip divergence at
/// all. Sum ≥ 72, the former whole-corpus floor.
const SKIP_FLOOR: [usize; SLICE_N] = [
    20, 36, 0, 0, 86, 71, 7, 0, 94, 58, 0, 0, 23, 38, 0, 0, 89, 74, 2, 0, 87, 57, 0, 0,
];

/// Build slice `k` of [`SLICE_N`] and assert its share of the corpus properties.
fn assert_geometry_slice(k: usize) {
    let outcomes = run_corpus_slice(Some((k, SLICE_N)));

    let compared = outcomes.iter().filter(|o| o.compared.is_some()).count();
    let converged = outcomes
        .iter()
        .filter(|o| o.compared.as_ref().is_some_and(|(d, r)| d == r))
        .count();
    let recommend_has_inv = outcomes
        .iter()
        .filter(|o| o.compared.as_ref().is_some_and(|(_, r)| r.contains("inv")))
        .count();
    let gated_diverged: Vec<&Outcome> = outcomes
        .iter()
        .filter(|o| {
            o.drop_compared
                .as_ref()
                .is_some_and(|(full, _skip, gated, _fired)| full != gated)
        })
        .collect();
    let skip_diverged = outcomes
        .iter()
        .filter(|o| {
            o.drop_compared
                .as_ref()
                .is_some_and(|(full, skip, _gated, _fired)| full != skip)
        })
        .count();
    let drop_compared = outcomes
        .iter()
        .filter(|o| o.drop_compared.is_some())
        .count();

    // Exact: the generator still emits this slice's geometry unchanged.
    assert_eq!(
        compared, EXPECTED_COMPARED[k],
        "slice {k}/{SLICE_N}: compared {compared}, pinned {} — the corpus geometry moved",
        EXPECTED_COMPARED[k],
    );
    // Non-vacuity: the geometry really produces inversions `normalize` types.
    assert!(
        recommend_has_inv >= INV_FLOOR[k],
        "slice {k}/{SLICE_N}: only {recommend_has_inv} carry an inv, floor {}",
        INV_FLOOR[k],
    );
    // Regression floor: convergence.
    assert!(
        converged >= CONVERGED_FLOOR[k],
        "slice {k}/{SLICE_N}: convergence regressed, {converged} converge, floor {}",
        CONVERGED_FLOOR[k],
    );
    // Non-vacuity: the raw unconditional skip still diverges on the shapes the
    // gate exists to catch, so `gated_diverged.is_empty()` below is not vacuous.
    assert!(
        skip_diverged >= SKIP_FLOOR[k],
        "slice {k}/{SLICE_N}: raw skip diverged on only {skip_diverged}, floor {} — \
         this slice's gate coverage went vacuous",
        SKIP_FLOOR[k],
    );
    // Non-vacuity for the gate below: every row must have produced a drop
    // comparison at all, or `gated_diverged.is_empty()` passes on an empty
    // population. This is the guarantee the deleted whole-corpus
    // `the_gated_drop_changes_no_output` carried (it pinned this count at
    // `GEOMETRY_CORPUS_SIZE`), restored per slice — `SKIP_FLOOR` is 0 on 10 of
    // the 24 slices, so the skip floor alone cannot stand in for it.
    assert_eq!(
        drop_compared, EXPECTED_COMPARED[k],
        "slice {k}/{SLICE_N}: only {drop_compared} rows produced a drop comparison, \
         pinned {} — gated_diverged.is_empty() would be vacuous",
        EXPECTED_COMPARED[k],
    );
    // The SHIPPED gate's safety property, strictly per row: it moved nothing.
    for o in gated_diverged.iter().take(20) {
        let (full, _skip, gated, _fired) = o.drop_compared.as_ref().unwrap();
        eprintln!(
            "GATED-DIVERGE\t{}\t{:?}\tfull={}\tgated={}",
            o.input, o.direction, full, gated
        );
    }
    assert!(
        gated_diverged.is_empty(),
        "slice {k}/{SLICE_N}: the SHIPPED gated Drop moved {} outputs (see GATED-DIVERGE \
         lines) — the gate let a divergence through",
        gated_diverged.len(),
    );
}

macro_rules! geometry_slice_tests {
    ($($name:ident => $k:expr),+ $(,)?) => {$(
        #[test]
        fn $name() {
            assert_geometry_slice($k);
        }
    )+};
}

geometry_slice_tests! {
    the_geometry_corpus_slice_00_converges_and_the_gate_holds => 0,
    the_geometry_corpus_slice_01_converges_and_the_gate_holds => 1,
    the_geometry_corpus_slice_02_converges_and_the_gate_holds => 2,
    the_geometry_corpus_slice_03_converges_and_the_gate_holds => 3,
    the_geometry_corpus_slice_04_converges_and_the_gate_holds => 4,
    the_geometry_corpus_slice_05_converges_and_the_gate_holds => 5,
    the_geometry_corpus_slice_06_converges_and_the_gate_holds => 6,
    the_geometry_corpus_slice_07_converges_and_the_gate_holds => 7,
    the_geometry_corpus_slice_08_converges_and_the_gate_holds => 8,
    the_geometry_corpus_slice_09_converges_and_the_gate_holds => 9,
    the_geometry_corpus_slice_10_converges_and_the_gate_holds => 10,
    the_geometry_corpus_slice_11_converges_and_the_gate_holds => 11,
    the_geometry_corpus_slice_12_converges_and_the_gate_holds => 12,
    the_geometry_corpus_slice_13_converges_and_the_gate_holds => 13,
    the_geometry_corpus_slice_14_converges_and_the_gate_holds => 14,
    the_geometry_corpus_slice_15_converges_and_the_gate_holds => 15,
    the_geometry_corpus_slice_16_converges_and_the_gate_holds => 16,
    the_geometry_corpus_slice_17_converges_and_the_gate_holds => 17,
    the_geometry_corpus_slice_18_converges_and_the_gate_holds => 18,
    the_geometry_corpus_slice_19_converges_and_the_gate_holds => 19,
    the_geometry_corpus_slice_20_converges_and_the_gate_holds => 20,
    the_geometry_corpus_slice_21_converges_and_the_gate_holds => 21,
    the_geometry_corpus_slice_22_converges_and_the_gate_holds => 22,
    the_geometry_corpus_slice_23_converges_and_the_gate_holds => 23,
}

/// Cheap meta-guard (no corpus build): the per-slice pins still encode the former
/// whole-corpus guarantees, so slicing did not silently weaken coverage. If a
/// slice floor is lowered, this fails unless the global guarantee still holds.
#[test]
fn slice_floors_preserve_the_global_guarantees() {
    let sum = |v: &[usize]| v.iter().sum::<usize>();
    assert_eq!(
        sum(&EXPECTED_COMPARED),
        GEOMETRY_CORPUS_SIZE,
        "the per-slice compared pins must tile the whole corpus exactly",
    );
    // The three former whole-corpus floors, restated as pin sums.
    assert!(
        sum(&CONVERGED_FLOOR) >= 60_731,
        "converged floors sum to {}, below the former whole-corpus floor 60_731",
        sum(&CONVERGED_FLOOR),
    );
    assert!(
        sum(&INV_FLOOR) >= 39_760,
        "inv floors sum to {}, below the former whole-corpus floor 39_760",
        sum(&INV_FLOOR),
    );
    assert!(
        sum(&SKIP_FLOOR) >= 72,
        "skip floors sum to {}, below the former whole-corpus floor 72",
        sum(&SKIP_FLOOR),
    );
}

/// Diagnostic (ignored): print the per-slice pin arrays for [`SLICE_N`], so the
/// slice guards can be re-tuned mechanically after an intended corpus change.
/// Run with `--run-ignored all --no-capture`, read the arrays, paste them above.
/// It also cross-checks that the slices tile the whole corpus: the summed slice
/// metrics must equal the whole-corpus (`None`) build's.
#[test]
#[ignore = "permanent diagnostic (#2161): prints per-slice pin arrays and verifies the slices \
            tile the corpus; a re-bless helper, not a CI gate, so it stays ignored"]
fn report_slice_pins() {
    let count = |os: &[Outcome], f: &dyn Fn(&Outcome) -> bool| os.iter().filter(|o| f(o)).count();
    let is_compared = |o: &Outcome| o.compared.is_some();
    let is_converged = |o: &Outcome| o.compared.as_ref().is_some_and(|(d, r)| d == r);
    let has_inv = |o: &Outcome| o.compared.as_ref().is_some_and(|(_, r)| r.contains("inv"));
    let skip_div = |o: &Outcome| o.drop_compared.as_ref().is_some_and(|(f, s, ..)| f != s);
    let gated_div = |o: &Outcome| o.drop_compared.as_ref().is_some_and(|(f, _, g, _)| f != g);

    let mut compared = vec![0usize; SLICE_N];
    let mut converged = vec![0usize; SLICE_N];
    let mut inv = vec![0usize; SLICE_N];
    let mut skip_d = vec![0usize; SLICE_N];
    let mut gated_d = vec![0usize; SLICE_N];
    for k in 0..SLICE_N {
        let os = run_corpus_slice(Some((k, SLICE_N)));
        compared[k] = count(&os, &is_compared);
        converged[k] = count(&os, &is_converged);
        inv[k] = count(&os, &has_inv);
        skip_d[k] = count(&os, &skip_div);
        gated_d[k] = count(&os, &gated_div);
    }
    let sum = |v: &[usize]| v.iter().sum::<usize>();

    // Tiling proof: the summed slices must equal the whole-corpus build.
    let full = run_corpus_slice(None);
    assert_eq!(
        sum(&compared),
        count(&full, &is_compared),
        "compared does not tile"
    );
    assert_eq!(
        sum(&converged),
        count(&full, &is_converged),
        "converged does not tile"
    );
    assert_eq!(sum(&inv), count(&full, &has_inv), "inv does not tile");
    assert_eq!(
        sum(&skip_d),
        count(&full, &skip_div),
        "skip_diverged does not tile"
    );
    assert_eq!(
        sum(&gated_d),
        count(&full, &gated_div),
        "gated_diverged does not tile"
    );

    eprintln!("SLICE_N={SLICE_N}");
    eprintln!("EXPECTED_COMPARED = {compared:?}  sum={}", sum(&compared));
    eprintln!("converged(raw)    = {converged:?}  sum={}", sum(&converged));
    eprintln!("inv(raw)          = {inv:?}  sum={}", sum(&inv));
    eprintln!("skip_diverged(raw)= {skip_d:?}  sum={}", sum(&skip_d));
    eprintln!(
        "gated_diverged    = {gated_d:?}  sum={} (MUST be 0)",
        sum(&gated_d)
    );
}

/// The fix, stated as pins: a flanked inversion the derivation surface used to
/// fragment now types as an `inv`, and `rederive(false)` (derive) agrees with
/// `rederive(true)` (derive + `normalize`) on it — in BOTH shuffle directions,
/// across three contig structures, for a deletion 5' of the inversion, a
/// deletion 3' of it, and a multi-position case.
///
/// Each row asserts CONVERGENCE (`derive == recommend`) rather than one frozen
/// string, because the property increment 1 establishes is that the two
/// surfaces agree on the inversion — a later shift-reconciliation increment may
/// move both together and must not fail this. The headline row additionally
/// pins the exact canonical string, since it is the documented reproducer.
#[test]
fn from_sequences_types_a_flanked_inversion_like_normalize() {
    let general = "GCTAGCATGCATGCGTACAGTCGATCGATCAAAAAATGCAGTCAGTGGATCCGATTACGATCAGCT";
    let inv_rep = "AACCGGTTAATCGATCGATTGCACGTACGTGCAATCGATCGATTAACCGGTTAACCGGTTAACCGG";
    let tandem = "GGCCAATTGGCCATATATATATATATGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATT";

    // (contig, input, the `inv` span every convergent form must carry).
    let rows = [
        (general, "NC_TEST.1:g.[14del;21_23inv]", "21_23inv"),
        (general, "NC_TEST.1:g.[16_18inv;24del]", "16_18inv"),
        (inv_rep, "NC_TEST.1:g.[16del;24_28inv]", "24_28inv"),
        (inv_rep, "NC_TEST.1:g.[20_24inv;30del]", "20_24inv"),
        (tandem, "NC_TEST.1:g.[15del;20_24inv]", "20_24inv"),
    ];

    for (seq, input, inv_span) in rows {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            // Direction-matched normalizer, as in `run_corpus`: the final
            // normalize in `rederive(true)` must shuffle the SAME direction the
            // derivation used, or the FivePrime iteration would put a 5'-derived
            // form through a 3' normalizer and stop testing the intended path.
            let nz = Normalizer::with_config(
                provider(seq),
                NormalizeConfig::default().with_direction(direction),
            );
            let derive =
                rederive(&nz, input, direction, false).expect("derivation must not decline");
            let recommend =
                rederive(&nz, input, direction, true).expect("recommend must not decline");
            assert_eq!(
                derive, recommend,
                "{input} ({direction:?}): derive must converge with normalize",
            );
            assert!(
                derive.contains(inv_span),
                "{input} ({direction:?}): the inversion {inv_span} must be typed, got {derive}",
            );
        }
    }

    // The documented headline reproducer, exact.
    let nz = Normalizer::new(provider(general));
    assert_eq!(
        rederive(
            &nz,
            "NC_TEST.1:g.[14del;21_23inv]",
            ShuffleDirection::ThreePrime,
            false,
        )
        .unwrap(),
        "NC_TEST.1:g.[14del;21_23inv]",
    );
}

/// `derive_block_members` runs a final `shrink_pieces_to_differences` AFTER
/// `coalesce_inversion_runs` (to re-close the separation-zero splits that pass
/// opens on a flush-adjacent net-deletion `delins`). A concern was raised that
/// this final shrink could instead damage a *genuine* flanked inversion;
/// measured, it does not. This pins the property with the shape that stresses it
/// most — an inversion sandwiched between two deletions, so the re-close runs at
/// BOTH of its boundaries: the inversion still survives the derivation as an
/// `inv`, and `rederive(false)` (derive) converges with `rederive(true)`
/// (derive + `normalize`) in both shuffle directions.
#[test]
fn a_flanked_inversion_survives_the_final_shrink() {
    // pos 14=C and 24=A are the flanking deletions; 18_20 = TGG inverts to CCA
    // (all three bases change), and the inversion is separated from each deletion
    // by unchanged bases — so it is a genuine flanked inv, not a flush-adjacent
    // one the re-close is meant to collapse.
    let seq = "ACGTACGTACGTACGGTTGGGACACGTACGT";
    let input = "NC_TEST.1:g.[14del;18_20inv;24del]";
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        // Direction-matched normalizer (see `run_corpus`).
        let nz = Normalizer::with_config(
            provider(seq),
            NormalizeConfig::default().with_direction(direction),
        );
        let derive = rederive(&nz, input, direction, false).expect("derive must not decline");
        let recommend = rederive(&nz, input, direction, true).expect("recommend must not decline");
        assert_eq!(
            derive, recommend,
            "{input} ({direction:?}): derive must converge with normalize",
        );
        assert!(
            derive.contains("inv"),
            "{input} ({direction:?}): the flanked inversion must survive the final \
             shrink as an inv, got {derive}",
        );
    }
}

/// The over-merge increment, stated as pins: an equal-length `sub`+`inv` block
/// the place chain merged into ONE spanning `delins` is now split back into its
/// canonical `[sub; inv]` members by `split_canonical_delins`, the way
/// `normalize`'s `apply_canonical_split` does. Before this, `from_sequences`
/// emitted `g.11_15delinsCTCGC` where `normalize` emitted `g.[11A>C;13_15inv]` —
/// the largest ThreePrime divergence class after increment 1.
///
/// These converge in ThreePrime, which is the direction the drop of the second
/// partition depends on; the FivePrime placement of the surrounding subs is a
/// separate (direction-semantics) question, so unlike the flanked-inversion pins
/// above these assert ThreePrime convergence only.
#[test]
fn from_sequences_decomposes_a_merged_sub_flanked_inversion() {
    let general = "GCTAGCATGCATGCGTACAGTCGATCGATCAAAAAATGCAGTCAGTGGATCCGATTACGATCAGCT";

    // sub 5' of the inv, sub 3' of it, a longer inv, and a sub/inv/sub triple —
    // every one a net-substitution block `decompose_delins` splits.
    let rows = [
        "NC_TEST.1:g.[11A>C;13_15inv]",
        "NC_TEST.1:g.[13_15inv;17A>C]",
        "NC_TEST.1:g.[11A>C;13_17inv]",
        "NC_TEST.1:g.[11A>C;13_15inv;17A>C]",
    ];
    let nz = Normalizer::new(provider(general));
    for input in rows {
        let derive = rederive(&nz, input, ShuffleDirection::ThreePrime, false)
            .expect("derivation must not decline");
        let recommend = rederive(&nz, input, ShuffleDirection::ThreePrime, true)
            .expect("recommend must not decline");
        assert_eq!(
            derive, recommend,
            "{input}: derive must decompose the merged delins to match normalize",
        );
        assert!(
            derive.contains("inv") && derive.contains('>'),
            "{input}: the split must carry both the substitution and the inversion, got {derive}",
        );
    }

    // Exact, on the isolated block: the mechanism, not just convergence.
    assert_eq!(
        rederive(
            &nz,
            "NC_TEST.1:g.[11A>C;13_15inv]",
            ShuffleDirection::ThreePrime,
            false,
        )
        .unwrap(),
        "NC_TEST.1:g.[11A>C;13_15inv]",
    );
}

/// Report mode: print the divergence landscape. Not an assertion — run with
/// `--no-capture` while developing the corpus or a later increment, then read the
/// floors off it.
#[test]
#[ignore = "permanent diagnostic (#2161): prints the divergence landscape used to set \
            the floors in the geometry_corpus_slice tests; it never asserts, so it \
            stays ignored"]
fn report_inversion_convergence_landscape() {
    let outcomes = run_corpus();
    let compared = outcomes.iter().filter(|o| o.compared.is_some()).count();
    let converged = outcomes
        .iter()
        .filter(|o| o.compared.as_ref().is_some_and(|(d, r)| d == r))
        .count();
    let recommend_has_inv = outcomes
        .iter()
        .filter(|o| o.compared.as_ref().is_some_and(|(_, r)| r.contains("inv")))
        .count();
    let skipped = outcomes.len() - compared;
    let mut diverged = 0usize;
    for o in &outcomes {
        if let Some((d, r)) = &o.compared {
            if d != r {
                diverged += 1;
                eprintln!(
                    "DDET\t{}\t{:?}\tderive={}\trecommend={}",
                    o.input, o.direction, d, r
                );
            }
        }
    }
    eprintln!(
        "SUMMARY total={} compared={} skipped={} diverged={} converged={} recommend_has_inv={}",
        outcomes.len(),
        compared,
        skipped,
        diverged,
        converged,
        recommend_has_inv,
    );

    // The landscape: how many outputs the raw Drop (unconditional skip) MOVES,
    // how many the SHIPPED gated path moves (must be 0), and how often the gate
    // fired (the fallback count — the cost the gate pays for correctness).
    let drop_compared = outcomes
        .iter()
        .filter(|o| o.drop_compared.is_some())
        .count();
    let mut skip_diverged = 0usize;
    let mut gated_diverged = 0usize;
    let mut gate_fired = 0usize;
    for o in &outcomes {
        if let Some((full, skip, gated, fired)) = &o.drop_compared {
            if full != skip {
                skip_diverged += 1;
            }
            if full != gated {
                gated_diverged += 1;
                eprintln!(
                    "GATEDDIV\t{}\t{:?}\tfull={}\tgated={}",
                    o.input, o.direction, full, gated
                );
            }
            // The ACTUAL gate decision on `normalize_core`'s output.
            if *fired {
                gate_fired += 1;
            }
        }
    }
    eprintln!(
        "DROPSUMMARY drop_compared={} skip_diverged={} gated_diverged={} gate_fired={} \
         gate_fallback_rate={:.3}%",
        drop_compared,
        skip_diverged,
        gated_diverged,
        gate_fired,
        100.0 * gate_fired as f64 / drop_compared.max(1) as f64,
    );
}

// ===========================================================================
// #2161 dense-member cis corpus: gate CORRECTNESS + FALLBACK rate
// ===========================================================================
//
// The geometry corpus above is inversion-heavy by construction. This one mirrors
// a distribution class the geometry corpus does not: multi-member
// cis alleles authored from short (2-4 nt) del / ins / sub / dup members at small
// separations (0-2), 2-4 members — the shapes drawn from the public issue reproducers
// (`[del;ins]`, `[del;dup]`, `[dup;sub]`, `[sub;sub]`, `[del;del]`, `[ins;del]`,
// `[del;sub]`, homopolymer del runs, 3-member mixes). Crucially, NO member is
// authored as an inv or delins — those only ARISE from normalization, exactly as
// in the real workload. So this measures two things the geometry corpus cannot:
//
//   1. gate CORRECTNESS on the target distribution: `gated == full` for every
//      generated allele (asserted).
//   2. gate FALLBACK rate: what fraction of rederive(true) calls the gate sends to
//      the full re-partition (printed). That is the Drop's speed cost on this
//      distribution — the number that decides whether the gate is worth it.

/// An insertion of `bases` between 1-based positions `pos` and `pos+1`.
fn ins(pos: u64, bases: &str) -> Member {
    Member {
        text: format!("{pos}_{}ins{bases}", pos + 1),
        start: pos,
        end: pos,
    }
}

/// A tandem duplication of the 1-based span `[a, b]`.
fn dup(a: u64, b: u64) -> Member {
    let text = if a == b {
        format!("{a}dup")
    } else {
        format!("{a}_{b}dup")
    };
    Member {
        text,
        start: a,
        end: b,
    }
}

#[derive(Clone, Copy, Debug)]
enum Ek {
    Sub,
    Del,
    Ins,
    Dup,
}

/// Build one member of kind `kind` starting at 1-based `pos`, returning it and the
/// number of reference bases it occupies (for non-overlapping placement).
fn build_member(kind: Ek, seq: &[u8], pos: u64) -> Option<(Member, u64)> {
    match kind {
        Ek::Sub => sub(seq, pos).map(|m| (m, 1)),
        Ek::Del => Some((del(pos), 1)),
        // An insertion consumes no reference base, but advance one slot so the
        // next member never lands on the same interbase.
        Ek::Ins => Some((ins(pos, "T"), 1)),
        Ek::Dup => Some((dup(pos, pos + 1), 2)),
    }
}

/// Contigs with real bases, including the 125 nt `TEMPLATE` the issue reproducers
/// use, plus structured runs/palindromes so normalization can produce delins/inv.
fn dense_contigs() -> Vec<&'static str> {
    vec![
        // TEMPLATE (125 nt), verbatim from the #1229-#1421 reproducers.
        "ATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCGTACGTACGATCGGCATGCATGCTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGA",
        // General mixed + homopolymer + repeat + palindrome (from the geometry set).
        "GCTAGCATGCATGCGTACAGTCGATCGATCAAAAAATGCAGTCAGTGGATCCGATTACGATCAGCT",
        "GGCCAATTGGCCATATATATATATATGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATT",
        "ATATATATATGGGCGCGCCCATATATATATGGGCGCGCCCATATATATATGGGCGCGCCCATATAT",
    ]
}

/// Every authored member pattern, drawn from the public issue reproducers.
const DENSE_MEMBER_PATTERNS: &[&[Ek]] = &[
    &[Ek::Del, Ek::Ins],
    &[Ek::Del, Ek::Dup],
    &[Ek::Dup, Ek::Sub],
    &[Ek::Sub, Ek::Sub],
    &[Ek::Del, Ek::Del],
    &[Ek::Ins, Ek::Del],
    &[Ek::Del, Ek::Sub],
    &[Ek::Dup, Ek::Del],
    &[Ek::Sub, Ek::Del],
    &[Ek::Ins, Ek::Sub],
    &[Ek::Sub, Ek::Del, Ek::Dup],
    &[Ek::Del, Ek::Del, Ek::Del],
    &[Ek::Sub, Ek::Ins, Ek::Dup],
    &[Ek::Del, Ek::Sub, Ek::Ins],
    &[Ek::Sub, Ek::Sub, Ek::Sub],
    &[Ek::Del, Ek::Ins, Ek::Del],
];

struct GateOutcome {
    full: String,
    skip: String,
    gated: String,
    /// The ACTUAL shipped-gate decision on `normalize_core`'s output — whether it
    /// ran the full re-partition. This is the true fallback signal.
    gate_fired: bool,
}

/// Generate the dense-member corpus and run the three normalization paths on each.
/// Returns the collected outcomes and the number of inputs on which a
/// normalization **panicked** — tracked, not silently dropped, so a gate that
/// aborts on some shape cannot hide as a skipped row.
fn run_dense_member_corpus() -> (Vec<GateOutcome>, usize) {
    let mut outcomes = Vec::new();
    let mut panics = 0usize;
    for seq in dense_contigs() {
        let bytes = seq.as_bytes();
        let len = seq.len() as u64;
        let nz = Normalizer::new(provider(seq));
        for pattern in DENSE_MEMBER_PATTERNS {
            for sep in [0u64, 1, 2] {
                // Walk start positions across the contig interior.
                let mut start = 12u64;
                while start + 8 < len {
                    // Build the members left to right, non-overlapping.
                    let mut pos = start;
                    let mut members = Vec::new();
                    let mut ok = true;
                    for &kind in pattern.iter() {
                        match build_member(kind, bytes, pos) {
                            Some((m, width)) => {
                                members.push(m);
                                pos += width + sep;
                            }
                            None => {
                                ok = false;
                                break;
                            }
                        }
                    }
                    start += 3;
                    if !ok || pos as usize >= seq.len() - 4 {
                        continue;
                    }
                    let Some(input) = assemble(&members) else {
                        continue;
                    };
                    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
                        let Some(variant) = parse_hgvs(&input).ok() else {
                            continue;
                        };
                        let options = FromSequencesOptions::default().with_direction(direction);
                        let outcome =
                            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                                let derived = nz.rederive(&variant, &options, false).ok()?;
                                // `normalize_core`'s output — what the gate actually
                                // inspects — is the skip path's variant.
                                let skip_variant =
                                    nz.normalize_skipping_repartition(&derived).ok()?;
                                let gate_fired = nz.repartition_gate_fires(&skip_variant);
                                let full = nz.normalize(&derived).ok()?.to_string();
                                let gated =
                                    nz.normalize_gated_repartition(&derived).ok()?.to_string();
                                Some((full, skip_variant.to_string(), gated, gate_fired))
                            }));
                        match outcome {
                            Ok(Some((full, skip, gated, gate_fired))) => {
                                outcomes.push(GateOutcome {
                                    full,
                                    skip,
                                    gated,
                                    gate_fired,
                                });
                            }
                            // A clean decline (parse/derive/normalize returned `None`).
                            Ok(None) => {}
                            // A panic — count it so a gate that aborts is not silent.
                            Err(_) => panics += 1,
                        }
                    }
                }
            }
        }
    }
    (outcomes, panics)
}

/// Gate CORRECTNESS on the dense-member distribution class: the shipped gated Drop must be
/// byte-identical to the full normalization on every generated allele, and the raw
/// unconditional skip must still diverge somewhere (or the gate is untested here).
#[test]
fn the_gate_is_correct_on_dense_member_alleles() {
    let (outcomes, panics) = run_dense_member_corpus();
    let compared = outcomes.len();
    let skip_diverged = outcomes.iter().filter(|o| o.full != o.skip).count();
    let gated_diverged: Vec<&GateOutcome> = outcomes.iter().filter(|o| o.full != o.gated).collect();

    // No input may panic on any of the three paths — a gate that aborts on a
    // shape must fail loudly, not vanish as a dropped row.
    assert_eq!(
        panics, 0,
        "{panics} dense-member inputs panicked during normalization",
    );
    assert!(
        compared >= 5_000,
        "the dense-member corpus compared only {compared} alleles — generator regressed",
    );
    assert!(
        skip_diverged >= 1,
        "the raw skip never diverged over {compared} dense-member alleles, so this \
         guard is vacuous — the generator stopped producing the divergence shapes",
    );
    for o in gated_diverged.iter().take(20) {
        eprintln!("DENSE-GATED-DIVERGE\tfull={}\tgated={}", o.full, o.gated);
    }
    assert!(
        gated_diverged.is_empty(),
        "the gated Drop moved {} of {compared} dense-member outputs (raw skip moved \
         {skip_diverged}) — the gate let a divergence through",
        gated_diverged.len(),
    );
}

/// Report mode: the gate FALLBACK rate on the dense-member distribution — the speed cost.
#[test]
#[ignore = "permanent diagnostic (#2161): prints the gate fallback rate over the \
            dense-member corpus; never asserts"]
fn report_gate_fallback_on_dense_member_alleles() {
    let (outcomes, _panics) = run_dense_member_corpus();
    let compared = outcomes.len();
    // The REAL gate decision: how often the shipped gate ran the full re-partition.
    let fired = outcomes.iter().filter(|o| o.gate_fired).count();
    // How many actually NEEDED it (raw skip would have diverged).
    let needed = outcomes.iter().filter(|o| o.full != o.skip).count();
    // Wasted fallbacks: gate fired but the skip would have been correct anyway.
    let wasted = outcomes
        .iter()
        .filter(|o| o.gate_fired && o.full == o.skip)
        .count();
    // Member count of the gated input (normalize_core's output = skip).
    let members_of = |s: &str| -> usize {
        if s.contains('[') {
            s.matches(';').count() + 1
        } else {
            1
        }
    };
    let wasted_single = outcomes
        .iter()
        .filter(|o| o.gate_fired && o.full == o.skip && members_of(&o.skip) == 1)
        .count();
    let wasted_multi = outcomes
        .iter()
        .filter(|o| o.gate_fired && o.full == o.skip && members_of(&o.skip) >= 2)
        .count();
    for o in outcomes.iter().filter(|o| o.full != o.skip) {
        eprintln!("NEEDED\tskip={}\tfull={}", o.skip, o.full);
    }
    let mut seen = std::collections::BTreeSet::new();
    for o in outcomes
        .iter()
        .filter(|o| o.gate_fired && o.full == o.skip && members_of(&o.skip) >= 2)
    {
        if seen.insert(o.skip.clone()) && seen.len() <= 30 {
            eprintln!("WASTED-MULTI\tskip={}", o.skip);
        }
    }
    eprintln!(
        "WASTEDBREAKDOWN wasted_single_member={wasted_single} wasted_multi_member={wasted_multi}"
    );
    eprintln!(
        "DENSEFALLBACK compared={compared} gate_fired={fired} ({:.2}%) actually_needed={needed} \
         ({:.3}%) wasted_fallbacks={wasted} ({:.2}%)",
        100.0 * fired as f64 / compared.max(1) as f64,
        100.0 * needed as f64 / compared.max(1) as f64,
        100.0 * wasted as f64 / compared.max(1) as f64,
    );
}

/// Is a LONE (single-member) `delins` output ever changed by the second block
/// partition? If not, the gate need not fire on it — which removes ~80% of the
/// gate fallbacks on the dense-member distribution. The risk case is a `delins` whose
/// payload is the reverse complement of its reference span (which the
/// re-partition would re-type to `inv`), so this probe generates exactly those,
/// plus unchanged-interior and arbitrary payloads, over real bases, and asserts
/// the unconditional skip is byte-identical to the full re-partition on every one.
#[test]
fn a_lone_delins_is_a_fixed_point_of_the_repartition() {
    fn revcomp(s: &[u8]) -> String {
        s.iter()
            .rev()
            .map(|&b| match b {
                b'A' => 'T',
                b'T' => 'A',
                b'C' => 'G',
                b'G' => 'C',
                _ => 'N',
            })
            .collect()
    }
    let mut compared = 0usize;
    let mut panics = 0usize;
    let mut diverged: Vec<(String, String, String)> = Vec::new();
    for seq in dense_contigs() {
        let bytes = seq.as_bytes();
        let nz = Normalizer::new(provider(seq));
        // Every span [a, b] with 2..=6 reference bases, interior enough not to hit
        // the contig edge.
        for a in 12u64..(seq.len() as u64 - 12) {
            for width in 2u64..=6 {
                let b = a + width - 1;
                if b as usize + 8 >= seq.len() {
                    continue;
                }
                let span = &bytes[(a - 1) as usize..b as usize];
                // Payloads that stress each re-typing route.
                let payloads = [
                    revcomp(span),              // -> should become inv
                    "A".repeat(width as usize), // homopolymer replacement
                    {
                        // unchanged-interior: keep first/last ref base, change middle
                        let mut p: Vec<u8> = span.to_vec();
                        if p.len() >= 3 {
                            let mid = p.len() / 2;
                            p[mid] = if p[mid] == b'A' { b'C' } else { b'A' };
                        }
                        String::from_utf8(p).unwrap()
                    },
                ];
                for payload in payloads {
                    let input = format!("NC_TEST.1:g.{a}_{b}delins{payload}");
                    let Ok(variant) = parse_hgvs(&input) else {
                        continue;
                    };
                    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
                        let options = FromSequencesOptions::default().with_direction(direction);
                        let outcome =
                            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                                let derived = nz.rederive(&variant, &options, false).ok()?;
                                // Only lone (single-member) outputs are in scope.
                                if matches!(derived, ferro_hgvs::HgvsVariant::Allele(_)) {
                                    return None;
                                }
                                let skip = nz
                                    .normalize_skipping_repartition(&derived)
                                    .ok()?
                                    .to_string();
                                if skip.contains('[') {
                                    return None; // became multi-member; out of scope
                                }
                                // Only lone delins outputs are the question here.
                                if !skip.contains("delins") {
                                    return None;
                                }
                                let full = nz.normalize(&derived).ok()?.to_string();
                                Some((full, skip))
                            }));
                        match outcome {
                            Ok(Some((full, skip))) => {
                                compared += 1;
                                if full != skip {
                                    diverged.push((input.clone(), skip, full));
                                }
                            }
                            // A clean decline (out of scope: multi-member, or not a
                            // lone delins). Legitimate, not a fixed-point failure.
                            Ok(None) => {}
                            // A panic — count it so a normalization that aborts is
                            // not silently dropped from `compared`, which would make
                            // the "0 divergences" floor read better than the truth.
                            Err(_) => panics += 1,
                        }
                    }
                }
            }
        }
    }
    for (input, skip, full) in diverged.iter().take(30) {
        eprintln!("LONE-DELINS-DIVERGE\tin={input}\tskip={skip}\tfull={full}");
    }
    eprintln!(
        "LONEDELINS compared={compared} diverged={} panics={panics}",
        diverged.len()
    );
    assert_eq!(
        panics, 0,
        "{panics} lone-delins probe inputs panicked during normalization — a gate \
         measured over a corpus that swallows panics is not measured",
    );
    assert!(
        compared >= 200,
        "the lone-delins probe compared only {compared} outputs — generator too narrow",
    );
    assert!(
        diverged.is_empty(),
        "{} lone delins outputs were CHANGED by the re-partition (see \
         LONE-DELINS-DIVERGE) — the gate MUST fire on lone delins; not firing on \
         them is unsafe",
        diverged.len(),
    );
}

/// Is a LONE (single-member) `inv` output ever changed by the second block
/// partition? A lone inv has no sibling to be re-placed against and is already a
/// whole-span reverse complement (the canonical `inv`), so it should be a fixed
/// point. Confirming it lets the gate fire on inv only in MULTI-member alleles.
#[test]
fn a_lone_inv_is_a_fixed_point_of_the_repartition() {
    let mut compared = 0usize;
    let mut panics = 0usize;
    let mut diverged: Vec<(String, String, String)> = Vec::new();
    for seq in dense_contigs() {
        let bytes = seq.as_bytes();
        let nz = Normalizer::new(provider(seq));
        for a in 12u64..(seq.len() as u64 - 12) {
            for width in 2u64..=8 {
                let b = a + width - 1;
                if b as usize + 8 >= seq.len() {
                    continue;
                }
                let Some(member) = inv(bytes, a, b) else {
                    continue; // palindrome: not an inversion
                };
                let input = format!("NC_TEST.1:g.{}", member.text);
                let Ok(variant) = parse_hgvs(&input) else {
                    continue;
                };
                for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
                    let options = FromSequencesOptions::default().with_direction(direction);
                    let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        let derived = nz.rederive(&variant, &options, false).ok()?;
                        if matches!(derived, ferro_hgvs::HgvsVariant::Allele(_)) {
                            return None;
                        }
                        let skip = nz
                            .normalize_skipping_repartition(&derived)
                            .ok()?
                            .to_string();
                        if skip.contains('[') || !skip.contains("inv") {
                            return None; // became multi-member or not a lone inv
                        }
                        let full = nz.normalize(&derived).ok()?.to_string();
                        Some((full, skip))
                    }));
                    match outcome {
                        Ok(Some((full, skip))) => {
                            compared += 1;
                            if full != skip {
                                diverged.push((input.clone(), skip, full));
                            }
                        }
                        // A clean decline (out of scope: multi-member, or not a
                        // lone inv). Legitimate, not a fixed-point failure.
                        Ok(None) => {}
                        // A panic — count it so a normalization that aborts is not
                        // silently dropped from `compared`, which would make the
                        // "0 divergences" floor read better than the truth.
                        Err(_) => panics += 1,
                    }
                }
            }
        }
    }
    for (input, skip, full) in diverged.iter().take(30) {
        eprintln!("LONE-INV-DIVERGE\tin={input}\tskip={skip}\tfull={full}");
    }
    eprintln!(
        "LONEINV compared={compared} diverged={} panics={panics}",
        diverged.len()
    );
    assert_eq!(
        panics, 0,
        "{panics} lone-inv probe inputs panicked during normalization — a gate \
         measured over a corpus that swallows panics is not measured",
    );
    assert!(
        compared >= 100,
        "the lone-inv probe compared only {compared} outputs — generator too narrow",
    );
    assert!(
        diverged.is_empty(),
        "{} lone inv outputs were CHANGED by the re-partition (see LONE-INV-DIVERGE) \
         — the gate must fire on lone inv",
        diverged.len(),
    );
}
