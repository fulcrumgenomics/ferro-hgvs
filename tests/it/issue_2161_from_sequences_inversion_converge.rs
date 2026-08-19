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
//! - [`the_inversion_geometry_corpus_converges_without_regression`] — the corpus
//!   is fully built (a fixed count) and the number of convergent variants holds a
//!   floor: a regression lowers it, a later increment only raises it.
//!
//! The guard is proven able to fail: on the pre-fix tree every curated row above
//! diverges (the `NC_TEST.1:g.[14del;21_23inv]` reproducer is one of them), and
//! the corpus carried 26,481 divergences against the shipped 9,170.

use ferro_hgvs::{parse_hgvs, FromSequencesOptions, JsonProvider, Normalizer, ShuffleDirection};
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
}

/// Generate the whole corpus and run every variant through both surfaces.
///
/// Every generated variant is a flanked inversion — a genuine `inv` span beside
/// at least one `del`/`sub` sibling — so the count of outcomes is the geometry
/// this file guards, floored below.
fn run_corpus() -> Vec<Outcome> {
    let mut outcomes = Vec::new();

    for (_name, seq) in contigs() {
        let bytes = seq.as_bytes();
        let len = seq.len() as u64;
        let nz = Normalizer::new(provider(seq));

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
                            if let Some(p) = five {
                                match sibling(p) {
                                    Some(m) => members.push(m),
                                    None => continue,
                                }
                            }
                        }
                        // The inversion always participates.
                        members.push(Member {
                            text: inv_member.text.clone(),
                            start: inv_member.start,
                            end: inv_member.end,
                        });
                        if use3 {
                            if let Some(p) = three {
                                match sibling(p) {
                                    Some(m) => members.push(m),
                                    None => continue,
                                }
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
                            let derive = rederive(&nz, &input, direction, false);
                            let recommend = rederive(&nz, &input, direction, true);
                            outcomes.push(Outcome {
                                input: input.clone(),
                                direction,
                                compared: derive.zip(recommend),
                            });
                        }
                    }
                }
            }
        }
    }

    outcomes
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
        let nz = Normalizer::new(provider(seq));
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
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

/// The geometry corpus, as a regression floor. Increments 1 (the flanked/interior
/// inversion port) and its over-merge follow-up (the `decompose_delins` split)
/// took the pre-port tree's divergences on this corpus from **26,481** down to
/// **8,100** and opened **none** (measured baseline-vs-fix by `input`+`direction`:
/// 0 rows that converged before now diverge); what remains is
/// the shift and over-merge convergence that increments 2–3 address. So this
/// asserts floors, not full convergence:
///
/// - the corpus is fully built (a fixed, deterministic count — a drop means the
///   generator stopped emitting the geometry, the vacuity this file guards
///   against);
/// - `normalize` actually types an inversion on a large share of the corpus (the
///   non-vacuity that matters: a generator that stopped producing genuine flanked
///   inversions — every span palindromic, say — would still build 64,352 rows but
///   type almost no `inv`, and the convergence floor alone could not tell);
/// - the number of variants on which `derive` converges with `normalize` does not
///   fall below what the ports achieved. A future regression lowers it; a later
///   increment that closes the remaining gaps only raises it, so the floor is
///   `>=`, not `==`.
///
/// Measured after the over-merge increment: 64,352 comparisons, 0 skipped, 56,252
/// convergent, 40,900 carrying an `inv` in the `normalize` output. Of the
/// remaining 8,100 divergences, all but **252** are FivePrime — where `normalize`
/// emits its direction-invariant 3'-canonical form and the derivation respects
/// the caller's direction, which is direction semantics, not a convergence gap.
/// The 252 ThreePrime residue is entirely repeat notation (`X[n]`), the next
/// increment; on ThreePrime — the direction the second-partition drop depends on —
/// inversion typing and over-merge are both closed.
#[test]
fn the_inversion_geometry_corpus_converges_without_regression() {
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

    // The corpus is hermetic and deterministic, so this is exact: every generated
    // variant put both surfaces through cleanly (0 skipped on the port).
    assert_eq!(
        compared, 64_352,
        "the geometry corpus changed size — the generator is no longer emitting \
         what this guard measures",
    );
    // Non-vacuity: the geometry really does produce inversions `normalize` types,
    // so the convergence floor below is a claim about flanked inversions and not
    // about a corpus that degenerated to palindromes or lone indels.
    assert!(
        recommend_has_inv >= 40_900,
        "the corpus stopped producing genuine inversions: only {recommend_has_inv} \
         of {compared} carry an inv in the normalize output, floor is 40,900",
    );
    assert!(
        converged >= 56_252,
        "convergence regressed: {converged} of {compared} converge, floor is 56,252",
    );
}

/// Report mode: print the divergence landscape. Not an assertion — run with
/// `--no-capture` while developing the corpus or a later increment, then read the
/// floors off it.
#[test]
#[ignore = "permanent diagnostic (#2161): prints the divergence landscape used to set \
            the floor in the_inversion_geometry_corpus_converges_without_regression; it \
            never asserts, so it stays ignored"]
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
}
