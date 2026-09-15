//! Critical-pair census over the ruled-cut fixpoint rule set (design §4.6,
//! migration step 6).
//!
//! Two rules form a **critical pair** on a cut when both fire and, after each
//! branch is closed with the full engine, the two normal forms render
//! differently. Because the Φ-decreasing engine terminates, Newman's lemma makes
//! "every such local diamond rejoins" equivalent to "the engine's normal form is
//! order-independent" — so a surviving critical pair is a place where the
//! registration order is load-bearing, and design §4.6 requires each one to carry
//! a **decided priority record**, cited by the higher-priority rule's authority.
//!
//! **The census now finds NO critical pair on the sweep — the engine is confluent
//! over it.** Before the engine fixes it found exactly two, both a **whole-span
//! reverse complement** where `run-inv` competed with a coincidence coarsen
//! (`compensating-gaps`) or a tandem-dup recut (`tandem-dup-run`):
//!
//! ```text
//!   compensating-gaps X run-inv        GACA -> TGTC
//!   run-inv           X tandem-dup-run GCA  -> TGC
//! ```
//!
//! Both are **pre-empted at the seed** by F1: `Seed::canonical` types a whole-span
//! reverse complement as ONE locked `Inv` before any rule runs
//! (`whole-span-reverse-complement-types-as-inv`), so the DAG never cuts it into
//! the competing substitutions/delins that `run-inv` and the coarsen/dup used to
//! race over — the diamond never opens. F4 removes the other way a diamond could
//! open here: `run-inv` declines its sep-0 flanked route (a `[del; inv; del]` that
//! `sep-zero` owns as a delins under
//! `delins-adjacent-members-when-both-consume-reference`), so no `Inv` × coarsen
//! divergence survives. The two former witnesses are kept in `FORMER_WHOLE_SPAN`
//! and each is pinned to now resolve to a single locked `Inv`
//! (`the_former_whole_span_pairs_now_resolve_to_a_locked_inv`), which is F1's
//! realization of the ruling on the ruled arm.
//!
//! **This changes the ruled arm, not the shipping default.** The engine fixes are
//! reached only via `partition_ruled`, so the DEFAULT partition stays
//! byte-identical delins; the ruled arm realizing `inv` here is the ruling's own
//! owed, separately-disclosed representation change, measured pre-flip.
//!
//! **The empty census is the confluence pin.** A future rule change that reopens a
//! diamond (breaks confluence without a record) adds a pair and fails
//! `the_critical_pair_census_is_empty_so_the_engine_is_confluent`; a non-vacuity
//! floor keeps "empty" from meaning "nothing ran".
//!
//! **Instrument choice — full-engine (Newman) closure, not pairwise 2-rule
//! commutation.** Design §4.6 phrases the check as "A∘B and B∘A reach the same
//! normal form". This census closes each branch with the *whole* rule set rather
//! than the pair alone, so a pair that diverges under the pair in isolation but
//! reconverges through a third rule is not flagged — the question that matters for
//! the fixpoint driver is whether the ENGINE's output is order-independent, and
//! Newman's lemma makes that exactly "every local diamond rejoins under the full
//! system". Two consequences worth stating, because they differ from what the
//! design's R6 anticipated under the stricter 2-rule reading: the codon-precedence
//! pairs (`codon-exception-vs-coincidence-carve-out-precedence`,
//! `codon-carve-out-excludes-the-compensated-zero-width-gap-two-pair`) reconverge
//! and are not flagged — their priority is already carried by registration order —
//! and the `Inv` × `sep-zero` question (#2007) likewise reconverges, because step
//! 5's lock makes `sep-zero` decline a locked `Inv`, so the engine reaches one
//! normal form regardless of order. #2007 is therefore answered in code here; its
//! ledger record (step 4) documents a question this instrument shows the engine
//! does not actually leave open.

use ferro_hgvs::partition::adapters::{
    CodonException, CompensatingGaps, PayloadCoincidence, RunInv, TandemDupRun,
};
use ferro_hgvs::partition::block_ctx::{BlockCtx, FrameContext, Molecule, Provenance};
use ferro_hgvs::partition::driver::partition_ruled;
use ferro_hgvs::partition::ruled::{run_engine, Authority, Cut, Rule};
use ferro_hgvs::partition::rules::SepZero;
use ferro_hgvs::partition::seed::Seed;
use ferro_hgvs::ShuffleDirection;
use std::collections::{BTreeMap, BTreeSet};

use crate::common::fixture_gen::fixture_path;

/// The six fixpoint rules, in the driver's registration (priority) order — the
/// same set and order `driver::partition_ruled` runs to a fixed point.
fn fixpoint_rules() -> [&'static dyn Rule; 6] {
    [
        &SepZero,
        &PayloadCoincidence,
        &CompensatingGaps,
        &TandemDupRun,
        &RunInv,
        &CodonException,
    ]
}

/// A pinned critical pair: the two rule ids (in `pair_key` order), the decided
/// ledger record that governs it, and a witness block that reproduces the
/// divergence deterministically (independent of the sweep depth).
struct ExpectedPair {
    a: &'static str,
    b: &'static str,
    ruling: &'static str,
    reference: &'static [u8],
    resulting: &'static [u8],
}

/// The measured critical set (design §4.6). **Empty**: F1's whole-span `Inv` seed
/// pre-empts both former pairs at the seed, and F4's flanked-route decline removes
/// the other way an `Inv` diamond could open, so the engine is confluent over the
/// sweep. See the module doc.
const EXPECTED: &[ExpectedPair] = &[];

/// The two pairs the census found BEFORE the engine fixes, both a whole-span
/// reverse complement. They are no longer critical (F1 seeds them as one locked
/// `Inv`), and `the_former_whole_span_pairs_now_resolve_to_a_locked_inv` pins that
/// each now resolves to a single `Inv` on the ruled arm — the realization of
/// `run-inv`'s cited authority.
const FORMER_WHOLE_SPAN: &[ExpectedPair] = &[
    ExpectedPair {
        a: "compensating-gaps",
        b: "run-inv",
        ruling: "whole-span-reverse-complement-types-as-inv",
        reference: b"GACA",
        resulting: b"TGTC",
    },
    ExpectedPair {
        a: "run-inv",
        b: "tandem-dup-run",
        ruling: "whole-span-reverse-complement-types-as-inv",
        reference: b"GCA",
        resulting: b"TGC",
    },
];

/// Sweep depth: all ACGT blocks of length `1..=SWEEP_LEN`. Default 4 (reproduces
/// both pinned pairs — `GACA->TGTC` is a 4-mer); `FERRO_CRITICAL_PAIR_FULL` widens
/// it to 6 as a deeper net for a new pair the default depth would miss.
fn sweep_len() -> usize {
    if std::env::var("FERRO_CRITICAL_PAIR_FULL").is_ok() {
        6
    } else {
        4
    }
}

/// All cuts reachable from `seed` by single rule applications (BFS to closure).
/// Finite: every application strictly decreases Φ.
fn reachable(seed: &Cut, ctx: &BlockCtx, rules: &[&dyn Rule]) -> Vec<Cut> {
    let mut seen: Vec<Cut> = vec![seed.clone()];
    let mut frontier = vec![seed.clone()];
    while let Some(cut) = frontier.pop() {
        for r in rules {
            if !r.scope().admits(ctx) {
                continue;
            }
            if let Some(next) = r.apply(&cut, ctx) {
                if !seen.contains(&next) {
                    seen.push(next.clone());
                    frontier.push(next);
                }
            }
        }
    }
    seen
}

/// At `cut`, if both `a` and `b` are admitted and both fire, close each branch
/// with the full engine and return `true` when the two normal forms **render**
/// differently. `to_partition` drops the partition-stage label/lock bookkeeping,
/// so a difference in it is a genuine output divergence — the inv/dup/delins
/// distinction survives, but two forms reached with different internal labels do
/// not count.
fn diverges(cut: &Cut, ctx: &BlockCtx, a: &dyn Rule, b: &dyn Rule, all: &[&dyn Rule]) -> bool {
    if !a.scope().admits(ctx) || !b.scope().admits(ctx) {
        return false;
    }
    let (Some(ca), Some(cb)) = (a.apply(cut, ctx), b.apply(cut, ctx)) else {
        return false;
    };
    run_engine(ca, ctx, all).to_partition() != run_engine(cb, ctx, all).to_partition()
}

/// The frames swept: genomic/noncoding, the three coding grid phases, and a
/// coding frame carrying a CDS-end seam mid-window (the geometry step 5's moves
/// live on).
fn sweep_frames(len: usize) -> Vec<FrameContext> {
    let mut frames = vec![
        FrameContext::NonCoding,
        FrameContext::coding(0, None, None),
        FrameContext::coding(2, Some(1), None),
        FrameContext::coding(1, Some(2), None),
    ];
    if len >= 2 {
        frames.push(FrameContext::coding(0, None, Some(len / 2)));
    }
    frames
}

/// Every ACGT string of length `1..=max_len`.
fn sequences(max_len: usize) -> Vec<Vec<u8>> {
    const ALPHABET: &[u8] = b"ACGT";
    let mut out = Vec::new();
    for len in 1..=max_len {
        for mut n in 0..ALPHABET.len().pow(len as u32) {
            let mut s = vec![0u8; len];
            for slot in s.iter_mut() {
                *slot = ALPHABET[n % ALPHABET.len()];
                n /= ALPHABET.len();
            }
            out.push(s);
        }
    }
    out
}

/// Canonical order-independent pair key by rule id.
fn pair_key(a: &dyn Rule, b: &dyn Rule) -> (&'static str, &'static str) {
    if a.id() <= b.id() {
        (a.id(), b.id())
    } else {
        (b.id(), a.id())
    }
}

/// One `(reference, resulting)` witness per critical pair found on the sweep, plus
/// the counts the non-vacuity floor reads.
struct SweepResult {
    critical: BTreeMap<(&'static str, &'static str), (String, String)>,
    seeds_examined: usize,
    multi_applicable_visits: usize,
}

/// Run the exhaustive small-block sweep and return the critical set it finds.
fn sweep(max_len: usize) -> SweepResult {
    let rules = fixpoint_rules();
    let rule_slice: Vec<&dyn Rule> = rules.to_vec();
    let mut critical = BTreeMap::new();
    let mut seeds_examined = 0usize;
    let mut multi_applicable_visits = 0usize;

    let seqs = sequences(max_len);
    for reference in &seqs {
        for resulting in &seqs {
            if reference == resulting {
                continue;
            }
            for frame in sweep_frames(reference.len().max(resulting.len())) {
                let prov = Provenance::none();
                let ctx = BlockCtx {
                    reference,
                    resulting,
                    frame: &frame,
                    molecule: Molecule::Dna,
                    provenance: &prov,
                };
                let Ok(seed) = Seed::canonical(&ctx) else {
                    continue;
                };
                seeds_examined += 1;
                for cut in reachable(&seed, &ctx, &rule_slice) {
                    let applicable = rules
                        .iter()
                        .filter(|r| r.scope().admits(&ctx) && r.apply(&cut, &ctx).is_some())
                        .count();
                    if applicable >= 2 {
                        multi_applicable_visits += 1;
                    }
                    for i in 0..rules.len() {
                        for j in (i + 1)..rules.len() {
                            let (a, b) = (rules[i], rules[j]);
                            if diverges(&cut, &ctx, a, b, &rule_slice) {
                                critical.entry(pair_key(a, b)).or_insert_with(|| {
                                    (
                                        String::from_utf8_lossy(reference).into_owned(),
                                        String::from_utf8_lossy(resulting).into_owned(),
                                    )
                                });
                            }
                        }
                    }
                }
            }
        }
    }
    SweepResult {
        critical,
        seeds_examined,
        multi_applicable_visits,
    }
}

/// The decided ruling ids in the committed ledger.
fn decided_ruling_ids() -> BTreeSet<String> {
    #[derive(serde::Deserialize)]
    struct LedgerRuling {
        id: String,
        status: String,
    }
    #[derive(serde::Deserialize)]
    struct Ledger {
        #[serde(default)]
        rulings: Vec<LedgerRuling>,
    }
    let path = fixture_path("tests/fixtures/grammar/hgvs_spec_normalization_overrides.json");
    let text = std::fs::read_to_string(&path).expect("read ledger");
    let ledger: Ledger = serde_json::from_str(&text).expect("parse ledger");
    ledger
        .rulings
        .into_iter()
        .filter(|r| r.status == "decided")
        .map(|r| r.id)
        .collect()
}

/// The ledger ruling a rule cites, if it cites one (rather than a clause or a
/// house choice).
fn cited_ruling(rule: &dyn Rule) -> Option<&'static str> {
    match rule.authority() {
        Authority::Ruling(id) => Some(id),
        _ => None,
    }
}

fn rule_by_id(id: &str) -> &'static dyn Rule {
    *fixpoint_rules()
        .iter()
        .find(|r| r.id() == id)
        .unwrap_or_else(|| panic!("no fixpoint rule named {id}"))
}

fn ctx_noncoding<'a>(
    reference: &'a [u8],
    resulting: &'a [u8],
    frame: &'a FrameContext,
    prov: &'a Provenance,
) -> BlockCtx<'a> {
    BlockCtx {
        reference,
        resulting,
        frame,
        molecule: Molecule::Dna,
        provenance: prov,
    }
}

#[test]
fn the_critical_pair_census_is_empty_so_the_engine_is_confluent() {
    let result = sweep(sweep_len());

    // Non-vacuity: the sweep really enumerated blocks and reached cuts where more
    // than one rule was live, so an empty critical set means "confluent", not
    // "nothing ran".
    assert!(
        result.seeds_examined > 100_000,
        "the sweep barely ran: {} seeds",
        result.seeds_examined
    );
    assert!(
        result.multi_applicable_visits > 0,
        "no cut had two rules applicable — the census cannot see a critical pair"
    );

    // The measured set is exactly the pinned set (empty): the engine is confluent
    // over the sweep. A NEW critical pair (a rule change that reopens a diamond
    // without a record) adds a key and fails here.
    let measured: BTreeSet<(&str, &str)> = result.critical.keys().copied().collect();
    let expected: BTreeSet<(&str, &str)> = EXPECTED.iter().map(|e| (e.a, e.b)).collect();
    assert_eq!(
        measured, expected,
        "critical-pair set drifted; witnesses = {:?}",
        result.critical
    );

    // The two pairs the census found before the engine fixes no longer diverge:
    // F1 seeds each whole-span reverse complement as ONE run, so no diamond opens.
    // (Sabotage: delete the F1 branch in `Seed::canonical` and each witness's seed
    // splits into competing members again, resurrecting the divergence.) Each still
    // carries a decided governing record, cited by run-inv's authority — F1 makes
    // that seed decision eagerly rather than leaving it to the fixpoint.
    let all: Vec<&dyn Rule> = fixpoint_rules().to_vec();
    let decided = decided_ruling_ids();
    for e in FORMER_WHOLE_SPAN {
        let frame = FrameContext::NonCoding;
        let prov = Provenance::none();
        let ctx = ctx_noncoding(e.reference, e.resulting, &frame, &prov);
        let seed = Seed::canonical(&ctx).expect("witness seeds");
        assert_eq!(
            seed.runs().len(),
            1,
            "F1 must seed {}->{} as one run, not the competing partition",
            String::from_utf8_lossy(e.reference),
            String::from_utf8_lossy(e.resulting),
        );
        let (a, b) = (rule_by_id(e.a), rule_by_id(e.b));
        let hit = reachable(&seed, &ctx, &all)
            .iter()
            .any(|cut| diverges(cut, &ctx, a, b, &all));
        assert!(
            !hit,
            "{} X {} still diverges on {}->{} — the diamond was not pre-empted",
            e.a,
            e.b,
            String::from_utf8_lossy(e.reference),
            String::from_utf8_lossy(e.resulting)
        );

        // The governing record is decided and is the authority run-inv cites.
        assert!(
            decided.contains(e.ruling),
            "governing record {} for {} X {} is not decided in the ledger",
            e.ruling,
            e.a,
            e.b
        );
        assert_eq!(
            cited_ruling(&RunInv),
            Some(e.ruling),
            "the governing record must be the authority run-inv cites",
        );
    }
}

fn render(cut: &Cut) -> String {
    cut.to_partition()
        .members
        .iter()
        .map(|m| {
            format!(
                "{:?}[{}..{}]={}",
                m.kind,
                m.ref_start,
                m.ref_end,
                String::from_utf8_lossy(&m.inserted)
            )
        })
        .collect::<Vec<_>>()
        .join(" ")
}

#[test]
fn the_former_whole_span_pairs_now_resolve_to_a_locked_inv() {
    // Pins the post-F1 resolution: each former critical witness is a whole-span
    // reverse complement, so `Seed::canonical` types it as ONE locked `Inv` and the
    // engine's normal form is that single `Inv` — never the spanning `delins` the
    // pre-F1 engine produced. This realizes `whole-span-reverse-complement-types-as-inv`
    // on the RULED arm; the shipping default (reached without `partition_ruled`) is
    // unchanged and still delins — a separately-disclosed representation change.
    use ferro_hgvs::partition::ruled::Label;
    let all: Vec<&dyn Rule> = fixpoint_rules().to_vec();
    let frame = FrameContext::NonCoding;
    let prov = Provenance::none();
    for e in FORMER_WHOLE_SPAN {
        let ctx = ctx_noncoding(e.reference, e.resulting, &frame, &prov);
        let seed = Seed::canonical(&ctx).expect("seed");
        let nf = run_engine(seed, &ctx, &all);
        assert_eq!(
            nf.runs().len(),
            1,
            "{}->{} resolves to a single run, got {}",
            String::from_utf8_lossy(e.reference),
            String::from_utf8_lossy(e.resulting),
            render(&nf)
        );
        assert_eq!(
            nf.runs()[0].label,
            Label::Inv,
            "{}->{} resolves to an inv, got {}",
            String::from_utf8_lossy(e.reference),
            String::from_utf8_lossy(e.resulting),
            render(&nf)
        );
        assert_eq!(
            nf.to_partition().members[0].kind,
            ferro_hgvs::partition::output::EditKind::Inv,
            "{}->{} renders as an inversion, got {}",
            String::from_utf8_lossy(e.reference),
            String::from_utf8_lossy(e.resulting),
            render(&nf)
        );
    }
}

#[test]
fn the_ruled_arm_member_count_is_direction_independent_step_7() {
    // Step 7 retired the direction mirror; the 5' arm is now the 3' partition
    // re-placed 5' (`bridge::place_only`), so the member count cannot depend on the
    // shuffle direction — by construction, over the exhaustive small-block sweep.
    // This is the ruled-arm analogue of `merge.rs`'s `direction_symmetry` module,
    // which sweeps only `PARTITION_RULE_NAMES` (the ruled arm is not among them). A
    // raw count guards against a `place_only` that silently drops or merges a run.
    //
    // Swept over EVERY frame — noncoding, the three coding phases, and a CDS-end
    // seam. Coding frames also exercise the codon exception's lateral, idempotent
    // Coarsen (`AAAA->CACC` re-partitions at a constant run count), which used to
    // trip the engine's strict-Φ `debug_assert` and forced this sweep to noncoding.
    // The engine now checks PROGRESS (`Φ` non-increasing AND the cut changed) rather
    // than strict decrease (`run_engine`'s comment), so the coding sweep runs — and
    // the direction-symmetry property still holds by construction (`place_only` is a
    // pure coordinate shift, frame-independent).
    let prov = Provenance::none();
    let seqs = sequences(4);
    let mut compared = 0usize;
    for reference in &seqs {
        for resulting in &seqs {
            if reference == resulting {
                continue;
            }
            for frame in sweep_frames(reference.len().max(resulting.len())) {
                let ctx = BlockCtx {
                    reference,
                    resulting,
                    frame: &frame,
                    molecule: Molecule::Dna,
                    provenance: &prov,
                };
                let (Ok(three), Ok(five)) = (
                    partition_ruled(&ctx, ShuffleDirection::ThreePrime),
                    partition_ruled(&ctx, ShuffleDirection::FivePrime),
                ) else {
                    continue; // a declined block declines in both directions
                };
                assert_eq!(
                    three.runs().len(),
                    five.runs().len(),
                    "member count depends on direction: {}->{} [{frame:?}]  3'={}  5'={}",
                    String::from_utf8_lossy(reference),
                    String::from_utf8_lossy(resulting),
                    render(&three),
                    render(&five),
                );
                compared += 1;
            }
        }
    }
    assert!(
        compared > 10_000,
        "the sweep barely ran: {compared} comparisons"
    );
}
