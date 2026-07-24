//! Property-based normalization idempotency test over reference-backed
//! synthetic fixtures (issue #1157 follow-up campaign).
//!
//! The pre-existing idempotency coverage is reference-free (empty provider, so
//! nothing shifts) or substitution-only (subs never shift). This fuzzes random
//! **repeat-biased** cores (homopolymers + short tandems, where del/dup/ins
//! actually 3'-shift) across the genomic and CDS (both strands) coordinate
//! systems, and every edit type, asserting `norm(norm(x)) == norm(x)`.
//!
//! Two properties are fuzzed:
//!   * `norm(norm(x)) == norm(x)` (idempotency), covering the full edit/coord
//!     matrix in the **3' (default, HGVS-canonical)** direction;
//!   * multiple spellings of ONE homopolymer-boundary haplotype normalize to ONE
//!     form (confluence), in **both** directions. Idempotency is necessary but
//!     NOT sufficient — a mode can be idempotent yet non-confluent (each spelling
//!     its own fixed point), which the confluence property catches. The 5'
//!     boundary/homopolymer shuffle was reworked to a coordinate-identity clamp
//!     keyed on the post-shuffle edit (see `src/normalize/mod.rs` CDS-start/-end
//!     clamps + the dup-routing recursion in `normalize_na_edit`) so those
//!     spellings now converge.
//!
//! Scope note: a *general* 5' idempotency fuzz is intentionally NOT run here yet.
//! The 5' direction still has separate, pre-existing bugs outside the boundary
//! rework — a `convert-to-dup` step that selects the wrong tract for
//! heterogeneous content, and a tandem-repeat `dup` whose repeat unit rotates on
//! re-normalization — each tracked as its own follow-up. The confluence property
//! plus the explicit `five_prime_boundary_delins_unification` regression tests
//! cover the 5' classes this change actually fixes.
//!
//! Cases are generated **valid by construction** (edit endpoints derived from
//! the core length), so there is no skip path at all: an unparseable render or a
//! failing first-pass normalize is a bug in the generator or the normalizer and
//! fails the test rather than silently shrinking the exercised set. Failure
//! seeds are recorded by proptest under
//! `tests/it/proptest-regressions/normalize_idempotency_proptest.txt`.

use crate::common::synthetic::{SyntheticBuilder, PAD_OFFSET};
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};
use proptest::prelude::*;
use proptest::test_runner::{Config as ProptestConfig, TestCaseError};

fn normalizer_for(core: &str, sys: Sys, dir: ShuffleDirection) -> Normalizer<MockProvider> {
    Normalizer::with_config(
        provider_for(core, sys),
        NormalizeConfig::default().with_direction(dir),
    )
}

fn both_directions() -> impl Strategy<Value = ShuffleDirection> {
    prop::sample::select(vec![
        ShuffleDirection::ThreePrime,
        ShuffleDirection::FivePrime,
    ])
}

/// Coordinate system under test. All use `cds_start == 1` / `cds_end == len`
/// for the CDS cases so `c.k` maps directly to core position `k` and both
/// transcript boundaries are exercised.
#[derive(Debug, Clone, Copy)]
enum Sys {
    Genomic,
    CdsPlus,
    CdsMinus,
}

#[derive(Debug, Clone)]
struct Case {
    core: String,
    sys: Sys,
    hgvs: String,
}

fn dna(n: std::ops::RangeInclusive<usize>) -> impl Strategy<Value = String> {
    prop::collection::vec(prop::sample::select(vec!['A', 'C', 'G', 'T']), n)
        .prop_map(|v| v.into_iter().collect())
}

/// Repeat-biased core: 2-5 runs, each a 1-3bp unit repeated 1-5 times, so
/// homopolymers and short tandems dominate. Length (of the body) clamped to
/// 6..=40.
///
/// The body is wrapped in fixed `G…C` breaker flanks. `SyntheticBuilder` pads
/// the genomic contig with a perfect `ACGT…` tandem, and its base immediately
/// 5' of the core is `T`, immediately 3' is `A`. A body that starts with `A` or
/// ends with `T` (common for `ACG`-ish cores) would *continue* that pad tandem,
/// creating a 256bp+ repeat an edit can 3'-shift into — but only ~100bp per pass
/// (the normalizer window), which is a window-limited-tandem artifact of the
/// harness, not a shuffle-logic bug. `G` (!= `A`, and `T`+`G` is no run) and `C`
/// (!= `T`, and `C`+`A` is no run) break the pad period on both sides, so all
/// shuffling stays inside the <=42bp core, well within the window.
fn core_strategy() -> impl Strategy<Value = String> {
    prop::collection::vec((dna(1..=3), 1usize..=5), 2..=5)
        .prop_map(|runs| {
            runs.into_iter()
                .flat_map(|(u, n)| u.repeat(n).chars().collect::<Vec<_>>())
                .collect::<String>()
        })
        .prop_filter("body length 6..=40", |s| s.len() >= 6 && s.len() <= 40)
        .prop_map(|body| format!("G{body}C"))
}

fn case_strategy() -> impl Strategy<Value = Case> {
    core_strategy().prop_flat_map(move |core| {
        let len = core.len();
        let sys = prop::sample::select(vec![Sys::Genomic, Sys::CdsPlus, Sys::CdsMinus]);
        // kind: 0=del 1=dup 2=ins 3=delins
        (sys, 0u8..=3, 1usize..=len, 0usize..=3, dna(1..=3)).prop_map(
            move |(sys, kind, start, span, ins)| {
                let end = (start + span).min(len);
                let (s, e) = (start.min(end), start.max(end));
                let render_pos = |p: usize| -> u64 {
                    match sys {
                        Sys::Genomic => (PAD_OFFSET as usize + p) as u64,
                        _ => p as u64, // cds_start == 1 => c.p maps to core p
                    }
                };
                let (prefix, acc) = match sys {
                    Sys::Genomic => ("g.", "NC_TEST.1"),
                    _ => ("c.", "NM_TEST.1"),
                };
                let range = |s: usize, e: usize| {
                    if s == e {
                        format!("{}", render_pos(s))
                    } else {
                        format!("{}_{}", render_pos(s), render_pos(e))
                    }
                };
                let body = match kind {
                    0 => format!("{}del", range(s, e)),
                    1 => format!("{}dup", range(s, e)),
                    2 => {
                        // insertion sits between p and p+1; keep p in 1..len
                        let p = s.min(len - 1);
                        format!("{}_{}ins{}", render_pos(p), render_pos(p + 1), ins)
                    }
                    _ => format!("{}delins{}", range(s, e), ins),
                };
                Case {
                    core: core.clone(),
                    sys,
                    hgvs: format!("{acc}:{prefix}{body}"),
                }
            },
        )
    })
}

fn provider_for(core: &str, sys: Sys) -> MockProvider {
    let len = core.len() as u64;
    match sys {
        Sys::Genomic => SyntheticBuilder::genomic(core).build(),
        Sys::CdsPlus => SyntheticBuilder::cds(core, 1, len, Strand::Plus).build(),
        Sys::CdsMinus => SyntheticBuilder::cds(core, 1, len, Strand::Minus).build(),
    }
}

/// `norm(norm(x)) == norm(x)` for one case in one direction, and the normalized
/// output must always re-parse — a normalizer that emits invalid HGVS (e.g. a
/// degenerate single-position `c.1insAA`) fails the re-parse step.
///
/// The render and first normalize fail LOUDLY, not by silent skip: cases are
/// valid by construction and the provider supplies the very core the coordinates
/// were derived from, so a parse/normalize error is a generator or normalizer
/// bug — and a skip branch here would let the whole property go vacuous unnoticed
/// (measured: 0 parse errors, 0 first-pass normalize errors over 4000 cases).
fn check_idempotent(case: &Case, dir: ShuffleDirection) -> Result<(), TestCaseError> {
    let norm = normalizer_for(&case.core, case.sys, dir);
    let v = parse_hgvs(&case.hgvs)
        .unwrap_or_else(|e| panic!("generator produced unparseable HGVS {}: {e}", case.hgvs));
    let n1 = norm
        .normalize(&v)
        .unwrap_or_else(|e| panic!("first normalize failed: {}: {e}", case.hgvs));
    let f1 = n1.to_string();
    // The normalized output MUST re-parse and re-normalize cleanly; an error
    // here (parse or normalize) is a bug.
    let reparsed = parse_hgvs(&f1).unwrap_or_else(|e| {
        panic!(
            "normalized output does not re-parse: {} -> {f1}: {e}\n  (dir={dir:?}, core={})",
            case.hgvs, case.core
        )
    });
    let n2 = norm
        .normalize(&reparsed)
        .unwrap_or_else(|e| panic!("re-normalize failed: {} -> {f1}: {e}", case.hgvs));
    prop_assert_eq!(
        &f1,
        &n2.to_string(),
        "NOT IDEMPOTENT\n  core={} sys={:?} dir={:?}\n  input={}\n  once ={}\n  twice={}",
        case.core,
        case.sys,
        dir,
        case.hgvs,
        f1,
        n2
    );
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(4000))]

    /// 3' (default, HGVS-canonical) direction, fully general edit content.
    /// (5' general idempotency is deferred — see the module-level scope note.)
    #[test]
    fn normalization_is_idempotent_three_prime(case in case_strategy()) {
        check_idempotent(&case, ShuffleDirection::ThreePrime)?;
    }
}

/// A homopolymer confluence case: a run of `base` of length `run_len`, embedded
/// in the CDS, with two distinct HGVS spellings of the SAME haplotype (inserting
/// `unit_len` more copies of `base` into the run at two different positions).
/// Both spellings describe an identical resulting sequence, so a confluent
/// normalizer must map them to the same normal form.
#[derive(Debug, Clone)]
struct ConfluenceCase {
    core: String,
    sys: Sys,
    a: String,
    b: String,
}

/// Build a pure-homopolymer or boundary-adjacent run and two equivalent
/// insertion spellings within it. Using `cds_start == 1`/`cds_end == len`
/// exercises the CDS-start and CDS-end boundary clamps (the confluence-critical
/// region), where all spellings must collapse to one `c.<boundary>delins…`.
fn confluence_case_strategy() -> impl Strategy<Value = ConfluenceCase> {
    (
        prop::sample::select(vec!['A', 'C', 'G', 'T']),
        4usize..=12, // run length
        1usize..=3,  // number of extra copies inserted
        prop::sample::select(vec![Sys::CdsPlus, Sys::CdsMinus]),
    )
        .prop_flat_map(|(base, run_len, extra, sys)| {
            // The whole CDS is the homopolymer run, so c.k maps to core k and
            // both transcript boundaries participate.
            let core: String = std::iter::repeat_n(base, run_len).collect();
            let unit: String = std::iter::repeat_n(base, extra).collect();
            // Insertion sits between positions p and p+1 (1-based), p in 1..run_len.
            // Any two distinct p give the same haplotype (run grows by `extra`).
            let core_c = core.clone();
            let unit_c = unit.clone();
            // `p2` is drawn as a non-zero offset from `p1` modulo the position
            // count, so the two spellings are always DISTINCT. Sampling both
            // independently would let `p1 == p2` through, and comparing a
            // normalized form against itself is trivially equal — a silent
            // weakening of the confluence property on those cases.
            let slots = run_len - 1; // positions are 1..run_len
            (1usize..run_len, 1usize..slots).prop_map(move |(p1, offset)| {
                let p2 = 1 + ((p1 - 1 + offset) % slots);
                let acc = "NM_TEST.1";
                let spell = |p: usize| format!("{acc}:c.{}_{}ins{}", p, p + 1, unit_c);
                ConfluenceCase {
                    core: core_c.clone(),
                    sys,
                    a: spell(p1),
                    b: spell(p2),
                }
            })
        })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(4000))]

    /// Two spellings of one haplotype (insert the same unit at two positions in
    /// a homopolymer run) must normalize to the SAME form, in both directions.
    /// This catches non-confluence that idempotency alone cannot: before the
    /// boundary-clamp rework the 5' direction mapped `c.1_2insAA`,
    /// `c.5_6insAA`, … to distinct (each individually idempotent)
    /// `c.1delinsAAA` / `c.1_3delinsAAAAA` / … forms.
    #[test]
    fn normalization_is_confluent(case in confluence_case_strategy(), dir in both_directions()) {
        let norm = normalizer_for(&case.core, case.sys, dir);
        // No skip branches: the generator only emits well-formed spellings, so a
        // parse or normalize failure is a BUG, not a case to pass over. Returning
        // `Ok(())` here would let the whole property go vacuous unnoticed — the
        // exact hazard `check_idempotent` documents and avoids.
        let va = parse_hgvs(&case.a)
            .unwrap_or_else(|e| panic!("generator produced unparseable HGVS {}: {e}", case.a));
        let vb = parse_hgvs(&case.b)
            .unwrap_or_else(|e| panic!("generator produced unparseable HGVS {}: {e}", case.b));
        let na = norm
            .normalize(&va)
            .unwrap_or_else(|e| panic!("normalize failed: {}: {e}", case.a));
        let nb = norm
            .normalize(&vb)
            .unwrap_or_else(|e| panic!("normalize failed: {}: {e}", case.b));
        prop_assert_eq!(
            na.to_string(),
            nb.to_string(),
            "NOT CONFLUENT\n  core={} sys={:?} dir={:?}\n  a={} b={}",
            case.core,
            case.sys,
            dir,
            case.a,
            case.b
        );
    }
}
