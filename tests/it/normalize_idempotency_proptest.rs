//! Property-based normalization idempotency test over reference-backed
//! synthetic fixtures (issue #1157 follow-up campaign).
//!
//! The pre-existing idempotency coverage is reference-free (empty provider, so
//! nothing shifts) or substitution-only (subs never shift). This fuzzes random
//! **repeat-biased** cores (homopolymers + short tandems, where del/dup/ins
//! actually 3'-shift) across the genomic and CDS (both strands) coordinate
//! systems, and every edit type, asserting `norm(norm(x)) == norm(x)`.
//!
//! Direction: this test covers the **3' (default, HGVS-canonical) direction**
//! only. The 5' direction is a separate, systemically under-baked mode with
//! multiple independent non-idempotency / non-confluence bugs in the shuffle +
//! boundary-clamp interaction; it is being reworked as its own dedicated effort
//! and is intentionally out of scope here (see the campaign handoff).
//!
//! Cases are generated **valid by construction** (edit endpoints derived from
//! the core length), so there is no skip path at all: an unparseable render or a
//! failing first-pass normalize is a bug in the generator or the normalizer and
//! fails the test rather than silently shrinking the exercised set. Failure
//! seeds are recorded by proptest under
//! `tests/it/proptest-regressions/normalize_idempotency_proptest.txt`.

use crate::common::synthetic::{SyntheticBuilder, PAD_OFFSET};
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};
use proptest::prelude::*;
use proptest::test_runner::Config as ProptestConfig;

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
    core_strategy().prop_flat_map(|core| {
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

proptest! {
    #![proptest_config(ProptestConfig::with_cases(4000))]

    /// `norm(norm(x)) == norm(x)` for every reference-backed synthetic case in
    /// the default 3' direction.
    #[test]
    fn normalization_is_idempotent_three_prime(case in case_strategy()) {
        let norm = Normalizer::new(provider_for(&case.core, case.sys));
        // Both of these are generator/normalizer bugs, not "cases not under
        // test", so they fail loudly rather than silently skipping — a skip
        // branch here would let the whole property go vacuous unnoticed. The
        // renders are valid by construction and the provider supplies the very
        // core the coordinates were derived from, so neither fires today
        // (measured: 0 parse errors, 0 first-pass normalize errors over 4000
        // cases). Should a future change reintroduce either, that is the signal.
        let v = parse_hgvs(&case.hgvs)
            .unwrap_or_else(|e| panic!("generator produced unparseable HGVS {}: {e}", case.hgvs));
        let n1 = norm
            .normalize(&v)
            .unwrap_or_else(|e| panic!("first normalize failed: {}: {e}", case.hgvs));
        let f1 = n1.to_string();
        // The normalized output MUST re-normalize cleanly; an error here is a bug.
        let n2 = norm
            .normalize(&n1)
            .unwrap_or_else(|e| panic!("re-normalize failed: {} -> {f1}: {e}", case.hgvs));
        prop_assert_eq!(
            &f1,
            &n2.to_string(),
            "NOT IDEMPOTENT\n  core={} sys={:?}\n  input={}\n  once ={}\n  twice={}",
            case.core,
            case.sys,
            case.hgvs,
            f1,
            n2
        );
    }
}
