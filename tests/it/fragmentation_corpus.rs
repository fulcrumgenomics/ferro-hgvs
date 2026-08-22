//! Fragmentation campaign corpus — every reprex the #2155 / #2174 / #2175 family
//! asked for, in one place, measured on BOTH derivation surfaces.
//!
//! Each case pins ferro's current output on the two surfaces that can produce a
//! description of a `(reference, resulting)` pair:
//!
//! - **`from_seq`** — [`from_sequences`], which derives a description from the
//!   two sequences directly. This is where the #2174/#2175 fixes landed.
//! - **`normalized`** — [`normalize`] fed the recommended form (`wanted`). A
//!   correct normalizer holds the recommended form (it denotes the same
//!   sequence, so idempotency requires it), so `normalized != wanted` is a bug
//!   on the normalize surface even when `from_seq` already reaches it.
//!
//! `wanted` is the form the issue asked for; it is the shared target both
//! surfaces are measured against. A case where a surface's output equals
//! `wanted` is CLOSED on that surface; a case where they differ is OPEN WORK,
//! deliberate and visible rather than hidden in a skipped test.
//!
//! Two censuses count the closed cases per surface, so nothing closes or
//! regresses in silence — the same discipline as `reported_partition_verdicts`.
//! Both surfaces now reach every case: the `from_sequences` peel/coalesce guards
//! and the `collapse_overlapping_cis_edits` reference-tandem dup guard (#2175)
//! keep a dup-that-extends-a-reference-tandem out of a spanning `delins` on each
//! path. The censuses stay per-surface so a regression on either is caught alone.

use crate::common::cis_apply_oracle::{apply, normalize, normalize_in, provider};
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{from_sequences, parse_hgvs, FromSequencesOptions, HgvsVariant, ShuffleDirection};
use proptest::prelude::*;

struct Case {
    issue: &'static str,
    label: &'static str,
    reference: &'static str,
    resulting: &'static str,
    /// The form the issue asked for — the shared target for both surfaces.
    wanted: &'static str,
    /// What [`from_sequences`] prints for `(reference, resulting)` today.
    from_seq: &'static str,
    /// What [`normalize`] prints for `wanted` today. Equal to `wanted` for every
    /// case now that the `collapse_overlapping_cis_edits` reference-tandem dup
    /// guard (#2175) stops the normalize surface collapsing a dup-beside-a-change
    /// back to a `delins`; a future regression that reopened that gap would show
    /// here as `normalized != wanted`.
    normalized: &'static str,
}

const CORPUS: &[Case] = &[
    Case {
        issue: "2155",
        label: "contiguous CTTAGTTA->AAACAAAC",
        reference: "AGAACCCCCCTTAGTTAAGAACAAAAGCAACAATCTTCGTGGTCCTGG",
        resulting: "AGAACCCCCAAACAAACAGAACAAAAGCAACAATCTTCGTGGTCCTGG",
        wanted: "TEMPLATE:g.10_17delinsAAACAAAC",
        from_seq: "TEMPLATE:g.10_17delinsAAACAAAC",
        normalized: "TEMPLATE:g.10_17delinsAAACAAAC",
    },
    Case {
        issue: "2174",
        label: "11_14 GACT→ACTG",
        reference: "ACGTTCAGGTGACTTTAGCTAGCTAG",
        resulting: "ACGTTCAGGTACTGTTAGCTAGCTAG",
        wanted: "TEMPLATE:g.11_14delinsACTG",
        from_seq: "TEMPLATE:g.11_14delinsACTG",
        normalized: "TEMPLATE:g.11_14delinsACTG",
    },
    Case {
        issue: "2174",
        label: "11_13 ATA→TAC",
        reference: "ACGTTCAGGTATATTAGCTAGCTAG",
        resulting: "ACGTTCAGGTTACTTAGCTAGCTAG",
        wanted: "TEMPLATE:g.11_13delinsTAC",
        from_seq: "TEMPLATE:g.11_13delinsTAC",
        normalized: "TEMPLATE:g.11_13delinsTAC",
    },
    Case {
        issue: "2174",
        label: "11_13 CGA→GAG",
        reference: "ACGTTCAGGTCGATTAGCTAGCTAG",
        resulting: "ACGTTCAGGTGAGTTAGCTAGCTAG",
        wanted: "TEMPLATE:g.11_13delinsGAG",
        from_seq: "TEMPLATE:g.11_13delinsGAG",
        normalized: "TEMPLATE:g.11_13delinsGAG",
    },
    Case {
        issue: "2174",
        label: "11_13 GAC→ACT",
        reference: "ACGTTCAGGTGACTTAGCTAGCTAG",
        resulting: "ACGTTCAGGTACTTTAGCTAGCTAG",
        wanted: "TEMPLATE:g.11_13delinsACT",
        from_seq: "TEMPLATE:g.11_13delinsACT",
        normalized: "TEMPLATE:g.11_13delinsACT",
    },
    Case {
        issue: "2174",
        label: "11_13 CGC→GCT",
        reference: "ACGTTCAGGTCGCTTAGCTAGCTAG",
        resulting: "ACGTTCAGGTGCTTTAGCTAGCTAG",
        wanted: "TEMPLATE:g.11_13delinsGCT",
        from_seq: "TEMPLATE:g.11_13delinsGCT",
        normalized: "TEMPLATE:g.11_13delinsGCT",
    },
    Case {
        issue: "2174",
        label: "11_14 GTCA→TCAT",
        reference: "ACGTTCAGGTGTCATTAGCTAGCTAG",
        resulting: "ACGTTCAGGTTCATTTAGCTAGCTAG",
        wanted: "TEMPLATE:g.11_14delinsTCAT",
        from_seq: "TEMPLATE:g.11_14delinsTCAT",
        normalized: "TEMPLATE:g.11_14delinsTCAT",
    },
    Case {
        issue: "2175",
        label: "CA[2]->CA[3] +sub (dup surfaced)",
        reference: "ACGTTCAGGTCACAATTAGCTAGCTAG",
        resulting: "ACGTTCAGGTCACACACTTAGCTAGCTAG",
        wanted: "TEMPLATE:g.[13_14dup;15A>C]",
        from_seq: "TEMPLATE:g.[13_14dup;15A>C]",
        normalized: "TEMPLATE:g.[13_14dup;15A>C]",
    },
    Case {
        issue: "2175",
        label: "AG[2]->AG[3] +sub",
        reference: "ACGTTCAGGTAGAGCTTAGCTAGCTAG",
        resulting: "ACGTTCAGGTAGAGAGATTAGCTAGCTAG",
        wanted: "TEMPLATE:g.[13_14dup;15C>A]",
        from_seq: "TEMPLATE:g.[13_14dup;15C>A]",
        normalized: "TEMPLATE:g.[13_14dup;15C>A]",
    },
    Case {
        issue: "2175",
        label: "GT[3]->GT[4] +sub",
        reference: "ACGTTCAGGTGTGTGTCTTAGCTAGCTAG",
        resulting: "ACGTTCAGGTGTGTGTGTATTAGCTAGCTAG",
        wanted: "TEMPLATE:g.[15_16dup;17C>A]",
        from_seq: "TEMPLATE:g.[15_16dup;17C>A]",
        normalized: "TEMPLATE:g.[15_16dup;17C>A]",
    },
    Case {
        issue: "2175",
        label: "isolated CACA->CACACA",
        reference: "TAGTAAACCATTTTACGGAGGATCACAAATTCCTCCTTAT",
        resulting: "TAGTAAACCATTTTACGGAGGATCACACAAATTCCTCCTTAT",
        wanted: "TEMPLATE:g.26_27dup",
        from_seq: "TEMPLATE:g.26_27dup",
        normalized: "TEMPLATE:g.26_27dup",
    },
    Case {
        issue: "2175",
        label: "CACA->CACACA +28A>C",
        reference: "TAGTAAACCATTTTACGGAGGATCACAAATTCCTCCTTAT",
        resulting: "TAGTAAACCATTTTACGGAGGATCACACACATTCCTCCTTAT",
        wanted: "TEMPLATE:g.[26_27dup;28A>C]",
        from_seq: "TEMPLATE:g.[26_27dup;28A>C]",
        normalized: "TEMPLATE:g.[26_27dup;28A>C]",
    },
];

/// How many cases the `from_sequences` surface already prints as the wanted form.
const FROM_SEQ_REACHED_CENSUS: usize = 12;

/// How many cases the `normalize` surface already prints as the wanted form.
///
/// Now equal to [`FROM_SEQ_REACHED_CENSUS`]: the `collapse_overlapping_cis_edits`
/// reference-tandem dup guard (#2175) stops the normalize surface folding a
/// dup-that-extends-a-reference-tandem back into a `delins`, so both surfaces
/// reach every case.
const NORMALIZE_REACHED_CENSUS: usize = 12;

/// Derive against the `TEMPLATE` accession, so the `"TEMPLATE"`/anchor/`ERR:`
/// convention lives in one place and the two surfaces cannot drift apart.
fn from_seq_of(reference: &str, resulting: &str) -> String {
    from_sequences(
        "TEMPLATE",
        1,
        reference,
        resulting,
        &FromSequencesOptions::default(),
    )
    .map(|v| v.to_string())
    .unwrap_or_else(|e| format!("ERR:{e}"))
}

fn from_seq(c: &Case) -> String {
    from_seq_of(c.reference, c.resulting)
}

/// Every case prints exactly its pinned `from_seq` form — the regression floor
/// for the surface the #2174/#2175 fixes landed on.
#[test]
fn every_case_derives_its_pinned_from_sequences_form() {
    for c in CORPUS {
        assert_eq!(
            &from_seq(c),
            c.from_seq,
            "{} (#{}) moved off its from_sequences pin\n  ref:     {}\n  result:  {}\n  wanted:  {}",
            c.label,
            c.issue,
            c.reference,
            c.resulting,
            c.wanted
        );
    }
}

/// Every case's `normalize` output equals its pinned `normalized` form.
///
/// `wanted` is fed as the input: a correct normalizer holds the recommended
/// form, so a case whose `normalized` differs from `wanted` is the open
/// normalize-surface gap, pinned so it cannot widen or narrow silently.
#[test]
fn every_case_normalizes_to_its_pinned_form() {
    for c in CORPUS {
        assert_eq!(
            normalize(c.reference, c.wanted),
            c.normalized,
            "{} (#{}) moved off its normalize pin\n  ref:     {}\n  wanted:  {}",
            c.label,
            c.issue,
            c.reference,
            c.wanted
        );
    }
}

/// The per-surface reached-count census: closing an open case on either surface
/// must bump the matching constant in the same commit.
#[test]
fn the_reached_censuses_hold() {
    let from_seq_reached = CORPUS.iter().filter(|c| c.from_seq == c.wanted).count();
    assert_eq!(
        from_seq_reached,
        FROM_SEQ_REACHED_CENSUS,
        "{}/{} cases reach `wanted` via from_sequences; census pins {}.",
        from_seq_reached,
        CORPUS.len(),
        FROM_SEQ_REACHED_CENSUS
    );

    let normalize_reached = CORPUS.iter().filter(|c| c.normalized == c.wanted).count();
    assert_eq!(
        normalize_reached,
        NORMALIZE_REACHED_CENSUS,
        "{}/{} cases reach `wanted` via normalize; census pins {}. An open case \
         closing (or a closed one regressing) must move this number in the same commit.",
        normalize_reached,
        CORPUS.len(),
        NORMALIZE_REACHED_CENSUS
    );
}

/// Cross-surface agreement on the shapes the #2175 dup guard touches.
///
/// `normalize(input)` and `from_sequences(reference, apply(input))` describe one
/// variant two ways; a confluent normalizer settles both on one form. These pin
/// that the two derivation surfaces stay consistent — the scoped
/// `collapse_overlapping_cis_edits` dup guard keeps a `[dup;sub]` on both, and a
/// `[dup;del]` merges to one delins on BOTH rather than fragmenting on either.
/// Guards the scope against re-widening: the `[dup;del]` row would go red the
/// moment the guard started declining that shape again (the measured regression).
///
/// Only shapes whose two surfaces already agree are pinned here. A trailing
/// repeat member still renders `dup` on one surface and `unit[N]` on the other
/// (`open-issues.md:88`, repeat-vs-dup), a pre-existing divergence unrelated to
/// this guard and deliberately out of scope.
const CROSS_SURFACE: &[(&str, &str, &str)] = &[
    // (reference, input, agreed canonical form)
    // dup beside a substitution — the #2175 shape; the dup is kept on both.
    (
        "ACGTTCAGGTAGAGCTTAGCTAGCTAG",
        "TEMPLATE:g.[13_14dup;15C>A]",
        "TEMPLATE:g.[13_14dup;15C>A]",
    ),
    // dup beside a deletion — merges to one delins on both (NOT fragmented). This
    // is the row that catches the broad-guard regression the scope excludes.
    (
        "ACGTTCAGGTCACACAGTTAGCTAGCTAG",
        "TEMPLATE:g.[11_12dup;17_18del]",
        "TEMPLATE:g.17_18delinsCA",
    ),
    // dup beside an insertion — the dup is kept on both, but NOT by the #2175
    // guard: a group of two insertion-like edits fails
    // `collapse_overlapping_cis_edits`' `has_repl` requirement and is returned
    // untouched long before the guard is consulted. Pinned here because #2175's
    // scope paragraph names this shape as deliberately undecided
    // (`delins-adjacent-members-when-both-consume-reference`, 2026-08-14), so
    // the row is a tripwire on that boundary rather than on the guard's arity.
    (
        "ACGTTCAGGTAGAGCTTAGCTAGCTAG",
        "TEMPLATE:g.[13_14dup;16_17insGG]",
        "TEMPLATE:g.[13_14dup;16_17insGG]",
    ),
];

/// The two derivation surfaces agree on every cross-surface shape, and on the
/// pinned form.
#[test]
fn the_two_surfaces_agree_on_the_shapes_the_guard_touches() {
    for (reference, input, expected) in CROSS_SURFACE {
        let via_normalize = normalize(reference, input);
        let resulting = apply(reference, input).expect("input applies to reference");
        let via_from_seq = from_seq_of(reference, &resulting);
        assert_eq!(
            via_normalize, *expected,
            "normalize moved for {input}\n  got:      {via_normalize}\n  expected: {expected}"
        );
        assert_eq!(
            via_from_seq, *expected,
            "from_sequences moved for {input}\n  got:      {via_from_seq}\n  expected: {expected}"
        );
        assert_eq!(
            via_normalize, via_from_seq,
            "the two surfaces disagree for {input}\n  normalize:      {via_normalize}\n  from_sequences: {via_from_seq}"
        );
    }
}

/// Spelling-convergence (confluence) guard: several spellings of ONE variant all
/// normalize to one pinned form. The #2175 dup guard MOVES the attractor for the
/// dup-beside-a-sub shape from the delins to the dup; this pins that the move did
/// not SPLIT the attractor. Before the guard the delins was the single attractor
/// (the dup spelling collapsed to it); after, the dup is — and the delins
/// spelling must still converge onto it, or normalization has two fixed points
/// for one variant and is no longer confluent. The `[dup;del]` row is the
/// control: its two spellings converge on the merged delins (the guard's scope
/// leaves that shape merged), so both directions of the scope are pinned.
const CONVERGENCE: &[(&str, &str, &[&str])] = &[
    // (reference, the one form every spelling must reach, spellings)
    (
        "ACGTTCAGGTAGAGCTTAGCTAGCTAG",
        "TEMPLATE:g.[13_14dup;15C>A]",
        &["TEMPLATE:g.[13_14dup;15C>A]", "TEMPLATE:g.15delinsAGA"],
    ),
    (
        "TAGTAAACCATTTTACGGAGGATCACAAATTCCTCCTTAT",
        "TEMPLATE:g.[26_27dup;28A>C]",
        &["TEMPLATE:g.[26_27dup;28A>C]", "TEMPLATE:g.28delinsCAC"],
    ),
    // control: a dup beside a DELETION stays merged, so both spellings converge
    // on the delins rather than on a dup (the scope's other direction).
    (
        "ACGTTCAGGTCACACAGTTAGCTAGCTAG",
        "TEMPLATE:g.17_18delinsCA",
        &["TEMPLATE:g.[11_12dup;17_18del]", "TEMPLATE:g.17_18delinsCA"],
    ),
];

// ============================================================================
// #2192 — the anti-narrowness cross-product corpus (generated, no fixtures)
//
// The #2174 corpus above pins single-member runs one shape at a time; that is
// exactly the blindness #2192 lived in. This section GENERATES a cross-product
// of {run shape} x {run count} x {member count} x {both shuffle directions} on
// BOTH surfaces, and asserts each cell reaches the same recommended form its
// runs reach in ISOLATION. A fix scoped to single-run or single-member fails
// loudly on the multi cells (on `origin/main` every multi-run cell fragments);
// this is the signal the #2174 corpus lacked.
// ============================================================================

/// A run template: a local `ref_seg -> alt_seg` edit, plus the member count its
/// recommended (non-fragmented) single-run form has. The counts are the #2174 /
/// #2175 results (a balanced equal-length run is one spanning `delins`; a
/// reference-tandem expansion abutting a substitution is a `dup` plus that
/// `sub`) and are asserted independently below, so a run that fragments in
/// isolation fails here rather than passing this corpus vacuously.
struct RunTemplate {
    name: &'static str,
    ref_seg: &'static str,
    alt_seg: &'static str,
    members: usize,
}

/// Net-length-preserving run templates: a balanced equal-length rotation is one
/// spanning `delins` (#2174), a point substitution is one `sub`. Each is one
/// member, and each returns the running length delta to zero, so the per-run
/// grouping isolates any number of them in any order (`coalesce_by_run`'s
/// zero-shift boundary). These fill the cross-product's run slots.
const BASE_TEMPLATES: &[RunTemplate] = &[
    RunTemplate {
        name: "balanced4",
        ref_seg: "GACT",
        alt_seg: "ACTG",
        members: 1,
    },
    RunTemplate {
        name: "balanced3",
        ref_seg: "ATA",
        alt_seg: "TAC",
        members: 1,
    },
    RunTemplate {
        name: "sub",
        ref_seg: "C",
        alt_seg: "A",
        members: 1,
    },
];

/// A reference tandem (`AGAG`) expanded by one unit beside a substitution ->
/// `[dup; sub]`, TWO members (#2175). This is the member-count (arity) axis: one
/// run that is two members. A `dup` is net-imbalanced and content-defined, so the
/// zero-shift grouping isolates it only when nothing follows it (see
/// `coalesce_by_run`'s deferral note); the cross-product therefore places it only
/// as the FINAL run. `dup`-before-another-run is deferred to a segment-first
/// partition and is asserted UNCHANGED (still fragmenting) by
/// `a_dup_before_another_run_is_left_for_category_one` below, so the boundary is
/// pinned rather than silent.
const TRAILING_DUP: RunTemplate = RunTemplate {
    name: "dupsub",
    ref_seg: "AGAGC",
    alt_seg: "AGAGAGA",
    members: 2,
};

/// Flanks and separators. The pads give the shuffle room to move at both ends;
/// the spacer is long enough (>= 2 unchanged bases, `general.md:33`) that the
/// runs it separates never interact, so a cell's runs are genuinely independent.
const PAD: &str = "CACGTACGTC";
const SPACER: &str = "TCGGATCAG";

/// One generated cell: the template filling each run slot, left to right.
struct Cell {
    templates: Vec<&'static RunTemplate>,
}

impl Cell {
    /// The reference, and the 0-based offset at which each run's `ref_seg` sits.
    fn reference(&self) -> (String, Vec<usize>) {
        let mut reference = String::from(PAD);
        let mut offsets = Vec::new();
        for (i, template) in self.templates.iter().enumerate() {
            if i > 0 {
                reference.push_str(SPACER);
            }
            offsets.push(reference.len());
            reference.push_str(template.ref_seg);
        }
        reference.push_str(PAD);
        (reference, offsets)
    }

    /// The reference with `which` run slots applied (all of them for the multi
    /// case, exactly one for an isolated per-run derivation). Built left to
    /// right so the offsets stay valid whatever the length changes.
    fn applied(&self, reference: &str, offsets: &[usize], which: impl Fn(usize) -> bool) -> String {
        let bytes = reference.as_bytes();
        let mut out = String::with_capacity(bytes.len());
        let mut cursor = 0usize;
        for (i, template) in self.templates.iter().enumerate() {
            let start = offsets[i];
            out.push_str(&reference[cursor..start]);
            if which(i) {
                out.push_str(template.alt_seg);
            } else {
                out.push_str(template.ref_seg);
            }
            cursor = start + template.ref_seg.len();
        }
        out.push_str(&reference[cursor..]);
        out
    }
}

/// `from_sequences` for the `TEMPLATE` accession in one shuffle direction.
fn from_seq_dir(reference: &str, resulting: &str, direction: ShuffleDirection) -> HgvsVariant {
    from_sequences(
        "TEMPLATE",
        1,
        reference,
        resulting,
        &FromSequencesOptions::default().with_direction(direction),
    )
    .unwrap_or_else(|e| panic!("from_sequences({reference}, {resulting}): {e}"))
}

/// The members of a rendered description, each as `(interbase_start, body)`
/// where `body` is the spelling without the `TEMPLATE:g.` prefix, sorted into
/// the ascending interbase order both surfaces render in (#1098). Comparing
/// these makes a single-member `TEMPLATE:g.X` and a one-member allele equal, and
/// is independent of how the members happened to be grouped.
fn member_bodies(reference: &str, desc: &str) -> Vec<String> {
    let template = provider(reference);
    let parsed = parse_hgvs(desc).unwrap_or_else(|e| panic!("parse {desc}: {e:?}"));
    let members: Vec<HgvsVariant> = match parsed {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    let mut keyed: Vec<(u64, String)> = members
        .iter()
        .map(|member| {
            let start = hgvs_to_spdi(member, &template)
                .unwrap_or_else(|e| panic!("{desc}: member {member} has no SPDI: {e}"))
                .position;
            let rendered = member.to_string();
            let body = rendered
                .strip_prefix("TEMPLATE:g.")
                .unwrap_or(&rendered)
                .to_string();
            (start, body)
        })
        .collect();
    keyed.sort_by_key(|(start, _)| *start);
    keyed.into_iter().map(|(_, body)| body).collect()
}

/// Every tuple of `len` base templates, in base-`n` counting order (so every
/// shape mix, with repetition, appears).
fn base_tuples(len: u32) -> Vec<Vec<&'static RunTemplate>> {
    let n = BASE_TEMPLATES.len();
    let count = n.pow(len);
    (0..count)
        .map(|mut code| {
            (0..len)
                .map(|_| {
                    let idx = code % n;
                    code /= n;
                    &BASE_TEMPLATES[idx]
                })
                .collect()
        })
        .collect()
}

/// The largest run count the exhaustive cross-product enumerates. Five net-zero
/// runs (and a sixth when a trailing `dup` is appended) exercises the four- and
/// five-member alleles #2192 is about at full arity, in every ordering — the whole
/// point being that the per-run coalesce must not depend on how many members
/// precede or follow a run.
const MAX_RUNS: u32 = 5;

/// Enumerate every cell of the cross-product:
///
/// - all base-template tuples of length `1..=MAX_RUNS` (run count and member count
///   both `1..=5`, every shape mix and ordering), and
/// - all base-template tuples of length `0..MAX_RUNS` with a `TRAILING_DUP`
///   appended (the member-count axis: the two-member `dup` run raises member count
///   as high as six while run count stays `1..=5`, `dup` always last so the
///   zero-shift grouping isolates it).
fn cells() -> Vec<Cell> {
    let mut cells = Vec::new();
    for len in 1..=MAX_RUNS {
        for templates in base_tuples(len) {
            cells.push(Cell { templates });
        }
    }
    for len in 0..MAX_RUNS {
        for mut templates in base_tuples(len) {
            templates.push(&TRAILING_DUP);
            cells.push(Cell { templates });
        }
    }
    cells
}

/// Sort member bodies into the ascending interbase order both surfaces render
/// in, so a per-run combination can be compared to a whole-allele output.
fn sort_bodies(reference: &str, mut bodies: Vec<String>) -> Vec<String> {
    // Build the provider once, and cache each member's key so it is parsed and
    // resolved once per body rather than once per comparison (`sort_by_key` may
    // call its closure `O(n log n)` times).
    let provider = provider(reference);
    bodies.sort_by_cached_key(|body| {
        hgvs_to_spdi(
            &parse_hgvs(&format!("TEMPLATE:g.{body}")).expect("parse member"),
            &provider,
        )
        .expect("member spdi")
        .position
    });
    bodies
}

/// Wrap member bodies back into a `TEMPLATE:g.` description (an allele when there
/// is more than one member).
fn wrap_members(bodies: &[String]) -> String {
    if bodies.len() == 1 {
        format!("TEMPLATE:g.{}", bodies[0])
    } else {
        format!("TEMPLATE:g.[{}]", bodies.join(";"))
    }
}

/// The cross-product's core assertion, run per shuffle direction and judged
/// PER SURFACE: each surface's multi-run output must equal that same surface's
/// per-run ISOLATED outputs combined. That is #2192's property — every run
/// derived independently of its neighbours — stated without coupling the two
/// surfaces to each other, since `from_sequences` and `normalize` may
/// legitimately place a `dup` at different shuffle anchors (an orthogonal
/// question, visible only on the 5' surface). A fix scoped to single-run or
/// single-member fails here: the multi output fragments while the isolated
/// per-run outputs do not, so the two stop being equal.
fn assert_cross_product(direction: ShuffleDirection) {
    for cell in cells() {
        assert_cell(&cell, direction);
    }
}

/// One cell's per-surface, per-run assertion — factored out of
/// [`assert_cross_product`] so the exhaustive enumeration and the proptest fuzzer
/// (`fuzz_the_cross_product`) check the identical property on their respective
/// corpora.
fn assert_cell(cell: &Cell, direction: ShuffleDirection) {
    let (reference, offsets) = cell.reference();
    let label: Vec<&str> = cell.templates.iter().map(|t| t.name).collect();
    let context = format!("cell {label:?} {direction:?}");
    let resulting = cell.applied(&reference, &offsets, |_| true);

    // Per-run isolated forms, on each surface. `from_sequences` derives from
    // the single-run sequence; `normalize` re-derives that run's own
    // recommended form. The per-template member-count floor is an independent
    // anti-fragmentation check: a run that fragments in isolation fails here
    // rather than letting the equality below pass vacuously.
    let mut expected_from_seq: Vec<String> = Vec::new();
    let mut expected_normalize: Vec<String> = Vec::new();
    for (i, template) in cell.templates.iter().enumerate() {
        let single = cell.applied(&reference, &offsets, |j| j == i);
        let fs = member_bodies(
            &reference,
            &from_seq_dir(&reference, &single, direction).to_string(),
        );
        assert_eq!(
            fs.len(),
            template.members,
            "{context}: run {i} ({}) fragmented in ISOLATION into {fs:?} — the per-template \
                 floor is wrong or #2174/#2175 regressed",
            template.name
        );
        let run = member_bodies(
            &reference,
            &normalize_in(&reference, &wrap_members(&fs), direction),
        );
        expected_from_seq.extend(fs);
        expected_normalize.extend(run);
    }
    let expected_from_seq = sort_bodies(&reference, expected_from_seq);
    let expected_normalize = sort_bodies(&reference, expected_normalize);

    // from_sequences surface: the multi-run derivation equals the per-run one.
    let multi = from_seq_dir(&reference, &resulting, direction).to_string();
    assert_eq!(
            member_bodies(&reference, &multi),
            expected_from_seq,
            "{context}: from_sequences fragmented the multi-run allele\n  got:      {multi}\n  expected: {expected_from_seq:?}"
        );

    // normalize surface: normalizing the whole recommended allele equals
    // normalizing each run's recommended form on its own.
    let recommended = wrap_members(&expected_from_seq);
    let normalized = normalize_in(&reference, &recommended, direction);
    assert_eq!(
            member_bodies(&reference, &normalized),
            expected_normalize,
            "{context}: normalize fragmented the multi-run allele\n  input:    {recommended}\n  expected: {expected_normalize:?}"
        );

    // meaning preservation on ALL THREE surfaces — the member-body equality above
    // only proves the spellings match, not that they denote the resulting sequence.
    // A per-run coalesce that drops or duplicates a base the SAME way in the
    // multi-run allele and in each isolated run passes the equality yet corrupts
    // the sequence, so normalize's own output must reach the apply oracle too.
    assert_eq!(
        apply(&reference, &recommended).as_deref(),
        Some(resulting.as_str()),
        "{context}: recommended form does not denote the resulting sequence\n  form: {recommended}"
    );
    assert_eq!(
        apply(&reference, &multi).as_deref(),
        Some(resulting.as_str()),
        "{context}: from_sequences output does not denote the resulting sequence\n  out: {multi}"
    );
    assert_eq!(
        apply(&reference, &normalized).as_deref(),
        Some(resulting.as_str()),
        "{context}: normalize output does not denote the resulting sequence\n  out: {normalized}"
    );
}

/// The cross-product reaches its recommended form on the 3' surface.
#[test]
fn the_cross_product_reaches_the_recommended_form_three_prime() {
    assert_cross_product(ShuffleDirection::ThreePrime);
}

/// The cross-product reaches its recommended form on the 5' surface.
#[test]
fn the_cross_product_reaches_the_recommended_form_five_prime() {
    assert_cross_product(ShuffleDirection::FivePrime);
}

/// The three reproductions from issue #2192, verbatim — accession `CONTIG1`,
/// sequences and coordinates exactly as the issue states them — pinned on BOTH
/// surfaces the issue names. Each is a contiguous run beside a non-contiguous
/// member that fragmented before this fix: Ex1 a `del+ins` rotation, Ex2 a
/// `del+dup` (the dup smeared 3' to `g.15`), Ex3 a run with its sibling on the 5'
/// side. This is the guard a reviewer looks for — the exact strings from the
/// report, not a paraphrase of their geometry.
#[test]
fn the_issue_2192_reprexes_reach_the_recommended_form() {
    struct Reprex {
        name: &'static str,
        reference: &'static str,
        resulting: &'static str,
        recommended: &'static str,
        // The fragmented 3'-shuffled output the issue records for `main`, fed back
        // through the normalize surface to prove it coalesces (not merely that the
        // recommended form is a fixed point).
        fragmented: &'static str,
    }
    let reprexes = [
        Reprex {
            name: "Ex1 del+ins",
            reference: "ACGTTCAGGTGACTTTAGCTAGCTAGCTAGCTAG",
            resulting: "ACGTTCAGGTACTGTTAGCTAGCTAGATAGCTAG",
            recommended: "g.[11_14delinsACTG;27C>A]",
            fragmented: "g.[11del;14_15insG;27C>A]",
        },
        Reprex {
            name: "Ex2 del+dup",
            reference: "ACGTTCAGGTGACTTAGCTAGCTAGCTAGCTAG",
            resulting: "ACGTTCAGGTACTTTAGCTAGCTAGCGAGCTAG",
            recommended: "g.[11_13delinsACT;27T>G]",
            fragmented: "g.[11del;15dup;27T>G]",
        },
        Reprex {
            name: "Ex3 5'-side change",
            reference: "ACGTTCAGGTATATTAGCTAGCTAGCTAGCTAG",
            resulting: "ACGTGCAGGTTACTTAGCTAGCTAGCTAGCTAG",
            recommended: "g.[5T>G;11_13delinsTAC]",
            fragmented: "g.[5T>G;11del;13_14insC]",
        },
    ];
    for r in reprexes {
        // Surface 1 — from_sequences, the issue's own reproducer, under CONTIG1 and
        // in BOTH shuffle directions. The recommended form is direction-independent
        // here (no run 3'-shifts across its sibling), so both must reach it.
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            let got = from_sequences(
                "CONTIG1",
                1,
                r.reference,
                r.resulting,
                &FromSequencesOptions::default().with_direction(direction),
            )
            .unwrap_or_else(|e| panic!("{}: from_sequences: {e}", r.name))
            .to_string();
            assert_eq!(
                got,
                format!("CONTIG1:{}", r.recommended),
                "{}: from_sequences {direction:?} fragmented",
                r.name
            );
        }
        // Surface 2 — normalize/rederive (the issue's `rederive(recommended_form=…)`
        // path). Driven through the TEMPLATE-keyed helper on the identical bytes, so
        // the coordinates match. The fragmented `main` output must coalesce, and the
        // recommended form must be a fixed point.
        let templated = |body: &str| format!("TEMPLATE:{body}");
        assert_eq!(
            normalize(r.reference, &templated(r.fragmented)),
            templated(r.recommended),
            "{}: normalize did not de-fragment the main output",
            r.name
        );
        assert_eq!(
            normalize(r.reference, &templated(r.recommended)),
            templated(r.recommended),
            "{}: recommended form is not a normalize fixed point",
            r.name
        );
        // Every spelling denotes the same bases as the resulting sequence.
        for body in [r.recommended, r.fragmented] {
            assert_eq!(
                apply(r.reference, &templated(body)).as_deref(),
                Some(r.resulting),
                "{}: {body} does not denote the resulting sequence",
                r.name
            );
        }
    }
}

proptest! {
    // Randomised extension of the exhaustive cross-product: a cell of 1..=6 net-zero
    // runs in any order, with an optional trailing `dup`, in either shuffle
    // direction. It samples the same shape `cells()` enumerates but out to higher
    // arity and beyond what an exhaustive `MAX_RUNS` can reach, so a fragmentation
    // that only appears at some member count or ordering the fixed set happens to
    // miss still fails. `assert_cell` is the identical property the two exhaustive
    // tests check. The strategy draws only from known-good templates, so every cell
    // is well-formed — a failure is a real fragmentation, never a malformed input.
    #![proptest_config(ProptestConfig::with_cases(512))]
    #[test]
    fn fuzz_the_cross_product(
        base_idx in proptest::collection::vec(0usize..BASE_TEMPLATES.len(), 1..=6),
        trailing_dup in any::<bool>(),
        five_prime in any::<bool>(),
    ) {
        let mut templates: Vec<&'static RunTemplate> =
            base_idx.iter().map(|&i| &BASE_TEMPLATES[i]).collect();
        if trailing_dup {
            templates.push(&TRAILING_DUP);
        }
        let direction = if five_prime {
            ShuffleDirection::FivePrime
        } else {
            ShuffleDirection::ThreePrime
        };
        assert_cell(&Cell { templates }, direction);
    }
}

/// The one shape the per-run grouping deliberately does NOT reach: a
/// net-imbalanced `dup` sitting 5' of another run. The `dup`'s length change is
/// banked but content-defined, so the zero-shift grouping cannot isolate it from
/// what follows (see `coalesce_by_run`'s deferral note), and both the `dup` and
/// the trailing run fragment. This is the Track-2 / segment-first-partition work
/// tracked by #2194. This pins that boundary — asserting the DEFECT, in the manner
/// of `issue_1610`'s `a_lone_unequal_length_delins_is_not_split` — so a future
/// segment-first partition that closes it flips this test loudly rather than
/// moving the row in silence. The output still denotes the right sequence.
#[test]
fn a_dup_before_another_run_is_left_for_category_one() {
    // A `dup` run (`AGAGC -> AGAGAGA`, positions 11-15) 5' of a balanced4 run
    // (`GACT -> ACTG`, positions 25-28).
    let reference = "CACGTACGTCAGAGCTCGGATCAGGACTCACGTACGTC";
    let resulting = "CACGTACGTCAGAGAGATCGGATCAGACTGCACGTACGTC";
    let recommended = "TEMPLATE:g.[13_14dup;15C>A;25_28delinsACTG]";
    // Meaning-preserving sanity: the recommended form denotes what the runs make.
    assert_eq!(
        apply(reference, recommended).as_deref(),
        Some(resulting),
        "recommended form must denote the resulting sequence"
    );
    for (direction, deferred) in [
        (
            ShuffleDirection::ThreePrime,
            "TEMPLATE:g.[15delinsAGA;25del;28_29insG]",
        ),
        (
            ShuffleDirection::FivePrime,
            "TEMPLATE:g.[15delinsAGA;24del;28_29insG]",
        ),
    ] {
        let got = from_seq_dir(reference, resulting, direction).to_string();
        assert_ne!(
            got, recommended,
            "{direction:?}: dup-before-run unexpectedly reached the recommended form — \
             if a segment-first partition now closes this, move the cell into the \
             cross-product above"
        );
        assert_eq!(
            got, deferred,
            "{direction:?}: deferred fragmentation form moved"
        );
        // Whatever it emits, it still denotes the resulting sequence.
        assert_eq!(
            apply(reference, &got).as_deref(),
            Some(resulting),
            "{direction:?}: the deferred output changed the sequence"
        );
    }
}

/// Every spelling of each variant converges on the one pinned form.
///
/// First establishes — through the independent `apply` oracle, not the
/// normalizer — that the spellings really denote ONE variant (they and the
/// `converged` form all apply to the same resulting sequence). Otherwise a typo
/// in `CONVERGENCE` that named two genuinely different variants which happened
/// to normalize alike would read as a convergence success.
#[test]
fn every_spelling_of_a_variant_converges() {
    for (reference, converged, spellings) in CONVERGENCE {
        let denoted =
            |desc: &str| apply(reference, desc).unwrap_or_else(|| panic!("{desc} applies"));
        let one_sequence = denoted(converged);
        for spelling in *spellings {
            assert_eq!(
                denoted(spelling),
                one_sequence,
                "{spelling} denotes a different variant than {converged} — not one equivalence class"
            );
            let out = normalize(reference, spelling);
            assert_eq!(
                &out, converged,
                "spelling did not converge\n  spelling:  {spelling}\n  got:       {out}\n  converged: {converged}"
            );
        }
    }
}
