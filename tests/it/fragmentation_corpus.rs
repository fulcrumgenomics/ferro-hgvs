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

use crate::common::cis_apply_oracle::{apply, normalize};
use ferro_hgvs::{from_sequences, FromSequencesOptions};

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
