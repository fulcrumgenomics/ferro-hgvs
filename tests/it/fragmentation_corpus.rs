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
//! Today the two disagree: `from_seq` reaches every case, while the normalize
//! surface still collapses four #2175 dup-beside-a-change cases to a `delins`.
//! That gap is the open work, and it is a census here rather than a comment.

use crate::common::cis_apply_oracle::normalize;
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
    /// What [`normalize`] prints for `wanted` today. Equal to `wanted` when the
    /// normalize surface holds the recommended form; unequal where it collapses
    /// a dup-beside-a-change back to a `delins` (the open #2175 gap).
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
        // OPEN: the normalize surface collapses the dup back to a delins.
        normalized: "TEMPLATE:g.15delinsCAC",
    },
    Case {
        issue: "2175",
        label: "AG[2]->AG[3] +sub",
        reference: "ACGTTCAGGTAGAGCTTAGCTAGCTAG",
        resulting: "ACGTTCAGGTAGAGAGATTAGCTAGCTAG",
        wanted: "TEMPLATE:g.[13_14dup;15C>A]",
        from_seq: "TEMPLATE:g.[13_14dup;15C>A]",
        // OPEN: the normalize surface collapses the dup back to a delins.
        normalized: "TEMPLATE:g.15delinsAGA",
    },
    Case {
        issue: "2175",
        label: "GT[3]->GT[4] +sub",
        reference: "ACGTTCAGGTGTGTGTCTTAGCTAGCTAG",
        resulting: "ACGTTCAGGTGTGTGTGTATTAGCTAGCTAG",
        wanted: "TEMPLATE:g.[15_16dup;17C>A]",
        from_seq: "TEMPLATE:g.[15_16dup;17C>A]",
        // OPEN: the normalize surface collapses the dup back to a delins.
        normalized: "TEMPLATE:g.17delinsGTA",
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
        // OPEN: the normalize surface collapses the dup back to a delins.
        normalized: "TEMPLATE:g.28delinsCAC",
    },
];

/// How many cases the `from_sequences` surface already prints as the wanted form.
const FROM_SEQ_REACHED_CENSUS: usize = 12;

/// How many cases the `normalize` surface already prints as the wanted form.
///
/// Four short of [`FROM_SEQ_REACHED_CENSUS`]: the normalize surface still
/// collapses the four #2175 dup-beside-a-change cases to a `delins`. Closing one
/// is one edit to its `normalized` field plus a bump here.
const NORMALIZE_REACHED_CENSUS: usize = 8;

fn from_seq(c: &Case) -> String {
    from_sequences(
        "TEMPLATE",
        1,
        c.reference,
        c.resulting,
        &FromSequencesOptions::default(),
    )
    .map(|v| v.to_string())
    .unwrap_or_else(|e| format!("ERR:{e}"))
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
