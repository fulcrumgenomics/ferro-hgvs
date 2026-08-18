//! Coding-axis rederive prototype (issue #2155 follow-up).
//!
//! Exercises the experimental `FERRO_CODING_REDERIVE` path: a CDS-body `c.`
//! allele re-derived from its resulting sequence through the confluent
//! sequence-first core, so the `c.[..]` spellings of one variant converge on one
//! description. Like `proto_genomic_rederive`, these assert **convergence +
//! idempotence**, never *which* form the core picks — the form is what the
//! operator re-adjudicates.
//!
//! Convergence alone is a self-consistency check: two spellings could converge on
//! one *wrong* form and still pass. So each case also checks
//! **sequence-preservation against an independent oracle** —
//! `cis_apply_oracle::apply_with`, which reaches the bases through `hgvs_to_spdi`
//! and an SPDI splice rather than through the normalizer — so a merge onto the
//! wrong bases fails here even though both spellings would still "converge".
//!
//! The path is **off by default**; each test enables it via the env flag. Nextest
//! runs each test in its own process, so the flag (read once into a process-wide
//! `OnceLock`) is scoped to the test that sets it before its first normalize.

use crate::common::cis_apply_oracle::apply_with;
use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};
use ferro_hgvs::reference::transcript::Strand;

/// Enable the experimental coding rederive for this test process. Must run
/// before the first `normalize`, which is why every test calls it first.
fn enable_coding_rederive() {
    // `set_var` is safe on edition 2021; nextest isolates each test in its own
    // process, so this never races another test's reading of the flag.
    std::env::set_var("FERRO_CODING_REDERIVE", "1");
}

/// Assert that two `c.` spellings of one variant converge under coding rederive,
/// and that the convergence is **sequence-preserving** — checked against an oracle
/// independent of the normalizer, so a merge onto the wrong bases fails here.
///
/// `core` is the transcript sequence (`cds_start = 1`, so `c.N` is transcript
/// position `N`); `a`/`b` are the two spellings. Returns the converged form.
fn assert_converges(core: &str, cds_end: u64, a: &str, b: &str) -> String {
    let p = SyntheticBuilder::cds(core, 1, cds_end, Strand::Plus).build();

    let na = normalize_to_string(p.clone(), a);
    let nb = normalize_to_string(p.clone(), b);
    assert_eq!(
        na, nb,
        "the two c. spellings must converge under coding rederive:\n  {a}\n  {b}"
    );

    // Independent oracle: apply each input and the converged output to the
    // transcript sequence via SPDI (NOT the normalizer, and not the core's own
    // `verify_round_trip`, which shares the applier). Without this, a merge onto
    // the wrong sequence would still pass the equality above.
    let bases_a = apply_with(&p, core, a).expect("input A denotes a sequence");
    let bases_b = apply_with(&p, core, b).expect("input B denotes a sequence");
    let bases_out = apply_with(&p, core, &na).expect("converged output denotes a sequence");
    assert_eq!(
        bases_a, bases_b,
        "fixture sanity: the two inputs must denote the same bases"
    );
    assert_eq!(
        bases_out, bases_a,
        "coding rederive must preserve the denoted bases (independent SPDI oracle)"
    );

    // Idempotent: normalizing the converged form again is a no-op.
    assert_eq!(
        normalize_to_string(p, &na),
        na,
        "coding rederive must be idempotent"
    );
    na
}

#[test]
fn coding_rederive_converges_the_2155_pair() {
    enable_coding_rederive();
    // The #2155 window: c.10_17 = CTTAGTTA, replaced by AAACAAAC. Two spellings of
    // one variant — the minimal-alignment split and the spanning delins.
    let core = "ATGCCCGGGCTTAGTTAGGCCAATTCCGGAT";
    assert_eq!(&core[9..17], "CTTAGTTA", "fixture: c.10_17 == CTTAGTTA");
    assert_converges(
        core,
        30,
        "NM_TEST.1:c.[10_12delinsAA;14_16delinsCAA;17_18insC]",
        "NM_TEST.1:c.10_17delinsAAACAAAC",
    );
}

#[test]
fn coding_rederive_merges_the_one_codon_exception() {
    enable_coding_rederive();
    // codon 4 = c.10..c.12 = CTT; the one-amino-acid exception (`general.md:35`):
    // `[10C>T;12T>A]` and `c.10_12delinsTTA` both denote CTT -> TTA.
    let core = "ATGCCCGGGCTTAGTTAGGCCAATTCCGGAT";
    assert_eq!(&core[9..12], "CTT", "fixture: c.10_12 == CTT");
    assert_converges(
        core,
        30,
        "NM_TEST.1:c.[10C>T;12T>A]",
        "NM_TEST.1:c.10_12delinsTTA",
    );
}
