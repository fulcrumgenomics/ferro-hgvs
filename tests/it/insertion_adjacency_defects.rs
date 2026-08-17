//! Defects found while auditing the insertion-adjacency and copy-range
//! families. **Every test in this module is expected to fail** until the
//! behaviour it describes is fixed or the policy question it raises is decided.
//!
//! # Why they are live tests and not `#[ignore]`
//!
//! An ignored guard is how this tree has lost coverage before: the pinned
//! proptest seeds in `cis_allele_confluence_proptest` guarded a bug for months
//! while the test that owned them was ignored, so nothing went red when the
//! guard stopped guarding. A red test is noisy and honest; an ignored one is
//! quiet and worthless. If these need to stop failing CI before they are fixed,
//! the right move is to decide that deliberately, not to hide them here.
//!
//! # Three different kinds of red, do not conflate them
//!
//! - **DEFECT** — ferro produces a wrong or unusable answer. Not contingent on
//!   any ruling: `emits_a_description_its_own_parser_rejects` and the two
//!   past-cds-end tests are defects on any reading.
//! - **POLICY** — ferro's behaviour contradicts a *published spec example*, but
//!   the spec contradicts itself nearby, so "correct" needs adjudication before
//!   the expectation here can be trusted.
//! - **UNRANKED** — deliberately absent. Ferro never emits a copy-range payload
//!   for a large insert, and it is tempting to call that a gap; the spec ranks
//!   neither form (`open-issues.md:77`), so there is nothing to fail against and
//!   no test is written. Recorded here so its absence reads as a decision.
//!
//! # Relationship to the non-idempotence class already pinned
//!
//! Two of these live next door to defects this tree already records, and neither
//! should be read as an unrelated discovery:
//!
//! - `defect_seam_allele_emits_a_description_its_own_parser_rejects` is at the
//!   **same CDS/3'UTR seam** as `defect_non_idempotent_outputs`' "mechanism 1"
//!   (the four `cds-end` families, explicitly untouched at both 3' and 5') and as
//!   `spec_corpus_regressions::an_insertion_at_the_cds_end_is_not_a_fixed_point`,
//!   which pins `c.*1delinsCTT` -> `c.72_*1insCT` -> `c.72delinsCCT`. Plausibly
//!   one root cause. The **symptom** is new: those pin an output that moves on a
//!   second pass, this one an output that `parse_hgvs` rejects outright, which
//!   the idempotency oracle structurally cannot see (it verifies by
//!   re-normalizing, and an unparseable output never gets that far — the gap
//!   `FERRO_ASSERT_REPARSE` exists to cover).
//! - `defect_emitted_form_is_refused_by_its_own_normalizer` is **genomic**, with
//!   no region boundary in play, whereas the pinned class is coding-axis and
//!   CDS-end. I did not find it in that class; if it turns out to share the
//!   mechanism, this pin should move there rather than be duplicated.
//!
//! Neither is caught by an existing armed-oracle run, for the same reason: no
//! test in the tree normalizes either input, and an oracle only fires on
//! normalizations some test performs.
//!
//! # Not covered here
//!
//! `DNA/duplication.md:90`'s own coordinates cannot be exercised on a synthetic
//! fixture — they are intronic offsets against a genomic parent, and ferro dies
//! on them before reaching any behaviour this module is about
//! (`Could not convert intronic position 1598+-703 to genomic`, itself malformed).
//! That needs a cdot-backed transcript, so it is left to the reference-aware
//! lane rather than faked here.

use std::panic::{catch_unwind, AssertUnwindSafe};

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, NormalizeConfig, Normalizer};

use crate::common::cis_apply_oracle::apply_with;
use crate::common::hg38_window::{base_at, local_desc, provider, HG38_WINDOW, LOCAL_CONTIG};

/// Normalize a local-coordinate genomic description in strict mode.
fn normalize_g(body: &str) -> Result<String, String> {
    let input = local_desc(body);
    let variant: HgvsVariant = parse_hgvs(&input).map_err(|e| format!("parse: {e}"))?;
    Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&variant)
    .map(|v| v.to_string())
    .map_err(|e| format!("{e}"))
}

/// A coding transcript with a real CDS/3'UTR seam.
///
/// The CDS ends at `c.21`, so `c.21` is the last coding base and `c.*1` is the
/// first 3'UTR base — the seam every test below is anchored on. The 3'UTR opens
/// with `G` so `c.*1G>A` is a well-formed substitution against it.
fn seam_transcript() -> MockProvider {
    const UTR5: &str = "CCCCCCCCCC";
    const CDS: &str = "ATGCCCGGGCATGACACCTAA";
    const UTR3: &str = "GTGTGTGTGTGTGTGTGTGT";
    assert_eq!(CDS.len(), 21);
    assert_eq!(UTR3.as_bytes()[0], b'G', "c.*1 must be G");

    let sequence = format!("{UTR5}{CDS}{UTR3}");
    let cds_start = (UTR5.len() + 1) as u64;
    let cds_end = (UTR5.len() + CDS.len()) as u64;
    let tx_len = sequence.len() as u64;

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_SEAM.1".to_string(),
        Some("SEAM_TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(cds_start),
        Some(cds_end),
        vec![Exon::new(1, 1, tx_len)],
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

fn seam_normalizer() -> Normalizer<MockProvider> {
    Normalizer::with_config(
        seam_transcript(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
}

/// Normalize, catching a panic so the failure reports the defect rather than an
/// opaque unwind. `Err(String)` is a refusal; a panic is reported as such.
fn seam_normalize(input: &str) -> Result<String, String> {
    let outcome = catch_unwind(AssertUnwindSafe(|| {
        let variant: HgvsVariant = parse_hgvs(input).map_err(|e| format!("parse: {e}"))?;
        seam_normalizer()
            .normalize(&variant)
            .map(|v| v.to_string())
            .map_err(|e| format!("{e}"))
    }));
    match outcome {
        Ok(r) => r,
        Err(_) => Err("PANIC".to_string()),
    }
}

// ---------------------------------------------------------------------------
// DEFECT — the normalizer emits a description its own parser rejects.
// ---------------------------------------------------------------------------

/// **PINNED FINDING** — a well-defined allele at the CDS/3'UTR seam is refused.
///
/// `c.[21_*1del;21_*1insAC]` is well-formed HGVS: `c.21` is the last CDS base,
/// `c.*1` the first 3'UTR base, the anchor names two adjacent positions, and the
/// input round-trips through `parse_hgvs` unchanged. It also denotes a single
/// unambiguous sequence — delete both bases, insert `AC` at the junction between
/// them — which is exactly `c.21_*1delinsAC`.
///
/// Ferro refuses it as `OverlapConflictingEdits / W5002`, reporting two
/// coincident edits at `c.21`. Correct behaviour is to accept and emit
/// `c.21_*1delinsAC`; the members do not compete for a base, because the
/// deletion removes the bases the insertion's junction is positioned against.
/// That is the `#1406` deletion-exemption family.
///
/// Pinned rather than fixed because the fix belongs to whichever detector owns
/// the exemption, and #1749 records that overlap is defined four times on three
/// geometries. The refusal is measured on both this branch and its merge base,
/// so it is not a consequence of the dup+ins ruling implemented here.
///
/// An **earlier version of this test asserted the pair was accepted and emitted
/// an unparseable string.** That premise was never true on either arm — it was
/// my error, not a behaviour change — and it is recorded here because a test
/// whose premise is false reads exactly like a defect finding.
#[test]
fn pinned_seam_allele_is_refused_though_it_denotes_one_sequence() {
    let input = "NM_SEAM.1:c.[21_*1del;21_*1insAC]";
    parse_hgvs(input).expect("the input itself is valid HGVS");

    let err = seam_normalize(input).expect_err("PINNED FINDING: currently refused");
    assert!(
        err.contains("W5002") && err.contains("OverlapConflictingEdits"),
        "PINNED FINDING — expected the W5002 overlap refusal, got: {err}"
    );
}

/// The genomic control for the test above.
///
/// Included so a future reader does not have to re-derive that the seam is the
/// variable. This one passes today; if it ever fails, the defect above is not
/// seam-specific and its diagnosis needs redoing.
#[test]
fn control_the_genomic_analogue_of_the_seam_shape_is_clean() {
    let input = format!("{LOCAL_CONTIG}:g.[302_303del;302_303insAC]");
    let variant: HgvsVariant = parse_hgvs(&input).expect("parse");
    let out = Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&variant)
    .expect("genomic analogue normalizes")
    .to_string();
    parse_hgvs(&out).unwrap_or_else(|e| panic!("genomic control emitted unparseable {out}: {e}"));
}

// ---------------------------------------------------------------------------
// DEFECT — an out-of-range coding position behaves differently inside an allele.
// ---------------------------------------------------------------------------

/// A position past `cds_end` is refused alone and is now refused in an allele.
///
/// `c.22` is past this transcript's `cds_end` (`c.21`), so it names no base on
/// the `c.` axis — `background/numbering.md:21` starts `c.1` at the `A` of the
/// initiator and ends at "the last nucleotide of the translation termination
/// (stop) codon"; past that the axis is `*N`. Ferro refuses `c.22del` outright
/// in strict mode.
///
/// This used to be a PINNED FINDING: inside a cis allele the same coordinate was
/// silently accepted and remapped to `*1` (`c.[22_23insG;*1G>A]` -> `c.21dup`),
/// which changed the description's meaning without saying so. #2018 reconciled
/// the two paths — the allele path now applies the same mode-dependent rule as
/// the lone path (see `issue_2018_past_cds_end_mode_dependent`). This test now
/// asserts the reconciled behaviour: in the strict `seam_normalize` mode both
/// the lone position and the allele refuse with the same `W4004`.
#[test]
fn pinned_past_cds_end_position_is_refused_alone_and_in_an_allele() {
    let alone = seam_normalize("NM_SEAM.1:c.22del")
        .expect_err("c.22 is past cds_end and is refused alone in strict mode");
    assert!(
        alone.contains("W4004") && alone.contains("PositionPastEnd"),
        "expected the lone position to refuse with W4004, got {alone:?}"
    );

    // RECONCILED (#2018): the allele path now refuses the past-cds-end
    // coordinate with the same W4004 the lone description gets, rather than
    // silently remapping it across the seam and normalizing to `c.21dup`.
    let in_allele = seam_normalize("NM_SEAM.1:c.[22_23insG;*1G>A]")
        .expect_err("c.22 in a cis allele is now refused with W4004, as it is alone");
    assert!(
        in_allele.contains("W4004") && in_allele.contains("PositionPastEnd"),
        "expected the allele to refuse with W4004, as the lone position does, got {in_allele:?}"
    );
}

/// A straddling anchor past `cds_end` must be refused, not crash.
///
/// `c.22_*1ins` mixes a coordinate past `cds_end` with a valid 3'UTR
/// coordinate. The anchor is nonsense, and the correct response is the same
/// refusal `c.22del` gets. Instead the pair reaches the 3'-shuffler, where
/// `new_end - start` underflows on `u64` (`src/normalize/shuffle.rs:94`).
///
/// The consequence depends on the build profile, and the shipped one is the
/// worse half: `[profile.release]` does not set `overflow-checks`, so a release
/// binary does not panic — it wraps and emits invalid HGVS (a single-position
/// insertion, or a reversed range), neither of which re-parses.
///
/// Sibling of the open **#1917**, which is the same "nothing establishes
/// `start <= end`" underflow at a different call site
/// (`src/normalize/mod.rs`), and of the closed #472 / #488.
///
/// **PINNED FINDING.** This pins the underflow as it stands today, so the
/// tripwire fires the moment the anchor is validated. A debug build panics on
/// the subtraction; the correct behaviour is the same `W4004` refusal
/// `c.22del` receives, and a clean refusal will fail this assertion.
#[test]
fn pinned_straddling_past_cds_end_anchor_underflows_the_shuffler() {
    match seam_normalize("NM_SEAM.1:c.[22_*1insG;*1G>A]") {
        Err(e) if e == "PANIC" => { /* the pinned defect */ }
        other => panic!(
            "PINNED FINDING — `shuffle.rs:94` underflows on a past-cds-end \
             anchor instead of refusing it, and this pins that. Correct \
             behaviour is the W4004 refusal `c.22del` gets; a release build \
             wraps instead and emits invalid HGVS. Sibling of open #1917. \
             Got {other:?} — if that is a clean refusal, the defect is fixed \
             and this pin should be deleted."
        ),
    }
}

// ---------------------------------------------------------------------------
// POLICY — ferro contradicts a published spec example.
// ---------------------------------------------------------------------------

/// A `dup` and an insertion at the junction it writes to.
///
/// `DNA/duplication.md:90` publishes
/// `NC_000001.11(NM_206933.2):c.[675-542_1211-703dup;1211-703_1211-702insGTAAA]`
/// as a **correct** description, glossed "a duplication of the sequence from
/// nucleotide position `c.675-542` to `c.1211-703`, **followed by** the
/// insertion of the sequence `GTAAA`". The insertion's junction is immediately
/// 3' of the duplication's last base — the junction the `dup` writes its copy
/// into — and the only thing the attached NOTE rejects is the `dupins` spelling.
///
/// Ferro refuses this geometry (`#1446`) on the ground that two writers at one
/// slot have no order between them. The spec's gloss supplies exactly that
/// order.
///
/// **POLICY, not DEFECT**: the spec contradicts itself in this family —
/// `delins.md:86-89` mandates merging an abutting pair while `duplication.md:90`
/// and `complex.md:113,117` publish zero-separation pairs split — and nothing
/// scopes either against the other. The expectation below follows the published
/// example; if a ruling goes the other way, this test is what should change.
#[test]
fn policy_duplication_md_90_publishes_a_dup_with_an_insertion_at_its_junction() {
    let input = format!("{LOCAL_CONTIG}:g.[302_304dup;304_305insG]");
    let variant: HgvsVariant = parse_hgvs(&input).expect("parse");
    let result = Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&variant);
    assert!(
        result.is_ok(),
        "DNA/duplication.md:90 publishes this geometry as correct, but it was refused: {:?}",
        result.err().map(|e| e.to_string())
    );
}

/// The normalizer must not emit a form its own **normalizer** then refuses.
///
/// **This was a live defect and the junction ruling fixed it.** Measured by A/B:
/// with `origin/main`'s `overlap.rs` and `merge.rs` substituted in, this test
/// fails with
///
/// ```text
/// strict mode emitted NC_TESTWIN.1:g.[200_203T[6];203_204insT] and then
/// refused it: ... 2 coincident cis-allele edits (repeat, ins) ...
/// ```
///
/// The refusal was the `(repeat, ins)` coincident-junction conflict, which is
/// exactly what `duplication.md:90` orders and this ruling now accepts. Kept as
/// a guard rather than deleted: the emitted form is *parseable*, so a re-parse
/// check cannot see this class and only an idempotence check catches it.
///
/// Found by `adjacency_confluence_proptest::normalization_is_idempotent`;
/// pinned here deterministically so it survives any later narrowing of that
/// generator, and reproduced against the assembly at
/// `NC_000001.11:g.[1001501_1001502insTT;1001502dup]`.
///
/// `[202_203insTT;203dup]` is accepted in strict mode and normalizes to
/// `[200_203T[6];203_204insT]`. That output **re-parses** — so the round-trip
/// check in `insertion_adjacency_corpus` does not see it — but feeding it back
/// to the same strict normalizer is refused as a `(repeat, ins)` overlap
/// conflict.
///
/// So strict mode emits a description strict mode rejects. Any consumer that
/// normalizes defensively, or that re-normalizes stored output after a version
/// bump, turns a stored value into a hard error.
///
/// Distinct from the seam defect above in two ways worth keeping: this one is on
/// the **genomic** axis with no region boundary involved, and the output is
/// parseable, so it needs an idempotence check rather than a re-parse check to
/// catch at all.
#[test]
fn an_emitted_form_is_accepted_by_its_own_normalizer() {
    let once = normalize_g("[202_203insTT;203dup]").expect("this allele is accepted");
    let stripped = once
        .strip_prefix(&format!("{LOCAL_CONTIG}:g."))
        .expect("accession survives rendering");

    parse_hgvs(&once).expect("precondition: the emitted form does re-parse");

    let twice = normalize_g(stripped);
    assert!(
        twice.is_ok(),
        "strict mode emitted {once} and then refused it: {}",
        twice.unwrap_err()
    );
    assert_eq!(
        twice.as_deref().unwrap(),
        once,
        "normalization is not idempotent"
    );
}

/// An abutting `inv` + `sub` is re-partitioned rather than kept as published.
///
/// `[302_310inv;311G>A]` becomes `[302del;304_309inv;311delinsAA]`. The
/// **sequence is preserved** — asserted below through the independent oracle, so
/// this is a form question and not a wrong answer — but the shape is neither of
/// the two the spec offers for this geometry: not the published split
/// (`DNA/complex.md:113`, an inversion with a substitution at its break point
/// kept as two members) and not a single merged `delins`.
///
/// **POLICY.** Re-derivation from the resulting sequence is the ledger's decided
/// rule (`canonical-form-choice-when-both-legal`, 2026-08-07: emit the form that
/// falls out, do not preserve the input's spelling), so this is that rule
/// working as decided. What it collides with is a published worked example, and
/// the spec does not scope one against the other — `complex.md:113` and `:117`
/// publish zero-separation pairs split while `delins.md:86-89` mandates merging
/// one. Ferro picks a third option.
///
/// The expectation asserted here is the published one. If a ruling endorses
/// re-partitioning at this scale, this test is what should change — and the
/// ruling should say what `complex.md:113` is then taken to mean.
#[test]
#[ignore = "#2007: policy — complex.md:113 publishes the split; ferro re-partitions. Sequence is preserved, so this is a form question awaiting a ruling."]
fn policy_an_abutting_inv_and_sub_are_repartitioned() {
    assert_eq!(
        base_at(311),
        b'G',
        "case assumes the base 3' of the inversion"
    );
    let body = "[302_310inv;311G>A]";
    let got = normalize_g(body).expect("accepted");

    // Sequence first: without this, a form complaint could be masking a wrong
    // answer, and the two deserve very different responses.
    let p = provider();
    let before = apply_with(&p, HG38_WINDOW, &local_desc(body)).expect("input denotes a sequence");
    let after = apply_with(&p, HG38_WINDOW, &got).expect("output denotes a sequence");
    assert_eq!(before, after, "re-partitioning changed the sequence");

    assert_eq!(
        got,
        local_desc(body),
        "DNA/complex.md:113 publishes an inv with an abutting substitution as two members; \
         ferro re-partitioned it into a third shape (sequence is preserved)"
    );
}

/// A `dup` and an `ins` at one junction settle identically in either authored
/// order — the defect this started as, closed by the ruling in this branch.
///
/// Found by
/// `adjacency_confluence_proptest::member_order_is_not_part_of_what_an_allele_denotes`.
/// Before the `duplication.md:90` ruling both spellings were *refused*, and the
/// two refusals disagreed on the location and on the member order:
///
/// ```text
/// g.200_201 has 2 coincident cis-allele edits (ins, dup)
/// g.200     has 2 coincident cis-allele edits (dup, ins)
/// ```
///
/// Measured on the merge base, so the divergence was real. Accepting the pair
/// removes the diagnostic entirely and the two orders now converge on one
/// description, which is the stronger property — an order-dependent *message*
/// was only ever a hint that two paths reported one conflict.
///
/// Kept as a guard rather than deleted: it is the regression test for that
/// convergence, and it fails loudly if either order starts refusing again.
#[test]
fn a_dup_and_an_ins_at_one_junction_settle_identically_in_either_order() {
    let a = normalize_g("[200_201insT;200dup]").expect("the dup+ins pair is accepted");
    let b = normalize_g("[200dup;200_201insT]").expect("member order is not part of the allele");
    assert_eq!(a, b, "the two authored orders settled apart: {a} vs {b}");
    assert_eq!(
        a, "NC_TESTWIN.1:g.200_203T[6]",
        "both orders settle on the tract spelling"
    );
}
