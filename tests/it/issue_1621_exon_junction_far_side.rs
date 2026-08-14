//! The exon-junction clamp binds from the FAR side too (#1621).
//!
//! `general.md:44` exempts deletions and duplications around exon/exon
//! junctions from the 3' rule on `c.`/`r.`/`n.` references. Ferro applied that
//! exemption in one direction only: a description approaching a junction was
//! stopped at it, but one already spelled *past* it was left standing, so the
//! two spellings of one transcript sequence were two fixed points.
//!
//! The operator ruling `exon-junction-dup-converge-from-the-far-side` (in
//! `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`, `decided`)
//! reads the exemption as a clamp binding in **both** directions: the canonical
//! position is the most 3' position that does not cross a junction, *reached
//! from either side*. `DNA/duplication.md` states the outcome three times —
//! `:26`, `:60` and `:148` — the last giving the reason, that `c.3922dup`
//! translated back to a genomic position lands "at the wrong nucleotide, in the
//! wrong exon".
//!
//! The ruling's own row, `LRG_199t1:c.3922dup` -> `c.3921dup`, is asserted by
//! `spec_worked_examples::the_decided_target_is_convergence_on_the_near_side`
//! against the committed reference slice. This module carries the rest of the
//! scope the ruling names — **deletions as well as duplications, and the `n.`
//! axis as well as `c.`** — on a synthetic transcript built here, and it
//! asserts **both directions** of the clamp on every one of them. A fix that
//! only pulls the far side back would leave the near side free to shift across
//! the junction, and a fix that only holds the near side is what shipped
//! before; each case below therefore pins the near-side spelling as a fixed
//! point *and* the far-side spelling converging onto it.
//!
//! **The fixture must be multi-exon, and that is the whole point.** A clamp at
//! an exon/exon junction is unreachable on a single-exon transcript, so a
//! generator that emits one reports a structural zero for this entire class —
//! the blindness recorded on `dump_normalized_corpus` as #1478. The transcript
//! below has two exons with the run straddling their junction, and
//! `the_fixture_really_straddles_a_junction` asserts that geometry directly
//! rather than trusting it.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Transcript bases, 40 of them, with a two-base `T` run at transcript 20-21.
///
/// The exon boundary below falls **between** those two `T`s, so they are one
/// contiguous run in spliced transcript space and two different exons in
/// genomic space — the reduced form of `LRG_199t1`'s `c.3921`/`c.3922` pair,
/// which sits at the exon 27/28 junction of `NM_004006.2` and projects 2,790 bp
/// apart for exactly this reason.
///
/// The flanks are chosen so the run is exactly two bases: transcript 19 is `A`
/// and transcript 22 is `C`, so neither end of the run extends. Without that
/// the 3' shuffle inside exon 2 would move the far-side spelling for an
/// unrelated reason and the test would pass while proving nothing.
const TX_SEQUENCE: &str = "ACGACGACGACGACGACGATTCGACGACGACGACGACGAC";

/// Transcript position of the last base of exon 1 — the near side of the clamp.
const NEAR_TX: u64 = 20;
/// Transcript position of the first base of exon 2 — the far side.
const FAR_TX: u64 = 21;

/// `cds_start` is deliberately **not** 1.
///
/// With `cds_start = 1` the `c.` axis is the identity on transcript
/// coordinates, which is the second half of #1478's blindness: a bug in the
/// axis translation is invisible, and the variant can never sit in a 5'UTR. At
/// 5 the translation is a real offset (`c.N` is transcript `N + 4`) and
/// transcript 1-4 form a genuine `c.-4`..`c.-1` 5'UTR.
const CDS_START: u64 = 5;
const CDS_END: u64 = 40;

/// `c.` position of a transcript position, under [`CDS_START`].
fn cds_of(tx: u64) -> u64 {
    tx - CDS_START + 1
}

/// Exon 1 is transcript 1-20 at genomic 101-120; exon 2 is transcript 21-40 at
/// genomic 201-220. The 80-base gap between them is a real intron, so the
/// junction is a junction in genomic space and not merely an annotation.
fn exons() -> Vec<Exon> {
    vec![
        Exon::with_genomic(1, 1, NEAR_TX, 101, 120),
        Exon::with_genomic(2, FAR_TX, 40, 201, 220),
    ]
}

/// The genomic contig the two exons are placed on.
///
/// The intron's first base is `A`, never `T`. That is load-bearing: a del/dup
/// resting on exon 1's last base triggers the #670 junction-crossing
/// continuation, which re-runs the 3' shuffle in genomic space where the intron
/// *is* visible. With a `T` there the near-side answer would legitimately shift
/// into the intron as `c.16+1dup`, and this module would be measuring that
/// instead of the clamp.
fn contig() -> String {
    let filler = "ACGT".repeat(25); // 100 bases
    let intron = "AC".repeat(40); // 80 bases, first base `A`
    let mut s = String::new();
    s.push_str(&filler); // g.1-100
    s.push_str(&TX_SEQUENCE[..NEAR_TX as usize]); // g.101-120 = tx 1-20
    s.push_str(&intron); // g.121-200
    s.push_str(&TX_SEQUENCE[NEAR_TX as usize..]); // g.201-220 = tx 21-40
    s.push_str(&filler); // g.221-320
    s
}

fn transcript(id: &str, cds: Option<(u64, u64)>) -> Transcript {
    Transcript::new(
        id.to_string(),
        None,
        Strand::Plus,
        TX_SEQUENCE.to_string(),
        cds.map(|(s, _)| s),
        cds.map(|(_, e)| e),
        exons(),
        Some("NC_SYNTH.1".to_string()),
        Some(101),
        Some(220),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    )
}

/// A provider carrying one coding and one non-coding transcript over the same
/// two-exon geometry, so the `c.` and `n.` axes are compared on identical bases.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_transcript(transcript("NM_SYNTH.1", Some((CDS_START, CDS_END))));
    provider.add_transcript(transcript("NR_SYNTH.1", None));
    provider.add_genomic_sequence("NC_SYNTH.1", contig());
    provider
}

fn normalize(input: &str) -> String {
    let normalizer = Normalizer::new(provider());
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize {input}: {e}"));
    format!("{normalized}")
}

/// The premise, asserted rather than assumed.
///
/// If any of these three drifts — the run stops being two bases, or the exon
/// boundary stops falling inside it — every assertion below would still pass
/// while testing nothing about an exon junction. This is the structural check
/// the corpus generators kept lacking (#1456, #1460, #1478): name the property
/// the change keys on, and show the fixture can carry it.
#[test]
fn the_fixture_really_straddles_a_junction() {
    let bases = TX_SEQUENCE.as_bytes();
    assert_eq!(bases.len(), 40, "transcript length");
    assert_eq!(
        (
            bases[NEAR_TX as usize - 2] as char,
            bases[NEAR_TX as usize - 1] as char,
            bases[FAR_TX as usize - 1] as char,
            bases[FAR_TX as usize] as char,
        ),
        ('A', 'T', 'T', 'C'),
        "transcript {NEAR_TX} and {FAR_TX} must be the same base, and the run must stop there — \
         otherwise the two spellings are not two spellings of one sequence"
    );

    let exons = exons();
    assert_eq!(
        exons.len(),
        2,
        "a single-exon fixture cannot test a junction"
    );
    assert_eq!(
        (exons[0].end, exons[1].start),
        (NEAR_TX, FAR_TX),
        "the exon/exon junction must fall between the two bases of the run"
    );
    assert_eq!(
        (exons[0].genomic_end, exons[1].genomic_start),
        (Some(120), Some(201)),
        "the two exons must be genuinely far apart in genomic space, so crossing the junction \
         is the ~2,790 bp class of error `DNA/duplication.md:148` describes and not an off-by-one"
    );
    assert_eq!(contig().len(), 320, "contig length");
}

/// A duplication spelled past the junction is pulled back to the clamp; one
/// already at the clamp stays.
///
/// This is the ruling's own shape, reduced: `c.17dup` is `LRG_199t1:c.3922dup`
/// and `c.16dup` is `c.3921dup`.
#[test]
fn a_duplication_converges_on_the_near_side_of_the_junction() {
    let near = format!("NM_SYNTH.1:c.{}dup", cds_of(NEAR_TX));
    let far = format!("NM_SYNTH.1:c.{}dup", cds_of(FAR_TX));

    assert_eq!(
        normalize(&near),
        near,
        "the near side is at the clamp and must not shift across the junction — this is the \
         direction ferro already had, and a pull-back that broke it would trade one half of the \
         exception for the other"
    );
    assert_eq!(
        normalize(&far),
        near,
        "the far side is past the clamp and must be pulled back to it \
         (`exon-junction-dup-converge-from-the-far-side`, `DNA/duplication.md:26`/`:60`/`:148`)"
    );
}

/// The same, for a deletion. The ruling's scope is "deletions and
/// duplications", and `general.md:44` names both; a fix keyed on `dup` alone
/// would satisfy the ruling's worked example and miss half its scope.
#[test]
fn a_deletion_converges_on_the_near_side_of_the_junction() {
    let near = format!("NM_SYNTH.1:c.{}del", cds_of(NEAR_TX));
    let far = format!("NM_SYNTH.1:c.{}del", cds_of(FAR_TX));

    assert_eq!(normalize(&near), near, "the near side stays at the clamp");
    assert_eq!(normalize(&far), near, "the far side is pulled back to it");
}

/// The `n.` axis, on the same bases and the same geometry.
///
/// `general.md:44` names `c.`, `r.` and `n.` together, and ferro's forward
/// clamp is already applied on all three. The far-side half must not be `c.`
/// only, or the two axes disagree about which spelling of one sequence is
/// canonical.
#[test]
fn the_noncoding_axis_converges_from_the_far_side_too() {
    for edit in ["dup", "del"] {
        let near = format!("NR_SYNTH.1:n.{NEAR_TX}{edit}");
        let far = format!("NR_SYNTH.1:n.{FAR_TX}{edit}");

        assert_eq!(
            normalize(&near),
            near,
            "n.{NEAR_TX}{edit}: the near side stays at the clamp"
        );
        assert_eq!(
            normalize(&far),
            near,
            "n.{FAR_TX}{edit}: the far side is pulled back to it"
        );
    }
}

/// The pulled-back answer is a fixed point.
///
/// A pull-back that re-fires on its own output would be a non-idempotency of
/// the kind `FERRO_ASSERT_IDEMPOTENT` exists to catch, and it is the obvious
/// failure mode of implementing this as "step back over the junction" rather
/// than as "find where the run begins": the answer's own run begins inside its
/// own exon, so nothing pulls it further.
#[test]
fn convergence_is_a_fixed_point() {
    for input in [
        format!("NM_SYNTH.1:c.{}dup", cds_of(FAR_TX)),
        format!("NM_SYNTH.1:c.{}del", cds_of(FAR_TX)),
        format!("NR_SYNTH.1:n.{FAR_TX}dup"),
        format!("NR_SYNTH.1:n.{FAR_TX}del"),
    ] {
        let once = normalize(&input);
        let twice = normalize(&once);
        assert_eq!(twice, once, "{input} normalized twice");
    }
}

/// A del/dup whose run does not straddle a junction is untouched.
///
/// The negative control. The pull-back runs for every `c.`/`n.` del/dup, so
/// without this the module could not distinguish "the clamp binds from both
/// sides" from "everything drifts 5'". Transcript 8 sits mid-exon in a
/// period-3 `ACG` tract, and transcript 22 is the `C` immediately 3' of the run
/// — the first position past the junction that the pull-back must **not** reach.
#[test]
fn a_run_inside_one_exon_is_not_pulled_anywhere() {
    for tx in [8_u64, 22] {
        let input = format!("NR_SYNTH.1:n.{tx}del");
        let out = normalize(&input);
        let position: u64 = out
            .trim_start_matches("NR_SYNTH.1:n.")
            .trim_end_matches("del")
            .parse()
            .unwrap_or_else(|_| panic!("unexpected shape for {input}: {out}"));
        assert!(
            position >= tx,
            "{input} -> {out}: a 3'-direction normalization must never move a variant 5', and \
             this run does not reach a junction"
        );
    }
}

/// The pull-back must not swallow a warning the input earned.
///
/// The clamp re-normalizes a description derived from its own *output*, so the
/// second pass never sees the input: a `REFSEQ_MISMATCH` raised because the
/// input stated a base the reference contradicts is not re-derivable there
/// (`normalize_na_edit` reads a `dup` payload from the reference, #219). If the
/// second pass's warnings simply replaced the first's, the disclosure would
/// vanish.
///
/// The near-side spelling is the control, and it is what makes this a real
/// assertion rather than a change detector: the two inputs differ only in which
/// side of the junction they are spelled on, and they produce the **same**
/// output, so any disagreement about the warning is the clamp's doing and
/// nothing else's.
#[test]
fn the_pull_back_keeps_the_warnings_the_input_earned() {
    use ferro_hgvs::normalize::NormalizationWarning;

    fn mismatch_reported(input: &str) -> (String, bool) {
        let normalizer = Normalizer::new(provider());
        let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
        let result = normalizer
            .normalize_with_diagnostics(&variant)
            .unwrap_or_else(|e| panic!("normalize {input}: {e}"));
        let reported = result
            .warnings
            .iter()
            .any(|w| matches!(w, NormalizationWarning::RefSeqMismatch { .. }));
        (format!("{}", result.result), reported)
    }

    // Transcript 20 and 21 are both `T`; stating `G` earns a corrected
    // REFSEQ_MISMATCH on either spelling.
    let far = format!("NM_SYNTH.1:c.{}dupG", cds_of(FAR_TX));
    let near = format!("NM_SYNTH.1:c.{}dupG", cds_of(NEAR_TX));
    let expected = format!("NM_SYNTH.1:c.{}dup", cds_of(NEAR_TX));

    let (near_out, near_reported) = mismatch_reported(&near);
    let (far_out, far_reported) = mismatch_reported(&far);

    assert_eq!(
        (near_out.as_str(), far_out.as_str()),
        (expected.as_str(), expected.as_str()),
        "the premise: both spellings converge on the near side, so the only thing that can \
         differ between them is whether the warning survived"
    );
    assert!(
        near_reported,
        "the control must report the mismatch — if it does not, this test is measuring \
         something other than the clamp"
    );
    assert!(
        far_reported,
        "the pulled-back spelling dropped the REFSEQ_MISMATCH its input earned: the clamp's \
         re-normalization must ADD to the first pass's warnings, never replace them"
    );
}

/// An `n.*N` position is refused by the clamp rather than read as a plain
/// transcript offset.
///
/// `TxPos::is_intronic` reads only `offset`, and the `*` marker is a spelling
/// flag on the position rather than a transcript region — so on the `n.` axis
/// there is no `boundary::axis_region_of` check standing behind the clamp the
/// way there is on `c.`. Read as a bare offset, `n.*21` would be pulled to the
/// junction's near side and rebuilt **without** the flag, silently re-zoning the
/// description into a zone `background/numbering.md:52` does not give this axis.
///
/// Built through the AST rather than parsed, deliberately: #1748 refuses the
/// spelling at parse in every mode, while `TxPos::downstream` stays public API.
/// That is exactly the reachability the re-parse oracle's own exemption list
/// records for this shape — no string entry point mints one, and a Rust caller
/// can.
#[test]
fn a_downstream_noncoding_position_is_not_pulled_across_the_junction() {
    use ferro_hgvs::hgvs::edit::NaEdit;
    use ferro_hgvs::hgvs::interval::TxInterval;
    use ferro_hgvs::hgvs::location::TxPos;
    use ferro_hgvs::hgvs::variant::{Accession, HgvsVariant, LocEdit, TxVariant};

    let mut position = TxPos::new(FAR_TX as i64);
    position.downstream = true;
    let variant = HgvsVariant::Tx(TxVariant {
        accession: Accession::new("NR", "SYNTH", Some(1)),
        gene_symbol: None,
        loc_edit: LocEdit::new(
            TxInterval::point(position),
            NaEdit::Deletion {
                sequence: None,
                length: None,
            },
        ),
    });
    assert_eq!(
        format!("{variant}"),
        format!("NR_SYNTH.1:n.*{FAR_TX}del"),
        "the premise: the constructed position really carries the downstream marker"
    );

    let normalizer = Normalizer::new(provider());
    let normalized = normalizer.normalize(&variant).expect("normalize");
    assert_eq!(
        format!("{normalized}"),
        format!("NR_SYNTH.1:n.*{FAR_TX}del"),
        "a downstream position names no transcript offset the clamp can read, so it must be \
         left alone — pulling it back drops the `*` and states a different variant"
    );
}
