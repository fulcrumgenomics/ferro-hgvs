//! Two in-band markers that a boundary was dropping or leaking.
//!
//! Both are the same class — data the AST carries in-band, silently discarded
//! at one seam or emitted raw at another — and both are pinned here so the two
//! halves are read together.
//!
//! Both classes had further instances beyond the first two found, and the
//! sweeps that missed them are worth naming. The drop was fixed at
//! `CoordinateMapper::tx_to_cds` while **two service handlers flattened a
//! `TxPos` to a bare `u64` themselves** — they go through `CdotMapper`, whose
//! entry points take the flattened number and so never see the flag, which is
//! why a fix at the `CoordinateMapper` seam could not reach them. Those two are
//! pinned beside the code, in `service::handlers::tx_position`,
//! `service::handlers::convert` and `service::handlers::vcf_convert`, because
//! the handler functions are private. The leak was fixed across the four
//! **axis** `Display` impls while a **fifth** `Display` renders an offset —
//! `IvsPos` — and is pinned at the bottom of this file.
//!
//! The lesson both halves share: the sweep was over the wrong population. "Every
//! axis" is not "every renderer", and "the conversion seam" is not "every
//! caller that flattens a position".
//!
//! # A THIRD drop route is KNOWN, MEASURED, AND DELIBERATELY LEFT OPEN
//!
//! `CoordinateMapper::tx_to_genomic` — the direct sibling of the `tx_to_cds`
//! this file's first section fixes, on the same type — reads `tx_pos.base` and
//! screens only `base < 1`. It never consults `downstream`. Measured on a
//! single exon at tx 1..10 / genome 1000..1009:
//!
//! ```text
//! tx_to_genomic(n.5)  = Ok(Some(1004))
//! tx_to_genomic(n.*5) = Ok(Some(1004))     <- collapse, same as the two fixed here
//! ```
//!
//! It is **not fixed here** because its blast radius is not the service layer's:
//! it has production callers in `vcf/from_hgvs.rs` (`convert_tx`, the library's
//! own n.→VCF path), `normalize/mod.rs`, `equivalence/checker.rs` and three more
//! inside `mapper.rs` itself, and its return type (`Result<Option<u64>, _>`)
//! makes "decline" and "error" two different answers that want their own
//! decision. That is a change owing its own measurement and disclosure.
//!
//! Recorded here rather than left to be re-found: an earlier census of this
//! defect enumerated **13** other marker-consulting sites and concluded the two
//! service handlers were the only unguarded ones. All 13 do hold — but the
//! census was of sites that *already consult* the marker, which cannot discover
//! one that never mentions it. `tx_to_genomic` is invisible to a `grep` for
//! `downstream`, which is exactly why it survived.
//!
//! # Reachability: real defect, unreachable spelling
//!
//! All three routes are reached only by a **Rust caller constructing a
//! `TxPos` directly**. `n.*N` is refused at parse in *every* mode by #1748 —
//! measured on this head, `parse_hgvs` and `parse_hgvs_lenient` both return
//! `E1003 InvalidPosition` for `NM_TEST.1:n.*5del` — so no string entry point,
//! the web service's included, can deliver the marker. `TxPos::downstream` is
//! nonetheless public API, and these guards are the same defence-in-depth the
//! PR's own `tx_to_cds` refusal is: the contract belongs with the function
//! whether or not today's grammar can reach it.
//!
//! **`TxPos.downstream` (`n.*N`).** `CoordinateMapper::tx_to_cds` classified a
//! position purely by comparing `base` against the CDS bounds, so `n.*5` was
//! answered as if it were in-transcript position 5. Per
//! `background/numbering.md:52` `n.` numbering runs "`n.1`, `n.2`, `n.3`, ...,
//! etc., from the first to the last nucleotide of the reference sequence" — it
//! admits no `*` zone at all — and `:54` states that "it is **not** allowed to
//! describe variants in nucleotides beyond the boundaries of a transcript
//! reference sequence, using that transcript reference sequence". So `n.*N`
//! names a nucleotide the `n.` axis cannot number, and there is no `c.` position
//! to map it onto: the conversion must refuse rather than answer for a
//! different base. This is not a new judgement — `project_*_direct` already
//! refused `*N` up front *because* `tx_to_cds` ignored the flag; the refusal has
//! moved to the function whose contract it belongs to.
//!
//! **The counter-evidence, recorded rather than left to be re-discovered.** The
//! competing shape was to *honour* the flag: `Transcript::sequence_length` is
//! available, so `n.*N` is mechanically `c.*(sequence_length - cds_end + N)` and
//! the mapping is exact. Two things decided against it. The derived `c.*M` would
//! itself be past the transcript's last base, so honouring converts one
//! description `numbering.md:54` forbids into another rather than into a legal
//! one. And the inverse would have to follow — `cds_to_tx` mapping an
//! over-length `c.*M` back to `n.*N` — which moves outputs on the coding axis,
//! where every real description lives, to serve an axis on which the shape does
//! not occur. What is genuinely on the other side of the ledger is ferro's own
//! poly-A carve-out (#797), which deliberately entertains a `c.*` position in
//! the post-transcriptional tail: ferro is not uniformly literal about `:54`.
//! That is a ferro policy rather than a clause, so it does not outrank `:52`
//! here, but it is the reason this is a decision and not a foregone conclusion.
//!
//! **The unknown-offset sentinels (`+?` / `-?`).** The parser stores them
//! in-band as `i64::MAX` / `i64::MIN`. `GenomePos` and `CdsPos` render them back
//! as `+?` / `-?`; `TxPos`, `RnaPos` and `IvsPos` did not, so a parsed `n.5+?`
//! printed as `5+9223372036854775807` and a `c.100+?` carried into IVS notation
//! printed as `IVS7+9223372036854775807`. Nothing adjudicates what these should
//! print: the spelling is the one the parser consumed, and `CdsPos::Display` has
//! rendered it that way since the sentinels were introduced.
//!
//! Four spellings move, not three: `parse_offset` is shared by every axis, so
//! `r.<b>-?` is accepted and moves exactly as `r.<b>+?` does.

use ferro_hgvs::convert::CoordinateMapper;
use ferro_hgvs::hgvs::location::{CdsPos, IvsNotation, IvsPos, RnaPos, TxPos};
use ferro_hgvs::hgvs::parser::position::{OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE};
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};

/// 5bp 5'UTR + 30bp CDS + 3bp 3'UTR over a single 38-base exon, matching the
/// shape `convert_tests` uses so the two read against the same geometry.
fn coding_transcript() -> Transcript {
    Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        "AAAAATGCCCAAAGGGTTTAGGCCCAAAGGGTTATAAA".to_string(),
        Some(6),
        Some(35),
        vec![Exon::new(1, 1, 38)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    )
}

// ---------------------------------------------------------------------------
// `TxPos.downstream` must not be silently dropped by tx -> CDS conversion
// ---------------------------------------------------------------------------

/// The defect in its bare form: `n.*5` and `n.5` are different nucleotides, so
/// the conversion must not answer the same thing for both.
#[test]
fn tx_to_cds_does_not_answer_a_downstream_position_as_its_in_transcript_twin() {
    let tx = coding_transcript();
    let mapper = CoordinateMapper::new(&tx);

    let plain = mapper.tx_to_cds(&TxPos::new(5));
    let downstream = mapper.tx_to_cds(&TxPos::downstream(5));

    assert!(
        plain.is_ok(),
        "control: a plain in-transcript n.5 must still convert"
    );
    assert!(
        downstream.is_err(),
        "n.*5 names a nucleotide past the transcript's last base \
         (background/numbering.md:52, :54) and has no c. position; \
         tx_to_cds answered {:?} instead of refusing",
        downstream.map(|p| p.to_string())
    );
}

/// The refusal is on the flag, not on the base — every `*N` is refused, and the
/// message names the flag so a caller can tell this decline from a bounds one.
#[test]
fn tx_to_cds_refuses_every_downstream_position_by_name() {
    let tx = coding_transcript();
    let mapper = CoordinateMapper::new(&tx);

    for base in [1, 5, 31, 4000] {
        let err = mapper
            .tx_to_cds(&TxPos::downstream(base))
            .expect_err("n.*{base} must be refused");
        let msg = err.to_string();
        assert!(
            msg.contains("n.*"),
            "the decline must name the notation it refused, got: {msg}"
        );
    }
}

/// An offset on a downstream position is refused for the same reason — the
/// anchor it is measured from is itself unnumberable.
#[test]
fn tx_to_cds_refuses_a_downstream_position_carrying_an_offset() {
    let tx = coding_transcript();
    let mapper = CoordinateMapper::new(&tx);
    assert!(
        mapper
            .tx_to_cds(&TxPos::downstream_with_offset(5, 10))
            .is_err(),
        "n.*5+10 must be refused like n.*5"
    );
}

/// The sibling direction is deliberately unchanged: `c.*N` is the 3'UTR of a
/// coding transcript, which is *inside* the transcript, so its `n.` form is a
/// plain base and `downstream` is correctly false. Pinned so the M1 fix is not
/// later "completed" by making `cds_to_tx` set the flag.
#[test]
fn cds_to_tx_keeps_a_coding_three_prime_utr_position_in_transcript() {
    let tx = coding_transcript();
    let mapper = CoordinateMapper::new(&tx);

    let n = mapper
        .cds_to_tx(&CdsPos::utr3(1))
        .expect("c.*1 is the first 3'UTR base and converts");
    assert_eq!(n.base, 36, "c.*1 is transcript base 36 in this fixture");
    assert!(
        !n.downstream,
        "c.*1 is inside the transcript, so its n. form must NOT carry the \
         past-the-end `*` marker; got {n}"
    );
}

// ---------------------------------------------------------------------------
// The unknown-offset sentinels must not escape as 19-digit integers
// ---------------------------------------------------------------------------

#[test]
fn tx_pos_renders_the_unknown_offset_sentinels_as_question_marks() {
    assert_eq!(
        TxPos::with_offset(5, OFFSET_UNKNOWN_POSITIVE).to_string(),
        "5+?"
    );
    assert_eq!(
        TxPos::with_offset(5, OFFSET_UNKNOWN_NEGATIVE).to_string(),
        "5-?"
    );
    assert_eq!(
        TxPos::downstream_with_offset(5, OFFSET_UNKNOWN_POSITIVE).to_string(),
        "*5+?"
    );
}

#[test]
fn rna_pos_renders_the_unknown_offset_sentinels_as_question_marks() {
    assert_eq!(
        RnaPos::with_offset(5, OFFSET_UNKNOWN_POSITIVE).to_string(),
        "5+?"
    );
    assert_eq!(
        RnaPos::with_offset(5, OFFSET_UNKNOWN_NEGATIVE).to_string(),
        "5-?"
    );
}

/// The two axes that already rendered the sentinel keep doing so — this is the
/// form the fix is being made consistent *with*, not a new decision.
#[test]
fn cds_pos_still_renders_the_unknown_offset_sentinels_as_question_marks() {
    assert_eq!(
        CdsPos::with_offset(5, OFFSET_UNKNOWN_POSITIVE).to_string(),
        "5+?"
    );
    assert_eq!(
        CdsPos::with_offset(5, OFFSET_UNKNOWN_NEGATIVE).to_string(),
        "5-?"
    );
}

/// End to end, through the public entry point: the sentinel is reachable from a
/// parsed description, so the leak was a round-trip failure and not a latent
/// constructor-only shape.
///
/// All **four** moving shapes are here. `r.<b>-?` is the one the first pass
/// missed: `parse_offset` is shared by every axis, so the `r.` parser accepts
/// `-?` exactly as it accepts `+?`, and `RnaPos::Display` had no arm for either.
#[test]
fn a_parsed_unknown_offset_round_trips_on_the_n_and_r_axes() {
    for input in [
        "NM_TEST.1:n.5+?del",
        "NM_TEST.1:n.5-?del",
        "NM_TEST.1:r.5+?del",
        "NM_TEST.1:r.5-?del",
    ] {
        let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
        assert_eq!(
            parsed.to_string(),
            input,
            "an unknown offset must render back as `?`, not as the i64 sentinel"
        );
    }
}

// ---------------------------------------------------------------------------
// The FIFTH `Display` that has to ask the same question
// ---------------------------------------------------------------------------

/// `IvsPos` renders an offset too, and was the arm the four-axis sweep above
/// missed.
///
/// It is reachable from a parsed description through public API and nothing
/// else: [`IvsNotation::to_ivs`] maps a `CdsPos`/`TxPos` offset **straight in**
/// (`self.offset.map(|offset| IvsPos::new(intron_number, offset))`), so a
/// parsed `c.100+?` arrives here carrying the sentinel and printed as
/// `IVS3+9223372036854775807` — a 19-digit intronic distance no reader could
/// tell from a measured one.
///
/// Asserted as a *divergence* from the sentinel's numeric spelling, not only as
/// equality with `IVS3+?`, so the guard still fails if a future arm renders
/// some other digits.
#[test]
fn ivs_pos_renders_the_unknown_offset_sentinels_as_question_marks() {
    let positive = IvsPos::new(3, OFFSET_UNKNOWN_POSITIVE).to_string();
    let negative = IvsPos::new(3, OFFSET_UNKNOWN_NEGATIVE).to_string();

    assert_eq!(positive, "IVS3+?");
    assert_eq!(negative, "IVS3-?");
    for rendered in [&positive, &negative] {
        assert!(
            !rendered.contains(|c: char| c.is_ascii_digit() && c != '3'),
            "the sentinel must not escape as an integer, got: {rendered}"
        );
    }

    // A measured offset is untouched — this is a marker guard, not a rewrite of
    // how distances print.
    assert_eq!(IvsPos::new(3, 9).to_string(), "IVS3+9");
    assert_eq!(IvsPos::new(3, -9).to_string(), "IVS3-9");
}

/// The reachable route, rather than the constructor: a parsed `c.100+?` /
/// `c.100-?` converted to IVS notation must not leak the sentinel either.
#[test]
fn a_parsed_unknown_offset_converted_to_ivs_does_not_leak_the_sentinel() {
    for (input, expected) in [
        ("NM_TEST.1:c.100+?del", "IVS7+?"),
        ("NM_TEST.1:c.100-?del", "IVS7-?"),
    ] {
        let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
        let HgvsVariant::Cds(v) = &parsed else {
            panic!("{input} must parse as a c. variant");
        };
        let pos = v
            .loc_edit
            .location
            .start
            .inner()
            .expect("a point c. position");
        let ivs = pos
            .to_ivs(7)
            .expect("an intronic position has IVS notation");
        assert_eq!(
            ivs.to_string(),
            expected,
            "an unknown offset carried into IVS notation must render as `?`, \
             not as the i64 sentinel"
        );
    }
}
