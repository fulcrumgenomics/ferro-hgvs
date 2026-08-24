//! #2193 — a lone contiguous run where EVERY base differs (equal-length) must be
//! the single spanning `delins` `DNA/delins.md:16`/`:44-47` recommend, not a
//! `[dup;del]` whose `dup` is anchored 5' *outside* the changed span.
//!
//! `ATGCGAG -> GGTATGC` (7 bases each) used to derive `g.[8_10dup;15_17del]`: a
//! coincidental `GGT` tandem 5' of the run let the minimal-edit alignment park a
//! copy as a boundary insertion, so `coalesce_solid_run` saw a "multi-base tandem
//! dup" and declined to collapse the equal-length hull. Per
//! `duplication-must-ranks-the-label-not-the-partition` (narrowed for #2193 to
//! net-longer duplications only), `DNA/duplication.md:18`'s MUST ranks the
//! *label* of a change already identified, never a manufactured partition, and
//! `general.md:57` forbids the `[del;dup]` shape by name. So the run collapses to
//! one `delins`.
//!
//! The carve-out survives only for a genuine, net-LONGER duplication (#2175); the
//! net-SHORTER direction was probed (contiguous runs beginning with the 5' flank)
//! and already renders as a spanning `delins`, so there is no net-shorter twin to
//! guard against — the consistency cases below pin that.

use ferro_hgvs::{from_sequences, FromSequencesOptions, ShuffleDirection};

fn derive(refe: &str, alt: &str) -> String {
    from_sequences("CONTIG1", 1, refe, alt, &FromSequencesOptions::default())
        .expect("from_sequences must derive this block")
        .to_string()
}

fn derive_5p(refe: &str, alt: &str) -> String {
    from_sequences(
        "CONTIG1",
        1,
        refe,
        alt,
        &FromSequencesOptions::default().with_direction(ShuffleDirection::FivePrime),
    )
    .expect("from_sequences must derive this block")
    .to_string()
}

/// The #2193 reproductions: a lone equal-length all-differing run at `g.11_17`
/// whose alt begins with a copy of the 5' flank `GGT` (`g.8_10`). Each used to
/// render as an outside-anchored `g.[8_10dup;…del]`; each must now be the single
/// spanning `delins`, in both shuffle directions. All three rows are the issue's
/// own reproduction table (verified there on 0.17.0 via `from_sequences`).
#[test]
fn the_equal_length_run_is_a_spanning_delins() {
    // g.11_17 ATGCGAG -> GGTATGC  (was g.[8_10dup;15_17del])
    let (refe, alt) = (
        "ACGTTCAGGTATGCGAGTTAGCTAGCTAG",
        "ACGTTCAGGTGGTATGCTTAGCTAGCTAG",
    );
    assert_eq!(derive(refe, alt), "CONTIG1:g.11_17delinsGGTATGC");
    assert_eq!(derive_5p(refe, alt), "CONTIG1:g.11_17delinsGGTATGC");

    // g.11_17 TCCCGGT -> GGTTCCC  (was g.[8_10dup;15_17del])
    let (refe, alt) = (
        "ACGTTCAGGTTCCCGGTTTAGCTAGCTAG",
        "ACGTTCAGGTGGTTCCCTTAGCTAGCTAG",
    );
    assert_eq!(derive(refe, alt), "CONTIG1:g.11_17delinsGGTTCCC");
    assert_eq!(derive_5p(refe, alt), "CONTIG1:g.11_17delinsGGTTCCC");

    // g.11_17 CAAGTGA -> GGTCAAG  (was g.[8_10dup;16_18del])
    let (refe, alt) = (
        "ACGTTCAGGTCAAGTGATTAGCTAGCTAG",
        "ACGTTCAGGTGGTCAAGTTAGCTAGCTAG",
    );
    assert_eq!(derive(refe, alt), "CONTIG1:g.11_17delinsGGTCAAG");
    assert_eq!(derive_5p(refe, alt), "CONTIG1:g.11_17delinsGGTCAAG");
}

/// The #422 shape re-derived directly from sequences: the equal-length run
/// `CGTACG -> TACGTA` at `g.7_12`. A coincidental `TA` tandem used to spell it
/// `g.[..del;..dup]`; it is now the spanning `delins`. (The cross-reference
/// expansion path pins the same shape in `issue_422_cross_reference_ins`.)
#[test]
fn the_issue_422_shape_is_a_spanning_delins() {
    let refe = "AAAAAACGTACGTTTTTT";
    let alt = "AAAAAATACGTATTTTTT";
    assert_eq!(derive(refe, alt), "CONTIG1:g.7_12delinsTACGTA");
    assert_eq!(derive_5p(refe, alt), "CONTIG1:g.7_12delinsTACGTA");
}

/// The net-LONGER carve-out is untouched: a genuine tandem expansion stays a
/// `dup`, whether isolated or abutting a substitution (#2175). This is what the
/// #2193 fix must NOT regress — the equal-length gate, not the removed guard, is
/// what keeps these out of `coalesce_solid_run`.
#[test]
fn a_net_longer_tandem_expansion_stays_a_dup() {
    let refe = "TAGTAAACCATTTTACGGAGGATCACAAATTCCTCCTTAT";
    let isolated = "TAGTAAACCATTTTACGGAGGATCACACAAATTCCTCCTTAT"; // CACA -> CACACA (+2)
    assert_eq!(derive(refe, isolated), "CONTIG1:g.26_27dup");
    let abutting = "TAGTAAACCATTTTACGGAGGATCACACACATTCCTCCTTAT"; // +2 then g.28A>C
    assert_eq!(derive(refe, abutting), "CONTIG1:g.[26_27dup;28A>C]");
}

/// Net-SHORTER consistency: a contiguous all-differing run whose alt payload
/// BEGINS with a copy of the 5' flank (`GG` / `GGT`) still renders as a single
/// spanning `delins` — the fix does not leak an outside-anchored `dup` into the
/// net-shorter direction. `ATGCGAG` (`g.11_17`, 7 bases) shortened to a 3- or
/// 4-base payload beginning with the flank, plus a flank-free baseline.
#[test]
fn a_net_shorter_run_beginning_with_the_flank_is_a_spanning_delins() {
    let refe = "ACGTTCAGGTATGCGAGTTAGCTAGCTAG";
    // payload GGCA (net -3) begins with the flank copy `GG`
    assert_eq!(
        derive(refe, "ACGTTCAGGTGGCATTAGCTAGCTAG"),
        "CONTIG1:g.11_17delinsGGCA"
    );
    // payload GGT (net -4) is the exact 5' flank
    assert_eq!(
        derive(refe, "ACGTTCAGGTGGTTTAGCTAGCTAG"),
        "CONTIG1:g.11_17delinsGGT"
    );
    // flank-free baseline
    assert_eq!(
        derive(refe, "ACGTTCAGGTCCAGCTTAGCTAGCTAG"),
        "CONTIG1:g.11_17delinsCCAGC"
    );
}
