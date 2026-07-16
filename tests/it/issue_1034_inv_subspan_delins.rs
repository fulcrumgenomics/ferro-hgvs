//! Issue #1034: ferro must not fragment a contiguous `delins` to expose an
//! internal reverse-complement sub-run as an `inv`.
//!
//! Per HGVS an inversion describes a *maximal contiguous run* whose alt equals
//! the reverse complement of the reference. A reverse-complement **sub-run**
//! may not be carved out of a longer contiguous change — that change must stay
//! a single `delins` (`DNA/inversion.md`, `DNA/delins.md`,
//! `general.md:34` + `:56`).
//!
//! v0.8.0 regressed here: the `#160`/`#166` sub-span inv typer matched
//! reverse-complement identity on sub-runs of a contiguous mismatch run rather
//! than requiring the whole run to be a reverse complement, so a contiguous
//! `delins` fragmented into `[…; inv; …]`. This suite pins the spec-correct
//! single-`delins` output while confirming genuine full-run inversions (whole
//! maximal run is a reverse complement) still normalize to `inv`.

use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};

const SEQID: &str = "NC_TEST.1";

// ---------------------------------------------------------------------------
// The exact minimal example from issue #1034.
// ---------------------------------------------------------------------------

#[test]
fn contiguous_delins_with_internal_revcomp_subrun_stays_delins() {
    // Core "CTG" at 257-259 (ref 257C 258T 259G).
    //   raw:  [257C>A; 258T>C; 259G>A]
    //   window 257_259: ref CTG -> alt ACA (length-preserving).
    //   revcomp(CTG) = CAG != ACA  =>  the whole run is NOT an inversion.
    //   The 2-nt sub-run 258_259 (TG -> CA) satisfies revcomp(TG)=CA, but it
    //   may NOT be carved out of the contiguous change.
    // Expected: a single delins, NOT `[257C>A; 258_259inv]`.
    let p = SyntheticBuilder::genomic("CTG").build();
    let result = normalize_to_string(p, &format!("{}:g.[257C>A;258T>C;259G>A]", SEQID));
    assert_eq!(result, format!("{}:g.257_259delinsACA", SEQID));
}

#[test]
fn contiguous_delins_with_leading_revcomp_subrun_stays_delins() {
    // Core "TCC" at 257-259.
    //   raw:  [257T>G; 258C>A; 259C>G]  ->  window ref TCC / alt GAG.
    //   revcomp(TCC) = GGA != GAG  =>  not a full inversion.
    //   The 2-nt prefix 257_258 (TC -> GA) satisfies revcomp(TC)=GA but must
    //   not be split out. This is the exact `#160` CHANGELOG example, now
    //   spec-corrected.
    let p = SyntheticBuilder::genomic("TCC").build();
    let result = normalize_to_string(p, &format!("{}:g.[257T>G;258C>A;259C>G]", SEQID));
    assert_eq!(result, format!("{}:g.257_259delinsGAG", SEQID));
}

// ---------------------------------------------------------------------------
// Guard: a genuine full-run inversion (whole maximal run is a reverse
// complement) must STILL normalize to `inv`.
// ---------------------------------------------------------------------------

#[test]
fn full_run_reverse_complement_still_emits_inv() {
    // Core "GCT" at 257-259.
    //   raw:  [257G>A; 258C>G; 259T>C]  ->  window ref GCT / alt AGC.
    //   revcomp(GCT) = AGC == alt  =>  the whole run IS an inversion.
    let p = SyntheticBuilder::genomic("GCT").build();
    let result = normalize_to_string(p, &format!("{}:g.[257G>A;258C>G;259T>C]", SEQID));
    assert_eq!(result, format!("{}:g.257_259inv", SEQID));
}

#[test]
fn two_base_full_run_reverse_complement_still_emits_inv() {
    // Core "GG" at 257-258.
    //   raw:  [257G>C; 258G>C]  ->  window ref GG / alt CC.
    //   revcomp(GG) = CC == alt  =>  the whole run IS an inversion.
    let p = SyntheticBuilder::genomic("GG").build();
    let result = normalize_to_string(p, &format!("{}:g.[257G>C;258G>C]", SEQID));
    assert_eq!(result, format!("{}:g.257_258inv", SEQID));
}
