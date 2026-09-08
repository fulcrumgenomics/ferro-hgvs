//! A clean single-copy tandem duplication is never re-partitioned into
//! `[sub; dup]` (`pure-tandem-expansion-insertion-is-not-repartitioned`).
//!
//! `DNA/duplication.md:17-18`: a change that is a copy of the sequence directly
//! 5' of it must be described as a duplication. When the resulting sequence's
//! minimal derivation is already one clean dup, the ruled partitioner's Part-B
//! peel ranks its placements by the substitutions each would invent, so the
//! zero-substitution placement (the dup itself) wins. It used to keep only
//! one-mismatch placements and return the 3'-most, which minted a substitution
//! on a base the variant never changes and moved the dup along the tract.
//!
//! Each case is a real reference window, taken from the prepared reference and
//! padded so `window[i]` sits at `TEMPLATE:g.(257 + i)`. The windows reproduce
//! corpus rows that moved under the ruled partitioner:
//!
//! - `NG_034246.1:g.29133_29134dup` -> `g.[29135A>T;29184_29185dup]`
//!   (`CA` tract, the dup 49 nt away);
//! - `NM_002247.3:c.28_30dup` -> `c.[31A>G;34_36dup]` (`GGC` tract), here on
//!   the genomic axis over the same bases.
//!
//! The negative control is an insertion with NO clean placement, which keeps its
//! `ins -> [dup; sub]` re-derivation.

use crate::common::cis_apply_oracle::normalize;

/// `ACGT` x 64 on both sides, the padding `duplication_label_not_partition`
/// uses, so `window[i]` is at 1-based `257 + i`.
fn padded(window: &str) -> String {
    format!("{}{}{}", "ACGT".repeat(64), window, "ACGT".repeat(64))
}

/// `NG_034246.1` `g.29081`-`g.29220`. `g.29133_29134` is `TEMPLATE:g.309_310`.
const NG_034246_WINDOW: &str = "ATCTAGAAGAAAAAATAGGTTCCTAGAGGTAAACTGCACCTAATCTTCCCAGTCACACACACACACACACACACACACACACACACACACACACACACACACACATCCTTCACTAGAACACACTGCCCCACTGCCAAGAC";

/// `NM_002247.3` `c.-20`-`c.100`. `c.28_30` is `TEMPLATE:g.304_306`.
const NM_002247_WINDOW: &str = "CCGGCTGCCCGTTGCTAGCTATGGCAAATGGTGGCGGCGGCGGCGGCGGCAGCAGCGGCGGCGGCGGCGGCGGCGGAGGCAGCAGTCTTAGAATGAGTAGCAATATCCACGCGAACCATC";

/// `NC_000002.11` `g.73613031`-`g.73613090`. `g.73613061_73613062` is
/// `TEMPLATE:g.287_288`.
const NC_000002_WINDOW: &str = "TGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAAGAGGAGGAGGCTGCAGCGG";

#[test]
fn ng_034246_clean_dup_stays_one_dup() {
    let seq = padded(NG_034246_WINDOW);
    assert_eq!(
        normalize(&seq, "TEMPLATE:g.309_310dup"),
        "TEMPLATE:g.309_310dup",
        "a clean single-copy dup must not become [sub; dup] (duplication.md:17-18)",
    );
}

#[test]
fn nm_002247_clean_codon_dup_stays_one_dup() {
    let seq = padded(NM_002247_WINDOW);
    assert_eq!(
        normalize(&seq, "TEMPLATE:g.304_306dup"),
        "TEMPLATE:g.304_306dup",
        "a clean single-copy dup must not become [sub; dup] (duplication.md:17-18)",
    );
}

/// Every spelling of a variant converges (`canonical-form-choice-when-both-legal`
/// derives the form from the resulting sequence). The `[sub; dup]` spelling the
/// defect emitted denotes the same bases as the clean dup, so it must normalize
/// to the clean dup too, never be kept because it is the input's spelling.
#[test]
fn the_sub_dup_spelling_converges_to_the_clean_dup() {
    let ng = padded(NG_034246_WINDOW);
    assert_eq!(
        normalize(&ng, "TEMPLATE:g.[311A>T;360_361dup]"),
        "TEMPLATE:g.309_310dup",
    );
    let nm = padded(NM_002247_WINDOW);
    assert_eq!(
        normalize(&nm, "TEMPLATE:g.[307A>G;310_312dup]"),
        "TEMPLATE:g.304_306dup",
    );
}

/// Negative control. `NC_000002.11:g.73613061_73613062insAGAGGA` has no
/// zero-substitution placement, so the fewest-substitutions ranking still picks
/// the 3'-most one-substitution placement and the ruled partitioner keeps its
/// `g.[73613050_73613055dup;73613056G>A]` re-derivation.
#[test]
fn an_insertion_with_no_clean_placement_still_becomes_dup_plus_sub() {
    let seq = padded(NC_000002_WINDOW);
    assert_eq!(
        normalize(&seq, "TEMPLATE:g.287_288insAGAGGA"),
        "TEMPLATE:g.[276_281dup;282G>A]",
    );
}
