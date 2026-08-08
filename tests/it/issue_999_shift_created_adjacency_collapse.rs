//! Regression coverage for issues #999 and #1000 — one root cause.
//!
//! When a cis allele member is an insertion that 3'-shifts to become
//! directly adjacent (0-nt gap) to a downstream substitution, the two
//! must coalesce into a single `delins` (HGVS: adjacent changes with no
//! intervening nucleotide are one edit). ferro left them as two bracketed
//! members because the cross-member collapse
//! (`merge::collapse_overlapping_cis_edits`) ran on the *raw input*
//! members, before the per-member 3'-shift that creates the adjacency, and
//! was never re-run afterwards.
//!
//! Two observable symptoms, same defect:
//!
//!   #999  — `NC_000001.11:g.[305_306insC;307G>T]`
//!           was: `g.[306dup;307G>T]`   want: `g.307delinsCT`
//!
//!   #1000 — `NC_000001.11:g.[611_612insCG;613C>A]` was non-idempotent:
//!           pass 1 gave `g.[612_613insGC;613C>A]`, pass 2 gave
//!           `g.613delinsGCA`. A single pass must reach the fixed point.
//!
//! These pins use synthetic MockProvider sequences so future refactors of
//! the merge/collapse/shift pipeline can't silently regress the rule.

use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// Build a MockProvider whose genomic contig `NC_000001.11` carries
/// `substantive` starting at 1-based `start`, padded to 1000 bases of
/// filler `A` so the default normalize window always lands inside.
fn provider_with(start: usize, substantive: &str) -> MockProvider {
    let mut seq = vec![b'A'; 1000];
    for (i, b) in substantive.bytes().enumerate() {
        seq[start - 1 + i] = b;
    }
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_000001.11", String::from_utf8(seq).unwrap());
    provider
}

fn normalize(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    let variant = parse_hgvs(input).expect("parse");
    format!("{}", normalizer.normalize(&variant).expect("normalize"))
}

#[test]
fn cis_insertion_shifted_adjacent_to_substitution_coalesces_to_delins() {
    // #999. The inserted `C` matches ref[306]=C, so it 3'-shifts to a dup
    // at gap 306, landing directly 5' of the `307G>T` substitution. With no
    // intervening nucleotide the two describe one contiguous change: the
    // inserted `C` then the substituted `T` at 307 -> `307delinsCT`.
    let out = normalize(
        provider_with(300, "CATCCTCGCTCCT"),
        "NC_000001.11:g.[305_306insC;307G>T]",
    );
    assert_eq!(out, "NC_000001.11:g.307delinsCT");
}

#[test]
fn cis_shift_created_adjacency_reaches_fixed_point_in_one_pass() {
    // #1000. `insCG` at 611_612 legitimately 3'-shifts to `insGC` at
    // 612_613 (ref[612]=C lets the dinucleotide roll one position; 613
    // blocks a second roll), landing adjacent to `613C>A`.
    //
    // **Re-blessed back to `g.613delinsGCA` under the partition model**, which is
    // the form #1000 itself named as the fixed point ("pass 2 gave
    // `g.613delinsGCA`"; see the module header). The intervening
    // `g.[612_613insG;613_614insA]` came from the two-gap alignment #1260 chose
    // when the partition was re-derived from `C -> GCA`; that re-derivation is
    // what the partition model removes, so the question is decided by the input's
    // own members again.
    //
    // And they merge, unavoidably. The allele asserts two members. After the 3'
    // shift the insertion occupies the `612|613` junction and `613C>A` claims base
    // 613, so **their separation is zero** — no unchanged nucleotide lies between
    // them. `general.md:34` is a rule about members "separated by one or more
    // nucleotides" and does not reach separation 0; the clause that does is
    // `DNA/substitution.md:32`, "changes involving two or more consecutive
    // nucleotides are described as deletion-insertion (delins)", stated there as a
    // prohibition on the two-member spelling. `DNA/delins.md:86-89` is the worked
    // case: `NM_007294.3:c.[2077G>A;2077_2078insTA]` has "The correct
    // description" `c.2077delinsATA`, and `:89` records that the sentence
    // permitting the two-member form was removed by the committee. So an abutting
    // pair merges unconditionally, whatever the input spelled.
    //
    // What #1000 asks for is a **fixed point in one pass**, and that has held
    // through every form this row has carried.
    // `cis_shift_created_adjacency_is_idempotent` below asserts the invariant
    // directly and is unchanged.
    let out = normalize(
        provider_with(600, "TGACTTCAGTCACCTGACTGACTG"),
        "NC_000001.11:g.[611_612insCG;613C>A]",
    );
    assert_eq!(out, "NC_000001.11:g.613delinsGCA");
}

#[test]
fn cis_shift_created_adjacency_is_idempotent() {
    // #1000, stated as the invariant: normalize(normalize(x)) == normalize(x).
    let input = "NC_000001.11:g.[611_612insCG;613C>A]";
    let pass1 = normalize(provider_with(600, "TGACTTCAGTCACCTGACTGACTG"), input);
    let pass2 = normalize(provider_with(600, "TGACTTCAGTCACCTGACTGACTG"), &pass1);
    assert_eq!(pass1, pass2, "normalize must be idempotent; pass1 != pass2");
}

#[test]
fn cis_insertion_shifted_to_non_adjacent_stays_separate() {
    // Negative control: when the shift does NOT bring the insertion flush
    // against the substitution (an unchanged base remains between them),
    // the two must stay as separate members — the collapse is for
    // adjacent/overlapping edits only. Here ins `C` at 305_306 shifts to a
    // dup at 306 (gap 306), but the substitution is at 308, leaving ref[307]
    // untouched between them.
    let out = normalize(
        provider_with(300, "CATCCTCGCTCCT"),
        "NC_000001.11:g.[305_306insC;308C>A]",
    );
    assert_eq!(out, "NC_000001.11:g.[306dup;308C>A]");
}
