//! A repeat's tract-wide span must not swallow a cis sibling.
//!
//! `deletion_to_repeat` (`repeated.md` B2) re-expresses a deletion inside a
//! tandem tract as a copy count over the **whole tract**, not over the deleted
//! bases. On a nine-`T` run, `g.1_2del` becomes `g.1_9T[7]` — correct for a lone
//! variant, wrong the moment a cis sibling lives inside that tract:
//!
//! ```text
//! reference   TTTTTTTTTAATATATTTTA        (positions 1-9 are T)
//! input       TEMPLATE:g.[1_2del;4del]
//! was      -> TEMPLATE:g.[1_9T[7];9del]   `9del` sits inside `1_9`
//! now      -> TEMPLATE:g.1_9T[6]
//! ```
//!
//! The conversion runs per member, deep inside `normalize_na_edit`, with no
//! sibling context reaching it, so the widened span is undone in the allele loop
//! instead: the repeat is re-spelled as the equivalent 3'-most deletion, the
//! sibling-crossing clamp pulls it back, and the merge finishes the job. Repeat
//! notation is restored once the members have coalesced and there is no sibling
//! left to span.
//!
//! This is the residual left by #1254 and #1234, which fixed the *shift* into a
//! sibling but not the *span expansion* over one — `NaEdit::Repeat` claims no
//! reference bases as far as the clamp is concerned, so it skipped these. It is
//! pre-existing on `main` and independent of both.

use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{
    parse_hgvs, HgvsVariant, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection,
};

/// Nine `T` at positions 1-9, then `AA` breaking the run.
const TRACT: &str = "TTTTTTTTTAATATATTTTA";

fn provider(seq: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("TEMPLATE", seq.to_string());
    provider
}

fn normalize(seq: &str, input: &str) -> String {
    normalize_in(seq, input, ShuffleDirection::ThreePrime)
}

fn normalize_in(seq: &str, input: &str, direction: ShuffleDirection) -> String {
    let normalizer = Normalizer::with_config(
        provider(seq),
        NormalizeConfig::default().with_direction(direction),
    );
    let variant = parse_hgvs(input).expect("parse");
    format!("{}", normalizer.normalize(&variant).expect("normalize"))
}

/// Apply `descriptor` to `seq` **independently of the normalizer**, via SPDI
/// triples. Returns `None` for an unconvertible or overlapping description —
/// an overlap has no single well-defined resulting sequence, so it must fail a
/// test rather than be double-spliced into a sequence neither side denotes.
fn apply(seq: &str, descriptor: &str) -> Option<String> {
    let provider = provider(seq);
    let members: Vec<HgvsVariant> = match parse_hgvs(descriptor).ok()? {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    let mut triples = Vec::with_capacity(members.len());
    for member in &members {
        triples.push(hgvs_to_spdi(member, &provider).ok()?);
    }
    triples.sort_by_key(|t| std::cmp::Reverse(t.position));
    let reference = seq.as_bytes();
    let mut edited = reference.to_vec();
    let mut claimed_from = reference.len();
    for triple in &triples {
        let start = usize::try_from(triple.position).ok()?;
        let end = start.checked_add(triple.deletion.len())?;
        if end > reference.len() || end > claimed_from {
            return None;
        }
        if !reference[start..end].eq_ignore_ascii_case(triple.deletion.as_bytes()) {
            return None;
        }
        edited.splice(start..end, triple.insertion.bytes());
        claimed_from = start;
    }
    String::from_utf8(edited).ok()
}

/// Assert `input` normalizes to `expected` and that both denote one sequence.
fn assert_normalizes_preserving(seq: &str, input: &str, expected: &str) {
    let actual = normalize(seq, input);
    assert_eq!(actual, expected, "normalizing {input}");
    let from_input = apply(seq, input).expect("input applies");
    let from_output = apply(seq, &actual)
        .unwrap_or_else(|| panic!("{actual} has no well-defined resulting sequence"));
    assert_eq!(
        from_output, from_input,
        "{input} -> {actual} changed the sequence"
    );
}

#[test]
fn two_deletions_in_one_tract_reach_a_single_repeat() {
    // The reported case. Three bases leave the nine-`T` tract, so the whole
    // allele is one repeat of six copies — not a seven-copy repeat with a
    // stray deletion inside it.
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1_2del;4del]", "TEMPLATE:g.1_9T[6]");
}

#[test]
fn three_deletions_in_one_tract_reach_the_same_repeat() {
    // Same total, spelled as three single-base members.
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1del;5del;9del]", "TEMPLATE:g.1_9T[6]");
}

#[test]
fn a_substitution_inside_the_tract_blocks_the_repeat_form() {
    // A substituted base inside the tract cannot be described by a copy count,
    // so the deletion stays a deletion, clamps against the substitution, and
    // the two coalesce.
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1_2del;4T>A]", "TEMPLATE:g.2_4delinsA");
}

#[test]
fn a_substitution_at_the_tract_end_blocks_it_too() {
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1_2del;9T>A]", "TEMPLATE:g.7_9delinsA");
}

#[test]
fn the_repeat_and_deletion_spellings_reach_one_canonical_form() {
    // Confluence. `[1_2del;4T>A]` used to reach the overlapping
    // `[1_9T[7];4T>A]` while `[2_3del;4T>A]` — the same variant, already at the
    // clamp position — reached `2_4delinsA`. Both now land on the latter.
    for input in ["TEMPLATE:g.[1_2del;4T>A]", "TEMPLATE:g.[2_3del;4T>A]"] {
        assert_normalizes_preserving(TRACT, input, "TEMPLATE:g.2_4delinsA");
    }
}

#[test]
fn a_lone_deletion_still_becomes_a_repeat() {
    // Negative control: with no sibling there is nothing to span, so B2 applies
    // exactly as before. Guards against "fixing" this by disabling the repeat
    // form for deletions.
    assert_eq!(normalize(TRACT, "TEMPLATE:g.1_2del"), "TEMPLATE:g.1_9T[7]");
    assert_eq!(normalize(TRACT, "TEMPLATE:g.1_3del"), "TEMPLATE:g.1_9T[6]");
}

#[test]
fn a_sibling_outside_the_tract_leaves_the_repeat_alone() {
    // Negative control: `10A>G` is past the tract's last base (9), so the
    // repeat span claims nothing the sibling claims and must survive.
    assert_normalizes_preserving(
        TRACT,
        "TEMPLATE:g.[1_2del;10A>G]",
        "TEMPLATE:g.[1_9T[7];10A>G]",
    );
}

#[test]
fn a_short_tract_with_a_clear_sibling_keeps_its_repeat() {
    // Negative control on a four-base tract (`T` at 4-7): `5_6del` becomes
    // `4_7T[2]`, and the substitution at 8 is outside it.
    assert_normalizes_preserving(
        "ACGTTTTACGTACGTACGTA",
        "TEMPLATE:g.[5_6del;8A>G]",
        "TEMPLATE:g.[4_7T[2];8A>G]",
    );
}

#[test]
fn a_five_prime_shuffle_is_unaffected() {
    // B2 is gated on 3' shuffle, so under FivePrime no repeat is produced and
    // this pass has nothing to undo.
    let out = normalize_in(
        TRACT,
        "TEMPLATE:g.[1_2del;4del]",
        ShuffleDirection::FivePrime,
    );
    assert!(
        !out.contains('['),
        "5' shuffle should not produce repeat notation, got {out}"
    );
    assert_eq!(
        apply(TRACT, &out).expect("output applies"),
        apply(TRACT, "TEMPLATE:g.[1_2del;4del]").expect("input applies"),
    );
}

#[test]
fn no_two_member_allele_normalizes_to_overlapping_members() {
    // The class. Every two-member cis allele of a deletion plus a downstream
    // deletion or substitution, over deterministic 20-mers, in both shuffle
    // directions — checked for overlapping output *and* for sequence
    // preservation against an independently applied reference.
    //
    // Overlapping outputs across this sweep: 2,578 before the sibling-crossing
    // clamp, 526 after it (all of them this repeat-span class), 0 now.
    let mut checked = 0usize;
    let mut overlapping: Vec<String> = Vec::new();
    let mut changed: Vec<String> = Vec::new();

    for seed in 0..64u32 {
        for alphabet in [b"AT".as_slice(), b"ACGT".as_slice()] {
            let mut state = u64::from(seed).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
            let seq: String = (0..20)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    alphabet[(state % alphabet.len() as u64) as usize] as char
                })
                .collect();

            for first_start in 1..=14usize {
                for first_len in 1..=2usize {
                    let first_end = first_start + first_len - 1;
                    let first = if first_len == 1 {
                        format!("{first_start}del")
                    } else {
                        format!("{first_start}_{first_end}del")
                    };
                    for second_start in first_end + 2..=19usize {
                        let base = seq.as_bytes()[second_start - 1] as char;
                        let alt = if base == 'A' { 'G' } else { 'A' };
                        for second in [
                            format!("{second_start}del"),
                            format!("{second_start}{base}>{alt}"),
                        ] {
                            for direction in
                                [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime]
                            {
                                let input = format!("TEMPLATE:g.[{first};{second}]");
                                let output = normalize_in(&seq, &input, direction);
                                checked += 1;
                                let Some(want) = apply(&seq, &input) else {
                                    continue; // unconvertible input
                                };
                                match apply(&seq, &output) {
                                    // `apply` declines an overlapping output, so
                                    // `None` here is the overlap signal.
                                    None if overlapping.len() < 10 => {
                                        overlapping.push(format!("{seq}: {input} -> {output}"));
                                    }
                                    None => {}
                                    Some(got) if got != want && changed.len() < 10 => {
                                        changed.push(format!(
                                            "{seq}: {input} -> {output} (want {want}, got {got})"
                                        ));
                                    }
                                    Some(_) => {}
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    assert!(checked > 100_000, "sweep covered too little: {checked}");
    assert!(
        overlapping.is_empty(),
        "overlapping or unconvertible output in {checked} cases: {overlapping:#?}"
    );
    assert!(
        changed.is_empty(),
        "sequence-changing normalizations in {checked} cases: {changed:#?}"
    );
}
