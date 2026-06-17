//! Issue #180 / #181: 3' shifting of del/dup/ins produced by allele-bracket merge.
//!
//! When consecutive edits inside an allele bracket merge into a `del`, `dup`,
//! or `ins`, the HGVS 3' rule still applies to the merged form. Before this
//! fix, only merged pure deletions were re-shifted; merged insertions (and
//! merged delins that canonicalize into a shiftable form) skipped the shift.
//!
//! The pre-fix code also ran the per-element 3' shift *before* merge, which
//! could collapse adjacent edits onto the same shifted position and emit
//! outputs with duplicate or overlapping positions — invalid HGVS. Issue
//! #181 documented that failure mode; the fix here (merge first, then run
//! the per-variant pipeline on every merged result) resolves both #180 and
//! #181 with one ordering change.
//!
//! Examples below use synthetic genomic providers via `SyntheticBuilder`,
//! which pads the core with 256 bases so HGVS pos 257 maps to core base 1.

use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// HGVS sequence-id used by `SyntheticBuilder::genomic`.
const SEQID: &str = "NC_TEST.1";

/// Core from issue #180's deletion examples (positions 1337..1351 in the
/// issue map to HGVS pos 257..271 here): `CCCAAAAATAAACGC` — a 5-A
/// homopolymer at HGVS 260..264, a T at 265, then AAA at 266..268.
const CORE_DEL: &str = "CCCAAAAATAAACGC";

/// Core from issue #180's duplication example (positions 1048..1065 in the
/// issue map to HGVS pos 257..274 here): `ACCTATGATCTGAAGAGC`. The inserted
/// `TGATC` between 261 and 262 matches the reference at 262..266, then
/// shifts right through TGATC|TGAAG until 265..269.
const CORE_DUP: &str = "ACCTATGATCTGAAGAGC";

// ---------------------------------------------------------------------------
// Deletions
// ---------------------------------------------------------------------------

#[test]
fn allele_two_single_dels_shift_to_3prime_end_of_homopolymer() {
    // Issue example 1: [260delA;261delA] merges to 260_261del. The 3' shift
    // through the A-run then composes with the homopolymer-repeat
    // canonicalization (matching del_shift_matrix's "k>1 in homopolymer →
    // repeat" rule) to produce 260_264A[3] — the 5-A tract at HGVS 260..264
    // becomes A[3] (3 A's remaining after deletion of 2).
    let p = SyntheticBuilder::genomic(CORE_DEL).build();
    let result = normalize_to_string(p, &format!("{}:g.[260delA;261delA]", SEQID));
    assert_eq!(result, format!("{}:g.260_264A[3]", SEQID));
}

#[test]
fn allele_seven_single_dels_spanning_run_boundary_shift_3prime() {
    // Issue example 2: 7-base deletion spanning the A-run, the T, and one A
    // of the AAA run. Merges to 260_266del, then 3'-shifts to 262_268del.
    let p = SyntheticBuilder::genomic(CORE_DEL).build();
    let result = normalize_to_string(
        p,
        &format!(
            "{}:g.[260delA;261delA;262delA;263delA;264delA;265delT;266delA]",
            SEQID
        ),
    );
    assert_eq!(result, format!("{}:g.262_268del", SEQID));
}

#[test]
fn allele_single_del_in_bracket_shifts_to_3prime_end_of_run() {
    // Issue example 3: a singleton allele [260delA] currently emits 260del.
    // Per the 3' rule it must shift to 264del (last A of the 5-A run).
    let p = SyntheticBuilder::genomic(CORE_DEL).build();
    let result = normalize_to_string(p, &format!("{}:g.[260delA]", SEQID));
    assert_eq!(result, format!("{}:g.264del", SEQID));
}

// ---------------------------------------------------------------------------
// Duplications (canonicalized from insertion)
// ---------------------------------------------------------------------------

#[test]
fn allele_ins_canonicalized_to_dup_then_shifted_3prime() {
    // Issue dup example: [261_262insTGATC] canonicalizes to 262_266dup, then
    // 3'-shifts through TGATC|TGAAG until ref[265]=T != ref[270]=A, leaving
    // 265_269dup.
    let p = SyntheticBuilder::genomic(CORE_DUP).build();
    let result = normalize_to_string(p, &format!("{}:g.[261_262insTGATC]", SEQID));
    assert_eq!(result, format!("{}:g.265_269dup", SEQID));
}

// ---------------------------------------------------------------------------
// Issue #181: no duplicate / overlapping positions in cis output.
//
// Same root cause as #180 (per-element 3' shift running before merge), but
// surfaces as *invalid HGVS output* — duplicate positions, overlapping
// ranges — rather than as a missed canonical form. The patterns below
// reproduce #181's three examples at the literal positions cited there.
// Cores are padded so HGVS pos N maps to core idx (N - 257); we use
// non-homopolymer filler (alternating CGCT) so the only shiftable run is
// the homopolymer we plant at the target positions.
// ---------------------------------------------------------------------------

/// Build a core where `snippet` lands at the requested HGVS position.
///
/// Prefix and trailing filler use a non-shiftable alternating pattern
/// (CGCT) so the only run available to the 3' shift is the snippet itself.
/// `hgvs_pos` is the 1-based HGVS coordinate where the first base of
/// `snippet` should sit; the helper computes the prefix length as
/// `hgvs_pos - 257` (since PAD_OFFSET=256 places core base 1 at HGVS 257).
fn core_with_snippet_at(hgvs_pos: u64, snippet: &str, trail: usize) -> String {
    let prefix_len = (hgvs_pos - 257) as usize;
    let mut s = String::with_capacity(prefix_len + snippet.len() + trail);
    for i in 0..prefix_len {
        s.push(['C', 'G', 'C', 'T'][i % 4]);
    }
    s.push_str(snippet);
    for i in 0..trail {
        s.push(['C', 'G', 'C', 'T'][i % 4]);
    }
    s
}

/// Assert that no sub-variant repeats or overlaps another in a cis allele
/// rendering. Catches the #181 anti-pattern (`[1220del;1220del]`,
/// `[1054_1059del;1061del;1061del]`) regardless of which canonical form
/// the normalizer happens to emit.
fn assert_no_duplicate_or_overlap(rendered: &str) {
    if !rendered.contains(";") {
        return; // not a multi-element allele
    }
    // Extract the body between the first '[' and the last ']'.
    let Some(open) = rendered.find('[') else {
        return;
    };
    let Some(close) = rendered.rfind(']') else {
        return;
    };
    let body = &rendered[open + 1..close];
    // Parse `<start>(_<end>)?<op>...` from each `;`-separated element and
    // collect inclusive position ranges. Anything that doesn't parse as a
    // numeric range (e.g. `delins` strings) is ignored — sub-variants in
    // a cis allele always lead with a positional anchor.
    let mut ranges: Vec<(u64, u64)> = Vec::new();
    for elt in body.split(';') {
        let digits: String = elt.chars().take_while(|c| c.is_ascii_digit()).collect();
        if digits.is_empty() {
            continue;
        }
        let start: u64 = digits.parse().unwrap();
        let rest = &elt[digits.len()..];
        let end = if let Some(stripped) = rest.strip_prefix('_') {
            let d2: String = stripped
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .collect();
            d2.parse().unwrap_or(start)
        } else {
            start
        };
        ranges.push((start, end));
    }
    ranges.sort();
    for w in ranges.windows(2) {
        let (_, end_a) = w[0];
        let (start_b, _) = w[1];
        assert!(
            end_a < start_b,
            "overlap or duplicate in cis allele output: {} (ranges {:?})",
            rendered,
            ranges
        );
    }
}

#[test]
fn issue_181_example_1_adjacent_dels_no_duplicate_position() {
    // Before the fix, ferro emitted `g.[1220del;1220del]` — same position
    // twice, invalid HGVS. With merge-first ordering the two dels collapse
    // to `1216_1217del` over the A-run at 1216..1220 and then 3'-shift,
    // canonicalizing into the repeat form `1216_1220A[3]` (5 A's → 3 A's
    // after deleting 2). The validity invariant: no sub-variant repeats.
    let core = core_with_snippet_at(1216, "AAAAA", 100);
    let p = SyntheticBuilder::genomic(&core).build();
    let result = normalize_to_string(p, &format!("{}:g.[1216del;1217del]", SEQID));
    assert_eq!(result, format!("{}:g.1216_1220A[3]", SEQID));
    assert_no_duplicate_or_overlap(&result);
}

#[test]
fn issue_181_example_2_eight_adjacent_dels_no_overlap() {
    // Before the fix, ferro emitted `g.[1054_1059del;1061del;1061del]` —
    // both a duplicate (`1061del` twice) and overlap between the merged
    // 1054_1059 range and the subsequent 1061del. With merge-first
    // ordering all 8 dels collapse into `1054_1061del` over the 10-A run
    // and 3'-shift / repeat-canonicalize into `1054_1063A[2]`.
    let core = core_with_snippet_at(1054, "AAAAAAAAAA", 100);
    let p = SyntheticBuilder::genomic(&core).build();
    let result = normalize_to_string(
        p,
        &format!(
            "{}:g.[1054del;1055del;1056del;1057del;1058del;1059del;1060del;1061del]",
            SEQID
        ),
    );
    assert_eq!(result, format!("{}:g.1054_1063A[2]", SEQID));
    assert_no_duplicate_or_overlap(&result);
}

#[test]
fn issue_181_example_3_mixed_dels_and_ins_no_duplicate_position() {
    // Before the fix, ferro emitted `g.[1006del;1008del;1008del;1010C[4]]`
    // — `1008del` twice. With merge-first ordering the three adjacent
    // dels collapse to `1006_1008del`; the trailing `1009_1010insCCC`
    // does not fold into the del because the merge step's strict
    // adjacency rule (`prev.end + 1 == next.start`) requires the
    // insertion's anchored start (`p+1 = 1010`) to equal `prev.end + 1`,
    // not `prev.end + 2`. Fully collapsing del+ins chains into a single
    // delins is a separate improvement; here we pin the validity
    // invariant — no duplicate or overlapping positions — which is what
    // #181 was actually filed about.
    let core = core_with_snippet_at(1006, "AGGTACGT", 100);
    let p = SyntheticBuilder::genomic(&core).build();
    let result = normalize_to_string(
        p,
        &format!("{}:g.[1006delA;1007delG;1008delG;1009_1010insCCC]", SEQID),
    );
    assert_eq!(
        result,
        format!("{}:g.[1006_1008del;1009_1010insCCC]", SEQID)
    );
    assert_no_duplicate_or_overlap(&result);
}

// ---------------------------------------------------------------------------
// Idempotency: a second normalization round must be a no-op.
// ---------------------------------------------------------------------------

fn normalize_twice(provider: MockProvider, input: &str) -> (String, String) {
    let normalizer = Normalizer::new(provider);
    let parsed = parse_hgvs(input).expect("parse");
    let once = normalizer.normalize(&parsed).expect("normalize 1");
    let twice = normalizer.normalize(&once).expect("normalize 2");
    (once.to_string(), twice.to_string())
}

#[test]
fn allele_norm_3prime_shift_is_idempotent() {
    let cases: &[(String, String)] = &[
        (CORE_DEL.to_string(), "g.[260delA;261delA]".to_string()),
        (
            CORE_DEL.to_string(),
            "g.[260delA;261delA;262delA;263delA;264delA;265delT;266delA]".to_string(),
        ),
        (CORE_DEL.to_string(), "g.[260delA]".to_string()),
        (CORE_DUP.to_string(), "g.[261_262insTGATC]".to_string()),
        (
            core_with_snippet_at(1216, "AAAAA", 100),
            "g.[1216del;1217del]".to_string(),
        ),
        (
            core_with_snippet_at(1054, "AAAAAAAAAA", 100),
            "g.[1054del;1055del;1056del;1057del;1058del;1059del;1060del;1061del]".to_string(),
        ),
        (
            core_with_snippet_at(1006, "AGGTACGT", 100),
            "g.[1006delA;1007delG;1008delG;1009_1010insCCC]".to_string(),
        ),
    ];
    for (core, tail) in cases {
        let p = SyntheticBuilder::genomic(core).build();
        let input = format!("{}:{}", SEQID, tail);
        let (once, twice) = normalize_twice(p, &input);
        assert_eq!(once, twice, "not idempotent for {}", input);
    }
}
