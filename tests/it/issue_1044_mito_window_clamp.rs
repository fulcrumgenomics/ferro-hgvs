//! Issue #1044 — the mitochondrial / circular (`m.`/`o.`) normalization path
//! had the same unclamped genome-fetch window as the genomic path did before
//! #1041/#1042: a non-wraparound variant within `window_size` (100 bp) of the
//! contig's 3' end fetched past EOF, `get_sequence` errored, and the variant
//! dropped into `mt_fallback` — skipping 3' shifting AND the
//! `delins -> inv/sub/dup` canonicalization. A whole-span reverse-complement
//! `delins` near the 3' end therefore stayed `delins` instead of becoming `inv`.
//!
//! This mirrors `issue_1041_repro.rs` on the `m.` axis. The fix clamps
//! `fetch_end` to the contig length, gated on `end <= len` so a span that runs
//! *past* the contig end still passes through the safe fallback unchanged.

use ferro_hgvs::{parse_hgvs, JsonProvider, Normalizer};
use std::io::Write;

/// Mitochondrial-capable provider: one `contig` of `len` cyclic ACGT bytes with
/// `payload` written 1-based at `pos1` (mirrors the #1041 repro harness).
fn provider_len(contig: &str, len: usize, pos1: usize, payload: &str) -> JsonProvider {
    let mut bases: Vec<u8> = "ACGT".bytes().cycle().take(len).collect();
    for (i, b) in payload.bytes().enumerate() {
        bases[pos1 - 1 + i] = b;
    }
    let seq = String::from_utf8(bases).unwrap();
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { contig: seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

fn norm(p: &JsonProvider, input: &str) -> String {
    Normalizer::new(p.clone())
        .normalize(&parse_hgvs(input).unwrap())
        .unwrap()
        .to_string()
}

#[test]
fn issue_1044_mt_whole_run_revcomp_is_inv_regardless_of_3prime_proximity() {
    // ref[200..202] = GCA, revcomp = TGC. Same variant, two contig lengths.
    let contig = "NC_012920.1";

    // Far from the end (600 bp): window [100, 302] fits comfortably.
    let far = provider_len(contig, 600, 200, "GCA");
    assert_eq!(
        norm(&far, &format!("{contig}:m.200_202delinsTGC")),
        format!("{contig}:m.200_202inv"),
    );

    // Within window_size (100) of the 3' end (250 bp): the window upper bound
    // 202 + 100 = 302 exceeds the contig length. Before the clamp fix this
    // errored and left the variant as `delins`; now it still canonicalizes.
    let near = provider_len(contig, 250, 200, "GCA");
    assert_eq!(
        norm(&near, &format!("{contig}:m.200_202delinsTGC")),
        format!("{contig}:m.200_202inv"),
    );
}

#[test]
fn issue_1044_mt_root_cause_also_restores_3prime_shift_near_end() {
    // The clamp also restores ordinary 3' shifting near the contig end,
    // confirming the root cause is the window fetch, not anything inv-specific.
    // A homopolymer run of A's at 1-based 195..=205; delete one A at 196 on a
    // 250 bp contig (196 + window overruns EOF before the clamp).
    let contig = "NC_012920.1";
    let mut bases: Vec<u8> = "ACGT".bytes().cycle().take(250).collect();
    for b in bases.iter_mut().take(205).skip(194) {
        *b = b'A';
    }
    let seq = String::from_utf8(bases).unwrap();
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { contig: seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    let p = JsonProvider::from_json(f.path()).unwrap();
    // 3'-shifts to the rightmost equivalent position (m.205), not m.196.
    assert_eq!(
        norm(&p, &format!("{contig}:m.196del")),
        format!("{contig}:m.205del"),
    );
}

#[test]
fn issue_1044_mt_span_past_contig_end_passes_through_unchanged() {
    // The clamp must NOT truncate a variant whose span runs past the contig
    // end: fetching a window shorter than the span would feed a truncated
    // reference to normalize_na_edit and mis-normalize. The clamp is gated on
    // `end <= len`, so these fall through to the safe fallback instead.
    let contig = "NC_012920.1";
    let p = provider_len(contig, 100, 1, "A"); // 100 bp contig, payload irrelevant
    for suffix in ["m.99_103inv", "m.105del", "m.105A>T"] {
        let input = format!("{contig}:{suffix}");
        assert_eq!(
            norm(&p, &input),
            input,
            "past-contig-end mt variant must pass through",
        );
    }
}
