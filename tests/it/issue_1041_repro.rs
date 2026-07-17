//! Issue #1041 — a whole-run reverse-complement delins near a contig's 3' end
//! was left as `delins` instead of canonicalizing to `inv`.
//!
//! Root cause: `normalize_genome` fetches a `+/-window_size` window around the
//! variant. On the ordinary (non-special) path the upper bound was not clamped
//! to the contig length, and every provider ERRORS on a past-EOF read rather
//! than clamping. So an indel within `window_size` (100 bp) of the contig 3'
//! end fetched past EOF, `get_sequence` errored, and the whole variant dropped
//! into the minimal-notation fallback — skipping 3' shift AND the
//! delins->inv/sub/dup canonicalization. Clamping `fetch_end` to the contig
//! length fixes it.

use ferro_hgvs::{parse_hgvs, JsonProvider, Normalizer};
use std::io::Write;

/// Genome-capable provider: one `contig` of `len` cyclic ACGT bytes with
/// `payload` written 1-based at `pos1`.
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
    let v = parse_hgvs(input).unwrap();
    Normalizer::new(p.clone())
        .normalize(&v)
        .unwrap()
        .to_string()
}

#[test]
fn issue_1041_whole_run_revcomp_is_inv_regardless_of_3prime_proximity() {
    // ref[200..202] = GCA, revcomp = TGC. Same variant, two contig lengths.

    // Far from the end (600 bp): window [100, 302] fits comfortably.
    let far = provider_len("c", 600, 200, "GCA");
    assert_eq!(norm(&far, "c:g.200_202delinsTGC"), "c:g.200_202inv");

    // Within window_size (100) of the 3' end (250 bp): the window upper bound
    // 202 + 100 = 302 exceeds the contig length. Before the clamp fix this
    // errored and left the variant as `delins`; now it still canonicalizes.
    let near = provider_len("c", 250, 200, "GCA");
    assert_eq!(norm(&near, "c:g.200_202delinsTGC"), "c:g.200_202inv");
}

#[test]
fn issue_1041_compound_subs_near_3prime_end_merge_to_inv() {
    // The reported shape: a compound of adjacent substitutions whose merged
    // span is a reverse complement, sited near the 3' end.
    let near = provider_len("c", 250, 200, "GCA");
    assert_eq!(norm(&near, "c:g.[200G>T;201C>G;202A>C]"), "c:g.200_202inv");
}

#[test]
fn issue_1041_root_cause_also_restores_3prime_shift() {
    // The clamp also restores ordinary 3' shifting near the contig end,
    // confirming the root cause is the window fetch, not anything inv-specific.
    // A homopolymer run of A's at 1-based 195..=205; delete one A at 196.
    let mut bases: Vec<u8> = "ACGT".bytes().cycle().take(250).collect();
    for b in bases.iter_mut().take(205).skip(194) {
        *b = b'A';
    }
    let seq = String::from_utf8(bases).unwrap();
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { "c": seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    let p = JsonProvider::from_json(f.path()).unwrap();
    // 3'-shifts to the rightmost equivalent position (g.205), not g.196.
    assert_eq!(norm(&p, "c:g.196del"), "c:g.205del");
}

#[test]
fn issue_1041_span_past_contig_end_passes_through_unchanged() {
    // The clamp must NOT truncate a variant whose span runs past the contig
    // end: fetching a window shorter than the span would feed a truncated
    // reference to normalize_na_edit and mis-normalize (e.g. `g.99_103inv`
    // collapsing to `g.99_103=`). The clamp is gated on `end <= len`, so these
    // fall through to the safe minimal-notation pass-through instead.
    let p = provider_len("c", 100, 1, "A"); // 100 bp contig, payload irrelevant
    for input in [
        "c:g.99_103inv",
        "c:g.105del",
        "c:g.98_105delinsGGGGGGGG",
        "c:g.105A>T",
    ] {
        assert_eq!(
            norm(&p, input),
            input,
            "past-contig-end variant must pass through"
        );
    }
}
