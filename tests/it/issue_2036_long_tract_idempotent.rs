//! #2036 — normalizing a long tandem tract must be idempotent.
//!
//! `long_tract_window_provenance` (in `examples/dump_normalized_corpus.rs`)
//! found that minting a repeat from an insertion at a tract junction under-reports
//! the tract's 5' extent when the tract is longer than the normalizer's 100 bp
//! shuffle window: the mint walks the tract only within that window, so its
//! reference span stops short of the true 5' boundary. Re-normalizing the minted
//! repeat notation then walks the *whole* tract and finds it starts three bases
//! (one `ACG` copy) further 5', bumping the copy count. Each pass finds another
//! copy, so `norm(norm(x)) != norm(x)`:
//!
//! ```text
//! input: NC_TEST.1:g.558_559insACGACG
//! once:  NC_TEST.1:g.460_558ACG[35]
//! twice: NC_TEST.1:g.457_558ACG[36]
//! ```
//!
//! The fixed point is the one that names the maximal reference tract
//! (`g.457_558ACG[36]`), which is the form `normalize_repeat` already converges on
//! when handed a repeat directly — so this is a determinate idempotency defect,
//! not a canonical-form policy choice.

use crate::common::synthetic::{SyntheticBuilder, PAD_OFFSET};
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

/// Repeated unit and design, mirroring `long_tract_window_provenance`'s 34-copy
/// row: a 102-base `ACG` tract, longer than the 100 bp shuffle window.
const UNIT: &str = "ACG";
const COPIES: usize = 34;
/// 0-based offset of the tract within the core, matching the corpus family.
const TRACT_START: usize = 200;
const CORE_LEN: usize = 800;

/// A base for flank position `at` that will not let the tract extend through it,
/// and does not create a two-base run against its neighbour. Mirrors
/// `dump_normalized_corpus::flank_base`.
fn flank_base(core: &[u8], at: usize, forbidden: u8) -> u8 {
    let neighbour = if at == 0 { None } else { Some(core[at - 1]) };
    *b"GCTA"
        .iter()
        .find(|&&candidate| {
            candidate != forbidden
                && Some(candidate) != neighbour
                && !(at == 0 && candidate == b'A')
                && !(at == core.len() - 1 && candidate == b'T')
        })
        .expect("four candidates cannot all be excluded by three constraints")
}

/// Build the corpus family's `tract_core`: a `GCTA`-period core with `COPIES`
/// copies of `UNIT` at `TRACT_START`, flanked so the tract cannot leak.
fn tract_core() -> String {
    let mut core: Vec<u8> = (0..CORE_LEN).map(|i| b"GCTA"[i % 4]).collect();
    let tract = UNIT.repeat(COPIES);
    core[TRACT_START..TRACT_START + tract.len()].copy_from_slice(tract.as_bytes());
    let unit = UNIT.as_bytes();
    core[TRACT_START - 1] = flank_base(&core, TRACT_START - 1, unit[unit.len() - 1]);
    let after = TRACT_START + tract.len();
    if after < CORE_LEN {
        core[after] = flank_base(&core, after, unit[0]);
    }
    String::from_utf8(core).expect("valid UTF-8")
}

fn normalizer(dir: ShuffleDirection) -> Normalizer<ferro_hgvs::reference::mock::MockProvider> {
    // The defect only exists when the tract outruns the shuffle window, so the
    // guard is worthless if the window ever grows past the tract. Read the window
    // rather than hard-coding it (#1925): a config change that lifts it above the
    // 102-base tract would otherwise let this test pass without reaching the
    // window-boundary path it exists to cover.
    let window = NormalizeConfig::default().window_size as usize;
    assert!(
        UNIT.len() * COPIES > window,
        "the {}-base tract must exceed the {window}-base shuffle window",
        UNIT.len() * COPIES,
    );
    Normalizer::with_config(
        SyntheticBuilder::genomic(&tract_core()).build(),
        NormalizeConfig::default().with_direction(dir),
    )
}

fn norm(n: &Normalizer<ferro_hgvs::reference::mock::MockProvider>, input: &str) -> String {
    n.normalize(&parse_hgvs(input).expect("parse"))
        .expect("normalize")
        .to_string()
}

#[test]
fn long_tract_insertion_mint_is_idempotent() {
    // Genomic (1-based) coordinates of the tract in the padded reference.
    let tract_start = PAD_OFFSET as usize + TRACT_START + 1; // 457
    let tract_end = tract_start + UNIT.len() * COPIES - 1; // 558
    assert_eq!((tract_start, tract_end), (457, 558), "design matches #2036");

    let n = normalizer(ShuffleDirection::ThreePrime);
    // Insert two copies at the 3' junction of the tract.
    let input = format!(
        "NC_TEST.1:g.{tract_end}_{}ins{}",
        tract_end + 1,
        UNIT.repeat(2)
    );
    assert_eq!(input, "NC_TEST.1:g.558_559insACGACG");

    let once = norm(&n, &input);
    let twice = norm(&n, &once);
    assert_eq!(
        once, twice,
        "NOT IDEMPOTENT\n  input={input}\n  once ={once}\n  twice={twice}"
    );

    // The fixed point names the maximal reference tract: 457..558 = 34 copies,
    // plus the two inserted copies = 36.
    assert_eq!(
        once, "NC_TEST.1:g.457_558ACG[36]",
        "must reach the maximal-tract fixed point in one pass"
    );
}

/// The 5' junction insertion under the 5' shuffle direction is the mirror image:
/// `insertion_to_repeat` walks the tract *3'*-ward and the same window truncation
/// applied, minting a repeat whose 3' extent stopped short of the tract's end.
/// It must also reach the maximal-tract fixed point in one pass.
#[test]
fn long_tract_five_prime_insertion_is_idempotent() {
    let tract_start = PAD_OFFSET as usize + TRACT_START + 1; // 457
    let tract_end = tract_start + UNIT.len() * COPIES - 1; // 558
    let n = normalizer(ShuffleDirection::FivePrime);
    let input = format!(
        "NC_TEST.1:g.{}_{tract_start}ins{}",
        tract_start - 1,
        UNIT.repeat(2)
    );
    assert_eq!(input, "NC_TEST.1:g.456_457insACGACG");

    let once = norm(&n, &input);
    let twice = norm(&n, &once);
    assert_eq!(
        once, twice,
        "NOT IDEMPOTENT (5')\n  input={input}\n  once ={once}\n  twice={twice}"
    );
    assert_eq!(
        once,
        format!("NC_TEST.1:g.{tract_start}_{tract_end}ACG[{}]", COPIES + 2),
        "5' direction must reach the maximal-tract fixed point in one pass"
    );
}
