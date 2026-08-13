//! Harvest every multi-member cis allele out of the bulk corpora (#1443).
//!
//! # Why this exists
//!
//! `canonicalize_from_sequence` is gated on `members.len() > 1`
//! (`src/normalize/mod.rs`), so a single-member input never reaches the
//! partitioner, the changed-column bound, or the member-ordering rules. Across
//! the four bulk corpora — 9 949 738 rows — only **592** inputs are
//! multi-member cis alleles. Everything else exercises the per-member path.
//!
//! Those 592 are the only real-world evidence the repo has for that code path,
//! and they currently sit inside bulk corpora whose test modules *skip green*
//! when the fixtures are absent (`clinvar_hgvs_tests`, `cmrg_exhaustive_tests`,
//! `paraphase_exhaustive_tests` all `return` on a missing file). Those corpora
//! are not in the git tree — they are release assets fetched by
//! `scripts/fetch-test-fixtures.sh` — so on a checkout that has not fetched
//! them the evidence is silently invisible.
//!
//! Extracting them into one small committed fixture makes the axis
//! unconditional: `multi_member_cis_axis` reads the fixture, not the corpora,
//! and therefore runs everywhere.
//!
//! # Usage
//!
//! ```text
//! cargo run --features dev --example harvest_multi_member_cis
//! cargo run --features dev --example harvest_multi_member_cis -- --check
//! cargo run --features dev --example harvest_multi_member_cis -- --output <path>
//! ```
//!
//! `--check` is a real gate here, unlike for the generated spec artifacts: the
//! output **is** committed, so an absent or stale file is drift with a
//! baseline to drift from. Requires the bulk fixtures
//! (`scripts/fetch-test-fixtures.sh`); without them it reports which are missing
//! and exits non-zero rather than writing a truncated file.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use ferro_hgvs::hgvs::variant::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

/// The bulk corpora to harvest, in a fixed order so the output is stable.
const SOURCES: &[&str] = &[
    "tests/fixtures/bulk/clinvar_hgvs_500k.json.gz",
    "tests/fixtures/bulk/clinvar_hgvs_unique.json.gz",
    "tests/fixtures/validation/cmrg_genes_exhaustive.json.gz",
    "tests/fixtures/validation/paraphase_genes_exhaustive.json.gz",
];

const DEFAULT_OUTPUT: &str = "tests/fixtures/cis/multi_member_cis_alleles.json";

/// One harvested row, plus enough provenance to find it again upstream.
#[derive(serde::Serialize, serde::Deserialize, PartialEq, Eq, PartialOrd, Ord)]
struct Row {
    /// The HGVS string exactly as the corpus carried it.
    input: String,
    /// Corpus basename, so a surprising row can be traced back.
    source: String,
    /// The corpus's own identifier where it had one.
    #[serde(skip_serializing_if = "Option::is_none")]
    variation_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    gene: Option<String>,
    /// How many members the *parser* found in the allele, so a nested
    /// inserted-range payload's own `;` is not counted as a separator.
    members: usize,
}

#[derive(serde::Serialize, serde::Deserialize)]
struct Fixture {
    description: String,
    /// How this file is produced, in the file, so a reader need not guess.
    generator: String,
    sources: Vec<String>,
    /// Rows scanned across every source — the denominator that makes the
    /// harvest's rarity claim checkable rather than asserted.
    rows_scanned: usize,
    rows: Vec<Row>,
}

/// Whether `input` is a **cis** allele with more than one member, and how many.
///
/// Two stages, in this order because of cost. Cheap string filters reject
/// essentially all of the ten million corpus rows without allocating; whatever
/// survives — a few thousand strings — is handed to `parse_hgvs`, which is what
/// actually decides and where the member count comes from.
///
/// **The parse is not an optimisation, it is the rule.** The question "is this
/// one allele bracket holding several cis members?" is not answerable by
/// scanning for `[`, `]` and `;`, because an inserted-range payload carries its
/// own brackets *and* its own separator:
///
/// ```text
/// LRG_676:g.32086_32087ins[N[1090];32074_32086]
/// ```
///
/// That is a single insertion, but the outermost `[` … `]` is the payload, so a
/// bracket scan reads it as a two-member allele. The same payload nested inside
/// a genuine allele inflates that allele's count instead — `LRG_501:g.[…insN;…
/// ins[SEQ;a_b];…insAT]` scans as four members and is three. Both shapes are
/// present in the corpora, and a depth-aware scan fixes only the second: depth
/// alone cannot tell an allele bracket from a payload bracket.
///
/// Two shapes are out of scope even when they parse, and the cheap filters
/// encode them so the parse is never reached:
///
/// * `[a];[b]` — trans, a different path entirely.
/// * `[a(;)b]` — unphased. The `(;)` separator means "cis or trans unknown", so
///   it is not evidence about the cis partitioner.
///
/// The parse re-checks both through [`AllelePhase`], so the filters are a
/// speed-up rather than the definition.
fn multi_member_cis(input: &str) -> Option<usize> {
    let open = input.find('[')?;
    let close = input.rfind(']')?;
    if close <= open {
        return None;
    }
    let inner = &input[open + 1..close];
    // No separator anywhere between the brackets: a payload or a repeat, never
    // a multi-member allele. Checked first because it is what discards the bulk.
    if !inner.contains(';') {
        return None;
    }
    // `];[` is trans, and `(;)` is unphased. Both are out of scope.
    if input.contains("];[") || inner.contains("(;)") {
        return None;
    }
    match parse_hgvs(input) {
        Ok(HgvsVariant::Allele(allele))
            if allele.phase == AllelePhase::Cis && allele.variants.len() > 1 =>
        {
            Some(allele.variants.len())
        }
        _ => None,
    }
}

fn read_gz(path: &Path) -> Option<Vec<u8>> {
    let file = std::fs::File::open(path).ok()?;
    let mut out = Vec::new();
    flate2::read::GzDecoder::new(file)
        .read_to_end(&mut out)
        .ok()?;
    Some(out)
}

/// Pull the row array out of a corpus, which wraps it in a metadata object
/// under a key that differs between corpora.
fn rows_of(value: &serde_json::Value) -> Option<&Vec<serde_json::Value>> {
    if let Some(array) = value.as_array() {
        return Some(array);
    }
    value
        .as_object()?
        .values()
        .find_map(|candidate| candidate.as_array())
}

fn harvest() -> Result<Fixture, String> {
    let mut rows: Vec<Row> = Vec::new();
    let mut scanned = 0usize;
    let mut missing: Vec<&str> = Vec::new();

    for source in SOURCES {
        let Some(bytes) = read_gz(Path::new(source)) else {
            missing.push(source);
            continue;
        };
        let value: serde_json::Value =
            serde_json::from_slice(&bytes).map_err(|e| format!("{source}: {e}"))?;
        let corpus = rows_of(&value).ok_or_else(|| format!("{source}: no row array"))?;
        let basename = Path::new(source)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or(source)
            .to_string();

        for entry in corpus {
            scanned += 1;
            let Some(input) = entry.get("input").and_then(|v| v.as_str()) else {
                continue;
            };
            let Some(members) = multi_member_cis(input) else {
                continue;
            };
            rows.push(Row {
                input: input.to_string(),
                source: basename.clone(),
                variation_id: entry
                    .get("variation_id")
                    .and_then(|v| v.as_str())
                    .map(str::to_string),
                gene: entry
                    .get("gene")
                    .and_then(|v| v.as_str())
                    .map(str::to_string),
                members,
            });
        }
    }

    if !missing.is_empty() {
        return Err(format!(
            "missing bulk corpora, so the harvest would be silently short: {}\n\
             These are release assets, not git objects; run \
             `scripts/fetch-test-fixtures.sh` and retry.",
            missing.join(", ")
        ));
    }

    // Sorted and deduplicated: the same allele appears in more than one corpus,
    // and a fixture whose order depended on corpus iteration would churn its
    // diff for no reason.
    //
    // The dedup is load-bearing because the two ClinVar corpora genuinely
    // overlap — **not** because either contains the other. They are drawn
    // differently and neither is a subset of the other:
    //
    //   `clinvar_hgvs_500k`    a stratified sample of ClinVar HGVS *expressions*
    //                          (a variant may contribute several)
    //   `clinvar_hgvs_unique`  one expression per `VariationID`, across all of
    //                          ClinVar — an order of magnitude larger
    //
    // Measured on the committed fixtures, their intersection is a small
    // minority of either side, and each holds inputs the other does not; among
    // the bracketed multi-member candidates this harvest cares about, a handful
    // are shared and both corpora contribute rows the other lacks. So dropping
    // either source would lose alleles, and skipping the dedup would double-count
    // the shared ones.
    //
    // (This comment previously claimed `clinvar_hgvs_unique` was "a subset of
    // the 500k by construction". It is not, in either direction.)
    rows.sort();
    rows.dedup_by(|a, b| a.input == b.input);

    Ok(Fixture {
        description: "Multi-member cis alleles harvested from the bulk corpora. These are the \
                      only inputs in those corpora that reach `canonicalize_from_sequence`, \
                      which is gated on `members.len() > 1`."
            .to_string(),
        generator: "cargo run --features dev --example harvest_multi_member_cis".to_string(),
        sources: SOURCES.iter().map(|s| s.to_string()).collect(),
        rows_scanned: scanned,
        rows,
    })
}

/// Members-per-allele histogram, printed so a harvest run reports the shape of
/// what it found rather than only a count.
fn census(fixture: &Fixture) -> BTreeMap<usize, usize> {
    let mut counts = BTreeMap::new();
    for row in &fixture.rows {
        *counts.entry(row.members).or_insert(0) += 1;
    }
    counts
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let check = args.iter().any(|a| a == "--check");
    let output: PathBuf = args
        .iter()
        .position(|a| a == "--output")
        .and_then(|i| args.get(i + 1))
        .map_or_else(|| PathBuf::from(DEFAULT_OUTPUT), PathBuf::from);

    let fixture = match harvest() {
        Ok(f) => f,
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    };

    let mut rendered = serde_json::to_string_pretty(&fixture).expect("serialize");
    rendered.push('\n');

    println!(
        "scanned {} rows across {} corpora -> {} multi-member cis alleles ({:.4}%)",
        fixture.rows_scanned,
        SOURCES.len(),
        fixture.rows.len(),
        100.0 * fixture.rows.len() as f64 / fixture.rows_scanned.max(1) as f64,
    );
    for (members, count) in census(&fixture) {
        println!("  {members} members: {count}");
    }

    if check {
        match std::fs::read_to_string(&output) {
            Ok(existing) if existing == rendered => {
                println!("up to date: {}", output.display());
            }
            Ok(_) => {
                eprintln!("stale: {} differs from a fresh harvest", output.display());
                std::process::exit(1);
            }
            Err(e) => {
                eprintln!("missing: {} ({e})", output.display());
                std::process::exit(1);
            }
        }
        return;
    }

    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).expect("create output directory");
    }
    std::fs::write(&output, rendered).expect("write fixture");
    println!("wrote {}", output.display());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn selects_only_multi_member_cis() {
        assert_eq!(multi_member_cis("NM_000208.3:c.[2504G>T;2525C>T]"), Some(2));
        assert_eq!(multi_member_cis("T:g.[1del;5del;9del]"), Some(3));
    }

    /// Each exclusion is a distinct shape, and each would poison the axis in a
    /// different way if admitted.
    #[test]
    fn rejects_everything_that_is_not_a_cis_allele() {
        // Single member: takes the per-member path, so it is not evidence.
        assert_eq!(multi_member_cis("NM_000208.3:c.2504G>T"), None);
        // Trans: a different code path entirely.
        assert_eq!(multi_member_cis("NM_000208.3:c.[2504G>T];[2525C>T]"), None);
        // Unphased: `(;)` means cis-or-trans-unknown, so it says nothing about
        // the cis partitioner.
        assert_eq!(multi_member_cis("NM_015046.5:c.[59G>A(;)3809C>T]"), None);
        // Brackets that are not an allele at all.
        assert_eq!(multi_member_cis("T:g.10_11ins[20_30]"), None);
        assert_eq!(multi_member_cis("T:g.1_3A[12]"), None);
        assert_eq!(multi_member_cis("NM_000208.3:c.2504G>T"), None);
    }

    /// An inserted-range payload carrying its own `;` is a **single** variant,
    /// not an allele. All four shapes below are real corpus rows that the
    /// previous `find('[')..rfind(']')` scan admitted as two- and three-member
    /// alleles; 12 such rows sat in the committed fixture.
    #[test]
    fn a_semicolon_bearing_payload_is_not_an_allele() {
        // `N[1090]` is a repeat inside the payload, so the payload nests too.
        assert_eq!(
            multi_member_cis("LRG_676:g.32086_32087ins[N[1090];32074_32086]"),
            None
        );
        // A payload mixing literal sequence with a range on another accession.
        assert_eq!(
            multi_member_cis("NM_017635.5:c.438_439ins[TCTT;KT192064.1:1_310]"),
            None
        );
        // The same trap on a `delins`, with three payload pieces.
        assert_eq!(
            multi_member_cis(
                "NC_000013.10:g.114819939_qterdelins[96729864_114814234inv;\
                 96735632_104289803]"
            ),
            None
        );
        assert_eq!(
            multi_member_cis(
                "NC_000017.11:g.80114186_80114187ins[80114172_80114186;\
                 NC_000020.11:g.2823027_2826302;AAA]"
            ),
            None
        );
    }

    /// A genuine allele whose *middle* member carries a payload with its own
    /// `;`. This one is real too (`LRG_501`, ClinVar 2535): the old scan counted
    /// four members where the parser finds three.
    #[test]
    fn a_nested_payload_does_not_inflate_a_genuine_allele() {
        assert_eq!(
            multi_member_cis(
                "LRG_501:g.[83418_83419insN;83419_83420ins[TAGCTTTGCTAGGTGTTAAATAGTAATGT\
                 AATTATATTTCCTATAATACAGCAATATA;84466_84575];83423_83424insAT]"
            ),
            Some(3)
        );
    }

    /// A malformed row must be declined rather than panic the harvest — it runs
    /// over ten million strings that no parser has yet judged.
    #[test]
    fn malformed_brackets_are_declined() {
        assert_eq!(multi_member_cis("T:g.]1del;5del["), None);
        assert_eq!(multi_member_cis("T:g.[1del;5del"), None);
        assert_eq!(multi_member_cis(""), None);
    }
}
