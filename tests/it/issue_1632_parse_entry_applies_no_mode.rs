//! Issue #1632 — `parse_hgvs` documents itself as strict but applies no
//! `ErrorConfig`, so it accepts what strict rejects.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/1632>.
//!
//! `parse_hgvs` runs the fast path and then the grammar. It never constructs an
//! `InputPreprocessor`, so none of the input-hygiene ladder runs on it — not the
//! repairs lenient performs, and not the refusals that make strict strict. Only
//! `parse_hgvs_with_config` applies a mode.
//!
//! This file pins that as *measured behaviour* so the doc comment and the code
//! cannot drift apart again. It deliberately does **not** change what
//! `parse_hgvs` accepts: routing it through `ErrorConfig::strict()` is the
//! crate's headline entry point newly refusing inputs it accepts today, which is
//! a behaviour change for callers to weigh.
//!
//! It also sits directly downstream of `absolute-prohibition-enforcement-stage`,
//! which #1634 decided as **mode-dependent** — an absolute prohibition is refused
//! at strict-mode normalize rather than unconditionally at parse. A mode-less
//! entry point therefore has no defined enforcement stage under that ruling,
//! which is an argument for giving `parse_hgvs` a mode rather than against it;
//! what stage that mode should enforce at is the part the ruling settles and this
//! PR does not touch.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::parse_hgvs;

/// Parse under an explicit mode, flattening to `Ok(rendered)` / `Err(message)`
/// so the two entries can be compared on equal terms.
fn with_config(input: &str, config: ErrorConfig) -> Result<String, String> {
    parse_hgvs_with_config(input, config)
        .map(|r| r.result.to_string())
        .map_err(|e| e.to_string())
}

/// The bare entry, flattened the same way.
fn bare(input: &str) -> Result<String, String> {
    parse_hgvs(input)
        .map(|v| v.to_string())
        .map_err(|e| e.to_string())
}

// Every refusal below is asserted on its **message**, not on `ErrorCode`, and
// that choice is measured rather than stylistic. `ErrorCode` does not
// discriminate these rungs: the redundant-base-label refusal and the
// size-count-suffix refusal both report `ErrorCode::InvalidEdit`, and the bare
// entry's whitespace refusal carries `code: None` at all. Pinning the code
// would therefore be barely stronger than `is_err()` — it would pass on an
// unrelated `InvalidEdit` from anywhere in the grammar, which is the exact
// weakness these assertions exist to remove. The messages name the rung, so
// they are what the stated contract is pinned on.

/// Guards the comment above. Those `ErrorCode` values are the whole reason the
/// assertions in this file key on messages, and an unpinned claim about
/// behaviour is what #1632 is — so pin it here rather than let it rot into a
/// second wrong comment about the same parser.
///
/// If this test fails because a rung gained its own code, that is good news:
/// re-point the assertions in this file at the code and delete the comment.
#[test]
fn the_error_codes_do_not_discriminate_these_rungs() {
    use ferro_hgvs::error::{ErrorCode, FerroError};

    // Generic over the `Ok` type: the two entries succeed with different types
    // (`ParseResultWithWarnings<HgvsVariant>` vs `HgvsVariant`), and only the
    // error side is under test here.
    fn code_of<T>(r: Result<T, FerroError>) -> Option<String> {
        match r {
            Err(FerroError::Parse { diagnostic, .. }) => {
                diagnostic.and_then(|d| d.code).map(|c| format!("{c:?}"))
            }
            _ => panic!("expected a parse error"),
        }
    }

    let invalid_edit = Some(format!("{:?}", ErrorCode::InvalidEdit));

    assert_eq!(
        code_of(parse_hgvs_with_config(
            "NM_024312.4:r.-6_-3g[6]",
            ErrorConfig::strict()
        )),
        invalid_edit,
        "the redundant-base-label rung reports the generic InvalidEdit"
    );
    assert_eq!(
        code_of(parse_hgvs_with_config(
            "NC_TEST.1:g.266del3",
            ErrorConfig::strict()
        )),
        invalid_edit,
        "so does the size-count-suffix rung — one code, two rungs, which is why \
         a code assertion would not discriminate them"
    );
    assert_eq!(
        code_of(parse_hgvs("NM_TEST.1:c.10_12 del")),
        None,
        "and the bare entry's whitespace refusal carries no code at all"
    );
}

/// The plain finding: two inputs the *actual* strict config refuses are
/// accepted, unchanged, by the entry point whose doc comment claims strict.
///
/// `r.-6_-3g[6]` and `r.-125_-123cug[4]` state a base label that the range
/// already determines — the `RedundantRepeatLabel` rung (`W3013`) of the
/// hygiene ladder. The ladder never runs on `parse_hgvs`, so it has no opinion.
#[test]
fn the_bare_entry_accepts_two_inputs_strict_refuses() {
    for input in ["NM_024312.4:r.-6_-3g[6]", "NM_024312.4:r.-125_-123cug[4]"] {
        assert_eq!(
            bare(input).as_deref(),
            Ok(input),
            "parse_hgvs accepts the input unchanged"
        );
        let refusal = with_config(input, ErrorConfig::strict())
            .expect_err("strict must refuse a redundant repeat base label");
        assert!(
            refusal.contains("Repeat description has redundant base label"),
            "the refusal must name this rung, not merely fail; got: {refusal}"
        );
    }
}

/// The same seam viewed from the other side: inputs the grammar refuses for its
/// own reasons. `parse_hgvs` reports a grammar error; the mode-aware entry
/// reports the *mode-appropriate* answer — a named refusal under strict, and a
/// repair under lenient. So the bare entry is not merely "strict without the
/// repairs"; it is a third behaviour that is neither mode.
#[test]
fn the_bare_entry_gives_a_grammar_error_where_a_mode_gives_a_mode_answer() {
    // `del3` — a size-count suffix (`checklist.md:49`). Both entries refuse,
    // and the two refusals are not the same refusal: the grammar cites the
    // clause and offers the rewrite, strict names its hygiene rung.
    let bare_del3 =
        bare("NC_TEST.1:g.266del3").expect_err("the grammar refuses a size-count suffix");
    assert!(
        bare_del3.contains("checklist.md:49"),
        "the bare entry refuses in the grammar, citing the clause; got: {bare_del3}"
    );
    let strict_del3 = with_config("NC_TEST.1:g.266del3", ErrorConfig::strict())
        .expect_err("strict refuses a size-count suffix");
    assert!(
        strict_del3.contains("uses size-count suffix"),
        "strict refuses at its own hygiene rung; got: {strict_del3}"
    );
    assert_ne!(
        bare_del3, strict_del3,
        "the two entries must refuse for different stated reasons — that is what \
         makes the bare entry a third behaviour rather than strict-without-repairs"
    );
    assert_eq!(
        with_config("NC_TEST.1:g.266del3", ErrorConfig::lenient()).as_deref(),
        Ok("NC_TEST.1:g.266_268del"),
        "lenient repairs the size-count suffix into the range it denotes"
    );

    // An internal space (`general.md:96`).
    let bare_space =
        bare("NM_TEST.1:c.10_12 del").expect_err("the grammar refuses an internal space");
    assert!(
        bare_space.contains("Unexpected trailing characters"),
        "the bare entry sees only unparsed trailing input; got: {bare_space}"
    );
    let strict_space = with_config("NM_TEST.1:c.10_12 del", ErrorConfig::strict())
        .expect_err("strict refuses an internal space");
    assert!(
        strict_space.contains("Extra whitespace"),
        "strict names the whitespace rung; got: {strict_space}"
    );
    assert_eq!(
        with_config("NM_TEST.1:c.10_12 del", ErrorConfig::lenient()).as_deref(),
        Ok("NM_TEST.1:c.10_12del"),
        "lenient strips the internal whitespace"
    );
}

/// Leading and trailing whitespace is the one hygiene concern `parse_hgvs`
/// handles on its own, via `trim_hgvs`. Pinned so the doc's "grammar only"
/// wording is not read as "nothing at all".
#[test]
fn the_bare_entry_still_trims_surrounding_whitespace() {
    assert_eq!(
        bare("  NM_000088.3:c.459del  ").as_deref(),
        Ok("NM_000088.3:c.459del"),
    );
}

/// **Negative control for the intronic question, and the reason it is here.**
///
/// A bare-transcript intronic position (`checklist.md:20`) is accepted by
/// `parse_hgvs` — but it is *equally* accepted by
/// `parse_hgvs_with_config(_, strict())`. So the missing `ErrorConfig` is not
/// what lets that shape through: there is no parse-stage rung for it in any
/// mode. Its refusal lives at strict-mode **normalize**
/// (`IntronicOnBareTranscript` / `W4007`), which is where the
/// `bare-transcript-intronic-position` ruling put it — see
/// `corpus_prohibited_inputs.rs::a_bare_transcript_intronic_position_is_refused_in_strict_only`.
///
/// This pin exists so that closing #1632 by routing `parse_hgvs` through a mode
/// is not mistaken for a fix to the intronic shape, in either direction.
#[test]
fn the_missing_config_is_not_what_admits_a_bare_transcript_intronic_position() {
    let input = "NM_TEST.1:c.20+2del";
    assert_eq!(bare(input).as_deref(), Ok(input));
    assert_eq!(
        with_config(input, ErrorConfig::strict()).as_deref(),
        Ok(input),
        "strict input hygiene has no rung for checklist.md:20 either — the \
         refusal is a normalize-stage one"
    );
    assert_eq!(
        with_config(input, ErrorConfig::lenient()).as_deref(),
        Ok(input)
    );
}
