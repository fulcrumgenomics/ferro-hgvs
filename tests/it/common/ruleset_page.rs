//! The normalization rules page, read by the guards that check a document
//! points at the rules rather than copying them.
//!
//! The page is the single source of the rules, so these helpers read its rule
//! names from the page itself instead of keeping a copy of them here.

use std::path::PathBuf;

/// The page that states the normalization rules.
pub const RULESET_PAGE: &str = "docs/src/reference/normalization-rules.md";

/// The bolded opening words of each numbered rule on the page, with a trailing
/// `.` or `:` removed: `Conformant`, `Recommended form`, and so on.
///
/// Panics if the page is missing or lists no rules, so a check built on this
/// cannot pass by testing nothing.
pub fn rule_openers() -> Vec<String> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(RULESET_PAGE);
    let page = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
    let openers: Vec<String> = page
        .lines()
        .filter_map(|line| {
            let line = line.trim_start();
            let unnumbered = line.trim_start_matches(|c: char| c.is_ascii_digit());
            if unnumbered.len() == line.len() {
                return None;
            }
            let rest = unnumbered.strip_prefix(". **")?;
            let (opener, _) = rest.split_once("**")?;
            Some(opener.trim_end_matches(['.', ':']).to_string())
        })
        .collect();
    assert!(
        !openers.is_empty(),
        "found no numbered rules on {RULESET_PAGE}, so a check built on them would test \
         nothing; the page's list format changed and this helper no longer reads it"
    );
    openers
}

/// The rule openers that `text` repeats in bold.
///
/// Matches `**<opener>` followed by anything, so it catches the page's own form
/// (`**Conformant.**`) and a summary form (`**Conformant** — …`) alike.
pub fn restated_rule_openers(text: &str) -> Vec<String> {
    rule_openers()
        .into_iter()
        .filter(|opener| text.contains(&format!("**{opener}")))
        .collect()
}
