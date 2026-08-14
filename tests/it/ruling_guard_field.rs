//! Every record in the adjudication ledger must say what enforces its ruling,
//! and what it says must resolve.
//!
//! # The gap this closes
//!
//! A record and the test that enforces it were connected by nothing but prose,
//! and often by no prose at all. Measured over the ledger as it stood when this
//! file was written — 28 records, 25 `decided`, 3 `undecided`. (The ledger has
//! since gained #1882's `c-description-against-an-unresolvable-cds-is-refused`,
//! so the current census is 29 / 26 / 3; the table below is the *historical*
//! measurement of the gap and is deliberately not restated against it.)
//!
//! | how a record pointed at its enforcement | records |
//! |---|---:|
//! | a `<path>.rs::<function>` citation somewhere in its prose | 5 |
//! | a bare `tests/….rs` path, naming a file but no proposition | 12 |
//! | an `equivalence_classes` array — the one structured pointer | 4 |
//! | **nothing, by any mechanism** | **8** |
//!
//! (The first three overlap; the last does not.) So for eight records the
//! ledger read as enforced while pointing at nothing, and for twelve more it
//! pointed at a file — which cannot go stale in any way a reader would notice,
//! because a file that still exists still resolves however much its contents
//! have moved on.
//!
//! # Why a field and not a scan
//!
//! Prose can be scanned, and `ruling_citation_currency.rs` scans it. But a scan
//! is *opportunistic*: it judges the records that happen to cite a guard and is
//! structurally blind to the ones that cite nothing, which is exactly the
//! population that most needs judging. No cleverer matcher fixes that — there
//! is nothing in the text to match.
//!
//! So the pointer is declared instead, in the ledger, as a field. The model is
//! the `Representation-Change:` trailer: **declining is a first-class answer,
//! and what is rejected is silence.** A record may name its guards, or it may
//! say that nothing enforces it and why — and both of those are answers. What
//! it may not do is omit the field, because an absent one is indistinguishable
//! from an unconsidered one.
//!
//! The precedent inside the ledger is `equivalence_classes`, which is already
//! enforced end to end: the generator refuses an unknown class id, and the
//! guard that reads it asserts and checks its own vacuity. This field is that
//! same shape, for the much larger population whose evidence is an ordinary
//! test rather than a curated class.
//!
//! # It composes with `ruling_citation_currency.rs`; neither subsumes the other
//!
//! They ask different questions and both are worth asking:
//!
//! - **This file** asks whether every record *has* an enforcement story, and
//!   whether the story resolves. It is exhaustive over records by construction.
//! - **That file** asks whether the discursive prose *around* a record stays
//!   true — that a comment does not call a settled record open, and that an id
//!   named in a rationale still exists. A record can carry a perfectly good
//!   `guard` field and still be described wrongly three files away.
//!
//! The interaction worth knowing about: a record's prose may name a guard that
//! has been renamed out from under it while this field names the live one. That
//! is not a contradiction to be resolved here — the prose sentence is the
//! record's own argument and repairing it is an adjudication, not a rename.
//! This file deliberately reads **only** the field.
//!
//! # What is checked, and what is deliberately not
//!
//! Checked: the field is **present** (here, in
//! [`every_record_declares_a_guard`] — see that test for why the check lives at
//! this stage rather than at parse); that it is well-formed, once present (in
//! `common::rulings`, at parse); that every citation names a `#[test]` which
//! exists in the file it names; that at least one of them **runs** (a record
//! guarded only by `#[ignore]`d tests is guarded by nothing today); and that
//! every named test contains an assertion, because a test that runs and asserts
//! nothing passes (#1858).
//!
//! The presence check and the well-formedness checks sit at different stages on
//! purpose. A record that has not answered the question yet is the ordinary
//! state of a PR in flight, and it must produce a failure that names the fix; a
//! record whose `guard` contradicts itself is broken, and there is nothing
//! downstream that can use it.
//!
//! Not checked: whether the named test enforces *this* ruling rather than a
//! neighbouring one. That is a judgement a matcher cannot make, and pretending
//! otherwise would put a veneer of mechanism over the thing that actually
//! carries it, which is review. What the field buys is that the judgement is
//! **written down, once, in the ledger** — so it can be argued with — instead of
//! being reconstructed from prose by each reader.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use crate::common::rulings::{self, split_guard_citation, Guard, LEDGER_RELATIVE_PATH};

/// Floor for how many guard citations the ledger must carry.
///
/// Measured at **46** over the 29 records currently in the ledger (42 when this
/// file was written, plus the four #1882's record names). The floor sits well
/// below that so ordinary ledger edits cannot trip it, while a reader that
/// stopped seeing citations — the failure that would make every check here pass
/// by looking at an empty set — still would.
const GUARD_CITATION_FLOOR: usize = 20;

/// Floor for how many `#[test]` functions the tree scan must find.
///
/// Measured at over nine thousand. A scan that found a handful would resolve
/// nothing and report every citation as missing, which is a loud failure; a
/// scan that found *most* of them would resolve most citations and silently
/// stop judging the rest, which is not.
const SCANNED_TEST_FLOOR: usize = 5_000;

/// Roots scanned for test functions, relative to the crate root.
const SCAN_ROOTS: &[&str] = &["src", "tests"];

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// One `#[test]` function found in the tree.
#[derive(Clone, Debug)]
struct TestFunction {
    /// 1-based line of the `fn` keyword, for failure messages.
    line: usize,
    /// Whether the function carries `#[ignore]`.
    ignored: bool,
    /// Whether anything in its body — directly or through a helper defined in
    /// the same file — asserts. See [`asserts_directly`].
    asserts: bool,
}

/// Every `#[test]` in the tree, keyed by `(scan-relative path, function name)`.
///
/// Keyed by path as well as by name because two files may legitimately hold
/// tests of the same name, and a citation that resolved to "some test called
/// this, somewhere" would let a guard move between modules unnoticed.
fn test_functions() -> BTreeMap<(String, String), TestFunction> {
    let root = crate_root();
    let mut found = BTreeMap::new();
    for scan_root in SCAN_ROOTS {
        let dir = root.join(scan_root);
        assert!(
            dir.is_dir(),
            "scan root {} does not exist — every citation would report as missing",
            dir.display()
        );
        collect(&dir, scan_root, &mut found);
    }
    found
}

fn collect(dir: &Path, rel: &str, out: &mut BTreeMap<(String, String), TestFunction>) {
    let entries =
        std::fs::read_dir(dir).unwrap_or_else(|e| panic!("read_dir {}: {e}", dir.display()));
    let mut paths: Vec<PathBuf> = entries
        .map(|entry| {
            entry
                .unwrap_or_else(|err| panic!("dir entry under {}: {err}", dir.display()))
                .path()
        })
        .collect();
    paths.sort();
    for path in paths {
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or_default()
            .to_string();
        let child_rel = format!("{rel}/{name}");
        if path.is_dir() {
            collect(&path, &child_rel, out);
            continue;
        }
        if path.extension().and_then(|e| e.to_str()) != Some("rs") {
            continue;
        }
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        // Parsed once and shared with the assertion pass. Parsing per file was
        // previously done twice — once here and once inside
        // `asserts_transitively` — over every `.rs` file in the tree.
        let functions = functions_in(&text);
        for function in &functions {
            if function.is_test {
                out.insert(
                    (child_rel.clone(), function.name.clone()),
                    TestFunction {
                        line: function.line,
                        ignored: function.ignored,
                        asserts: false, // filled in below, once helpers are known
                    },
                );
            }
        }
        // Assertion detection is per file, because a helper is resolved within
        // the file that defines it — see `asserts_transitively`.
        for (name, asserts) in asserts_transitively(&functions) {
            if let Some(entry) = out.get_mut(&(child_rel.clone(), name)) {
                entry.asserts = asserts;
            }
        }
    }
}

/// One function found in a source file, test or not.
struct ParsedFunction {
    name: String,
    line: usize,
    is_test: bool,
    ignored: bool,
    /// Whether the function carries `#[should_panic]`, which is an assertion
    /// living in an attribute rather than in the body. See [`asserts_directly`].
    should_panic: bool,
    body: String,
}

/// Which lines of `text` sit inside a raw string literal.
///
/// The scan reads source text, and this repository embeds Rust source *as
/// data*: [`the_assertion_detector_follows_a_local_helper`] puts three
/// `#[test] fn`s inside an `r#"…"#` literal precisely so the detector has
/// something to read. Without this mask those declarations register as real
/// tests, and a ledger citation naming one **resolves** — the record then
/// passes every check in this file while nothing enforces it, which is the
/// exact state the `guard` field exists to abolish, reached through this
/// file's own fixtures. Measured before this mask: six such declarations in
/// the tree, five of them here.
///
/// The mask is **generic over the hash count** (`r"`, `r#"`, `r##"`, …) and
/// closes only on a `"` followed by at least as many `#`. A mask hardcoded to
/// one hash count is not a floor but a bug waiting for the first nested
/// fixture: [`a_test_declared_inside_a_raw_string_literal_does_not_register`]
/// has to wrap a `r#"…"#` fixture in a `r##"…"##` one, and a `r#"`-only mask
/// would leave that test's *own* inner declarations registering as phantoms.
///
/// A line that only *opens* a literal is reported as outside, because the code
/// before the delimiter is real. An ordinary `"…"` literal cannot span a line
/// break, so it can never hide a `fn` declaration line and is not tracked; the
/// `r` must not be preceded by an identifier character, so `for`, `char` and
/// friends do not open anything.
/// Scans bytes rather than `char`s, and skips any line with no `"` at all.
/// Both matter: this runs over every `.rs` file under [`SCAN_ROOTS`], and a
/// first cut that collected each line into a `Vec<char>` took the five checks
/// here from ~7 s to ~68 s each — past nextest's 60 s SLOW threshold. Every
/// delimiter byte is ASCII, so byte indexing cannot split a multi-byte
/// character in a way that changes an answer.
fn raw_string_lines(lines: &[&str]) -> Vec<bool> {
    let mut inside = Vec::with_capacity(lines.len());
    let mut open: Option<usize> = None;
    for line in lines {
        inside.push(open.is_some());
        let bytes = line.as_bytes();
        // Every raw-string delimiter contains a `"`, so a line without one can
        // neither open nor close a literal. Most lines in the tree are of that
        // shape, and this runs over ~9,800 tests' worth of source.
        if !bytes.contains(&b'"') {
            continue;
        }
        // An opener is `r"` or `<#>"`, so a line holding neither byte pair can
        // only ever *close* one. Skipping the scan when nothing is open is the
        // difference between reading every quoted line and reading the few that
        // can matter.
        if open.is_none() && !has_opener_bytes(bytes) {
            continue;
        }
        let mut index = 0;
        while index < bytes.len() {
            match open {
                // Closing: `"` then at least `hashes` `#`.
                Some(hashes) => {
                    if bytes[index] == b'"' && run_of_hashes(bytes, index + 1) >= hashes {
                        open = None;
                        index += 1 + hashes;
                        continue;
                    }
                    index += 1;
                }
                // Opening: a non-identifier `r`, then `#`*n, then `"`.
                None => {
                    let preceded_by_ident = index
                        .checked_sub(1)
                        .is_some_and(|p| bytes[p].is_ascii_alphanumeric() || bytes[p] == b'_');
                    if bytes[index] == b'r' && !preceded_by_ident {
                        let hashes = run_of_hashes(bytes, index + 1);
                        if bytes.get(index + 1 + hashes) == Some(&b'"') {
                            open = Some(hashes);
                            index += 2 + hashes;
                            continue;
                        }
                    }
                    index += 1;
                }
            }
        }
    }
    inside
}

/// Whether `bytes` could hold a raw-string *opener* at all.
///
/// A necessary condition, not a sufficient one: an opener is `r` then zero or
/// more `#` then `"`, so the two bytes before the quote are either `r"` or
/// `#"`. A line containing neither pair cannot open a literal, whatever else it
/// contains — so the full scan can be skipped for it. Being merely necessary is
/// what keeps this safe: a false positive costs a scan that finds nothing, and
/// a false negative is impossible.
fn has_opener_bytes(bytes: &[u8]) -> bool {
    bytes
        .windows(2)
        .any(|pair| pair[1] == b'"' && (pair[0] == b'r' || pair[0] == b'#'))
}

/// How many consecutive `#` start at `from`.
fn run_of_hashes(bytes: &[u8], from: usize) -> usize {
    bytes[from.min(bytes.len())..]
        .iter()
        .take_while(|b| **b == b'#')
        .count()
}

/// Every function in `text`, with its attributes and its body.
///
/// A hand-rolled scanner rather than a parser: the properties wanted are
/// coarse (does this carry `#[test]`, does its body mention an assertion), and
/// a syntax-tree dependency for that would be a large amount of machinery
/// pointed at a small question.
///
/// Declarations inside a raw string literal are skipped — see
/// [`raw_string_lines`] for why that is load-bearing rather than tidiness.
fn functions_in(text: &str) -> Vec<ParsedFunction> {
    let lines: Vec<&str> = text.lines().collect();
    let in_raw_string = raw_string_lines(&lines);
    let mut found = Vec::new();
    for (index, line) in lines.iter().enumerate() {
        if in_raw_string[index] {
            continue;
        }
        let Some(name) = function_name(line) else {
            continue;
        };
        // Walk back over the contiguous attribute/comment/blank block above.
        //
        // `unclosed` tracks `]` that have not yet met their `[` as the walk moves
        // upward, which is how a **multi-line attribute** is recognised from
        // below. Without it the walk stops on the continuation line, never
        // reaches the `#[test]` above it, and a real test becomes invisible:
        //
        //     #[test]
        //     #[ignore = "one row still open: … see #1627. The \
        //                 other three rows pass"]
        //     fn the_decided_target_is_a_mode_gated_refusal() {
        //
        // That is not hypothetical — it is how this scanner met `main`. #1872
        // reflowed exactly that attribute onto four lines, and the ledger's
        // citation of that test began reporting as unresolved on rebase, while
        // the test sat there running. Measured over the tree at that point:
        // **30** `#[test]` functions were invisible for this reason, and the
        // failure mode is the worst available here, because it is the one this
        // file exists to prevent wearing the opposite sign — a guard that is
        // real and reads as missing today, and (once a citation is "repaired"
        // to something else) a record that resolves while nothing enforces it.
        //
        // Comment lines are excluded from the bracket count deliberately: a doc
        // comment may carry unbalanced brackets in prose or in a Markdown link,
        // and letting those open a phantom attribute would walk the scan up
        // through the *previous* function's body and let its `#[test]` leak onto
        // a plain `fn`. An attribute continuation is never a comment line.
        let mut is_test = false;
        let mut ignored = false;
        let mut should_panic = false;
        let mut above = index;
        let mut unclosed = 0i32;
        while above > 0 {
            let raw = lines[above - 1];
            let candidate = raw.trim_start();
            let is_comment = candidate.starts_with("//");
            let delta = if is_comment || candidate.is_empty() {
                0
            } else {
                raw.matches(']').count() as i32 - raw.matches('[').count() as i32
            };
            // Inside a multi-line attribute either because one was already open
            // below, or because this very line closes one that must open above.
            let inside_attribute = unclosed > 0 || unclosed + delta > 0;
            let is_preamble = inside_attribute
                || candidate.starts_with("#[")
                || is_comment
                || candidate.is_empty();
            if !is_preamble {
                break;
            }
            if candidate.starts_with("#[test") || candidate.starts_with("#[rstest") {
                is_test = true;
            }
            if candidate.starts_with("#[ignore") {
                ignored = true;
            }
            if candidate.starts_with("#[should_panic") {
                should_panic = true;
            }
            unclosed += delta;
            above -= 1;
        }
        found.push(ParsedFunction {
            name,
            line: index + 1,
            is_test,
            ignored,
            should_panic,
            body: body_from(&lines, index),
        });
    }
    found
}

/// The function name on a `fn` declaration line, if it is one.
fn function_name(line: &str) -> Option<String> {
    let rest = line.trim_start();
    let rest = rest.strip_prefix("pub ").unwrap_or(rest);
    let rest = rest.strip_prefix("async ").unwrap_or(rest);
    let rest = rest.strip_prefix("fn ")?;
    let name: String = rest
        .chars()
        .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
        .collect();
    if name.is_empty() || !rest[name.len()..].starts_with(['(', '<']) {
        return None;
    }
    Some(name)
}

/// The text of the function body starting at `start`, by brace matching.
///
/// **The signature is dropped**, everything up to and including the opening
/// brace. That is not tidiness: the assertion scan below is a substring test
/// for `assert`, and a function *named* `asserts_nothing_at_all` would
/// otherwise match on its own declaration line and report itself as asserting —
/// which is the #1858 shape the scan exists to catch, passed by the scan
/// meant to catch it. Caught by
/// [`the_assertion_detector_follows_a_local_helper`]'s negative control before
/// this file was ever run against the ledger.
///
/// Braces inside string literals and comments are not excluded. That is a
/// deliberate floor: the worst it can do is end a body early or late, and both
/// scans over this text are substring tests, so a slightly wrong extent changes
/// an answer only for a function whose assertion sits in the disputed region.
fn body_from(lines: &[&str], start: usize) -> String {
    let mut depth = 0i32;
    let mut opened = false;
    let mut body = String::new();
    for line in &lines[start..] {
        body.push_str(line);
        body.push('\n');
        depth += line.matches('{').count() as i32;
        depth -= line.matches('}').count() as i32;
        if line.contains('{') {
            opened = true;
        }
        if opened && depth <= 0 {
            break;
        }
    }
    match body.find('{') {
        Some(brace) => body.split_off(brace + 1),
        // No brace at all: a trait method or an `fn` declaration with no body.
        // Nothing to read, which is the honest answer rather than the signature.
        None => String::new(),
    }
}

/// Whether `body` asserts anything, reading only the body itself.
///
/// `assert` covers `assert!`, `assert_eq!`, `assert_ne!` and `debug_assert*`,
/// and also catches an assertion made through a helper *named* for it
/// (`assert_normalizes_preserving_in`) — which is the shape a macro-only
/// matcher gets wrong. `expect_err` and `unwrap_err` are how this repo spells
/// "this must fail".
///
/// **`#[should_panic]` is deliberately not read here**, and listing it here is
/// what the earlier revision got wrong: it is an *attribute*, and [`body_from`]
/// strips everything above the opening brace, so this function can never see
/// one. The term matched nothing and the guarantee it advertised was absent —
/// a real `#[should_panic]` test whose body holds no `assert` (there are eight
/// in this tree, e.g. `src/coords/mod.rs::test_one_based_rejects_zero`) was
/// reported as asserting nothing, so citing one as a guard failed
/// [`every_declared_guard_asserts`]. It is read in [`functions_in`] alongside
/// `#[ignore]` and folded in by [`asserts_transitively`] instead, which also
/// lets it propagate through a local helper.
fn asserts_directly(body: &str) -> bool {
    body.contains("assert") || body.contains("expect_err") || body.contains("unwrap_err")
}

/// Whether each function in `text` asserts, following calls to helpers defined
/// in the same file.
///
/// The transitive step is not a refinement, it is required for correctness on
/// this tree. Three of the guards the ledger cites have a body that is one call
/// to a local `converges_to(…)` helper, which is where the assertions live —
/// and a direct-only matcher reads all three as assertionless, reporting real
/// guards as the #1858 shape. Pinned by
/// [`the_assertion_detector_follows_a_local_helper`].
///
/// Resolution is within one file, because that is where it can be done without
/// a name resolver, and it is where this repo's assertion helpers live. A
/// helper imported from `common/` is not followed; such a test is reported as
/// assertionless, and the answer is to name a different guard or to make the
/// assertion visible, not to widen the matcher until it says yes to everything.
fn asserts_transitively(functions: &[ParsedFunction]) -> BTreeMap<String, bool> {
    let mut asserts: BTreeMap<String, bool> = functions
        .iter()
        .map(|f| (f.name.clone(), asserts_directly(&f.body) || f.should_panic))
        .collect();
    let bodies: BTreeMap<&str, &str> = functions
        .iter()
        .map(|f| (f.name.as_str(), f.body.as_str()))
        .collect();

    // Fixpoint. Bounded by the function count, and in practice by one or two
    // passes; the bound is what stops a mutually recursive pair looping.
    for _ in 0..functions.len().min(8) {
        let mut moved = false;
        for (name, body) in &bodies {
            if asserts.get(*name) == Some(&true) {
                continue;
            }
            let reached = bodies.keys().any(|callee| {
                callee != name
                    && asserts.get(*callee) == Some(&true)
                    && body.contains(&format!("{callee}("))
            });
            if reached {
                asserts.insert((*name).to_string(), true);
                moved = true;
            }
        }
        if !moved {
            break;
        }
    }
    asserts
}

// ---------------------------------------------------------------------------
// The checks
// ---------------------------------------------------------------------------

/// The scan is not vacuous, and neither is the ledger side of it.
///
/// Without this, every assertion below passes trivially when the roots move,
/// the `.rs` filter breaks, or the ledger reader stops returning citations —
/// the structural zero this repository's `CLAUDE.md` warns about.
#[test]
fn the_guard_scan_reads_the_tree_and_the_ledger() {
    let records = rulings::records();
    assert!(
        records.len() >= 10,
        "only {} records read from {LEDGER_RELATIVE_PATH}",
        records.len()
    );

    let citations: usize = records.iter().map(|r| r.guard.tests().len()).sum();
    assert!(
        citations >= GUARD_CITATION_FLOOR,
        "only {citations} guard citations in {LEDGER_RELATIVE_PATH}, below the floor of \
         {GUARD_CITATION_FLOOR} — the reader is probably broken, and every check below would \
         pass by resolving nothing"
    );

    let tests = test_functions();
    assert!(
        tests.len() >= SCANNED_TEST_FLOOR,
        "only {} test functions found under {SCAN_ROOTS:?}, below the floor of \
         {SCANNED_TEST_FLOOR} — the walker or the `fn` matcher is broken",
        tests.len()
    );

    // Both arms of the field must be exercised by the ledger, or one of them is
    // dead code that nothing has ever run.
    assert!(
        records.iter().any(|r| matches!(r.guard, Guard::Tests(_))),
        "no record names a guard — the citing arm is untested"
    );
    assert!(
        records
            .iter()
            .any(|r| matches!(r.guard, Guard::Declined(_))),
        "no record declares an exemption — the declining arm is untested, and a first-class \
         answer nobody uses is one nobody has checked works"
    );
}

/// Every record must say what enforces it. **This is where silence is
/// refused.**
///
/// # Why it is here and not at parse
///
/// It used to be at parse: `common::rulings` panicked on a record with no
/// `guard`, so a single silent record took down every test that reads the ledger
/// *and* `generate_spec_fixture`, which meant the generated fixture never
/// existed and everything downstream of it failed too. Measured on this branch
/// rebased onto `origin/main` at `426a944b`, with the two records `main` had
/// just added carrying no `guard`: **50 failed** — 28 on the reader's panic, 20
/// on the generator cascade, 2 on the generator's own precondition tests. Five
/// of the fifty were about guards.
///
/// That is the wrong failure to hand someone. Seven open PRs write ledger
/// records at any one time, so every one of them met that wall on rebase, and
/// what it told them was that the spec fixture was broken. The fix is one line
/// in a JSON file, and nothing in those fifty failures said so.
///
/// So the rejection moved out here, where it is **one** failure that names the
/// records and the two shapes of answer. The rule is unchanged — a record with
/// no `guard` still cannot merge, because this test is not optional — only the
/// blast radius is. Everything that is *malformed* rather than *absent* (both
/// keys set, an empty list, a blank reason, an unparseable citation) still fails
/// at parse, because those are broken records rather than unanswered ones and
/// there is nothing downstream that can do anything useful with them.
///
/// The staging is only sound because absence is representable without being
/// confusable: see [`rulings::Guard::Absent`], which no JSON an author can write
/// produces.
#[test]
fn every_record_declares_a_guard() {
    let records = rulings::records();
    assert!(
        records.len() >= 10,
        "only {} records read from {LEDGER_RELATIVE_PATH} — this check would pass by reading \
         nothing",
        records.len()
    );

    let silent: Vec<&str> = records
        .iter()
        .filter(|record| record.guard.is_absent())
        .map(|record| record.id.as_str())
        .collect();

    assert!(
        silent.is_empty(),
        "these records in {LEDGER_RELATIVE_PATH} carry no `guard`, so the ledger reads as \
         though something enforces them while nothing is named:\n  {}\n\nAdd the field to each. \
         Either name the tests that fail when the ruling stops holding:\n    \
         \"guard\": {{ \"tests\": [\"tests/it/foo.rs::a_guard_name\"] }}\nor declare the \
         exemption, with the reason it is one:\n    \
         \"guard\": {{ \"none\": \"why nothing enforces this\" }}\nDeclining is a first-class \
         answer. Silence is not one: an absent field is indistinguishable from an unconsidered \
         one.",
        silent.join("\n  ")
    );
}

/// Every guard a record names must exist, in the file the record names.
///
/// This is the check the field exists for. A guard renamed out from under a
/// record leaves the ledger reading as enforced while enforcing nothing, and
/// until the field existed there was nothing to check — for most records there
/// was no citation at all.
#[test]
fn every_declared_guard_resolves() {
    let tests = test_functions();
    let mut unresolved: Vec<String> = Vec::new();

    for record in rulings::records() {
        for citation in record.guard.tests() {
            let (path, function) = split_guard_citation(citation)
                .unwrap_or_else(|| panic!("{citation} passed the parser but not the splitter"));
            if !tests.contains_key(&(path.to_string(), function.to_string())) {
                let elsewhere: Vec<&String> = tests
                    .keys()
                    .filter(|(_, name)| name == function)
                    .map(|(file, _)| file)
                    .collect();
                let hint = if elsewhere.is_empty() {
                    format!("no `#[test] fn {function}` exists anywhere under {SCAN_ROOTS:?}")
                } else {
                    format!("a test of that name exists in {elsewhere:?}, but not in {path}")
                };
                unresolved.push(format!("`{}` names `{citation}`: {hint}", record.id));
            }
        }
    }

    assert!(
        unresolved.is_empty(),
        "these guard citations in {LEDGER_RELATIVE_PATH} name a test that cannot enforce the \
         record. Repair the citation — or, if the guard was renamed in a way that changed what \
         it claims, that is an adjudication and belongs in the record:\n  {}",
        unresolved.join("\n  ")
    );
}

/// A record must be guarded by at least one test that actually runs.
///
/// An `#[ignore]`d guard is a legitimate thing to name — this repository uses
/// one to record a deferral, asserting the decided answer with an issue as its
/// un-ignore criterion. What is not legitimate is a record guarded *only* by
/// tests that never run: that is a declaration of enforcement with no
/// enforcement behind it, which is the state the whole field exists to make
/// unrepresentable.
#[test]
fn no_record_is_guarded_only_by_tests_that_never_run() {
    let tests = test_functions();
    let mut dormant: Vec<String> = Vec::new();

    for record in rulings::records() {
        let citations = record.guard.tests();
        if citations.is_empty() {
            continue;
        }
        let running = citations.iter().filter(|citation| {
            split_guard_citation(citation)
                .and_then(|(path, function)| tests.get(&(path.to_string(), function.to_string())))
                .is_some_and(|found| !found.ignored)
        });
        if running.count() == 0 {
            dormant.push(format!("`{}` names only {citations:?}", record.id));
        }
    }

    assert!(
        dormant.is_empty(),
        "these records name guards that are all `#[ignore]`d, so nothing enforces them today. \
         Name a running guard alongside the deferral, or declare `guard.none` with the \
         reason:\n  {}",
        dormant.join("\n  ")
    );
}

/// Every named guard must assert something.
///
/// A test that runs and asserts nothing passes (#1858), and reads as coverage
/// from every angle except the one that matters. Naming such a test as a guard
/// is worse than naming none, because it survives this file's resolution check.
#[test]
fn every_declared_guard_asserts() {
    let tests = test_functions();
    let mut silent: Vec<String> = Vec::new();

    for record in rulings::records() {
        for citation in record.guard.tests() {
            let Some((path, function)) = split_guard_citation(citation) else {
                continue;
            };
            let Some(found) = tests.get(&(path.to_string(), function.to_string())) else {
                continue; // reported by `every_declared_guard_resolves`
            };
            if !found.asserts {
                silent.push(format!(
                    "`{}` names `{citation}` ({path}:{}), whose body asserts nothing",
                    record.id, found.line
                ));
            }
        }
    }

    assert!(
        silent.is_empty(),
        "these guards run and check nothing, so they cannot fail when the ruling stops \
         holding:\n  {}",
        silent.join("\n  ")
    );
}

/// A declared exemption may not smuggle in a dangling citation.
///
/// The reason is prose, so a `path::function` inside it is checked by nothing —
/// "the status quo is pinned by `tests/it/foo.rs::bar`" would read as a
/// complete declaration while naming a test that does not exist. That is the
/// same defect this file exists to close, one level up, and it is the reason
/// `RETIRED_RECORD_IDS`'s notes are checked in `ruling_citation_currency.rs`.
///
/// Naming no test at all is legitimate: several exemptions here are for
/// records that rule on nothing a test could observe.
#[test]
fn an_exemption_reason_names_no_test_that_is_missing() {
    let tests = test_functions();
    let mut dangling: Vec<String> = Vec::new();
    let mut checked = 0usize;

    for record in rulings::records() {
        let Some(reason) = record.guard.declined_reason() else {
            continue;
        };
        for token in backticked_tokens(reason) {
            let Some((path, function)) = split_guard_citation(&token) else {
                continue;
            };
            checked += 1;
            if !tests.contains_key(&(path.to_string(), function.to_string())) {
                dangling.push(format!("`{}`'s exemption names `{token}`", record.id));
            }
        }
    }

    assert!(
        dangling.is_empty(),
        "a `guard.none` reason names a test that does not exist:\n  {}",
        dangling.join("\n  ")
    );

    // Non-vacuity, scoped honestly: this only fires if some exemption names a
    // test, which is a property of the current ledger rather than a rule. The
    // assertion is on the extractor, not on the ledger's content.
    if checked == 0 {
        assert!(
            backticked_tokens("see `tests/it/a.rs::b` for the pin").contains("tests/it/a.rs::b"),
            "no exemption named a test, and the extractor cannot find one in a string that \
             holds one — it is broken, and this check would pass by reading nothing"
        );
    }
}

/// Backticked runs in `prose`.
fn backticked_tokens(prose: &str) -> BTreeSet<String> {
    prose
        .split('`')
        .enumerate()
        .filter(|(index, _)| index % 2 == 1)
        .map(|(_, chunk)| chunk.to_string())
        .collect()
}

// ---------------------------------------------------------------------------
// The scanner's own tests
// ---------------------------------------------------------------------------

/// The assertion detector must follow a call to a helper defined beside the
/// test, or it reports real guards as assertionless.
///
/// Not hypothetical: three of the guards this ledger cites have a body that is
/// exactly one `converges_to(…)` call. A direct-only matcher called all three
/// the #1858 shape — a false positive on the very defect the check is for.
/// Both directions are pinned, because too loose is as bad as too tight: a
/// matcher that said "yes" for any call at all would stop distinguishing an
/// assertionless test from an asserting one.
#[test]
fn the_assertion_detector_follows_a_local_helper() {
    let source = r#"
fn converges_to(a: &str, b: &str) {
    assert_eq!(a, b);
}

fn just_logs(a: &str) {
    println!("{a}");
}

#[test]
fn asserts_through_a_helper() {
    converges_to("x", "x");
}

#[test]
fn asserts_directly_itself() {
    assert!(true);
}

#[test]
fn asserts_nothing_at_all() {
    just_logs("x");
}
"#;
    let verdicts = asserts_transitively(&functions_in(source));
    assert_eq!(verdicts.get("asserts_through_a_helper"), Some(&true));
    assert_eq!(verdicts.get("asserts_directly_itself"), Some(&true));
    assert_eq!(
        verdicts.get("asserts_nothing_at_all"),
        Some(&false),
        "a call to a helper that asserts nothing must not read as an assertion"
    );
    assert_eq!(verdicts.get("just_logs"), Some(&false));
}

/// The `#[test]` and `#[ignore]` attributes are read from the block above the
/// declaration, and a plain function is not mistaken for a test.
#[test]
fn the_scanner_reads_test_and_ignore_attributes() {
    let source = r#"
/// A doc comment, then the attributes.
#[test]
#[ignore = "waiting on #1"]
fn a_deferred_guard() {
    assert!(true);
}

#[test]
fn a_running_guard() {
    assert!(true);
}

fn not_a_test_at_all() {
    assert!(true);
}
"#;
    let parsed = functions_in(source);
    let by_name: BTreeMap<&str, &ParsedFunction> =
        parsed.iter().map(|f| (f.name.as_str(), f)).collect();

    let deferred = by_name["a_deferred_guard"];
    assert!(deferred.is_test && deferred.ignored);
    let running = by_name["a_running_guard"];
    assert!(running.is_test && !running.ignored);
    let plain = by_name["not_a_test_at_all"];
    assert!(!plain.is_test, "a bare fn must not read as a test");
}

/// An attribute spread over several lines must not hide the `#[test]` above it.
///
/// The walk-back that reads attributes runs upward from the `fn`, and a
/// continuation line — the second and later lines of `#[ignore = "… \` — starts
/// with none of `#[`, `//` or nothing at all. A walk that stops there never sees
/// the `#[test]`, so the function is not a test as far as this file is
/// concerned, and a record citing it fails
/// [`every_declared_guard_resolves`] while the test runs perfectly well.
///
/// Not hypothetical, and not caught by review: this scanner was written against
/// a tree where `corpus_prohibited_inputs.rs::the_decided_target_is_a_mode_gated_refusal`
/// carried a one-line `#[ignore]`. #1872 reflowed it onto four lines, and on the
/// next rebase the ledger's citation of it reported as naming a test that "does
/// not exist anywhere". **30** tests in the tree were invisible for this reason
/// at that point.
///
/// Both directions are pinned. The negative control is the one that matters: the
/// mechanism keys on an unmatched `]`, so it must not run away up the file and
/// collect the *previous* function's attributes onto a plain `fn`.
#[test]
fn a_multi_line_attribute_does_not_hide_the_test_attribute() {
    let source = r#"
#[test]
#[ignore = "one row still open: see #1627. The \
            other three rows pass (#1627's delinsX), \
            as does checklist.md:20"]
fn a_deferred_guard_with_a_reflowed_ignore() {
    assert!(true);
}

#[test]
#[should_panic(
    expected = "a message long enough that rustfmt \
                puts it on its own line"
)]
fn a_negative_guard_with_a_wrapped_should_panic() {
    reject("nonsense");
}

fn a_plain_helper_after_all_that() {
    assert!(true);
}
"#;
    let parsed = functions_in(source);
    let by_name: BTreeMap<&str, &ParsedFunction> =
        parsed.iter().map(|f| (f.name.as_str(), f)).collect();

    let deferred = by_name["a_deferred_guard_with_a_reflowed_ignore"];
    assert!(
        deferred.is_test,
        "a `#[test]` above a reflowed `#[ignore]` must still register"
    );
    assert!(
        deferred.ignored,
        "the `#[ignore]` itself must still be read, or a dormant guard reads as running"
    );

    let negative = by_name["a_negative_guard_with_a_wrapped_should_panic"];
    assert!(negative.is_test);
    assert!(
        negative.should_panic,
        "`#[should_panic]` spread over lines is still the assertion"
    );

    // The negative control: the bracket tracking must not walk past a real body
    // and hand this helper the `#[test]` belonging to the function above it.
    let plain = by_name["a_plain_helper_after_all_that"];
    assert!(
        !plain.is_test,
        "a bare `fn` below a test body must not inherit that test's attributes"
    );
}

/// A doc comment carrying an unbalanced bracket must not open a phantom
/// attribute.
///
/// This is the failure mode the comment exclusion in [`functions_in`] exists to
/// prevent, and it is worth its own control because it is invisible from the
/// positive test above: prose and Markdown links routinely carry `]` without
/// `[`, and treating one as an attribute continuation would let the walk-back
/// run up through the previous function's body.
#[test]
fn an_unbalanced_bracket_in_prose_does_not_open_an_attribute() {
    let source = r#"
#[test]
fn a_real_test_whose_body_has_brackets() {
    let v = vec![1, 2, 3];
    assert_eq!(v[0], 1);
}

/// A doc comment whose prose closes a bracket it never opened: see 3] above.
fn not_a_test_despite_the_bracket() {
    assert!(true);
}
"#;
    let parsed = functions_in(source);
    let by_name: BTreeMap<&str, &ParsedFunction> =
        parsed.iter().map(|f| (f.name.as_str(), f)).collect();

    assert!(by_name["a_real_test_whose_body_has_brackets"].is_test);
    assert!(
        !by_name["not_a_test_despite_the_bracket"].is_test,
        "an unbalanced `]` in a doc comment must not make the walk-back climb \
         into the previous function and inherit its `#[test]`"
    );
}

/// A `fn` inside a raw string literal is source *data*, and must not register
/// as a test.
///
/// This is the negative control for [`raw_string_lines`], and the defect it
/// pins was live: citing
/// `tests/it/ruling_guard_field.rs::asserts_directly_itself` — one of this
/// file's own fixtures, which exists only inside an `r#"…"#` literal — made
/// every check in this file pass on a record that nothing enforces. Resolution
/// succeeded because the scan could not tell code from a string.
#[test]
fn a_test_declared_inside_a_raw_string_literal_does_not_register() {
    let source = r##"
#[test]
fn a_real_test() {
    assert!(true);
}

fn holds_a_fixture() {
    let fixture = r#"
#[test]
fn a_phantom_test() {
    assert!(true);
}
"#;
    drop(fixture);
}
"##;
    let names: Vec<String> = functions_in(source)
        .iter()
        .filter(|f| f.is_test)
        .map(|f| f.name.clone())
        .collect();
    assert!(
        names.iter().any(|n| n == "a_real_test"),
        "a real `#[test]` must still register; found {names:?}"
    );
    assert!(
        !names.iter().any(|n| n == "a_phantom_test"),
        "a `#[test] fn` inside a raw string literal is data, not a test — citing one as a guard \
         would resolve while nothing enforced the record; found {names:?}"
    );
}

/// `#[should_panic]` is an assertion, and it lives where [`body_from`] cannot
/// see it.
///
/// The negative control matters as much as the positive one: if the attribute
/// scan were widened to "any attribute counts", an assertionless test would
/// stop being distinguishable from a guarded one, which is the #1858 shape this
/// file exists to catch.
#[test]
fn a_should_panic_test_counts_as_asserting() {
    let source = r#"
#[test]
#[should_panic(expected = "boom")]
fn a_negative_guard() {
    reject("nonsense");
}

#[test]
#[ignore = "waiting on #1"]
fn a_silent_test() {
    reject("nonsense");
}
"#;
    let verdicts = asserts_transitively(&functions_in(source));
    assert_eq!(
        verdicts.get("a_negative_guard"),
        Some(&true),
        "`#[should_panic]` IS the assertion; the body holds none"
    );
    assert_eq!(
        verdicts.get("a_silent_test"),
        Some(&false),
        "a test that neither asserts nor declares a panic must still read as silent"
    );
}

/// The generator keeps its own copy of the citation grammar, because the two
/// crates do not link.
///
/// **What this pins is narrower than "the copies agree", and saying otherwise
/// would be the defect this whole file is about.** It cannot execute the
/// generator's copy — that binary does not link here — so it checks that the
/// generator still *carries* a splitter, and separately exercises the reader's
/// copy on the cases that discriminate the grammar. Agreement itself is held by
/// the generator's own unit test
/// (`examples/generate_spec_fixture.rs::the_guard_citation_grammar_admits_a_path_and_nothing_else`),
/// which runs the same case list against the other copy because that target
/// sets `test = true`. Two tests over one case list, not one test over two
/// copies.
#[test]
fn the_generator_and_the_ledger_reader_agree_on_the_citation_grammar() {
    let generator = std::fs::read_to_string(crate_root().join("examples/generate_spec_fixture.rs"))
        .expect("read the generator");
    assert!(
        generator.contains("fn split_guard_citation"),
        "the generator no longer carries a citation splitter; if the grammar moved somewhere \
         shared, delete this test rather than leaving it asserting about a copy that is gone"
    );
    // The discriminating cases, spelled out so the assertion is about behaviour
    // rather than about the source text matching itself.
    for (token, accepted) in [
        ("tests/it/foo.rs::a_guard_name", true),
        ("merge::a_guard_name", false),
        ("tests/it/foo.rs", false),
        ("tests/it/foo.rs::AGuard", false),
    ] {
        assert_eq!(
            split_guard_citation(token).is_some(),
            accepted,
            "{token} must {} parse",
            if accepted { "" } else { "not" }
        );
    }
}
