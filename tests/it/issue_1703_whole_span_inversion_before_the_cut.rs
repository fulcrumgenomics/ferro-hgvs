//! Issue #1703: a whole-span reverse complement must be **typed** before the
//! block is **cut**.
//!
//! The defect
//! ----------
//! On every sequence-first partition arm (`FERRO_PARTITION=shadow`, `canonical`,
//! `canonical-coalesced`) a 5 nt span whose replacement is its own reverse
//! complement came back as a `[del;dup]` pair rather than an `inv`. Every
//! spelling reached the same wrong string, so the arm was self-consistent and
//! wrong together and no confluence check could see it.
//!
//! Real GRCh38 reproducer, which the synthetic core below mirrors base for base:
//! `NC_000017.11:g.43044327_43044331` is `GTTAA`, `revcomp(GTTAA) = TTAAC`, and
//! columns `43044328`/`43044330` are self-complementary — so the substitution
//! spelling is three *separated* substitutions.
//!
//! ```text
//! g.        43044327  43044328  43044329  43044330  43044331
//! ref          G         T         T         A         A
//!              X         |         X         |         X
//! alt          T         T         A         A         C
//! ```
//!
//! `FERRO_PARTITION=live` emitted `g.43044327_43044331inv` for all three
//! spellings; the sequence-first arms emitted `g.[43044327del;43044332dup]`.
//! Mutalyzer converges all three on the `inv`.
//!
//! Why the wrong form is not admissible
//! ------------------------------------
//! `DNA/inversion.md:5` defines an inversion by the content of the **whole**
//! span — "**more than one nucleotide** replacing the original sequence is the
//! reverse complement of the original sequence" — and this span is an exact
//! reverse complement. `DNA/delins.md:5` independently excludes the spanning
//! `delins` spelling by definition, "**and which is not** a substitution or
//! inversion". Both are prohibitions read from prose force, which is rule 1's
//! scope in `README.md`, so this is a bug and not a disclosure.
//!
//! What the guards below are for
//! -----------------------------
//! `issue_1040_inv_overrecognition_probes` already pins the same 5 nt core on
//! the **`live`** arm, in process, through `assert_padded_preserving`. It cannot
//! reach an arm: `partition_rule()` caches its read of `FERRO_PARTITION` in a
//! `OnceLock`, and `tests/it/common/cis_apply_oracle.rs` records why a test must
//! not `set_var` to get around that (sound under nextest, false under `cargo
//! test`, which runs tests as threads in one process). So the arm is reached the
//! way `partition_switch_wiring.rs` reaches it — through the `ferro` binary, with
//! the variable in the child's environment.
//!
//! The two negative controls are the point of the file as much as the positive
//! case is, because the fix widens an admission gate:
//!
//! * `AAGCTA -> TAGCTT` is **also** a whole-span reverse complement, but its
//!   members are two lone substitutions, and `general.md:56` ranks substitution
//!   (1) above inversion (3). It must stay split — this is the #1040 control.
//! * `GTTAT -> TTATC` has the **identical geometry** to the positive case (a
//!   1-base deletion at the hull start plus a 1-base insertion at the hull end,
//!   four unchanged bases between them) and is *not* a reverse complement. It
//!   must stay split, which is what makes the positive case a statement about
//!   the sequence rather than about the shape.

use crate::common::cis_apply_oracle::apply;
use std::collections::BTreeSet;
use std::io::Write;
use std::path::Path;
use std::process::Command;

/// The contig every fixture below serves, and the accession
/// [`crate::common::cis_apply_oracle`] builds its own provider under — so the
/// same descriptor can be applied in process to cross-check what the child
/// process printed.
const CONTIG: &str = "TEMPLATE";

/// 0-based offset of the core inside a fixture sequence, i.e. the core's first
/// 1-based HGVS position is `PREFIX_LEN + 1` = 81.
const PREFIX_LEN: usize = 80;

/// Every arm `FERRO_PARTITION` accepts. The defect was present on the three
/// sequence-first arms and absent from `live`; the fix has to leave `live` where
/// it is, so it is asserted here rather than merely believed.
const ARMS: [&str; 4] = ["live", "shadow", "canonical", "canonical-coalesced"];

/// The reference a fixture serves: `ACGT…` padding, the core, a single `C`, then
/// more padding.
///
/// The trailing `C` is not decoration. It is what makes the 3'-shift render the
/// derived insertion as `86dup` rather than `85_86insC`, so the synthetic output
/// is the same *shape* as the real locus's `g.[43044327del;43044332dup]` — where
/// `NC_000017.11:g.43044332` is likewise a `C`.
fn fixture_sequence(core: &str) -> String {
    format!("{pad}{core}C{pad}", pad = "ACGT".repeat(PREFIX_LEN / 4))
}

/// Write `core`'s fixture into `dir` as an indexed single-contig FASTA.
///
/// `MultiFastaProvider::from_directory` — which `ferro normalize --reference`
/// takes when the directory holds no `manifest.json` — requires a `.fai`, so one
/// is written alongside. The index is exact rather than approximate: a
/// single-record FASTA wrapped at a fixed width has a closed-form index, and
/// getting it wrong would surface as a wrong base rather than as an error.
fn write_fixture(dir: &Path, core: &str) -> String {
    const LINE_BASES: usize = 60;
    let sequence = fixture_sequence(core);
    let header = format!(">{CONTIG}\n");

    let mut fasta = std::fs::File::create(dir.join("t.fa")).expect("create fasta");
    fasta.write_all(header.as_bytes()).expect("write header");
    for line in sequence.as_bytes().chunks(LINE_BASES) {
        fasta.write_all(line).expect("write bases");
        fasta.write_all(b"\n").expect("write newline");
    }
    fasta.flush().expect("flush fasta");

    std::fs::write(
        dir.join("t.fa.fai"),
        format!(
            "{CONTIG}\t{}\t{}\t{LINE_BASES}\t{}\n",
            sequence.len(),
            header.len(),
            LINE_BASES + 1,
        ),
    )
    .expect("write fai");

    sequence
}

/// Normalize `variant` against `core`'s fixture with `FERRO_PARTITION=arm`, and
/// return the normalized description.
///
/// `ferro normalize` prints either the description (when nothing changed) or
/// `<input> -> <output>`, so the output is whatever follows the last arrow.
fn normalize_on_arm(core: &str, arm: &str, variant: &str) -> String {
    let dir = tempfile::tempdir().expect("temp dir");
    write_fixture(dir.path(), core);

    let run = Command::new(env!("CARGO_BIN_EXE_ferro"))
        .arg("normalize")
        .arg(variant)
        .arg("--reference")
        .arg(dir.path())
        .env("FERRO_PARTITION", arm)
        .output()
        .expect("run ferro normalize");
    assert!(
        run.status.success(),
        "`{variant}` on arm `{arm}` failed; stderr was:\n{}",
        String::from_utf8_lossy(&run.stderr),
    );

    let stdout = String::from_utf8(run.stdout).expect("utf-8 stdout");
    let last = stdout
        .lines()
        .next_back()
        .unwrap_or_else(|| panic!("`{variant}` on arm `{arm}` printed nothing"))
        .trim();
    last.rsplit(" -> ")
        .next()
        .expect("a description")
        .to_string()
}

/// Assert that every `spelling` of one variant reaches `expected` on every arm,
/// and that each output still denotes the bases its input did.
///
/// Convergence is a property of the *set* of outputs, so the set is built and
/// its cardinality asserted rather than each spelling being compared against a
/// constant in isolation — a single-input assertion cannot observe convergence,
/// and three of them agreeing with a constant is not the same statement as the
/// three agreeing with each other.
///
/// The denoted-sequence half goes through
/// [`crate::common::cis_apply_oracle::apply`], which reaches the bases through
/// `hgvs_to_spdi` and an SPDI splice rather than through the normalizer — so a
/// re-typing that quietly changed the sequence fails here loudly instead of
/// merely disagreeing with a pinned string. This is the same discipline
/// `assert_padded_preserving` applies on the in-process `live` path, expressed
/// against a child process's output.
fn assert_converges(core: &str, spellings: &[String], expected: &str) {
    let reference = fixture_sequence(core);
    for arm in ARMS {
        let mut outputs = BTreeSet::new();
        for spelling in spellings {
            let output = normalize_on_arm(core, arm, spelling);
            let denoted_in = apply(&reference, spelling)
                .unwrap_or_else(|| panic!("`{spelling}` denotes no sequence"));
            let denoted_out = apply(&reference, &output).unwrap_or_else(|| {
                panic!("`{spelling}` -> `{output}` on arm `{arm}` denotes no sequence")
            });
            assert_eq!(
                denoted_out, denoted_in,
                "`{spelling}` -> `{output}` on arm `{arm}` changed the sequence",
            );
            outputs.insert(output);
        }
        assert_eq!(
            outputs.len(),
            1,
            "arm `{arm}` gave {} answers for one variant: {outputs:?}",
            outputs.len(),
        );
        assert_eq!(
            outputs.iter().next().map(String::as_str),
            Some(expected),
            "arm `{arm}` converged on the wrong form",
        );
    }
}

/// The core the real locus carries, and the one
/// `issue_1040_inv_overrecognition_probes::
/// every_spelling_of_a_derived_whole_block_inversion_converges_on_inv` already
/// pins on `live`.
const INVERSION_CORE: &str = "GTTAA";

#[test]
fn every_spelling_of_a_whole_span_reverse_complement_converges_on_inv_on_every_arm() {
    // `NC_000017.11:g.43044327_43044331inv`, its three-substitution spelling and
    // its spanning `delins` spelling, re-sited at position 81 of the fixture.
    //
    // The `delins` spelling is included even though `DNA/delins.md:5` makes it an
    // incorrect *output*: ferro accepts it as input, and must land on the same
    // normal form as the spellings that are correct.
    assert_converges(
        INVERSION_CORE,
        &[
            format!("{CONTIG}:g.81_85inv"),
            format!("{CONTIG}:g.[81G>T;83T>A;85A>C]"),
            format!("{CONTIG}:g.81_85delinsTTAAC"),
        ],
        &format!("{CONTIG}:g.81_85inv"),
    );
}

#[test]
fn a_reverse_complement_whose_members_are_substitutions_still_splits_on_every_arm() {
    // The #1040 control, and the assertion that stops the fix from becoming "merge
    // any whole-span reverse complement". `AAGCTA -> TAGCTT` is an exact reverse
    // complement whose only changed columns are its first and last, four unchanged
    // bases apart, and whose members are therefore two **lone substitutions**.
    //
    // `general.md:56` ranks substitution (1) above inversion (3), so the split
    // wins here and the spanning `inv` must not be reached. That ranking is what
    // `whole_block_inversion`'s `no_piece_is_a_lone_substitution` route encodes,
    // and widening that route is exactly what this issue's fix does — so this
    // control is load-bearing rather than decorative.
    assert_converges(
        "AAGCTA",
        &[
            format!("{CONTIG}:g.81_86inv"),
            format!("{CONTIG}:g.[81A>T;86A>T]"),
            format!("{CONTIG}:g.81_86delinsTAGCTT"),
        ],
        &format!("{CONTIG}:g.[81A>T;86A>T]"),
    );
}

#[test]
fn a_del_ins_block_that_is_not_a_reverse_complement_still_splits_on_the_arms() {
    // The geometric control. `GTTAT -> TTATC` partitions on the canonical arms
    // into precisely the same shape as the positive case — a 1-base deletion at
    // the hull start and a 1-base insertion at the hull end, four unchanged bases
    // between them, rendering as `g.[81del;86dup]` — and it is **not** a reverse
    // complement (`revcomp(GTTAT) = ATAAC`, not `TTATC`).
    //
    // So the two cases are separated by the sequence relation and by nothing else.
    // Without this, the positive test above would pass just as well against a rule
    // that merged every compensating deletion/insertion pair.
    //
    // Asserted per arm rather than through `assert_converges`, because the answer
    // legitimately differs between the arms: `live`'s single-gap search never
    // proposes the compensating pair at all.
    let spelling = format!("{CONTIG}:g.81_85delinsTTATC");
    assert_eq!(
        normalize_on_arm("GTTAT", "live", &spelling),
        format!("{CONTIG}:g.[81G>T;83_85delinsATC]"),
    );
    for arm in ["shadow", "canonical", "canonical-coalesced"] {
        assert_eq!(
            normalize_on_arm("GTTAT", arm, &spelling),
            format!("{CONTIG}:g.[81del;86dup]"),
            "arm `{arm}` merged a block that is not a reverse complement",
        );
    }
}

#[test]
fn the_wrong_form_denoted_the_right_bases_which_is_why_no_oracle_reached_it() {
    // Why this needed an issue rather than being caught by CI. The form the arms
    // emitted is a *mis-typing*, not a mis-application: it denotes exactly the
    // bases the `inv` does, it re-parses, its coordinates are in bounds, and it is
    // a fixed point. So all four seam oracles — `FERRO_ASSERT_IDEMPOTENT`,
    // `FERRO_ASSERT_REPARSE`, `FERRO_ASSERT_IN_BOUNDS` and the strongest of them,
    // `FERRO_ASSERT_SEQUENCE` — pass on it, and every spelling reached it, so the
    // confluence checks passed too.
    //
    // Pinned by applying both descriptions independently of the normalizer, so
    // this assertion stays true (and stays the explanation) whatever the
    // normalizer later does.
    let reference = fixture_sequence(INVERSION_CORE);
    let inversion = apply(&reference, &format!("{CONTIG}:g.81_85inv")).expect("the inv applies");
    let del_dup =
        apply(&reference, &format!("{CONTIG}:g.[81del;86dup]")).expect("the pair applies");
    assert_eq!(
        inversion, del_dup,
        "the two forms must denote one sequence — otherwise this was never a \
         typing question",
    );
    assert_ne!(inversion, reference, "the fixture must actually change");
}
