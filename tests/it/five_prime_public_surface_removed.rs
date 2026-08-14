//! The 5' shuffling direction is not reachable from any ferro entry point.
//!
//! `README.md` rule 6 states that there are no user options for normalization
//! form, and then concedes in its next sentence that the 3'/5' knob is not
//! orthogonal — it selects the frame every other rule is evaluated in. The knob
//! was therefore a user option for normalization form sitting inside the rule
//! forbidding user options for normalization form. Removing it from the public
//! surface is what makes rule 6 true.
//!
//! [`ShuffleDirection::FivePrime`] itself is **kept**, as an internal
//! differential oracle. 3'/5' disagreement has twice found a defect in the
//! shipped 3' output that nothing else could see (#1542 / PR #1840, where 7 of
//! 8 `FERRO_PARTITION` x direction configurations agreed and only `live`/3'
//! diverged; and `cis_confluence_axis.rs`'s `enclosing_exon` off-by-one, a
//! 5'-only symptom of a bug in shared code). Deleting the arm would delete that
//! detector, so ~75 test modules keep driving it — this module is about
//! *reachability*, never about the behaviour.
//!
//! **The failure shape this module exists to prevent is a silent one.** A
//! removal that merely stopped honouring `--direction 5prime`, or that let PyO3
//! ignore an unknown keyword, would hand a caller who explicitly asked for 5' a
//! 3' answer with nothing to tell them — a wrong description delivered under an
//! argument that appears to have been accepted. Every assertion below is that
//! the removed surface *refuses* rather than defaults.

use std::process::Command;

use ferro_hgvs::normalize::ShuffleDirection;
use ferro_hgvs::{from_sequences, FromSequencesOptions};

/// Run the `ferro` binary cargo built for this test run.
///
/// `CARGO_BIN_EXE_<name>` is set by cargo for a `[[bin]]` target in the same
/// package, so this cannot pick up a stale `ferro` from `PATH` — which would
/// make the whole module a test of whatever is installed on the machine.
fn ferro() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ferro"))
}

/// A variant the bundled mock provider normalizes cleanly under the default
/// `--error-mode strict`, so a non-zero exit below is about `--direction` and
/// not about the variant. (`c.100delA`, the obvious choice, is refused for
/// spelling its own deleted base — a real refusal that would mask this one.)
const VARIANT: &str = "NM_000088.3:c.4C>G";

/// `--direction` is rejected, with a non-zero exit.
///
/// This is the single most important assertion in the file. Before the flag was
/// removed, `parse_shuffle_direction`'s fallback arm was `_ => ThreePrime`
/// (#1863), so an unrecognized value 3'-shifted and reported success. Deleting
/// the *value* handling without deleting the *flag* would have kept exactly
/// that behaviour for `--direction 5prime`. Deleting the flag makes clap refuse
/// the argument instead.
#[test]
fn the_direction_flag_is_rejected_rather_than_ignored() {
    for value in ["5prime", "3prime", "5'", "banana"] {
        let output = ferro()
            .args(["normalize", VARIANT, "--direction", value])
            .output()
            .expect("run ferro normalize");

        assert!(
            !output.status.success(),
            "`--direction {value}` must fail; a silent 3' run for a caller who \
             asked for 5' is the failure this test exists to prevent. \
             status={:?} stdout={}",
            output.status,
            String::from_utf8_lossy(&output.stdout),
        );

        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("--direction") || stderr.contains("unexpected argument"),
            "the refusal must name the offending argument, got: {stderr}"
        );
    }
}

/// `--direction` is gone from the help text as well as from the parser.
///
/// A flag that is refused but still advertised would send a reader straight
/// back to the argument they cannot use.
#[test]
fn the_help_text_does_not_advertise_a_direction_flag() {
    let output = ferro()
        .args(["normalize", "--help"])
        .output()
        .expect("run ferro normalize --help");
    assert!(output.status.success(), "--help must succeed");

    let help = String::from_utf8_lossy(&output.stdout);
    assert!(
        !help.contains("--direction"),
        "`ferro normalize --help` must not advertise --direction:\n{help}"
    );
    assert!(
        !help.contains("5prime"),
        "`ferro normalize --help` must not mention a 5' direction:\n{help}"
    );
}

/// `ferro normalize` still works, and shifts 3'.
///
/// The negative tests above would all pass against a binary that had lost the
/// subcommand entirely, so this is the control: removing the flag must not have
/// removed the run, and the direction it did not remove is 3'.
#[test]
fn normalize_still_runs_and_is_three_prime() {
    let output = ferro()
        .args(["normalize", VARIANT, "-f", "text"])
        .output()
        .expect("run ferro normalize");
    assert!(
        output.status.success(),
        "plain `ferro normalize` must succeed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&output.stdout).trim(),
        VARIANT,
        "the shipped 3' path must still answer"
    );
}

/// The internal direction still *acts* — the oracle has not been hollowed out.
///
/// Deleting a public knob is easy to do in a way that also stops the underlying
/// machinery ever seeing anything but the default; that would leave ~75 test
/// modules running a differential against themselves and passing. This drives
/// [`FromSequencesOptions::with_direction`] on a pair whose two answers differ
/// and asserts they still differ.
///
/// The rows are the ones the removed Python tests asserted, kept verbatim so
/// the measurement is relocated rather than dropped: a single-base deletion
/// inside an `AAAA` run in a window starting at 10, which 3' places at 14 and
/// 5' at 11 (`tests/python/test_from_sequences.py`), and the `SequencePair`
/// window that gave `g.15del` / `g.12del` for `derive`.
#[test]
fn the_internal_direction_still_moves_the_answer() {
    let three = from_sequences(
        "NC_1",
        10,
        "CAAAAG",
        "CAAAG",
        &FromSequencesOptions::default().with_direction(ShuffleDirection::ThreePrime),
    )
    .expect("3' derivation");
    let five = from_sequences(
        "NC_1",
        10,
        "CAAAAG",
        "CAAAG",
        &FromSequencesOptions::default().with_direction(ShuffleDirection::FivePrime),
    )
    .expect("5' derivation");

    assert_eq!(format!("{three}"), "NC_1:g.14del");
    assert_eq!(format!("{five}"), "NC_1:g.11del");
    assert_ne!(
        format!("{three}"),
        format!("{five}"),
        "the internal 3'/5' differential must still differ, or every census \
         driving it is comparing an arm against itself"
    );
}

/// The default is 3', so a caller who names no direction gets the shipped form.
#[test]
fn the_default_direction_is_three_prime() {
    assert_eq!(ShuffleDirection::default(), ShuffleDirection::ThreePrime);

    let defaulted = from_sequences(
        "NC_1",
        10,
        "CAAAAG",
        "CAAAG",
        &FromSequencesOptions::default(),
    )
    .expect("default derivation");
    assert_eq!(format!("{defaulted}"), "NC_1:g.14del");
}
