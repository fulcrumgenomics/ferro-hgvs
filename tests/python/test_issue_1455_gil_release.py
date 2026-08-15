"""The reference-backed bindings must not hold the GIL (#1455).

``Normalizer.to_spdi``, ``.canonical_spdi`` and ``.apply_to_reference`` read
reference sequence through a provider — real file I/O for a FASTA- or mmap-backed
one, and ``canonical_spdi`` issues several reads as its window doubles through a
repeat tract. Running that with the GIL held serialises every other Python thread
in the process for the duration. ``.normalize``, ``.normalize_variant`` and
``.normalize_with_warnings`` reach the reference through the same provider, and
``EquivalenceChecker.check`` / ``.all_equivalent`` settle a verdict by
normalizing — and at the strongest tier applying — both descriptions against it.
All eight held the GIL; all eight now release it.

Two neighbours are deliberately **not** changed, having been checked rather than
assumed. ``VariantProjector`` and ``BatchProcessor`` already detach on every
reference-backed entry point. ``CoordinateMapper``'s methods touch the provider
for a single ``get_transcript`` metadata lookup and then do arithmetic on the
returned record — no sequence is read, so there is no I/O to overlap.

How this is measured
--------------------

A pure-Python **observer** thread counts iterations of a tight loop, twice: once
running alone, and once alongside a worker thread that hammers the binding. The
observable is the ratio of the two counts — the share of its solo rate the
observer keeps while the binding is running.

* GIL **held** across the call: the observer cannot execute a single bytecode
  until the call returns, so it is starved down to the slivers between calls.
* GIL **released**: the observer runs concurrently and keeps ~all of its rate.

This replaced an earlier wall-clock design — total time on one thread over the
same total on two — and the replacement was not cosmetic. That version had a
2x signal at best (~1.0 held, ~0.5 released) and flaked 4 runs in 10 of the full
suite on a build already proven to detach. The observer share separates by ~5x
with far tighter spreads, because it is a ratio of two counts over one fixed
wall-clock window rather than a comparison of two separately-timed phases.

``sys.setswitchinterval`` matters and is set here on purpose. At CPython's 5 ms
default the switch interval dwarfs a ~100 µs call, so a GIL-holding worker only
gets a few percent of the time and the observer sails on at ~0.5 either way — the
measurement simply does not discriminate. Dropped to 10 µs, well under the call
duration, a held GIL starves the observer as intended.

**Every threshold is measured on the host, not hard-coded.** Two calibration
workers run alongside the real ones: a pure-Python loop, which shares the GIL
*fairly*, and a Rust call known to *hold* it (rendering a ~1 Mb ``SpdiVariant``
to a Python ``str``). Fair sharing is the midpoint the question turns on — a call
that released the GIL scores above it, one that held the GIL scores below — so
the assertion is "beat fair sharing by 1.5x" rather than any absolute number.

That matters, and not hypothetically. An earlier revision of this file compared
against an absolute figure calibrated on a 12-core laptop; on CI's 4-vCPU runner
it skipped all eight tests, so the job was green and guarding nothing.

Measured on a 12-core release wheel, Python 3.14, on a quiet host:

=====================================  ==========  ==========
worker                                    unfixed       fixed
=====================================  ==========  ==========
Rust call known to hold the GIL              0.23        0.23
pure-Python worker ("fair")                  0.36        0.36
the eight entry points                  0.13-0.30   1.21-1.49
=====================================  ==========  ==========

The populations sit either side of fair sharing, which is what the threshold
keys on.

**This is still a timing measurement, and it has a stated limit.** A contended
host deschedules the *worker* mid-call, handing the observer time a GIL-holding
call would never have given it, and the two populations converge — at load
average 30 on this 12-core host the held-GIL calibration rose to 0.311 while the
weakest entry point fell to 0.349. Two guards detect that and **skip**, because a
measurement that cannot discriminate should not be reported as a failure: see
``MAX_LOAD_PER_CORE`` and ``MAX_CONTROL_OVER_FAIR``.

What every run emits
--------------------

Two kinds of warning, so a run on a machine you cannot log into is diagnosable
from the CI log alone:

* one module-scoped ``#1455 GIL calibration:`` line, carrying ``fair`` and
  ``held``; and
* one ``#1455 GIL margin:`` row **per parametrized entry point** — ``MEASURED``
  with the observed ratios, or ``STOOD-DOWN`` naming the guard that fired.

The second is the point of the exercise, and the invariant it rests on is that
the two shapes are exhaustive: a case emits one row or the other, never neither.
So a parametrized case with no row did not run at all — a stand-down is never
silent, and cannot be mistaken for a healthy quiet run. ``_margin_line``'s own
docstring says which of its fields may be pooled across runs and which may not.
"""

import json
import os
import statistics
import sys
import tempfile
import threading
import time
import warnings
from pathlib import Path

import pytest

import ferro_hgvs

# Long enough to hold `WIDE_DESCRIPTOR` below. The capped entry points never read
# more than `MAX_APPLY_WINDOW` (100 000 bases) of it however long it is.
SEQUENCE_LENGTH = 1_200_000

# `(i * 7 + (i // 4) * 3) % 4` has period 16, so the 16-base period is built once
# and repeated rather than evaluated 1.2 million times — the same string, far
# cheaper at import. `test_sequence_matches_its_generator` pins the equivalence,
# so the shortcut cannot silently become a different fixture.
SEQUENCE_PERIOD = "".join("ACGT"[(i * 7 + (i // 4) * 3) % 4] for i in range(16))
SEQUENCE = SEQUENCE_PERIOD * (SEQUENCE_LENGTH // 16)

# A ~98 kb deletion, just inside `MAX_APPLY_WINDOW`. Every entry point except
# `to_spdi` is capped at that window however wide the description names, so this
# is as much reference work as they can be given.
DESCRIPTOR = "TEMPLATE:n.1000_99000del"

# `to_spdi` is a transliteration and carries no window cap, so unlike the others
# it is not pinned at ~100 kb of work — and on `DESCRIPTOR` it costs only ~7 µs,
# small enough that the GIL hand-off is a visible fraction of it. Given a ~1 Mb
# deletion it costs ~65 µs and behaves like the rest.
WIDE_DESCRIPTOR = "TEMPLATE:n.1000_999000del"

# The same deletion as `DESCRIPTOR` spelled as two adjacent members, so comparing
# the two makes `EquivalenceChecker` do reference work rather than answer from a
# cheap syntactic tier.
DECOMPOSED = "TEMPLATE:n.[1000_50000del;50001_99000del]"

# Must sit well below the duration of one call (~65-300 µs here). See the module
# docstring: at CPython's 5 ms default this measurement does not discriminate at
# all. Restored after the module by `fine_grained_gil`.
SWITCH_INTERVAL = 1e-5

# One observation window. The observer is timed for this long alone and this long
# alongside the worker, so a case costs 2x this per round.
WINDOW_SECONDS = 0.08

# Rounds per measurement and attempts per case. Three rounds is the smallest a
# median can reject an outlier from. A retry cannot manufacture a false pass — a
# GIL-holding binding sits at the control's own share, nowhere near the threshold
# — but it does turn a transient CPU spike into a slower success.
ROUNDS = 3
# Raised 2 -> 3 after a merge-group run ejected a PR on a 1.1% miss (share 0.671
# against a 0.678 threshold) while a PR-level run on the SAME head passed sixteen
# minutes earlier, with a diff that compiles into nothing in the wheel. That is
# the transient-spike case this constant exists for, and two attempts did not
# absorb it. The reasoning above is unchanged and bounds the cost: a retry cannot
# manufacture a false pass, so the only price is a slower success on a noisy host.
ATTEMPTS = 3

# The negative control is measured once for the whole module and every threshold
# is derived from it, so it gets more rounds than a single case does: one unlucky
# control skips all eight tests rather than one, which was observed once in ten
# full-suite runs at `ROUNDS`.
CONTROL_ROUNDS = 5

# The entry point must keep at least this multiple of a **pure-Python worker's**
# share. That reference point is the whole design, so it is worth saying why.
#
# There are three regimes, and a pure-Python worker sits exactly between the two
# that matter:
#
#   worker holds the GIL in Rust  -> observer frozen for each call    -> < fair
#   worker is pure Python         -> the two share the GIL fairly     -> "fair"
#   worker released the GIL       -> observer runs concurrently       -> > fair
#
# So "did it release the GIL" is "is it above or below fair sharing", and fair
# sharing is measured on the machine under test rather than assumed. Idle here
# that is 1.21-1.49 against ~0.36 (3.4-4.1x) after the fix and 0.13-0.30 (0.4-0.8x)
# before, so the populations sit either side of 1.0 and this threshold splits them
# with margin on both sides.
#
# The earlier version of this file compared against an absolute number calibrated
# on a 12-core laptop instead. It skipped all eight tests on CI's 4-vCPU runner —
# green, and guarding nothing — which is the failure mode this replaces.
MIN_SHARE_OVER_FAIR = 1.5

# If a Rust call known to hold the GIL is *not* measurably worse than fair
# sharing, the harness is not measuring what it thinks it is, so skip rather than
# report. This is a property of the measurement, not a tuned constant: it stays
# meaningful on any machine because both sides of it are measured there.
#
# Deliberately a *second* guard alongside `MAX_LOAD_PER_CORE`, not a duplicate.
# Load average counts only the guest's own processes, so on a CI VM with a noisy
# neighbour it can read near zero while the host is saturated; this one measures
# the contention instead of asking the OS about it.
MAX_CONTROL_OVER_FAIR = 1.0

# Two GIL-bound Python threads cannot each run at their solo rate, so a fair
# share at or near 1.0 is arithmetically impossible and means the calibration was
# corrupted — typically by load swinging between the solo window and the paired
# one, which inflates the ratio.
#
# This guard is not decoration. Since `fair` is the pivot every threshold derives
# from, an inflated one raises the bar instead of lowering it, so a corrupted
# calibration produces eight *false failures* rather than a skip. Observed
# directly: with the load guard disabled on a 12-core host at load average 20,
# one run measured fair=1.024, and all eight correct entry points failed against
# the 1.536 threshold it implied.
MAX_FAIR_SHARE = 0.8

# Per-core 1-minute load above which this measurement is abandoned rather than
# reported. It is not a hedge: a saturated machine descheduled the *worker*
# thread mid-call, which hands the observer time a GIL-holding call would never
# have given it, and the two populations converge. Measured on a 12-core host at
# load average 30 (other work on the box, not this test), the control rose to
# 0.311 and the weakest entry point fell to 0.349 — a ratio of 1.1 where an idle
# host shows 5. At load ~1 the same build was 10 for 10 across full-suite runs.
#
# A CI runner executing only this job sits near 0.25-0.5, well under the bar. A
# developer laptop compiling something else does not, and there the honest answer
# is "cannot tell", not a red test.
MAX_LOAD_PER_CORE = 0.7


@pytest.fixture(scope="module", autouse=True)
def fine_grained_gil():
    """Drop the GIL switch interval below one call's duration for this module,
    and put it back afterwards — it is process-global state."""
    previous = sys.getswitchinterval()
    sys.setswitchinterval(SWITCH_INTERVAL)
    yield
    sys.setswitchinterval(previous)


@pytest.fixture(scope="module")
def reference_path():
    reference = {
        "transcripts": [
            {
                "id": "TEMPLATE",
                "gene_symbol": "T",
                "strand": "+",
                "sequence": SEQUENCE,
                "exons": [{"number": 1, "start": 1, "end": SEQUENCE_LENGTH}],
            }
        ]
    }
    path = Path(tempfile.mkdtemp()) / "reference.json"
    path.write_text(json.dumps(reference))
    return str(path)


@pytest.fixture(scope="module")
def normalizer(reference_path):
    return ferro_hgvs.Normalizer(reference_json=reference_path)


@pytest.fixture(scope="module")
def checker(reference_path):
    return ferro_hgvs.EquivalenceChecker(reference_json=reference_path)


@pytest.fixture(scope="module")
def variant():
    return ferro_hgvs.parse(DESCRIPTOR)


@pytest.fixture(scope="module")
def wide_variant():
    return ferro_hgvs.parse(WIDE_DESCRIPTOR)


@pytest.fixture(scope="module")
def decomposed():
    return ferro_hgvs.parse(DECOMPOSED)


def _load_per_core():
    """1-minute load average per core, or ``None`` where the platform has no
    load average to report (Windows) — in which case the measurement proceeds
    rather than skipping on every run."""
    if not hasattr(os, "getloadavg"):
        return None
    return os.getloadavg()[0] / (os.cpu_count() or 1)


def _spin_until(deadline):
    """Pure-Python bytecode, so it runs only while this thread holds the GIL."""
    count = 0
    while time.perf_counter() < deadline:
        count += 1
    return count


def _observer_share(call):
    """The share of its solo iteration rate a Python observer thread keeps while
    a worker thread calls ``call`` in a loop.

    ~1 means the worker released the GIL and the two ran concurrently; a small
    fraction means the observer was frozen inside every call.
    """
    for _ in range(10):
        call()  # Warm up, so first-call costs land outside the window.

    alone = _spin_until(time.perf_counter() + WINDOW_SECONDS)

    counted = []
    failures = []
    # The paired deadline is set **at the rendezvous**, by the barrier's own
    # action, so the paired window is exactly `WINDOW_SECONDS` — the same span
    # `alone` was counted over.
    #
    # It used to be `perf_counter() + WINDOW_SECONDS + 0.05`, evaluated before
    # the threads were even started. Thread startup is well under a millisecond,
    # so that 50 ms was not compensation for anything: it simply ran the paired
    # phase for ~130 ms against an 80 ms solo baseline and inflated every share
    # by ~1.6x. That mostly cancels in the relative thresholds, since `fair` is
    # measured the same way — but `MAX_FAIR_SHARE` is an *absolute* guard, and
    # inflating a genuine 0.5 fair-share to ~0.8 is exactly what trips it and
    # skips all eight tests.
    window = {}

    def start_clock():
        window["deadline"] = time.perf_counter() + WINDOW_SECONDS

    barrier = threading.Barrier(2, action=start_clock)

    def observe():
        barrier.wait()
        counted.append(_spin_until(window["deadline"]))

    def work():
        # A raising `call()` must not read as a released GIL. Without this the
        # worker dies at the barrier or mid-loop, the observer runs unopposed,
        # and the share lands near 1.0 — a *passing* result produced by a broken
        # call.
        try:
            barrier.wait()
            while time.perf_counter() < window["deadline"]:
                call()
        except BaseException as exc:  # noqa: BLE001 - re-raised in the test thread
            failures.append(exc)
            barrier.abort()

    observer = threading.Thread(target=observe, daemon=True)
    worker = threading.Thread(target=work, daemon=True)
    observer.start()
    worker.start()
    # Bounded, so a deadlock fails the test in seconds instead of hanging the
    # job until its timeout.
    join_timeout = WINDOW_SECONDS * 10 + 5
    observer.join(join_timeout)
    worker.join(join_timeout)
    assert not observer.is_alive(), "observer thread did not finish"
    assert not worker.is_alive(), "worker thread did not finish"
    if failures:
        raise failures[0]
    assert counted, "observer produced no measurement"
    return counted[0] / alone


def _median_share(call, rounds=ROUNDS):
    return statistics.median(_observer_share(call) for _ in range(rounds))


def _starved_share(call, rounds=CONTROL_ROUNDS):
    """The **minimum** share across rounds, which is the right statistic for the
    negative control specifically.

    Contamination is one-directional and points opposite ways for the two sides.
    A GIL-holding worker that gets descheduled hands the observer time it would
    not otherwise have had, so noise only ever pushes the control *up* — its
    minimum is its least-contaminated sample. A detached entry point is the
    reverse: competition only takes rate away, so noise pushes it *down* and a
    median is the robust choice there. Using one statistic for both let the
    control's upper tail (0.23 typical, 0.44 worst) drive the threshold and skip
    all eight tests.
    """
    return min(_observer_share(call) for _ in range(rounds))


def _entry_points(normalizer, checker, variant, wide_variant, decomposed):
    """Every binding that reads reference sequence and was not already detached."""
    return {
        "Normalizer.to_spdi": lambda: normalizer.to_spdi(wide_variant),
        "Normalizer.canonical_spdi": lambda: normalizer.canonical_spdi(variant),
        "Normalizer.apply_to_reference": lambda: normalizer.apply_to_reference(variant),
        "Normalizer.normalize_variant": lambda: normalizer.normalize_variant(variant),
        "Normalizer.normalize": lambda: normalizer.normalize(DESCRIPTOR),
        "Normalizer.normalize_with_warnings": (
            lambda: normalizer.normalize_with_warnings(DESCRIPTOR)
        ),
        "EquivalenceChecker.check": lambda: checker.check(variant, decomposed),
        "EquivalenceChecker.all_equivalent": (
            lambda: checker.all_equivalent([variant, decomposed])
        ),
    }


def _python_worker():
    """A worker that is GIL-bound by construction and *shares* the GIL fairly
    rather than monopolising it — the midpoint between a held and a released
    call. See ``MIN_SHARE_OVER_FAIR``."""
    total = 0
    for i in range(3_000):
        total += i * i
    return total


@pytest.fixture(scope="module")
def calibration(normalizer, wide_variant):
    """Measure what "fair sharing" and "GIL held" look like *on this machine*, so
    every threshold below is relative to the host rather than to a number
    recorded on the author's laptop.

    ``fair`` is a pure-Python worker. ``held`` is a Rust call known to hold the
    GIL — rendering a ~1 Mb ``SpdiVariant`` to a Python ``str`` — which is the
    very thing the entry points no longer do.

    Emitted as a warning rather than printed: pytest captures stdout for passing
    tests, and when this measurement misbehaves on a machine you cannot log into,
    these three numbers are the entire diagnosis. The warnings summary survives
    into a CI log; a print does not.
    """
    spdi = normalizer.to_spdi(wide_variant)
    fair = _median_share(_python_worker, rounds=CONTROL_ROUNDS)
    held = _starved_share(lambda: str(spdi))
    warnings.warn(
        f"#1455 GIL calibration: fair-share={fair:.3f} gil-held={held:.3f} "
        f"ratio={held / fair:.3f} (skips if ratio >= {MAX_CONTROL_OVER_FAIR})",
        stacklevel=1,
    )
    return fair, held


#: Every row this module emits about a single entry point opens with this, whether it
#: measured or stood down, so one grep over a CI log returns exactly one row per
#: parametrized case per run. A *missing* row therefore means the case did not run,
#: which is the one thing a silent skip made unreadable.
MARGIN_PREFIX = "#1455 GIL margin:"


def _stand_down_reason(load, fair, held):
    """Why this host cannot support the measurement, or ``None`` if it can.

    All three guards in one place, in the order they are checked, as a pure function
    of three numbers — so their boundaries can be driven from a test with fixed
    values instead of being reachable only by finding a machine in the right state.
    An unexercised guard is how a guard quietly stops guarding, and this module's
    guards are the only thing standing between a contended runner and a false red.
    """
    if load is not None and load > MAX_LOAD_PER_CORE:
        return (
            f"load average is {load:.2f} per core, above {MAX_LOAD_PER_CORE} — a "
            "contended machine deschedules the worker mid-call and the measurement "
            "stops discriminating"
        )
    if fair > MAX_FAIR_SHARE:
        return (
            f"fair-sharing calibration measured {fair:.3f}, which two GIL-bound "
            "Python threads cannot reach — the calibration was corrupted by load "
            "moving between its two windows, and it is the pivot every threshold "
            "here derives from"
        )
    if held >= fair * MAX_CONTROL_OVER_FAIR:
        return (
            f"a Rust call known to hold the GIL kept {held:.3f} of the observer's "
            f"rate against {fair:.3f} for fair sharing, so this machine does not "
            "show starvation at all and the measurement cannot discriminate"
        )
    return None


def _load_field(load):
    """``os.getloadavg`` has no answer on Windows, and "no answer" must be a value in
    the row rather than a missing field, or the row's shape varies by platform."""
    return "n/a" if load is None else f"{load:.2f}"


def _margin_line(entry_point, *, shares, threshold, fair, load):
    """The row a case that took a measurement emits.

    ``shares`` is **every** attempt in order, not only the one the verdict used.
    Three ratios come out of it and they are three different statistics — pooling
    the wrong one across runs is precisely the mistake this row exists to prevent:

    ``first/fair``
        Attempt 1, taken unconditionally before any outcome is known. This is the
        field to mine: it is an unbiased sample of this host's headroom, and its
        distribution is what answers "does CI normally clear at 1.6x or at 1.51x".
    ``decided/fair``
        The attempt the assertion acted on, i.e. the last one. On a **pass** this
        exceeds ``MIN_SHARE_OVER_FAIR`` *by construction*, because the loop stops at
        the first attempt that clears — so its passing population is truncated from
        below and is not a margin distribution. It is reported because it is what
        the verdict used, and for no other purpose.
    ``all/fair``
        Every attempt. Elements past the first exist only because an earlier attempt
        missed, so they are conditioned on a miss and over-represent the low tail.
        Read one run with it; do not pool it.

    ``fair`` divides all three and is measured once for the whole module, so a
    corrupted calibration moves every row of a run together rather than one of them.
    That is what ``_stand_down_reason``'s second and third guards exist to catch, and
    why ``load/core`` rides along here instead of only in the stand-down row.
    """
    ratios = [share / fair for share in shares]
    return (
        f"{MARGIN_PREFIX} {entry_point} MEASURED "
        f"first/fair={ratios[0]:.3f} decided/fair={ratios[-1]:.3f} "
        f"share={shares[-1]:.3f} threshold={threshold:.3f} "
        f"(needs > {MIN_SHARE_OVER_FAIR}) attempts={len(shares)}/{ATTEMPTS} "
        f"all/fair={[round(ratio, 3) for ratio in ratios]} "
        f"load/core={_load_field(load)}"
    )


def _stand_down_line(entry_point, reason, load, fair, held):
    """The row a case that could not take a measurement emits.

    It carries the same three numbers on every platform and in every reason, so the
    decision to stand down is auditable from the log without the prose being parsed.
    """
    return (
        f"{MARGIN_PREFIX} {entry_point} STOOD-DOWN "
        f"load/core={_load_field(load)} fair={fair:.3f} held={held:.3f} — {reason}"
    )


def _stand_down(entry_point, reason, load, fair, held):
    """Report the stand-down on the margin series, **then** skip.

    The order is the whole point and is pinned by a test. ``pytest.skip`` raises, so
    anything after it never runs; a case that emitted nothing would be
    byte-indistinguishable in a CI log from a case that never ran, which is how a
    green suite comes to be guarding nothing.
    """
    warnings.warn(_stand_down_line(entry_point, reason, load, fair, held), stacklevel=2)
    pytest.skip(reason)


@pytest.mark.parametrize(
    "entry_point",
    [
        "Normalizer.to_spdi",
        "Normalizer.canonical_spdi",
        "Normalizer.apply_to_reference",
        "Normalizer.normalize_variant",
        "Normalizer.normalize",
        "Normalizer.normalize_with_warnings",
        "EquivalenceChecker.check",
        "EquivalenceChecker.all_equivalent",
    ],
)
def test_reference_backed_calls_release_the_gil(
    normalizer, checker, variant, wide_variant, decomposed, entry_point, calibration
):
    fair, held = calibration
    load = _load_per_core()

    # A guard may stand the measurement down; it may not do so silently. `_stand_down`
    # emits a row under the same `MARGIN_PREFIX` before it skips, so a grep over a CI
    # log gets one row per parametrized case whether or not a number was obtained.
    #
    # This is not a hypothetical tidiness: with the load guard the only one firing, a
    # CI-shaped run of the whole Python suite reported `982 passed, 8 skipped` and
    # emitted the margin ZERO times, while `-v` printed a bare `SKIPPED` with no
    # reason. The absence of a number and the absence of its explanation coincided.
    reason = _stand_down_reason(load, fair, held)
    if reason is not None:
        _stand_down(entry_point, reason, load, fair, held)

    call = _entry_points(normalizer, checker, variant, wide_variant, decomposed)[entry_point]
    threshold = fair * MIN_SHARE_OVER_FAIR

    # Keep every attempt, not only the one the loop stops on.
    #
    # The loop still breaks on the first attempt that clears — that is the retry
    # semantics `ATTEMPTS` exists for and it is unchanged. What changes is that the
    # attempts it discards are still recorded. Reporting only the deciding attempt
    # censors the passing population from below at `MIN_SHARE_OVER_FAIR`: a row could
    # never show a margin *tighter* than the threshold, which is the only regime a
    # reader of these rows cares about. Measured on this file's own base, one probe
    # run logged `apply_to_reference share/fair=1.539 attempts=3/3` — two
    # sub-threshold samples taken and thrown away — and `to_spdi share/fair=2.177
    # attempts=2/3`, which reads comfortable for a run that was not.
    #
    # See `_margin_line` for what each reported ratio means and which of them may be
    # pooled across runs.
    shares: list[float] = []
    for _ in range(ATTEMPTS):
        shares.append(_median_share(call))
        if shares[-1] > threshold:
            break
    share = shares[-1]

    # Why the row exists at all: the calibration warning reports `fair` and `held`, so
    # a CI log says whether the measurement could discriminate — but nothing reported
    # how much headroom a *passing* case had. The assertion message below carries
    # those numbers only when it fires, so "is this test chronically marginal on CI,
    # or was that one unlucky?" was unanswerable from history. A merge-group ejection
    # was measured at 1.483x fair against a 1.5x threshold, where the docstring's
    # quiet-host table has these entry points at 3.4-4.1x. This makes every run that
    # takes a measurement answer it, and every run that cannot say so out loud.
    warnings.warn(
        _margin_line(entry_point, shares=shares, threshold=threshold, fair=fair, load=load),
        stacklevel=1,
    )

    assert share > threshold, (
        f"{entry_point} starved a competing Python thread: it kept only "
        f"{share:.3f} of its solo rate while the call ran, against {fair:.3f} for "
        f"fair GIL sharing and {held:.3f} for a Rust call known to hold the GIL. "
        "A call that released the GIL scores above fair sharing, not below it. "
        "The binding is holding the GIL across reference I/O (#1455)."
    )


@pytest.mark.parametrize(
    ("load", "fair", "held", "expected"),
    [
        # A healthy host, and the two ways "healthy" can sit exactly on a boundary.
        (0.1, 0.4, 0.2, None),
        (MAX_LOAD_PER_CORE, 0.4, 0.2, None),
        (0.1, MAX_FAIR_SHARE, 0.2, None),
        # No load average to read (Windows) is not a reason to stand down.
        (None, 0.4, 0.2, None),
        (None, 2.0, 0.2, "fair-sharing calibration"),
        # One guard past its bar at a time, so each reason is attributable.
        (MAX_LOAD_PER_CORE + 0.01, 0.4, 0.2, "load average"),
        (0.1, MAX_FAIR_SHARE + 0.01, 0.2, "fair-sharing calibration"),
        # `held >= fair * MAX_CONTROL_OVER_FAIR`: equality stands the run down,
        # a hair under it does not.
        (0.1, 0.4, 0.4, "known to hold the GIL"),
        (0.1, 0.4, 0.39999, None),
    ],
)
def test_the_stand_down_guards_fire_at_their_documented_boundaries(load, fair, held, expected):
    """Each guard's bar is a claim about when this measurement stops meaning anything,
    and the three bars are the only thing that can turn a red run green here.

    Driving them with fixed numbers is the only way to exercise them at all: reaching
    a boundary by finding a machine in the right state is not reproducible, so before
    this the comparison operators were unexercised in both directions.
    """
    reason = _stand_down_reason(load, fair, held)
    if expected is None:
        assert reason is None, f"stood down on a host it should have measured: {reason}"
    else:
        assert reason is not None, "measured on a host it should have stood down on"
        assert expected in reason, f"the wrong guard fired: {reason}"


def test_the_margin_row_reports_the_attempts_the_verdict_discarded():
    """A retried case's earlier, sub-threshold attempts must reach the log.

    The loop stops at the first attempt that clears the threshold, so the attempt the
    assertion acts on is above `MIN_SHARE_OVER_FAIR` whenever the case passes. Were
    that the only attempt reported, the recorded margin could never be *below* 1.5x
    however many attempts missed first — the passing population would be truncated
    from below, censoring exactly the regime the row is collected to characterise.
    """
    fair = 0.400
    threshold = fair * MIN_SHARE_OVER_FAIR
    line = _margin_line(
        "Namespace.entry_point",
        shares=[0.480, 0.900],
        threshold=threshold,
        fair=fair,
        load=0.25,
    )

    assert line.startswith(f"{MARGIN_PREFIX} Namespace.entry_point MEASURED ")
    # 0.480 / 0.400 = 1.200, i.e. a miss — the attempt the loop discarded.
    assert "first/fair=1.200" in line
    assert "decided/fair=2.250" in line
    assert "all/fair=[1.2, 2.25]" in line
    assert f"attempts=2/{ATTEMPTS}" in line
    assert "threshold=0.600" in line
    assert "load/core=0.25" in line


def test_a_single_attempt_row_reports_that_attempt_as_both_ratios():
    """The common case, pinned so the two ratios cannot drift apart when they agree —
    and so `load/core` still has a value on a platform with no load average."""
    line = _margin_line(
        "Namespace.entry_point", shares=[1.200], threshold=0.600, fair=0.400, load=None
    )

    assert "first/fair=3.000" in line
    assert "decided/fair=3.000" in line
    assert "all/fair=[3.0]" in line
    assert f"attempts=1/{ATTEMPTS}" in line
    assert "load/core=n/a" in line


def test_a_stand_down_reports_a_row_before_it_skips():
    """A guard may stand the measurement down; it may not do so silently.

    `pytest.skip` raises, so a report placed after it never runs — and a case that
    emits nothing is indistinguishable in a CI log from a case that never ran, which
    is how the margin came to be emitted zero times on a run reporting `8 skipped`.
    `pytest.warns` fails unless the warning is raised inside the block, and the block
    cannot be left except through the skip, so the order is what is being pinned.
    """
    with (
        pytest.warns(UserWarning, match=r"#1455 GIL margin: Namespace\.ep STOOD-DOWN") as caught,
        pytest.raises(pytest.skip.Exception),
    ):
        _stand_down("Namespace.ep", "the reason it stood down", load=9.31, fair=0.815, held=0.684)

    (record,) = caught
    message = str(record.message)
    assert "load/core=9.31" in message
    assert "fair=0.815" in message
    assert "held=0.684" in message
    assert "the reason it stood down" in message


def test_the_measurement_test_reports_on_both_of_its_exits():
    """Both ways out of the measurement test must emit a row, and this is the only
    level at which that is visible.

    The tests above pin what `_margin_line` and `_stand_down` *do*. Nothing in them
    stops the measurement test from dropping the margin warning, or from calling
    `pytest.skip` directly again — and neither regression reddens anything, because a
    skip is green and a warning nobody asserts on is free to disappear. Measured:
    restoring the direct `pytest.skip` makes this module report `8 skipped`, exit 0,
    and emit zero margin rows. So the guard reads the source.
    """
    source = Path(__file__).read_text()
    tail = source.split("def test_reference_backed_calls_release_the_gil(", 1)
    assert len(tail) == 2, "the measurement test was renamed; update this guard"
    # Top-level definitions are separated by two blank lines under `ruff format`, so
    # this is the test body plus, at worst, more of the file — and "more" only makes
    # the prohibition below stricter.
    body = tail[1].split("\n\n\n", 1)[0]

    assert "_margin_line(" in body, "the measuring exit no longer emits a margin row"
    assert "_stand_down(" in body, "the standing-down exit no longer emits a margin row"
    assert "pytest.skip(" not in body, (
        "the measurement test calls `pytest.skip` directly, so that stand-down emits "
        "no margin row — the failure this module's reporting exists to remove"
    )


def test_sequence_matches_its_generator():
    """`SEQUENCE` is built by repeating a 16-base period instead of evaluating the
    generator 1.2 million times. Pin that the shortcut is exact, so an edit to the
    formula cannot silently change the fixture every other test measures against.
    """
    generated = "".join("ACGT"[(i * 7 + (i // 4) * 3) % 4] for i in range(SEQUENCE_LENGTH))
    assert generated == SEQUENCE
    assert len(SEQUENCE) == SEQUENCE_LENGTH


def test_concurrent_calls_agree_with_serial_ones(normalizer, variant):
    """Releasing the GIL means the reference work genuinely runs concurrently, so
    pin that it still produces one answer — a detached call racing on shared
    provider state would show up here and nowhere else."""
    expected = str(normalizer.canonical_spdi(variant))
    results: list[str] = []
    lock = threading.Lock()
    barrier = threading.Barrier(4)

    def run():
        barrier.wait()
        local = [str(normalizer.canonical_spdi(variant)) for _ in range(50)]
        with lock:
            results.extend(local)

    failures: list[BaseException] = []

    def guarded():
        # Captured and re-raised below: an exception here otherwise only reaches
        # stderr, and the failure surfaces as a bare count mismatch that says
        # nothing about what actually went wrong.
        try:
            run()
        except BaseException as exc:  # noqa: BLE001 - re-raised in the test thread
            failures.append(exc)
            barrier.abort()

    threads = [threading.Thread(target=guarded, daemon=True) for _ in range(4)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(60)
    assert all(not t.is_alive() for t in threads), "a worker thread did not finish"
    if failures:
        raise failures[0]

    assert len(results) == 200
    assert set(results) == {expected}
