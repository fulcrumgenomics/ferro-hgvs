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
``MAX_LOAD_PER_CORE`` and ``MAX_CONTROL_OVER_FAIR``. The calibration numbers are
emitted as a warning on every run so a skip on a machine you cannot log into is
still diagnosable from the CI log.
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
ATTEMPTS = 2

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
    if load is not None and load > MAX_LOAD_PER_CORE:
        pytest.skip(
            f"load average is {load:.2f} per core, above {MAX_LOAD_PER_CORE} — a "
            "contended machine deschedules the worker mid-call and the measurement "
            "stops discriminating"
        )

    if fair > MAX_FAIR_SHARE:
        pytest.skip(
            f"fair-sharing calibration measured {fair:.3f}, which two GIL-bound "
            "Python threads cannot reach — the calibration was corrupted by load "
            "moving between its two windows, and it is the pivot every threshold "
            "here derives from"
        )

    if held >= fair * MAX_CONTROL_OVER_FAIR:
        pytest.skip(
            f"a Rust call known to hold the GIL kept {held:.3f} of the observer's "
            f"rate against {fair:.3f} for fair sharing, so this machine does not "
            "show starvation at all and the measurement cannot discriminate"
        )

    call = _entry_points(normalizer, checker, variant, wide_variant, decomposed)[entry_point]
    threshold = fair * MIN_SHARE_OVER_FAIR
    for _ in range(ATTEMPTS):
        share = _median_share(call)
        if share > threshold:
            break

    assert share > threshold, (
        f"{entry_point} starved a competing Python thread: it kept only "
        f"{share:.3f} of its solo rate while the call ran, against {fair:.3f} for "
        f"fair GIL sharing and {held:.3f} for a Rust call known to hold the GIL. "
        "A call that released the GIL scores above fair sharing, not below it. "
        "The binding is holding the GIL across reference I/O (#1455)."
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
