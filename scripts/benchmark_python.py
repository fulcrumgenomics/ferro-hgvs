#!/usr/bin/env python3
"""
Benchmark script for ferro-hgvs Python bindings.

Compares:
1. Python API overhead (time to call from Python into Rust)
2. Python API vs native Rust CLI (same patterns through both)
3. Python API vs biocommons hgvs library (if installed)

Usage:
    python scripts/benchmark_python.py
    python scripts/benchmark_python.py --iterations 10000
    python scripts/benchmark_python.py --include-biocommons
"""

import argparse
import math
import subprocess
import time
from collections.abc import Callable
from dataclasses import dataclass

# Test variants covering different HGVS types
TEST_VARIANTS = [
    # Genomic variants
    ("g.sub", "NC_000001.11:g.12345A>G"),
    ("g.del", "NC_000001.11:g.100del"),
    ("g.del_range", "NC_000001.11:g.100_200del"),
    ("g.dup", "NC_000001.11:g.100dup"),
    ("g.ins", "NC_000001.11:g.100_101insATG"),
    ("g.delins", "NC_000001.11:g.100_200delinsATG"),
    # CDS variants
    ("c.sub", "NM_000088.3:c.459A>G"),
    ("c.del", "NM_000088.3:c.459del"),
    ("c.ins", "NM_000088.3:c.459_460insATG"),
    ("c.intronic+", "NM_000088.3:c.100+5G>A"),
    ("c.intronic-", "NM_000088.3:c.100-10A>G"),
    ("c.dup", "NM_000088.3:c.100_102dup"),
    ("c.utr3", "NM_000088.3:c.*100A>G"),
    ("c.utr5", "NM_000088.3:c.-50A>G"),
    # Protein variants
    ("p.sub", "NP_000079.2:p.Val600Glu"),
    ("p.del", "NP_000079.2:p.Val600del"),
    ("p.fs", "NP_000079.2:p.Val600fs"),
    ("p.fsTer", "NP_000079.2:p.Val600fsTer15"),
    ("p.ext", "NP_000079.2:p.Met1ext-5"),
    ("p.identity", "NP_000079.2:p.="),
    # Non-coding variants
    ("n.sub", "NR_000001.1:n.100A>G"),
    ("n.del", "NR_000001.1:n.100del"),
    # RNA variants
    ("r.sub", "NM_000088.3:r.100a>g"),
]


@dataclass
class BenchmarkResult:
    """Result of a benchmark run."""

    name: str
    variant: str
    time_ns: float
    iterations: int

    @property
    def time_us(self) -> float:
        return self.time_ns / 1000

    @property
    def ops_per_sec(self) -> float:
        return 1_000_000_000 / self.time_ns if self.time_ns > 0 else 0


def benchmark_function(
    func: Callable[[str], object],
    variant: str,
    iterations: int,
    warmup: int = 100,
) -> float:
    """Benchmark a function, returning average time in nanoseconds."""
    # Warmup
    for _ in range(warmup):
        func(variant)

    # Benchmark
    start = time.perf_counter_ns()
    for _ in range(iterations):
        func(variant)
    elapsed = time.perf_counter_ns() - start

    return elapsed / iterations


def benchmark_ferro_python(iterations: int) -> list[BenchmarkResult]:
    """Benchmark ferro-hgvs Python API."""
    import ferro_hgvs

    results = []
    for name, variant in TEST_VARIANTS:
        time_ns = benchmark_function(ferro_hgvs.parse, variant, iterations)
        results.append(BenchmarkResult(name, variant, time_ns, iterations))

    return results


def benchmark_ferro_rust_cli(iterations: int) -> list[BenchmarkResult] | None:
    """Benchmark ferro-hgvs Rust CLI (if available)."""
    # Check if ferro CLI is available
    try:
        result = subprocess.run(
            ["ferro", "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode != 0:
            return None
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None

    results = []
    # For CLI, we batch all variants and measure total time
    # This is less precise but avoids subprocess overhead per variant
    all_variants = [v for _, v in TEST_VARIANTS]

    # Warmup
    for _ in range(10):
        subprocess.run(
            ["ferro", "parse"] + all_variants[:5],
            capture_output=True,
            timeout=30,
        )

    # Benchmark - run multiple iterations of parsing all variants
    batch_iterations = max(1, iterations // 100)  # Reduce iterations for CLI
    start = time.perf_counter_ns()
    for _ in range(batch_iterations):
        subprocess.run(
            ["ferro", "parse"] + all_variants,
            capture_output=True,
            timeout=30,
        )
    elapsed = time.perf_counter_ns() - start

    # Calculate average time per variant
    time_per_variant = elapsed / (batch_iterations * len(all_variants))

    for name, variant in TEST_VARIANTS:
        results.append(BenchmarkResult(name, variant, time_per_variant, batch_iterations))

    return results


def benchmark_biocommons_hgvs(iterations: int) -> list[BenchmarkResult] | None:
    """Benchmark biocommons hgvs library (if installed)."""
    try:
        import hgvs.parser

        parser = hgvs.parser.Parser()
    except ImportError:
        return None

    results = []
    for name, variant in TEST_VARIANTS:
        try:
            time_ns = benchmark_function(parser.parse_hgvs_variant, variant, iterations)
            results.append(BenchmarkResult(name, variant, time_ns, iterations))
        except Exception:
            # Some variants may not be supported
            results.append(BenchmarkResult(name, variant, float("nan"), iterations))

    return results


def print_comparison_table(
    ferro_py: list[BenchmarkResult],
    ferro_rust: list[BenchmarkResult] | None,
    biocommons: list[BenchmarkResult] | None,
) -> None:
    """Print a comparison table of all benchmarks."""
    print("\n" + "=" * 80)
    print("BENCHMARK RESULTS: HGVS Parsing Performance")
    print("=" * 80)

    # Header
    headers = ["Variant", "ferro-py (us)"]
    if ferro_rust:
        headers.append("ferro-rust (us)")
        headers.append("py/rust")
    if biocommons:
        headers.append("biocommons (us)")
        headers.append("ferro/bio")

    header_fmt = "{:<15}" + "{:>15}" * (len(headers) - 1)
    print(header_fmt.format(*headers))
    print("-" * (15 + 15 * (len(headers) - 1)))

    # Data rows
    for i, result in enumerate(ferro_py):
        row = [result.name, f"{result.time_us:.2f}"]

        if ferro_rust:
            rust_result = ferro_rust[i]
            row.append(f"{rust_result.time_us:.2f}")
            if rust_result.time_ns > 0:
                ratio = result.time_ns / rust_result.time_ns
                row.append(f"{ratio:.1f}x")
            else:
                row.append("N/A")

        if biocommons:
            bio_result = biocommons[i]
            if math.isnan(bio_result.time_ns):
                row.append("N/A")
                row.append("N/A")
            else:
                row.append(f"{bio_result.time_us:.2f}")
                if bio_result.time_ns > 0:
                    speedup = bio_result.time_ns / result.time_ns
                    row.append(f"{speedup:.1f}x")
                else:
                    row.append("N/A")

        print(header_fmt.format(*row))

    # Summary statistics
    print("-" * (15 + 15 * (len(headers) - 1)))

    # Average times
    avg_ferro_py = sum(r.time_ns for r in ferro_py) / len(ferro_py)
    row = ["AVERAGE", f"{avg_ferro_py / 1000:.2f}"]

    if ferro_rust:
        avg_rust = sum(r.time_ns for r in ferro_rust) / len(ferro_rust)
        row.append(f"{avg_rust / 1000:.2f}")
        row.append(f"{avg_ferro_py / avg_rust:.1f}x")

    if biocommons:
        valid_bio = [r for r in biocommons if not math.isnan(r.time_ns)]
        if valid_bio:
            avg_bio = sum(r.time_ns for r in valid_bio) / len(valid_bio)
            row.append(f"{avg_bio / 1000:.2f}")
            row.append(f"{avg_bio / avg_ferro_py:.1f}x")
        else:
            row.append("N/A")
            row.append("N/A")

    print(header_fmt.format(*row))

    # Throughput
    print("\n" + "-" * 40)
    print("THROUGHPUT (variants/second):")
    print("-" * 40)
    print(f"  ferro-hgvs (Python): {1_000_000_000 / avg_ferro_py:,.0f} ops/sec")
    if ferro_rust:
        print(f"  ferro-hgvs (Rust):   {1_000_000_000 / avg_rust:,.0f} ops/sec")
    if biocommons and valid_bio:
        print(f"  biocommons hgvs:     {1_000_000_000 / avg_bio:,.0f} ops/sec")


def print_overhead_analysis(ferro_py: list[BenchmarkResult]) -> None:
    """Analyze and print Python API overhead."""
    print("\n" + "=" * 80)
    print("PYTHON API OVERHEAD ANALYSIS")
    print("=" * 80)

    avg_time_ns = sum(r.time_ns for r in ferro_py) / len(ferro_py)
    min_time_ns = min(r.time_ns for r in ferro_py)
    max_time_ns = max(r.time_ns for r in ferro_py)

    print(f"\nParsing statistics (ferro-hgvs Python API):")
    print(f"  Average: {avg_time_ns / 1000:.2f} us ({avg_time_ns:.0f} ns)")
    print(f"  Min:     {min_time_ns / 1000:.2f} us ({min_time_ns:.0f} ns)")
    print(f"  Max:     {max_time_ns / 1000:.2f} us ({max_time_ns:.0f} ns)")

    # Estimate PyO3 overhead (typically 50-200ns per call)
    # This is a rough estimate based on typical PyO3 call overhead
    estimated_pyo3_overhead_ns = 100  # Conservative estimate
    print(f"\nEstimated PyO3 call overhead: ~{estimated_pyo3_overhead_ns} ns")
    print(f"Overhead as % of average:     ~{estimated_pyo3_overhead_ns / avg_time_ns * 100:.1f}%")

    # Batch processing recommendation
    print("\nNote: For maximum performance with many variants, use BatchProcessor")
    print("which amortizes Python/Rust boundary crossing overhead.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Benchmark ferro-hgvs Python bindings",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python scripts/benchmark_python.py
    python scripts/benchmark_python.py --iterations 50000
    python scripts/benchmark_python.py --include-biocommons
        """,
    )
    parser.add_argument(
        "--iterations",
        "-n",
        type=int,
        default=10000,
        help="Number of iterations per variant (default: 10000)",
    )
    parser.add_argument(
        "--include-biocommons",
        action="store_true",
        help="Include biocommons hgvs comparison (requires 'pip install hgvs')",
    )
    parser.add_argument(
        "--skip-rust-cli",
        action="store_true",
        help="Skip Rust CLI benchmark",
    )

    args = parser.parse_args()

    print(f"Running benchmarks with {args.iterations} iterations per variant...")
    print(f"Test variants: {len(TEST_VARIANTS)}")

    # Run benchmarks
    print("\n[1/3] Benchmarking ferro-hgvs Python API...")
    ferro_py = benchmark_ferro_python(args.iterations)

    ferro_rust = None
    if not args.skip_rust_cli:
        print("[2/3] Benchmarking ferro-hgvs Rust CLI...")
        ferro_rust = benchmark_ferro_rust_cli(args.iterations)
        if ferro_rust is None:
            print("      (ferro CLI not found, skipping)")

    biocommons = None
    if args.include_biocommons:
        print("[3/3] Benchmarking biocommons hgvs...")
        biocommons = benchmark_biocommons_hgvs(args.iterations)
        if biocommons is None:
            print("      (biocommons hgvs not installed, skipping)")
            print("      Install with: pip install hgvs")
    else:
        print("[3/3] Skipping biocommons hgvs (use --include-biocommons to enable)")

    # Print results
    print_comparison_table(ferro_py, ferro_rust, biocommons)
    print_overhead_analysis(ferro_py)

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    avg_us = sum(r.time_us for r in ferro_py) / len(ferro_py)
    throughput = 1_000_000 / avg_us
    print(f"ferro-hgvs Python API: {avg_us:.2f} us/parse, {throughput:,.0f} parses/sec")

    if biocommons:
        valid_bio = [r for r in biocommons if not math.isnan(r.time_ns)]
        if valid_bio:
            bio_avg_us = sum(r.time_us for r in valid_bio) / len(valid_bio)
            speedup = bio_avg_us / avg_us
            print(f"Speedup vs biocommons hgvs: {speedup:.1f}x faster")


if __name__ == "__main__":
    main()
