#!/usr/bin/env python3
"""Extract HGVS test patterns from the mutalyzer-hgvs-parser GitHub repository.

This script parses the test files from the cloned mutalyzer-hgvs-parser repository
and extracts HGVS variant descriptions for use as test fixtures.

Repository: https://github.com/mutalyzer/mutalyzer-hgvs-parser

Usage:
    python scripts/extract_mutalyzer_github.py [--repo-path PATH] [--output OUTPUT_PATH]
"""

import argparse
import ast
import json
import re
import time
from pathlib import Path
from typing import Any


def extract_parametrize_values(source: str, decorator_name: str = "parametrize") -> list[str]:
    """Extract values from @pytest.mark.parametrize decorators.

    Args:
        source: Python source code as string
        decorator_name: Name of the decorator to look for

    Returns:
        List of test values extracted from decorators
    """
    values = []

    # Parse the source code
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return values

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            for decorator in node.decorator_list:
                if not isinstance(decorator, ast.Call):
                    continue
                func = decorator.func
                if (
                    isinstance(func, ast.Attribute)
                    and func.attr == decorator_name
                    and isinstance(func.value, ast.Attribute)
                    and func.value.attr == "mark"
                    and len(decorator.args) >= 2
                    and isinstance(decorator.args[1], ast.List)
                ):
                    # Found @pytest.mark.parametrize
                    for elt in decorator.args[1].elts:
                        if isinstance(elt, ast.Constant):
                            values.append(elt.value)
                        elif isinstance(elt, ast.Tuple):
                            # For tuples like ("GRCh38", "NM_...")
                            for sub in elt.elts:
                                if isinstance(sub, ast.Constant):
                                    values.append(sub.value)

    return values


def extract_dict_keys(source: str) -> list[str]:
    """Extract dictionary keys that look like HGVS patterns.

    Args:
        source: Python source code as string

    Returns:
        List of HGVS-like dictionary keys
    """
    keys = []

    # Find dictionary assignments with HGVS-like keys
    # Look for patterns like "NM_...:c...." or "NP_...:p...."
    hgvs_pattern = re.compile(r'"([A-Z]{2}_\d+\.\d+(?:\([^)]+\))?:[cgnpmr]\.[^"]+)"', re.MULTILINE)
    for match in hgvs_pattern.finditer(source):
        keys.append(match.group(1))

    # Also look for LRG patterns
    lrg_pattern = re.compile(r'"(LRG_\d+[tp]\d*:[cgnpmr]\.[^"]+)"', re.MULTILINE)
    for match in lrg_pattern.finditer(source):
        keys.append(match.group(1))

    # NC/NG patterns
    nc_pattern = re.compile(r'"(N[CG]_\d+\.\d+(?:\([^)]+\))?:[gm]\.[^"]+)"', re.MULTILINE)
    for match in nc_pattern.finditer(source):
        keys.append(match.group(1))

    # Simple REF patterns for grammar testing
    ref_pattern = re.compile(r'"(REF(?:_\d+)?(?:\([^)]+\))?:[cgnpmr]?\.[^"]+)"', re.MULTILINE)
    for match in ref_pattern.finditer(source):
        keys.append(match.group(1))

    # PREF patterns (protein reference placeholders)
    pref_pattern = re.compile(r'"(PREF:[p]\.[^"]+)"', re.MULTILINE)
    for match in pref_pattern.finditer(source):
        keys.append(match.group(1))

    return keys


def categorize_variant(hgvs: str) -> dict[str, Any]:
    """Categorize an HGVS variant by type.

    Args:
        hgvs: HGVS variant description

    Returns:
        Dictionary with variant info and category
    """
    info = {
        "hgvs": hgvs,
        "coordinate_system": None,
        "variant_type": None,
        "reference_type": None,
        "is_predicted": False,
    }

    # Determine reference type
    if hgvs.startswith("NM_"):
        info["reference_type"] = "transcript"
    elif hgvs.startswith("NP_"):
        info["reference_type"] = "protein"
    elif hgvs.startswith("NC_"):
        info["reference_type"] = "chromosome"
    elif hgvs.startswith("NG_"):
        info["reference_type"] = "genomic"
    elif hgvs.startswith("NR_"):
        info["reference_type"] = "non_coding_rna"
    elif hgvs.startswith("LRG_"):
        info["reference_type"] = "lrg"
    elif hgvs.startswith(("REF", "PREF")):
        info["reference_type"] = "test_placeholder"

    # Determine coordinate system
    coord_match = re.search(r":([cgnpmr])\.", hgvs)
    if coord_match:
        coord = coord_match.group(1)
        info["coordinate_system"] = {
            "c": "coding",
            "g": "genomic",
            "n": "non_coding",
            "p": "protein",
            "m": "mitochondrial",
            "r": "rna",
        }.get(coord, coord)

    # Determine variant type
    if ">" in hgvs:
        info["variant_type"] = "substitution"
    elif "del" in hgvs and "ins" in hgvs:
        info["variant_type"] = "delins"
    elif "del" in hgvs:
        info["variant_type"] = "deletion"
    elif "ins" in hgvs:
        info["variant_type"] = "insertion"
    elif "dup" in hgvs:
        info["variant_type"] = "duplication"
    elif "inv" in hgvs:
        info["variant_type"] = "inversion"
    elif "con" in hgvs:
        info["variant_type"] = "conversion"
    elif "fs" in hgvs or "fsTer" in hgvs:
        info["variant_type"] = "frameshift"
    elif "ext" in hgvs:
        info["variant_type"] = "extension"
    elif "[" in hgvs and "]" in hgvs:
        if ";" in hgvs:
            info["variant_type"] = "allele"
        else:
            info["variant_type"] = "repeat"
    elif "=" in hgvs:
        info["variant_type"] = "no_change"
    elif "?" in hgvs:
        info["variant_type"] = "uncertain"

    # Check if predicted
    if hgvs.count("(") > hgvs.count("(NM") + hgvs.count("(NP"):
        info["is_predicted"] = True

    return info


def extract_from_file(file_path: Path) -> list[dict[str, Any]]:
    """Extract HGVS patterns from a Python test file.

    Args:
        file_path: Path to Python test file

    Returns:
        List of variant info dictionaries
    """
    patterns = []

    try:
        source = file_path.read_text()
    except Exception as e:
        print(f"  Error reading {file_path}: {e}")
        return patterns

    # Extract from parametrize decorators
    parametrize_values = extract_parametrize_values(source)
    for val in parametrize_values:
        if isinstance(val, str) and ":" in val:
            patterns.append(categorize_variant(val))

    # Extract from dictionary keys
    dict_keys = extract_dict_keys(source)
    for key in dict_keys:
        patterns.append(categorize_variant(key))

    return patterns


def extract_from_repository(repo_path: Path) -> dict[str, list[dict[str, Any]]]:
    """Extract all HGVS patterns from the repository.

    Args:
        repo_path: Path to cloned mutalyzer-hgvs-parser repository

    Returns:
        Dictionary mapping file names to extracted patterns
    """
    results = {}

    tests_dir = repo_path / "tests"
    if not tests_dir.exists():
        print(f"Tests directory not found: {tests_dir}")
        return results

    test_files = list(tests_dir.glob("test_*.py"))
    print(f"Found {len(test_files)} test files")

    for test_file in test_files:
        print(f"Processing: {test_file.name}")
        patterns = extract_from_file(test_file)
        if patterns:
            results[test_file.name] = patterns
            print(f"  Extracted {len(patterns)} patterns")

    return results


def generate_fixture(results: dict[str, list[dict[str, Any]]], output_path: Path) -> None:
    """Generate test fixture JSON file.

    Args:
        results: Dictionary of extracted patterns by file
        output_path: Path to write fixture file
    """
    # Deduplicate patterns
    seen = set()
    all_patterns = []

    for file_name, patterns in results.items():
        for p in patterns:
            hgvs = p["hgvs"]
            if hgvs not in seen:
                seen.add(hgvs)
                p["source_file"] = file_name
                all_patterns.append(p)

    # Categorize by type
    by_coordinate = {}
    by_variant_type = {}
    by_reference = {}

    for p in all_patterns:
        coord = p.get("coordinate_system", "unknown")
        vtype = p.get("variant_type", "unknown")
        ref = p.get("reference_type", "unknown")

        by_coordinate.setdefault(coord, []).append(p["hgvs"])
        by_variant_type.setdefault(vtype, []).append(p["hgvs"])
        by_reference.setdefault(ref, []).append(p["hgvs"])

    fixture = {
        "source": "mutalyzer-hgvs-parser GitHub repository",
        "repository": "https://github.com/mutalyzer/mutalyzer-hgvs-parser",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "total_patterns": len(all_patterns),
        "summary": {
            "by_coordinate_system": {k: len(v) for k, v in by_coordinate.items()},
            "by_variant_type": {k: len(v) for k, v in by_variant_type.items()},
            "by_reference_type": {k: len(v) for k, v in by_reference.items()},
        },
        "patterns": all_patterns,
        "by_file": {
            file_name: [p["hgvs"] for p in patterns] for file_name, patterns in results.items()
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(fixture, f, indent=2)

    print(f"\nGenerated fixture: {output_path}")
    print(f"Total unique patterns: {len(all_patterns)}")
    print("\nBy coordinate system:")
    for coord, patterns in sorted(by_coordinate.items(), key=lambda x: (x[0] is None, x[0])):
        print(f"  {coord}: {len(patterns)}")
    print("\nBy variant type:")
    for vtype, patterns in sorted(by_variant_type.items(), key=lambda x: (x[0] is None, x[0])):
        print(f"  {vtype}: {len(patterns)}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract HGVS patterns from mutalyzer-hgvs-parser repository"
    )
    parser.add_argument(
        "--repo-path",
        type=Path,
        default=Path("external-repos/mutalyzer-hgvs-parser"),
        help="Path to cloned repository",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/fixtures/grammar/mutalyzer_github.json"),
        help="Output fixture path",
    )

    args = parser.parse_args()

    if not args.repo_path.exists():
        print(f"Repository not found: {args.repo_path}")
        print("Please clone the repository first:")
        print(
            "  git clone https://github.com/mutalyzer/mutalyzer-hgvs-parser external-repos/mutalyzer-hgvs-parser"
        )
        return

    print(f"Extracting patterns from: {args.repo_path}")
    results = extract_from_repository(args.repo_path)

    if results:
        generate_fixture(results, args.output)
    else:
        print("No patterns extracted")


if __name__ == "__main__":
    main()
