#!/usr/bin/env python3
"""Extract HGVS examples from the official HGVS nomenclature specification.

This script extracts example HGVS expressions from hgvs-nomenclature.org
and organizes them by variant type and coordinate system.

Source: https://hgvs-nomenclature.org

Usage:
    python scripts/extract_hgvs_spec_examples.py [--output OUTPUT_PATH]
"""

import argparse
import json
import re
import time
from pathlib import Path
from typing import Any

import requests
from bs4 import BeautifulSoup

HGVS_BASE_URL = "https://hgvs-nomenclature.org"

# Page structure for HGVS nomenclature site
PAGES = {
    "substitution": {
        "dna": "/recommendations/DNA/substitution/",
        "rna": "/recommendations/RNA/substitution/",
        "protein": "/recommendations/protein/substitution/",
    },
    "deletion": {
        "dna": "/recommendations/DNA/deletion/",
        "rna": "/recommendations/RNA/deletion/",
        "protein": "/recommendations/protein/deletion/",
    },
    "duplication": {
        "dna": "/recommendations/DNA/duplication/",
        "rna": "/recommendations/RNA/duplication/",
        "protein": "/recommendations/protein/duplication/",
    },
    "insertion": {
        "dna": "/recommendations/DNA/insertion/",
        "rna": "/recommendations/RNA/insertion/",
        "protein": "/recommendations/protein/insertion/",
    },
    "indel": {
        "dna": "/recommendations/DNA/indel/",
        "rna": "/recommendations/RNA/indel/",
        "protein": "/recommendations/protein/indel/",
    },
    "inversion": {
        "dna": "/recommendations/DNA/inversion/",
    },
    "conversion": {
        "dna": "/recommendations/DNA/conversion/",
    },
    "repeat": {
        "dna": "/recommendations/DNA/repeated/",
    },
    "frameshift": {
        "protein": "/recommendations/protein/frameshift/",
    },
    "extension": {
        "protein": "/recommendations/protein/extension/",
    },
}


def extract_hgvs_from_text(text: str) -> list[str]:
    """Extract HGVS patterns from text.

    Args:
        text: Text to search for HGVS patterns

    Returns:
        List of HGVS expressions found
    """
    patterns = []

    # Match common HGVS patterns
    # NC/NG/NM/NR/NP accessions with variants
    accession_pattern = re.compile(
        r'((?:NC|NG|NM|NR|NP|LRG)_[\d\.]+(?:\([^)]+\))?:[cgmnpr]\.[^\s<>"\']+)',
        re.IGNORECASE
    )

    for match in accession_pattern.finditer(text):
        hgvs = match.group(1)
        # Clean up trailing punctuation
        hgvs = re.sub(r'[,;.)\]]+$', '', hgvs)
        if len(hgvs) > 10:  # Skip very short matches
            patterns.append(hgvs)

    return patterns


def fetch_page(url: str) -> str | None:
    """Fetch a page from the HGVS nomenclature site.

    Args:
        url: URL to fetch

    Returns:
        Page HTML content or None if failed
    """
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.text
        else:
            print(f"  HTTP {response.status_code} for {url}")
            return None
    except Exception as e:
        print(f"  Error fetching {url}: {e}")
        return None


def extract_examples_from_page(html: str, variant_type: str, level: str) -> list[dict[str, Any]]:
    """Extract HGVS examples from a page.

    Args:
        html: HTML content
        variant_type: Type of variant (substitution, deletion, etc.)
        level: Molecular level (dna, rna, protein)

    Returns:
        List of example dictionaries
    """
    examples = []
    soup = BeautifulSoup(html, 'html.parser')

    # Extract text content
    text = soup.get_text()

    # Find all HGVS patterns
    hgvs_patterns = extract_hgvs_from_text(text)

    for hgvs in hgvs_patterns:
        # Determine coordinate system
        coord_match = re.search(r':([cgmnpr])\.', hgvs)
        coord_system = coord_match.group(1) if coord_match else None

        coord_name = {
            'c': 'coding',
            'g': 'genomic',
            'm': 'mitochondrial',
            'n': 'non_coding',
            'p': 'protein',
            'r': 'rna',
        }.get(coord_system, 'unknown')

        examples.append({
            "hgvs": hgvs,
            "variant_type": variant_type,
            "level": level,
            "coordinate_system": coord_name,
            "source": "hgvs-nomenclature.org",
        })

    return examples


def deduplicate_examples(examples: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Remove duplicate HGVS patterns.

    Args:
        examples: List of example dictionaries

    Returns:
        Deduplicated list
    """
    seen = set()
    unique = []
    for ex in examples:
        hgvs = ex["hgvs"]
        if hgvs not in seen:
            seen.add(hgvs)
            unique.append(ex)
    return unique


def generate_fixture(examples: list[dict[str, Any]], output_path: Path) -> None:
    """Generate JSON fixture file.

    Args:
        examples: List of example dictionaries
        output_path: Path to write fixture
    """
    # Organize by type
    by_type: dict[str, list[str]] = {}
    by_coord: dict[str, list[str]] = {}
    by_level: dict[str, list[str]] = {}

    for ex in examples:
        vtype = ex.get("variant_type", "unknown")
        coord = ex.get("coordinate_system", "unknown")
        level = ex.get("level", "unknown")

        by_type.setdefault(vtype, []).append(ex["hgvs"])
        by_coord.setdefault(coord, []).append(ex["hgvs"])
        by_level.setdefault(level, []).append(ex["hgvs"])

    fixture = {
        "description": "HGVS examples from official HGVS nomenclature specification",
        "source": "https://hgvs-nomenclature.org",
        "version": "2024.01",
        "generated": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "total_examples": len(examples),
        "summary": {
            "by_variant_type": {k: len(v) for k, v in by_type.items()},
            "by_coordinate_system": {k: len(v) for k, v in by_coord.items()},
            "by_level": {k: len(v) for k, v in by_level.items()},
        },
        "examples": examples,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(fixture, f, indent=2)

    print(f"\nGenerated fixture: {output_path}")
    print(f"Total examples: {len(examples)}")
    print("\nBy variant type:")
    for vtype, patterns in sorted(by_type.items()):
        print(f"  {vtype}: {len(patterns)}")
    print("\nBy coordinate system:")
    for coord, patterns in sorted(by_coord.items()):
        print(f"  {coord}: {len(patterns)}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract HGVS examples from hgvs-nomenclature.org"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/fixtures/grammar/hgvs_spec_examples.json"),
        help="Output fixture path",
    )

    args = parser.parse_args()

    all_examples = []

    print("Extracting HGVS examples from hgvs-nomenclature.org...")

    for variant_type, levels in PAGES.items():
        for level, path in levels.items():
            url = HGVS_BASE_URL + path
            print(f"Fetching {variant_type}/{level}: {url}")

            html = fetch_page(url)
            if html:
                examples = extract_examples_from_page(html, variant_type, level)
                all_examples.extend(examples)
                print(f"  Found {len(examples)} examples")

            time.sleep(0.5)  # Rate limiting

    # Deduplicate
    unique_examples = deduplicate_examples(all_examples)
    print(f"\nTotal unique examples: {len(unique_examples)}")

    # Generate fixture
    generate_fixture(unique_examples, args.output)


if __name__ == "__main__":
    main()
