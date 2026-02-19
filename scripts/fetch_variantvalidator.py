#!/usr/bin/env python3
"""Fetch validated HGVS variants from VariantValidator API.

VariantValidator provides comprehensive HGVS validation, normalization,
and format conversion services.

API Documentation: https://rest.variantvalidator.org

Key endpoints:
- /VariantValidator/variantvalidator/{genome}/{variant}/{transcripts} - Validate variant
- /VariantValidator/tools/hgvs2reference/{hgvs} - Convert to reference coordinates

Rate Limit: 2 requests/second

Usage:
    python scripts/fetch_variantvalidator.py [--output OUTPUT_PATH]
"""

import argparse
import json
import time
import urllib.parse
from pathlib import Path
from typing import Any

import requests

VV_API = "https://rest.variantvalidator.org"
RATE_LIMIT_DELAY = 1.0  # seconds between requests (conservative)


# Test variants for validation
TEST_VARIANTS = [
    # Clinically important substitutions
    ("GRCh38", "NM_000546.6:c.215C>G"),
    ("GRCh38", "NM_000546.6:c.524G>A"),  # TP53 R175H
    ("GRCh38", "NM_004333.4:c.1799T>A"),  # BRAF V600E
    ("GRCh38", "NM_001126112.3:c.35G>T"),  # KRAS G12V
    ("GRCh38", "NM_001904.4:c.2369C>T"),  # CTNNB1
    ("GRCh38", "NM_000314.8:c.388C>T"),  # PTEN
    # Deletions
    ("GRCh38", "NM_000492.4:c.1521_1523del"),  # CFTR F508del
    ("GRCh38", "NM_007294.4:c.68_69del"),  # BRCA1 185delAG
    ("GRCh38", "NM_000059.4:c.5946del"),  # BRCA2
    ("GRCh38", "NM_000546.6:c.586delC"),
    ("GRCh38", "NM_000546.6:c.586_597del"),
    # Insertions and duplications
    ("GRCh38", "NM_000546.6:c.375_376insA"),
    ("GRCh38", "NM_000518.4:c.27dupG"),
    ("GRCh38", "NM_007294.4:c.5266dupC"),  # BRCA1
    ("GRCh38", "NM_000546.6:c.722_723insGTACC"),
    # Intronic variants
    ("GRCh38", "NM_000546.6:c.782+1G>A"),
    ("GRCh38", "NM_000546.6:c.672+1G>A"),
    ("GRCh38", "NM_000546.6:c.375+1G>A"),
    ("GRCh38", "NM_000546.6:c.376-1G>A"),
    ("GRCh38", "NM_000546.6:c.920-2A>G"),
    ("GRCh38", "NM_000251.3:c.942+3A>T"),  # MSH2
    # Complex indels
    ("GRCh38", "NM_000546.6:c.375delCinsAA"),
    ("GRCh38", "NM_000546.6:c.100_102delinsTTT"),
    ("GRCh38", "NM_000546.6:c.743_744delinsTT"),
    # Genomic variants
    ("GRCh38", "NC_000017.11:g.7674220C>T"),
    ("GRCh38", "NC_000017.11:g.7673802G>A"),
    ("GRCh38", "NC_000007.14:g.140753336A>T"),  # BRAF region
    ("GRCh38", "NC_000013.11:g.32316461C>T"),  # BRCA2 region
    # GRCh37 variants for testing build support
    ("GRCh37", "NM_000546.5:c.215C>G"),
    ("GRCh37", "NM_000546.5:c.524G>A"),
    ("GRCh37", "NM_004333.4:c.1799T>A"),
    # UTR variants
    ("GRCh38", "NM_000546.6:c.-28G>C"),
    ("GRCh38", "NM_000546.6:c.*1171T>C"),
    # Additional CFTR variants
    ("GRCh38", "NM_000492.4:c.254G>A"),
    ("GRCh38", "NM_000492.4:c.1000C>T"),
    ("GRCh38", "NM_000492.4:c.3909C>G"),
    # PIK3CA variants
    ("GRCh38", "NM_006218.4:c.1633G>A"),
    ("GRCh38", "NM_006218.4:c.3140A>G"),
    # Nonsense variants
    ("GRCh38", "NM_000546.6:c.637C>T"),  # p.Arg213Ter
    ("GRCh38", "NM_000546.6:c.916C>T"),  # p.Arg306Ter
    # More splice variants
    ("GRCh38", "NM_000546.6:c.919+1G>T"),
    ("GRCh38", "NM_007294.4:c.4358-2A>G"),  # BRCA1
    # Additional cancer genes
    ("GRCh38", "NM_001128425.1:c.1624G>T"),  # EGFR
    ("GRCh38", "NM_002524.5:c.35G>A"),  # NRAS
    ("GRCh38", "NM_004360.5:c.7397C>T"),  # CDH1
    ("GRCh38", "NM_024675.4:c.1633G>A"),  # PALB2
    ("GRCh38", "NM_004448.4:c.2573T>G"),  # ERBB2
    ("GRCh38", "NM_000249.4:c.1799T>A"),  # MLH1
]


def validate_variant(
    variant: str, genome: str = "GRCh38", transcripts: str = "all"
) -> dict[str, Any] | None:
    """Validate a variant using VariantValidator API.

    Args:
        variant: HGVS variant description
        genome: Genome build (GRCh37 or GRCh38)
        transcripts: Transcript selection (all, select, raw, etc.)

    Returns:
        API response dict or None on error
    """
    encoded = urllib.parse.quote(variant, safe="")
    url = f"{VV_API}/VariantValidator/variantvalidator/{genome}/{encoded}/{transcripts}"

    try:
        response = requests.get(url, timeout=60)
        if response.status_code == 200:
            return response.json()
        else:
            return {"error": f"HTTP {response.status_code}", "input": variant}
    except Exception as e:
        return {"error": str(e), "input": variant}


def hgvs_to_reference(variant: str) -> dict[str, Any] | None:
    """Convert HGVS to reference coordinates.

    Args:
        variant: HGVS variant description

    Returns:
        API response dict or None on error
    """
    encoded = urllib.parse.quote(variant, safe="")
    url = f"{VV_API}/VariantValidator/tools/hgvs2reference/{encoded}"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.json()
        else:
            return {"error": f"HTTP {response.status_code}", "input": variant}
    except Exception as e:
        return {"error": str(e), "input": variant}


def fetch_all_variants(
    variants: list[tuple[str, str]],
) -> list[dict[str, Any]]:
    """Fetch validation data for all variants.

    Args:
        variants: List of (genome, variant) tuples

    Returns:
        List of result dictionaries
    """
    results = []

    for i, (genome, variant) in enumerate(variants):
        print(f"[{i + 1}/{len(variants)}] Processing: {variant} ({genome})")

        result = {
            "input": variant,
            "genome_build": genome,
            "validation": None,
            "reference": None,
        }

        # Fetch validation
        validation = validate_variant(variant, genome, "all")
        if validation:
            result["validation"] = validation
        time.sleep(RATE_LIMIT_DELAY)

        # Fetch reference coordinates (for coding variants)
        if ":c." in variant or ":n." in variant:
            ref = hgvs_to_reference(variant)
            if ref:
                result["reference"] = ref
            time.sleep(RATE_LIMIT_DELAY)

        results.append(result)

        # Progress indicator
        if (i + 1) % 10 == 0:
            print(f"  Progress: {i + 1}/{len(variants)} complete")

    return results


def extract_key_info(validation: dict[str, Any]) -> dict[str, Any]:
    """Extract key information from validation response.

    Args:
        validation: Full validation response

    Returns:
        Simplified dict with key fields
    """
    info = {
        "valid": False,
        "normalized": None,
        "genomic": None,
        "protein": None,
        "warnings": [],
        "errors": [],
    }

    if "error" in validation:
        info["errors"].append(validation["error"])
        return info

    # Look for the primary result
    for key, value in validation.items():
        if isinstance(value, dict):
            if "validation_warnings" in value:
                info["warnings"].extend(value.get("validation_warnings", []))

            if "hgvs_transcript_variant" in value:
                info["normalized"] = value.get("hgvs_transcript_variant")
                info["valid"] = True

            if "primary_assembly_loci" in value:
                loci = value.get("primary_assembly_loci", {})
                if "grch38" in loci:
                    grch38 = loci["grch38"]
                    if "hgvs_genomic_description" in grch38:
                        info["genomic"] = grch38["hgvs_genomic_description"]

            if "hgvs_predicted_protein_consequence" in value:
                protein_info = value.get("hgvs_predicted_protein_consequence", {})
                if "slr" in protein_info:
                    info["protein"] = protein_info["slr"]

    return info


def generate_fixture(results: list[dict[str, Any]], output_path: Path) -> None:
    """Generate test fixture JSON file.

    Args:
        results: List of API results
        output_path: Path to write fixture file
    """
    # Process results to extract key info
    processed = []
    for r in results:
        entry = {
            "input": r["input"],
            "genome_build": r["genome_build"],
            "raw_validation": r["validation"],
            "reference": r["reference"],
        }

        if r["validation"] and "error" not in r["validation"]:
            entry["summary"] = extract_key_info(r["validation"])
        else:
            entry["summary"] = {"valid": False, "errors": [r["validation"].get("error", "Unknown error")]}

        processed.append(entry)

    fixture = {
        "source": "VariantValidator API",
        "api_base": VV_API,
        "generated": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "total_variants": len(processed),
        "successful": sum(1 for p in processed if p["summary"].get("valid", False)),
        "variants": processed,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(fixture, f, indent=2)

    print(f"\nGenerated fixture: {output_path}")
    print(f"Total variants: {fixture['total_variants']}")
    print(f"Successful: {fixture['successful']}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch validated variants from VariantValidator API"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/fixtures/validation/variantvalidator_api.json"),
        help="Output fixture path",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Limit number of variants to process (0 = all)",
    )

    args = parser.parse_args()

    variants = TEST_VARIANTS
    if args.limit > 0:
        variants = variants[: args.limit]

    print(f"Fetching {len(variants)} variants from VariantValidator API...")
    print(f"Rate limit: {RATE_LIMIT_DELAY}s between requests\n")

    results = fetch_all_variants(variants)
    generate_fixture(results, args.output)


if __name__ == "__main__":
    main()
