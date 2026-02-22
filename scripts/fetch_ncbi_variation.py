#!/usr/bin/env python3
"""Fetch HGVS/SPDI conversions from NCBI Variation Services API.

NCBI Variation Services provides HGVS to SPDI conversion and vice versa,
as well as rsID lookups and contextual alignment.

API Documentation: https://api.ncbi.nlm.nih.gov/variation/v0/

Key endpoints:
- /hgvs/{hgvs}/contextuals - Convert HGVS to SPDI with contextual info
- /spdi/{spdi}/hgvs - Convert SPDI to HGVS expressions
- /spdi/{spdi}/vcf - Convert SPDI to VCF format
- /refsnp/{rsid} - Get variant info by rsID

Rate Limit: 3 requests/second

Usage:
    python scripts/fetch_ncbi_variation.py [--output OUTPUT_PATH]
"""

import argparse
import json
import time
import urllib.parse
from pathlib import Path
from typing import Any

import requests

NCBI_API = "https://api.ncbi.nlm.nih.gov/variation/v0"
RATE_LIMIT_DELAY = 0.5  # seconds between requests


# Test HGVS variants for SPDI conversion
TEST_HGVS_VARIANTS = [
    # Simple substitutions
    "NM_000546.6:c.215C>G",
    "NM_000546.6:c.524G>A",  # TP53 R175H
    "NM_004333.4:c.1799T>A",  # BRAF V600E
    "NM_000518.5:c.20A>T",  # HBB sickle cell
    "NM_000518.4:c.27dupG",
    # Deletions
    "NM_000492.4:c.1521_1523del",  # CFTR F508del
    "NM_007294.4:c.68_69del",  # BRCA1
    "NM_000546.6:c.586delC",
    # Insertions
    "NM_000546.6:c.375_376insA",
    "NM_000546.6:c.722_723insGTACC",
    # Genomic coordinates
    "NC_000017.11:g.7674220C>T",
    "NC_000017.11:g.7673802G>A",
    "NC_000007.14:g.140753336A>T",
    "NC_000013.11:g.32316461C>T",
    # Complex indels
    "NM_000546.6:c.100_102delinsTTT",
    "NM_000546.6:c.743_744delinsTT",
    # More clinically important
    "NM_001126112.3:c.35G>T",  # KRAS G12V
    "NM_000059.4:c.5946del",  # BRCA2
    "NM_000314.8:c.388C>T",  # PTEN
    "NM_006218.4:c.1633G>A",  # PIK3CA
    "NM_006218.4:c.3140A>G",  # PIK3CA
    # Additional
    "NM_000249.4:c.1799T>A",  # MLH1
    "NM_000251.3:c.942+3A>T",  # MSH2
    "NM_002524.5:c.35G>A",  # NRAS
    "NM_004448.4:c.2573T>G",  # ERBB2
]

# Known rsIDs for validation
TEST_RSIDS = [
    "rs121913529",  # BRAF V600E
    "rs80357906",  # BRCA1 185delAG
    "rs113993960",  # CFTR F508del
    "rs28934578",  # TP53 R175H
    "rs63750447",  # BRCA2 6174delT
    "rs121913527",  # BRAF V600K
    "rs121913530",  # BRAF K601E
    "rs121912651",  # KRAS G12D
    "rs121913535",  # KRAS G12V
    "rs121913237",  # EGFR L858R
    "rs121908120",  # TP53 R248Q
    "rs121912666",  # TP53 R273H
    "rs121913279",  # PIK3CA H1047R
    "rs121913273",  # PIK3CA E545K
    "rs587776544",  # BRCA1 C61G
    "rs80358550",  # BRCA2 S1982Rfs*22
    "rs11571833",  # BRCA2 K3326*
    "rs1042522",  # TP53 P72R (common polymorphism)
    "rs334",  # HBB sickle cell
    "rs6025",  # F5 Leiden
]


def hgvs_to_spdi(hgvs: str) -> dict[str, Any] | None:
    """Convert HGVS to SPDI using NCBI API.

    Args:
        hgvs: HGVS variant description

    Returns:
        API response dict or None on error
    """
    encoded = urllib.parse.quote(hgvs, safe="")
    url = f"{NCBI_API}/hgvs/{encoded}/contextuals"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.json()
        else:
            return {"error": f"HTTP {response.status_code}", "input": hgvs}
    except Exception as e:
        return {"error": str(e), "input": hgvs}


def spdi_to_hgvs(spdi: str) -> dict[str, Any] | None:
    """Convert SPDI to HGVS using NCBI API.

    Args:
        spdi: SPDI variant representation

    Returns:
        API response dict or None on error
    """
    encoded = urllib.parse.quote(spdi, safe="")
    url = f"{NCBI_API}/spdi/{encoded}/hgvs"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.json()
        else:
            return {"error": f"HTTP {response.status_code}", "input": spdi}
    except Exception as e:
        return {"error": str(e), "input": spdi}


def spdi_to_vcf(spdi: str) -> dict[str, Any] | None:
    """Convert SPDI to VCF using NCBI API.

    Args:
        spdi: SPDI variant representation

    Returns:
        API response dict or None on error
    """
    encoded = urllib.parse.quote(spdi, safe="")
    url = f"{NCBI_API}/spdi/{encoded}/vcf"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.json()
        else:
            return {"error": f"HTTP {response.status_code}", "input": spdi}
    except Exception as e:
        return {"error": str(e), "input": spdi}


def get_rsid_info(rsid: str) -> dict[str, Any] | None:
    """Get variant info by rsID.

    Args:
        rsid: dbSNP rsID (e.g., rs121913529)

    Returns:
        API response dict or None on error
    """
    # Remove 'rs' prefix if present
    rsid_num = rsid.replace("rs", "")
    url = f"{NCBI_API}/refsnp/{rsid_num}"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            return response.json()
        else:
            return {"error": f"HTTP {response.status_code}", "input": rsid}
    except Exception as e:
        return {"error": str(e), "input": rsid}


def fetch_hgvs_conversions(variants: list[str]) -> list[dict[str, Any]]:
    """Fetch SPDI conversions for HGVS variants.

    Args:
        variants: List of HGVS variant descriptions

    Returns:
        List of result dictionaries
    """
    results = []

    for i, variant in enumerate(variants):
        print(f"[{i + 1}/{len(variants)}] HGVS: {variant}")

        result = {
            "input_hgvs": variant,
            "spdi_result": None,
            "roundtrip_hgvs": None,
            "vcf": None,
        }

        # Convert HGVS to SPDI
        spdi_result = hgvs_to_spdi(variant)
        if spdi_result:
            result["spdi_result"] = spdi_result

            # If we got SPDI, try roundtrip back to HGVS
            if "data" in spdi_result and "spdis" in spdi_result.get("data", {}):
                spdis = spdi_result["data"]["spdis"]
                if spdis and len(spdis) > 0:
                    # Use the first SPDI for roundtrip
                    first_spdi = spdis[0]
                    spdi_str = f"{first_spdi['seq_id']}:{first_spdi['position']}:{first_spdi['deleted_sequence']}:{first_spdi['inserted_sequence']}"

                    time.sleep(RATE_LIMIT_DELAY)
                    roundtrip = spdi_to_hgvs(spdi_str)
                    if roundtrip:
                        result["roundtrip_hgvs"] = roundtrip

                    time.sleep(RATE_LIMIT_DELAY)
                    vcf = spdi_to_vcf(spdi_str)
                    if vcf:
                        result["vcf"] = vcf

        results.append(result)
        time.sleep(RATE_LIMIT_DELAY)

        if (i + 1) % 10 == 0:
            print(f"  Progress: {i + 1}/{len(variants)} complete")

    return results


def fetch_rsid_data(rsids: list[str]) -> list[dict[str, Any]]:
    """Fetch variant data for rsIDs.

    Args:
        rsids: List of rsIDs

    Returns:
        List of result dictionaries
    """
    results = []

    for i, rsid in enumerate(rsids):
        print(f"[{i + 1}/{len(rsids)}] rsID: {rsid}")

        result = {
            "rsid": rsid,
            "info": None,
            "hgvs_expressions": [],
            "spdi": None,
        }

        info = get_rsid_info(rsid)
        if info and "error" not in info:
            result["info"] = info

            # Extract HGVS expressions if available
            if "primary_snapshot_data" in info:
                snapshot = info["primary_snapshot_data"]
                if "placements_with_allele" in snapshot:
                    for placement in snapshot["placements_with_allele"]:
                        if "alleles" in placement:
                            for allele in placement["alleles"]:
                                if "hgvs" in allele:
                                    result["hgvs_expressions"].append(allele["hgvs"])
                                if "spdi" in allele:
                                    result["spdi"] = allele["spdi"]
        else:
            result["info"] = info

        results.append(result)
        time.sleep(RATE_LIMIT_DELAY)

        if (i + 1) % 10 == 0:
            print(f"  Progress: {i + 1}/{len(rsids)} complete")

    return results


def generate_fixture(
    hgvs_results: list[dict[str, Any]],
    rsid_results: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Generate test fixture JSON file.

    Args:
        hgvs_results: HGVS conversion results
        rsid_results: rsID lookup results
        output_path: Path to write fixture file
    """
    fixture = {
        "source": "NCBI Variation Services API",
        "api_base": NCBI_API,
        "generated": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "hgvs_conversions": {
            "total": len(hgvs_results),
            "successful": sum(
                1 for r in hgvs_results if r.get("spdi_result") and "error" not in r["spdi_result"]
            ),
            "variants": hgvs_results,
        },
        "rsid_lookups": {
            "total": len(rsid_results),
            "successful": sum(
                1 for r in rsid_results if r.get("info") and "error" not in r["info"]
            ),
            "variants": rsid_results,
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(fixture, f, indent=2)

    print(f"\nGenerated fixture: {output_path}")
    print(
        f"HGVS conversions: {fixture['hgvs_conversions']['successful']}/{fixture['hgvs_conversions']['total']}"
    )
    print(
        f"rsID lookups: {fixture['rsid_lookups']['successful']}/{fixture['rsid_lookups']['total']}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch HGVS/SPDI conversions from NCBI Variation Services"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/fixtures/validation/ncbi_variation.json"),
        help="Output fixture path",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Limit number of variants per category (0 = all)",
    )
    parser.add_argument(
        "--skip-rsids",
        action="store_true",
        help="Skip rsID lookups (faster)",
    )

    args = parser.parse_args()

    hgvs_variants = TEST_HGVS_VARIANTS
    rsids = TEST_RSIDS

    if args.limit > 0:
        hgvs_variants = hgvs_variants[: args.limit]
        rsids = rsids[: args.limit]

    print("NCBI Variation Services API")
    print(f"Rate limit: {RATE_LIMIT_DELAY}s between requests\n")

    print(f"=== HGVS to SPDI Conversions ({len(hgvs_variants)} variants) ===")
    hgvs_results = fetch_hgvs_conversions(hgvs_variants)

    if args.skip_rsids:
        rsid_results = []
        print("\n=== Skipping rsID lookups ===")
    else:
        print(f"\n=== rsID Lookups ({len(rsids)} IDs) ===")
        rsid_results = fetch_rsid_data(rsids)

    generate_fixture(hgvs_results, rsid_results, args.output)


if __name__ == "__main__":
    main()
