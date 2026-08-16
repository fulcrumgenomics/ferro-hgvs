#!/usr/bin/env python3
"""Regenerate the shadow-spec reference slice from a prepared ferro reference.

Extracts the full, verbatim transcript sequence + CDS bounds for each accession the
shadow-spec examples normalize against, and writes a `transcripts.json` a `JsonProvider`
can serve. Bases are REAL (GRCh38 RefSeq), so the committed slice reproduces ferro's
production normalization for these accessions without a multi-GB reference.

Usage:
    FERRO_MANIFEST=<ref>/manifest.json python3 build_fixture.py
"""

import glob
import gzip
import json
import os
import sys

MANIFEST = os.environ.get("FERRO_MANIFEST")
if not MANIFEST:
    sys.exit("set FERRO_MANIFEST to a prepared reference manifest.json")
ROOT = os.path.dirname(MANIFEST)
CDOT = os.path.join(ROOT, "cdot", "cdot-0.2.32.refseq.GRCh38.json")

# Accessions the shadow-spec pages normalize against. Grow this as pages are added.
WANTED = ["NM_004006.3"]


def transcript_sequence(accession):
    for f in glob.glob(os.path.join(ROOT, "transcripts", "*.fna.gz")):
        with gzip.open(f, "rt") as fh:
            grab, buf = False, []
            for line in fh:
                if line.startswith(">"):
                    if grab:
                        return "".join(buf)
                    grab = line[1:].split()[0] == accession
                    buf = []
                elif grab:
                    buf.append(line.strip())
            if grab:
                return "".join(buf)
    sys.exit(f"{accession}: sequence not found under transcripts/")


with open(CDOT) as fh:
    cdot = json.load(fh)["transcripts"]
records = []
for acc in WANTED:
    seq = transcript_sequence(acc)
    t = cdot[acc]
    # Fixture convention (verified against NM_000088.4): cds_start = start_codon + 1
    # (1-based first CDS base), cds_end = stop_codon.
    cds_start = t["start_codon"] + 1
    cds_end = t["stop_codon"]
    records.append(
        {
            "id": acc,
            "gene_symbol": t.get("gene_name", acc),
            "strand": "+",
            "sequence": seq,
            "cds_start": cds_start,
            "cds_end": cds_end,
            # Flat single exon: substitution normalization crosses no junction, and c./n.
            # positions are flat transcript offsets (ruling c-and-n-positions-are-flat-...).
            "exons": [{"number": 1, "start": 1, "end": len(seq)}],
        }
    )

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "transcripts.json")
with open(out, "w") as fh:
    json.dump(records, fh, indent=2)
    fh.write("\n")
print(
    f"wrote {out}: {[(r['id'], len(r['sequence']), r['cds_start'], r['cds_end']) for r in records]}"
)
