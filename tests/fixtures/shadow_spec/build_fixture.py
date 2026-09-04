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

# Proteins the shadow-spec `p.` pages normalize against, keyed by the accession the
# pages spell, mapped to the transcript whose CDS translates to it. The protein
# sequence is the *exact translation* of that transcript's CDS in the committed slice
# (there is no protein FASTA in the prepared reference — `protein_fastas` is empty),
# so it is derived here rather than fetched. Grow this as protein pages are added.
WANTED_PROTEINS = {"NP_003997.1": "NM_004006.3"}

# Standard genetic code (NCBI table 1), for translating a CDS to its protein.
_CODON_BASES = "TCAG"
_CODON_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
GENETIC_CODE = {
    a + b + c: _CODON_AAS[i]
    for i, (a, b, c) in enumerate(
        (a, b, c) for a in _CODON_BASES for b in _CODON_BASES for c in _CODON_BASES
    )
}


def translate_cds(accession, sequence, cds_start, cds_end):
    """Translate the CDS `[cds_start, cds_end]` (1-based, inclusive) to one-letter
    residues, dropping the trailing stop codon.

    Refuses a CDS that is not a clean protein-coding frame — a non-triplet length, an
    untranslatable codon, an internal stop, or a missing terminal stop — rather than
    silently coercing it (dropping trailing bases, mapping to `X`, keeping an internal
    `*`). The committed `proteins` map is the shadow-spec normalization reference, so a
    malformed translation must fail the build, not become the baseline."""
    cds = sequence[cds_start - 1 : cds_end].upper()
    if len(cds) % 3 != 0:
        sys.exit(f"{accession}: CDS length {len(cds)} is not a multiple of 3")
    codons = [cds[k : k + 3] for k in range(0, len(cds), 3)]
    residues = [GENETIC_CODE.get(codon, "X") for codon in codons]
    if "X" in residues:
        sys.exit(f"{accession}: untranslatable codon {codons[residues.index('X')]!r} in CDS")
    if residues[-1] != "*":
        sys.exit(f"{accession}: CDS does not end in a stop codon")
    if "*" in residues[:-1]:
        sys.exit(f"{accession}: internal stop codon in CDS")
    return "".join(residues[:-1])


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

# Protein slice: the exact translation of each source transcript's CDS. Built from the
# same `records` above so the CDS bounds cannot drift between the two halves.
by_id = {r["id"]: r for r in records}
proteins = {}
for protein_acc, source_tx in WANTED_PROTEINS.items():
    t = by_id.get(source_tx)
    if t is None:
        sys.exit(f"{protein_acc}: source transcript {source_tx} is not in WANTED")
    proteins[protein_acc] = translate_cds(protein_acc, t["sequence"], t["cds_start"], t["cds_end"])

# The `{ transcripts, proteins }` envelope a `JsonProvider` deserializes. Writing a bare
# transcript array here would delete `proteins` and make regeneration non-idempotent.
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "transcripts.json")
with open(out, "w") as fh:
    json.dump({"transcripts": records, "proteins": proteins}, fh, indent=2)
    fh.write("\n")
print(
    f"wrote {out}: transcripts={[(r['id'], len(r['sequence']), r['cds_start'], r['cds_end']) for r in records]} "
    f"proteins={[(k, len(v)) for k, v in proteins.items()]}"
)
