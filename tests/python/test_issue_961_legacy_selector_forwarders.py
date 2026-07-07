"""Issue #961: the Python provider must forward the two legacy-selector-resolution
methods (`resolve_legacy_gene_selector`, `sole_hosted_transcript`) to its inner
provider, so `NG_(GENE):c.` legacy-selector resolution (#500/#792) and bare-`NG_:c.`
selector synthesis (#923) work through the Python API.

Before the fix, `PyProvider` let both methods fall through to the trait defaults
(which return `None`), so neither rewrite happened on the Python path even when the
wrapped `MultiFastaProvider` supported them: the selector was left as the gene symbol
(or the input stayed bare). The observable contract is that `Normalizer.normalize`
rewrites the genomic-reference selector to the reference-standard `NM_`.

The fixture is generated programmatically (no committed reference data): a synthetic
single-transcript reference — an all-`A` `NM_003002.4` (SDHD) whose CDS is supplied via
supplemental metadata (the manifest carries no cdot) — plus a RefSeqGene summary (the
global gene->transcript map) and an `ng_hosted_transcripts` artifact (the parent-relative
map). This mirrors the Rust-side coverage in `tests/issue_500_legacy_gene_selector.rs`.
"""

from __future__ import annotations

import json
from pathlib import Path

import ferro_hgvs

# Synthetic reference constants. Real accessions/gene, synthetic all-`A` sequence:
# the rewrite is driven by the selector maps, and `c.92A>T` matches the all-`A` ref.
NG_PARENT = "NG_012337.3"
GENE = "SDHD"
TRANSCRIPT = "NM_003002.4"
SEQUENCE = "A" * 300  # CDS 1..300 (a multiple of 3); c.92 -> transcript position 92.


def _write_single_record_fasta(path: Path, name: str, sequence: str) -> None:
    """Write a one-record FASTA and its `.fai` index (single-line sequence)."""
    header = f">{name}\n"
    path.write_text(f"{header}{sequence}\n")
    # .fai columns: name, length, offset-of-first-base, linebases, linewidth(+newline).
    offset = len(header.encode())
    path.with_suffix(path.suffix + ".fai").write_text(
        f"{name}\t{len(sequence)}\t{offset}\t{len(sequence)}\t{len(sequence) + 1}\n"
    )


def _build_reference(root: Path) -> Path:
    """Materialize the manifest + all referenced artifacts under ``root``; return the manifest path."""
    tx_dir = root / "transcripts"
    tx_dir.mkdir(parents=True)
    fasta = tx_dir / "sdhd.fna"
    _write_single_record_fasta(fasta, TRANSCRIPT, SEQUENCE)

    # Supplemental CDS metadata sidecar ({basename}.metadata.json) — the manifest
    # has no cdot, so the transcript's CDS frame comes from here (1-based inclusive).
    fasta.with_suffix(".metadata.json").write_text(
        json.dumps(
            {
                "transcripts": {
                    TRANSCRIPT: {
                        "cds_start": 1,
                        "cds_end": 300,
                        "gene_symbol": GENE,
                        "sequence_length": len(SEQUENCE),
                    }
                }
            }
        )
    )

    # RefSeqGene summary: the global gene-symbol -> reference-standard NM_ map. Tab
    # columns are tax_id, GeneID, Symbol, RSG, LRG, RNA, t, Protein, p, Category.
    (root / "refseqgene_summary.tsv").write_text(
        "\t".join(
            [
                "9606",
                "6392",
                GENE,
                NG_PARENT,
                "-",
                TRANSCRIPT,
                "t1",
                "NP_002993.1",
                "p1",
                "reference standard",
            ]
        )
        + "\n"
    )

    # ng_hosted_transcripts artifact: NG_ parent -> gene(upper) -> hosted transcripts.
    (root / "ng_hosted.json").write_text(
        json.dumps({"schema_version": 1, "records": {NG_PARENT: {GENE: [TRANSCRIPT]}}})
    )

    manifest = root / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "prepared_at": "2026-07-07T00:00:00Z",
                "transcript_fastas": ["transcripts/sdhd.fna"],
                "supplemental_fasta": "transcripts/sdhd.fna",
                "refseqgene_summary": "refseqgene_summary.tsv",
                "ng_hosted_transcripts": "ng_hosted.json",
                "genome_fasta": None,
                "cdot_json": None,
                "transcript_count": 1,
                "available_prefixes": ["NM_"],
            }
        )
    )
    return manifest


def test_ng_gene_selector_is_forwarded_and_rewritten(tmp_path: Path) -> None:
    """`NG_(GENE):c.` resolves the legacy gene-symbol selector to the NM_ via the Python path."""
    manifest = _build_reference(tmp_path)
    norm = ferro_hgvs.Normalizer.from_manifest(str(manifest))
    result = norm.normalize(f"{NG_PARENT}({GENE}):c.92A>T")
    assert result == f"{NG_PARENT}({TRANSCRIPT}):c.92A>T", (
        f"legacy gene-symbol selector should resolve to {TRANSCRIPT}; got {result}"
    )


def test_bare_ng_selector_is_forwarded_and_synthesized(tmp_path: Path) -> None:
    """Bare `NG_:c.` synthesizes the uniquely-hosted transcript selector via the Python path."""
    manifest = _build_reference(tmp_path)
    norm = ferro_hgvs.Normalizer.from_manifest(str(manifest))
    result = norm.normalize(f"{NG_PARENT}:c.92A>T")
    assert result == f"{NG_PARENT}({TRANSCRIPT}):c.92A>T", (
        f"bare-NG_ input should synthesize the {TRANSCRIPT} selector; got {result}"
    )
