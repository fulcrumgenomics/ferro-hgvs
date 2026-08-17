# Project a variant to another axis

`ferro project` re-expresses a variant on a chosen output axis — genomic (`g`), coding (`c`), non-coding (`n`), protein (`p`), or RNA (`r`) — against a transcript. It is how you move a genomic call onto a transcript, or read the protein consequence of a coding one.

```bash
ferro project --axis c --transcript NM_000059.4 --reference <dir> \
  'NC_000013.11:g.32316466T>G'
# NC_000013.11(NM_000059.4):c.6T>G
```

## Coding-axis rules apply on the projected axis

The projected `c.`/`r.` axis is normalized on the transcript, so it applies the coding-axis rules a bare genomic axis cannot — and produces the same string `ferro normalize` gives for the coding-authored form of the variant, whichever axis you started from.

The clearest case is the coding one-amino-acid `delins` exception (`DNA/delins.md:18`): two substitutions one nucleotide apart **within a single codon** are written as one `delins`, not as two members.

```bash
# BRCA2 codon 2 (c.4_6 = CCT): c.4 C>A and c.6 T>G, with c.5 unchanged.
ferro project --axis c --transcript NM_000059.4 --reference <dir> \
  'NC_000013.11:g.[32316464C>A;32316466T>G]'
# NC_000013.11(NM_000059.4):c.4_6delinsACG
```

Because both members fall in one codon, the coding axis merges them into `c.4_6delinsACG` — identical to `ferro normalize` on `NM_000059.4:c.[4C>A;6T>G]`. The genomic axis has no reading frame, so it keeps the members individual, and a pair that straddles a codon boundary or an exon junction stays split on every axis.

## Pairs with deriving from sequences

If you hold bases rather than a description — a window from a BAM, a VCF row — derive the genomic description with `from_sequences` first (see [Deriving a description from sequences](deriving-from-sequences.md)), then project it onto the transcript:

```bash
# 1. derive the g. description from the bases
# 2. project it onto the transcript
ferro project --axis c --transcript NM_000059.4 --reference <dir> '<the g. description>'
```
