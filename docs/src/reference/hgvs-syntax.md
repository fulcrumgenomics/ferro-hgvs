# Supported HGVS Syntax

| Type | Prefix | Example |
|------|--------|---------|
| Genomic | `g.` | `NC_000001.11:g.12345A>G` |
| Coding DNA | `c.` | `NM_000088.3:c.459A>G` |
| Non-coding | `n.` | `NR_000001.1:n.100A>G` |
| RNA | `r.` | `NM_000088.3:r.459a>g` |
| Protein | `p.` | `NP_000079.2:p.Val600Glu` |
| Mitochondrial | `m.` | `NC_012920.1:m.3243A>G` |

## Edit Types

- Substitution: `A>G`, `Val600Glu`
- Deletion: `del`, `100_200del`
- Insertion: `100_101insATG`
- Deletion-Insertion: `100_102delinsATG`
- Duplication: `100_102dup`
- Inversion: `100_200inv`
- Repeat: `100CAG[20]`
