# AMR Topology Analysis

Stratified permutation test (AMR labels permuted within strain, 10,000 draws, two-sided). AMR determinants = curated manifest loci present in the gene-node universe of each strain.

| feature | AMR mean | background mean | diff | p (perm) | p_adj (BH) | Cohen's d | CI95 low | CI95 high | n_AMR | n_bg |
|---|---|---|---|---|---|---|---|---|---|---|
| knowledge_degree | 4.593 | 2.107 | +2.486 | 0.0071 | 0.0089 | +0.572 | +1.414 | +3.621 | 27 | 25086 |
| proximity_degree | 1.741 | 0.816 | +0.925 | 0.0001 | 0.0002 | +2.382 | +0.743 | +1.085 | 27 | 25086 |
| protein_go_degree | 3.111 | 1.102 | +2.009 | 0.0010 | 0.0017 | +0.697 | +1.120 | +2.963 | 27 | 25086 |
| protein_ppi_degree | 0.000 | 1.024 | -1.024 | 0.1972 | 0.1972 | -0.201 | -1.089 | -0.962 | 27 | 25086 |
| hyperedge_participation | 6.704 | 2.554 | +4.150 | 0.0001 | 0.0002 | +1.425 | +3.091 | +5.115 | 27 | 25086 |

`p (perm)` = two-sided stratified permutation test (AMR labels permuted within strain, 10,000 draws). `p_adj (BH)` = Benjamini-Hochberg FDR across the 5 features. `Cohen's d` = standardized mean difference (effect size). `CI95` = bootstrap 95% percentile interval (2,000 resamples, stratified by strain).

## Genome-only vs multi-omics-supported AMR determinants

| feature | genome-only mean | multi-omics mean | n_go | n_mo |
|---|---|---|---|---|
### B36
| knowledge_degree | 4.67 | 4.80 | 3 | 5 |
| proximity_degree | 2.00 | 1.60 | 3 | 5 |
| protein_go_degree | 3.33 | 4.20 | 3 | 5 |
| protein_ppi_degree | 0.00 | 0.00 | 3 | 5 |
| hyperedge_participation | 3.33 | 8.00 | 3 | 5 |
### MS_14386
| knowledge_degree | 2.50 | 4.83 | 2 | 12 |
| proximity_degree | 2.00 | 1.75 | 2 | 12 |
| protein_go_degree | 0.00 | 3.08 | 2 | 12 |
| protein_ppi_degree | 0.00 | 0.00 | 2 | 12 |
| hyperedge_participation | 1.50 | 7.75 | 2 | 12 |
### MS_14387
| knowledge_degree | 1.00 | 5.33 | 1 | 3 |
| proximity_degree | 1.00 | 1.67 | 1 | 3 |
| protein_go_degree | 0.00 | 4.00 | 1 | 3 |
| protein_ppi_degree | 0.00 | 0.00 | 1 | 3 |
| hyperedge_participation | 1.00 | 8.33 | 1 | 3 |

## Hyperedge-family participation of AMR determinants

| family | determinant-memberships |
|---|---|
| amr_mechanism | 27 |
| cog_category | 25 |
| go_term | 21 |
| multi_omics_triple | 21 |
| kegg_pathway | 13 |
| metabolic_pathway | 11 |
