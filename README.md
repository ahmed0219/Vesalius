# Multi-Omics Graph Architecture for Bacterial Systems Biology

Heterogeneous graph + hypergraph construction from matched bacterial multi-omics data (genomics, transcriptomics, proteomics, metabolomics), enriched with external knowledge bases (UniProt, STRING, KEGG, eggNOG). Designed for GNN-based systems biology and **antimicrobial resistance (AMR) analysis** — the graph is fused across omics layers and used to *characterize and analyze* resistance mechanisms using curated AMR knowledge as a biological annotation layer.

## Repository Structure

```
├── multiomics_graph/              # Full data→graph pipeline (main project)
│   ├── main.py                    # Pipeline orchestrator (StrainPipeline, CLI)
│   ├── preprocessing/             # GFF parsing, RNA-seq, proteomics processing
│   ├── annotation/                # UniProt, STRING, KEGG, eggNOG / COG annotations
│   ├── amr/                       # Curated AMR manifest + loader (amr_manifest.json, amr.py)
│   ├── graph/                     # Heterogeneous graph + hypergraph builders (incl. AMR hyperedges)
│   ├── visualization/             # Network plots (incl. AMR highlight), distribution charts
│   ├── scripts/                   # amr_report.py, amr_figures.py, eggNOG-mapper helpers
│   ├── tests/                     # pytest suite (71 passing)
│   ├── AMR_GRAPH_INTEGRATION.md   # AMR knowledge-layer documentation
│   ├── METABOLOMICS_INTEGRATION.md
│   └── TECHNICAL_DOCUMENTATION.md # Full implementation details
│
├── strains/                       # Per-strain data and configuration
│   ├── B36/                       # E. coli B36 (ST131)
│   ├── MS_14384/                  # E. coli MS_14384 (ST963)
│   ├── MS_14385/                  # E. coli MS_14385
│   ├── MS_14386/                  # E. coli MS_14386 (ST224)
│   ├── MS_14387/                  # E. coli MS_14387 (ST69)
│   └── TEMPLATE/                  # Template for new strains
│
├── metabolic/                     # MetaboLights MTBLS2015 metabolite MAF files
├── docs/                          # Architecture documentation
└── README.md
```

## Quick Start

### 1. Build the heterogeneous graph from real data

```bash
cd multiomics_graph
pip install -r requirements.txt

# List available strains
python main.py --list-strains

# Full pipeline for B36: genomics → transcriptomics → proteomics → metabolomics
# → integration → annotation → graph → hypergraph → visualization → export
python main.py --strain B36

# Process all configured strains
python main.py --all-strains

# Skip online annotation API calls (use cached / precomputed data)
python main.py --strain B36 --skip-annotation

# Re-download all cached data
python main.py --strain B36 --no-cache
```

### 2. AMR analysis (Phase 2 / 3)

```bash
# Generate AMR per-strain + cross-strain reports from saved outputs
python scripts/amr_report.py --output reports

# Regenerate AMR highlight figures + verify 100% paper-locus coverage
python scripts/amr_figures.py
```

### 3. eggNOG / COG annotation

Each strain carries a real `eggnog_annotations.tsv` generated with eggNOG-mapper
2.1.13 (`emapper.py --dmnd_iterate no`) and converted via
`scripts/convert_emapper_to_tsv.py`. See `AMR_GRAPH_INTEGRATION.md` §8.3 and
`scripts/run_eggnog_mapper.py` for regeneration.

## Adding a New Strain

Each strain is self-contained in `strains/<name>/`. To add a new strain:

1. **Create the folder structure:**
   ```
   strains/<name>/
   ├── config.json          # Sample IDs, settings
   ├── genome/              # GFF, FNA, FAA files
   ├── transcriptomic/      # RNA-seq counts file (optional)
   ├── proteomics/          # SWATH-MS Excel
   └── eggnog_annotations.tsv  # Real eggNOG/COG annotations (optional)
   ```

   Copy the template to start:
   ```bash
   cp -r strains/TEMPLATE strains/MyNewStrain
   ```

2. **Place your data** in the appropriate subdirectories. The pipeline auto-discovers files by convention — no path configuration needed if you follow the expected structure. Genomes are auto-resolved via `rglob('*.gff')`; store NCBI `ncbi_dataset` folders containing `protein.faa`, `cds_from_genomic.fna`, and `genomic.gff`.

3. **Edit `strains/<name>/config.json`** with your sample IDs:
   - `transcriptomics.rpmi_samples` / `sera_samples` — RNA-seq column header IDs
   - `proteomics.rpmi_ids` / `sera_ids` — proteomics sample numeric IDs
   - `proteomics.swath.rpmi_template` / `sera_template` — column name patterns (use `{sid}` placeholder)

4. **Verify and run:**
   ```bash
   cd multiomics_graph
   python main.py --list-strains      # should show your new strain
   python main.py --strain MyNewStrain
   ```

### File path resolution rules

Files inside `strains/<name>/genome/`, `transcriptomic/`, and `proteomics/` are auto-discovered. If your data lives elsewhere, add a `paths` section to your `config.json` with relative paths from the strain folder:

```json
{
  "paths": {
    "genome_gff": "../../genome/900622635.1/.../genomic.gff",
    "rna_counts": "../../transcriptomic/GSE152966_gene_counts.txt"
  }
}
```

The STRING species fallback (CFT073 → K-12) and KEGG organism code (`eco`) are pre-configured for *E. coli*.

## Data Sources

| Data | Source | Identifier |
|------|--------|------------|
| Genome (B36) | NCBI Assembly | GCF_900622635.1 |
| Genome (MS_14384) | NCBI Assembly | GCF_900622695.1 |
| Genome (MS_14385) | NCBI Assembly | **GCF_900622655.1** |
| Genome (MS_14386) | NCBI Assembly | GCF_900622665.1 |
| Genome (MS_14387) | NCBI Assembly | GCF_900622685.1 |
| Transcriptomics (B36) | GEO | GSE152966 |
| Transcriptomics (MS_14384) | GEO | GSE152967 |
| Proteomics | SWATH-MS | `strains/<name>/proteomics/` |
| Metabolomics | MetaboLights | MTBLS2015 (GC-MS / LC-MS MAF) |
| PPI (primary) | STRING | E. coli CFT073 (tax 199310) |
| PPI (fallback) | STRING | E. coli K-12 (tax 511145) |
| Pathways | KEGG | `eco` organism code |
| COG categories | eggNOG-mapper 2.1.13 | `strains/<name>/eggnog_annotations.tsv` |
| UniProt | UniProtKB | RefSeq → UniRef mapping |

## Pipeline Stages (`multiomics_graph/main.py`)

1. **Genomics** — Parse GFF → gene table with GeneID/ProteinID mapping, compute genomic proximity edges
2. **Transcriptomics** — RNA-seq log2 CPM → condition means, logFC, regulation calls
3. **Proteomics** — SWATH-MS → abundance means, logFC, regulation calls
4. **Metabolomics** — MetaboLights MAF → metabolite abundances + KEGG/EC enrichment
5. **Integration** — Cross-omics alignment on GeneID (central dogma merge; 2-layer genome+proteome fallback when no RNA)
6. **Annotation** — UniProt functional annotation, STRING PPI (with species fallback), KEGG pathway membership, eggNOG/COG assignment (real per-strain annotations)
7. **Graph Construction** — PyG `HeteroData` with gene/protein/metabolite/annotation node types
8. **Hypergraph Construction** — Multi-omics triples, KEGG pathway, COG, PPI cluster, GO, **AMR mechanism** hyperedges
9. **Visualization** — Network layouts, incidence heatmaps, distribution charts, AMR highlight
10. **Export** — CSV, NPY, PT files for downstream GNN training (`gene_features.csv` includes `is_amr_gene` / `amr_class`)

## AMR Knowledge Layer

See **`multiomics_graph/AMR_GRAPH_INTEGRATION.md`** — full Phase 1 (WP/COG entity
resolution), Phase 2 (AMR knowledge hyperedges + visualization), Phase 3
(cross-strain AMR analysis). Key points:

- Curated per-strain resistance determinants in `amr/amr_manifest.json` (traced into each strain's GFF).
- `amr_mechanism` hyperedges per class (β-lactam, aminoglycoside, sulfonamide, trimethoprim, tetracycline, macrolide, …) with `AMR:<class>` mechanism nodes.
- **Genome-only markers are real determinants and are never removed** just because they lack RNA/protein quantitation. In `amr_highlight.png` they are drawn **outlined/hatched**, while multi-omics-supported AMR genes are **solid**; status stays explicitly `genome_only` (no values invented).
- `scripts/amr_figures.py` regenerates the figures and verifies **100%** of the paper's AMR loci are drawn (B36 8/8, MS_14384 1/1, MS_14386 14/14, MS_14387 4/4).
- `scripts/amr_path_analysis.py` traces each determinant across 8 evidence layers (genome → transcriptome → proteome → metabolome → KEGG → COG → GO → PPI/AMR) with bounded path/BFS context, and writes `reports/amr_path_*` (trace, evidence, cross-strain, candidates, report). Candidate ranking picks non-AMR genes/proteins/metabolites recurring in AMR context, resolves ids to names, and orders specific candidates before generic pathway-shared "hub" elements (e.g. central-metabolism compounds such as glucose 6-phosphate, ATP) so informative hits like `tolC`, `mrdA`, `cysK`, `oppB` surface first.

## AMR Determinant Coverage Recovery

`scripts/convert_emapper_to_tsv.py` carries the eggNOG-mapper `KEGG_ko` /
`KEGG_Pathway` columns through to `strains/<name>/eggnog_annotations.tsv`
(previously dropped), and `annotation/eggnog.py` exposes them to the graph.
This closes a coverage gap where acquired/plasmid-borne resistance genes —
invisible to the `eco`-bridged UniProt KEGG cross-references — had no KEGG
layer. Two recovery helpers run offline against the saved outputs (no network):

```bash
# Per-ID UniProt GO lookup for AMR determinants (bypasses the pipeline's
# 500-hits-per-batch truncation)
python scripts/recover_amr_uniprot.py

# Merge eggNOG KEGG/COG + recovered GO + encodes + genomic-proximity edges for
# determinants into the saved graphs (idempotent)
python scripts/enrich_amr_graph.py

# Re-trace and regenerate the report after enrichment
python scripts/amr_path_analysis.py
```

The analysis also falls back to `transcriptomics_expression.csv` for RNA when a
determinant has no aligned multi-omics row, so genome-only markers report a real
transcriptome layer instead of a fabricated zero.

## Representative Output (B36)

```
outputs/B36/
├── aligned_multiomics.csv           # Aligned gene–RNA–protein table
├── genome_genes.csv                 # Gene table from GFF
├── transcriptomics_expression.csv   # RNA expression + logFC
├── proteomics_abundance.csv         # Protein abundance + logFC
├── metabolomics_abundance.csv       # Metabolite abundances
├── uniprot_annotations.csv          # UniProt functional annotations
├── string_ppi.csv                   # STRING PPI edges (merged CFT073 + K-12)
├── gene_features.csv                # Gene features + is_amr_gene / amr_class
├── heterodata.pt                    # PyG HeteroData (graph)
├── graph/
│   ├── heterodata_metadata.txt      # Node/edge type summary
│   ├── hyperedges.csv               # Hyperedge definitions
│   ├── hypergraph_nodes.csv         # Hypergraph node table
│   ├── incidence_matrix.npy         # Incidence matrix
│   └── edges_*.csv                  # Per-type edge lists (incl. covariate)
├── figures/                         # Visualization PNGs (incl. amr_highlight.png)
└── graph/ ...                       # Heterogeneous graph exports
```

## Graph Summary (verified outputs)

Real per-strain numbers (verified by `multiomics_graph/scripts/verify_counts.py`
→ `multiomics_graph/reports/verified_stats.md`, after a full pipeline re-run +
AMR enrichment): B36 has 5,118 gene nodes (genome-annotated universe, of which
2,283 are aligned), 2,283 protein nodes, 219 metabolite nodes, 2,833 annotation
nodes, 10 directed edge types (genomic proximity, encodes, PPI, gene-pathway,
COG, GO, AMR, cluster, metabolite-pathway, metabolite-enzyme) and 3,699
hyperedges across 7 families. MS_14384/MS_14385/MS_14386/MS_14387 produce
comparable graphs (2,246–2,348 aligned genes; MS_14385 runs a 2-layer
genome+proteome alignment, no RNA, and has no AMR determinants).

Results figures and the topology/AMR analyses are regenerated by
`scripts/results_figures.py` and `scripts/amr_topology_analysis.py` from the
saved outputs (see `reports/`).

## Requirements

- Python ≥ 3.10
- pandas, numpy, matplotlib, networkx, requests
- torch (PyTorch Geometric optional, for PyG operators)
- openpyxl (for Excel-based proteomics data)

## Further Reading

- `multiomics_graph/AMR_GRAPH_INTEGRATION.md` — AMR knowledge layer (Phase 1–3, corrections, usage)
- `multiomics_graph/TECHNICAL_DOCUMENTATION.md` — Full implementation details
- `multiomics_graph/METABOLOMICS_INTEGRATION.md` — Metabolomics layer documentation
- `docs/ARCHITECTURE_OVERVIEW.md` — High-level architecture rationale