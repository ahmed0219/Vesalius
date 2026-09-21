# Multi-Omics Heterogeneous Graph Architecture

## Technical Documentation — Implementation Chapter

---

# 1. Project Overview

## 1.1 Objective

Build a multi-omics heterogeneous graph + hypergraph for bacterial systems biology from matched genomics, transcriptomics, and proteomics data of *Escherichia coli* B36, enriched with external knowledge from UniProt, STRING, KEGG, and eggNOG databases. The architecture is designed for downstream Graph Neural Network (GNN) learning.

## 1.2 Biological Motivation

Antimicrobial resistance (AMR) in bacterial pathogens involves complex, multi-layered regulatory responses that span the central dogma — from genomic organization through transcription to protein abundance. No single omics layer fully captures this response. The architecture integrates all three layers (DNA, RNA, protein) into a unified graph representation where a GNN can learn from inter-omics dependencies.

## 1.3 Computational Motivation

Heterogeneous graphs (PyG `HeteroData`) support typed nodes and typed, directed edges — ideal for representing distinct biological entities (genes, proteins) and their relationships (encodes, interacts, co-expresses). Hypergraphs complement this by capturing higher-order set-level relationships (pathways, functional categories, PPI modules) that pairwise edges cannot represent.

The pipeline produces two complementary representations:
- **Heterogeneous graph** for node-level and edge-level tasks (gene classification, interaction prediction)
- **Hypergraph** for set-level tasks (pathway activity prediction, module identification)

## 1.4 Implementation Scope

The pipeline covers 10 stages executed sequentially in `StrainPipeline.run_all()`:

```
Genomics → Transcriptomics → Proteomics → Metabolomics → Integration →
Annotation → Graph Construction → Hypergraph Construction →
Visualization → Export
```

## 1.5 Project Status

| Module | Status |
|--------|--------|
| Genomics (GFF parsing) | **Completed** |
| Transcriptomics (RNA-seq) | **Completed** |
| Proteomics (SWATH-MS) | **Completed** |
| Cross-omics Integration | **Completed** |
| UniProt Annotation | **Completed** |
| STRING PPI | **Completed** |
| KEGG Pathway | **Completed** |
| eggNOG COG | **Completed** |
| Heterogeneous Graph Construction | **Completed** |
| Hypergraph Construction | **Completed** |
| Visualization | **Completed** |
| Export | **Completed** |
| GO Term Hyperedges | **Completed** |
| GNN Training / Prediction | **Not Yet Implemented** |

## 1.6 Execution Workflow

```
python main.py --strain B36 [--skip-annotation] [--no-cache] [--correlation 0.7]
```

The pipeline processes a single strain end-to-end. Five strains are configured and runnable (B36, MS_14384, MS_14385, MS_14386, MS_14387); use `python main.py --all-strains` to process all.

---

# 2. Repository Structure

```
multiomics_graph/
├── main.py                         # Pipeline orchestrator, CLI, StrainPipeline class
├── __init__.py                     # Empty
├── preprocessing/
│   ├── __init__.py                 # Exports GenomicsProcessor, TranscriptomicsProcessor, ProteomicsProcessor
│   ├── genomics.py                 # GFF parser, gene-protein mapping, genomic proximity edges
│   ├── transcriptomics.py          # RNA-seq expression → condition means, logFC, regulation
│   └── proteomics.py               # SWATH-MS → protein abundance, logFC, regulation
├── annotation/
│   ├── __init__.py                 # Empty
│   ├── uniprot.py                  # RefSeq → UniProt mapping + functional annotation retrieval
│   ├── string_api.py               # STRING PPI network retrieval with species fallback
│   ├── kegg.py                     # KEGG pathway membership via Gene → KO → Pathway
│   └── eggnog.py                   # eggNOG/COG functional category assignment
├── graph/
│   ├── __init__.py                 # Exports HeterogeneousGraphBuilder, BiologicalHypergraphBuilder
│   ├── heterogeneous_graph.py      # PyG HeteroData construction with typed nodes and edges
│   └── hypergraph.py               # Incidence-matrix hypergraph with biological module hyperedges
├── visualization/
│   └── network_plots.py            # Matplotlib visualizations (dashboard, incidence matrix, concordance)
├── tests/
│   ├── __init__.py                 # Empty
│   └── test_uniprot.py             # 50 tests for UniProt annotation module
└── outputs/
    ├── cache/                      # JSON API response cache (UniProt, STRING, KEGG)
    └── B36/                        # Strain-specific output directory
        ├── string_ppi.csv
        ├── heterodata.pt
        ├── hypergraph_laplacian.pt
        ├── aligned_multiomics.csv
        ├── genome_genes.csv
        ├── transcriptomics_expression.csv
        ├── proteomics_abundance.csv
        ├── uniprot_annotations.csv
        ├── graph/
        │   ├── heterodata.pt
        │   ├── incidence_matrix.npy
        │   ├── hyperedges.csv
        │   ├── hypergraph_nodes.csv
        │   ├── hyperedge_type_counts.csv
        │   └── hypergraph_laplacian.pt
        └── figures/
            ├── summary_dashboard.png
            ├── regulation_concordance.png
            ├── incidence_matrix.png
            ├── cog_distribution.png
            ├── hyperedge_types.png
            └── heterogeneous_graph.png
```

## 2.1 File Details

### `main.py` (904 lines)
- **Purpose**: Pipeline orchestrator, CLI entry point, stage coordination
- **Public classes**: `StrainPipeline`
- **Public functions**: `main()`, `load_config()`
- **Dependencies**: All processing and annotation modules, PyTorch, pandas, numpy
- **Inputs**: Command-line arguments + STRAIN_CONFIGS
- **Outputs**: All CSV, NPY, PT files in `outputs/{strain}/`
- **Interaction**: Calls every other module in sequence; central hub

### `preprocessing/genomics.py` (190 lines)
- **Purpose**: Parse bacterial GFF → gene table, compute genomic proximity edges
- **Public class**: `GenomicsProcessor`
- **Public methods**: `parse_gff()`, `compute_genomic_edges()`, `summary()`, `get_gene_by_protein()`, `get_protein_by_gene()`
- **Dependencies**: pandas, numpy
- **Inputs**: NCBI RefSeq GFF3 file
- **Outputs**: `gene_table` DataFrame, `genomic_edges` DataFrame

### `preprocessing/transcriptomics.py` (223 lines)
- **Purpose**: Load RNA-seq expression, compute condition means and regulation
- **Public class**: `TranscriptomicsProcessor`
- **Public methods**: `load_log2_cpm()`, `compute_condition_means()`, `get_regulated_genes()`, `summary()`
- **Dependencies**: pandas, numpy
- **Inputs**: Gene × sample expression matrix (CSV/TSV/Excel)
- **Outputs**: `expression_table` DataFrame, `raw_data` DataFrame

### `preprocessing/proteomics.py` (251 lines)
- **Purpose**: Load SWATH-MS proteomics, compute abundance and regulation
- **Public class**: `ProteomicsProcessor`
- **Public methods**: `load_protein_groups()`, `filter_bacterial_proteins()`, `compute_condition_means()`, `summary()`
- **Dependencies**: pandas, numpy
- **Inputs**: SWATH-MS Excel (`Area - proteins` sheet)
- **Outputs**: `protein_table` DataFrame

### `annotation/uniprot.py` (678 lines)
- **Purpose**: Map RefSeq WP_ → UniProt, fetch functional annotations with caching
- **Public class**: `UniProtAnnotator`
- **Public dataclasses**: `GoTerm`, `MappingStats`, `AnnotationStats`
- **Internal classes**: `_CacheManager`, `_AnnotationParser`
- **Public methods**: `map_refseq_to_uniprot()`, `fetch_annotations()`, `get_go_term_info()`
- **Dependencies**: requests, pandas, json
- **Inputs**: List of RefSeq WP_ accessions
- **Outputs**: `uniprot_annotations` DataFrame (31 columns)
- **API**: `https://rest.uniprot.org/uniprotkb/search`

### `annotation/string_api.py` (243 lines)
- **Purpose**: STRING PPI network retrieval with UniProt→STRING mapping and caching
- **Public class**: `STRINGClient`
- **Public methods**: `get_string_ids()`, `fetch_interactions()`, `summary()`
- **Dependencies**: requests, pandas, json
- **Inputs**: List of UniProt accessions
- **Outputs**: PPI edge DataFrame with `combined_score`, `experimental_score`, `database_score`, `coexpression_score`, `cooccurrence_score`
- **API**: `https://string-db.org/api/json`

### `annotation/kegg.py` (308 lines)
- **Purpose**: Map bacterial genes to KEGG pathways via Gene→KO→Pathway chain
- **Public class**: `KEGGAnnotator`
- **Public methods**: `list_pathways()`, `get_ko_for_genes()`, `map_ko_to_pathways()`, `build_pathway_membership()`, `get_pathway_name()`, `summary()`
- **Dependencies**: requests, pandas, json
- **Inputs**: List of GeneIDs + optional GeneID→KEGG gene ID map
- **Outputs**: `pathway_membership` dict `{pathway_id: [GeneIDs]}`
- **API**: `https://rest.kegg.jp`

### `annotation/eggnog.py` (195 lines)
- **Purpose**: Assign COG functional categories to bacterial genes
- **Public class**: `EggNOGAnnotator`
- **Public methods**: `load_precomputed()`, `annotate_by_locus_tag()`, `get_genes_by_cog()`, `as_locus_tag_map()`, `kegg_pathway_membership()`, `summary()`
- **Dependencies**: pandas, json
- **Inputs**: List of GeneIDs (or precomputed file path)
- **Outputs**: `annotations` dict `{GeneID: {COG_classes, COG_category, Description, KEGG_ko, KEGG_Pathway}}`
- **Note**: `load_precomputed()` reads real eggNOG-mapper 2.1.13 TSV output (all five strains committed) and now parses the emapper `KEGG_ko` / `KEGG_Pathway` columns carried by `scripts/convert_emapper_to_tsv.py`. `kegg_pathway_membership()` re-keys `{pathway_id: [locus_tag, ...]}` onto graph gene nodes so acquired/plasmid-borne AMR determinants (missed by the `eco`-bridged UniProt KEGG cross-references) still receive a KEGG pathway layer. `annotate_by_locus_tag()` serves as an 'S' (Function unknown) fallback for genes without an ortholog

### `graph/heterogeneous_graph.py` (625 lines)
- **Purpose**: Build PyTorch Geometric `HeteroData` object with typed nodes and edges
- **Public class**: `HeterogeneousGraphBuilder`
- **Public methods**: `add_gene_nodes()`, `add_protein_nodes()`, `add_genomic_proximity_edges()`, `add_transcriptional_correlation_edges()`, `add_encodes_edges()`, `add_protein_ppi_edges()`, `add_abundance_correlation_edges()`, `add_functional_edges()`, `build()`, `save()`, `summary()`
- **Dependencies**: torch, pandas, numpy, PyTorch Geometric (optional — falls back to dict)
- **Inputs**: Gene IDs, protein IDs, feature DataFrames, edge DataFrames
- **Outputs**: `HeteroData` object → serialized to `heterodata.pt`

### `graph/hypergraph.py` (420 lines)
- **Purpose**: Build biological-module hypergraph as incidence matrix
- **Public class**: `BiologicalHypergraphBuilder`
- **Public methods**: `add_multi_omics_triples()`, `add_kegg_pathway_hyperedges()`, `add_cog_hyperedges()`, `add_go_hyperedges()`, `add_ppi_clusters()`, `add_metabolite_hyperedges()`, `build_incidence_matrix()`, `build_hypergraph_laplacian()`, `save()`, `summary()`
- **Dependencies**: numpy, torch, pandas, networkx
- **Inputs**: Aligned multi-omics data, pathway membership, COG annotations, PPI edges
- **Outputs**: Incidence matrix `H`, Laplacian tensor, CSV exports

### `visualization/network_plots.py` (461 lines)
- **Purpose**: Generate summary figures for all pipeline outputs
- **Public class**: `NetworkVisualizer`
- **Public methods**: `plot_heterogeneous_graph()`, `plot_incidence_matrix()`, `plot_regulation_concordance()`, `plot_cog_distribution()`, `plot_hyperedge_type_distribution()`, `plot_summary_dashboard()`
- **Dependencies**: matplotlib, numpy, pandas, networkx
- **Inputs**: Pipeline summary dicts, DataFrames, numpy arrays
- **Outputs**: PNG figures in `outputs/{strain}/figures/`

### `tests/test_uniprot.py` (494 lines)
- **Purpose**: 50 unit tests for the UniProt annotation module
- **Test classes**: `TestIdentifierNormalization`, `TestMappingCandidateRanking`, `TestCacheManager`, `TestAnnotationParser`, `TestSearchMapping`, `TestAnnotationPipeline`, `TestCacheBehavior`, `TestGoTermInfo`, `TestEdgeCases`, `TestBackwardCompatibility`
- **Dependencies**: pytest, pandas, unittest.mock

---

# 3. External Datasets

## 3.1 Genomic Data — NCBI RefSeq Assembly GCF_900622635.1

| Property | Value |
|----------|-------|
| Dataset name | *Escherichia coli* B36 genome assembly |
| Accession | GCF_900622635.1 |
| Organism | *Escherichia coli* B36 (clinical UPEC isolate) |
| File used | `genomic.gff` (GFF3 format) |
| Biotype | Uropathogenic *E. coli* (UPEC) |
| Replicons | 1 chromosome + plasmids |

### Why Selected
B36 is the specific strain for which matched RNA-seq and proteomics data were generated. The GFF is the central identifier backbone — it provides the `locus_tag` (GeneID) ↔ `protein_id` (WP_ accession) mapping that bridges transcriptomics and proteomics.

### Files Used
| File | Purpose | Key columns |
|------|---------|-------------|
| `genomic.gff` | Gene annotation, CDS extraction | `seqid`, `source`, `type`, `start`, `end`, `strand`, `attributes` |
| `cds_from_genomic.fna` | CDS nucleotide sequences | Not currently used |
| `protein.faa` | Protein amino acid sequences | Not currently used |

### Information Extracted
- **Gene ID** (`locus_tag` attribute, e.g., `EW036_RS00010`)
- **Protein ID** (`protein_id` attribute, e.g., `WP_000673464.1`)
- **Gene symbol** (`gene` attribute, when available)
- **Product description** (`product` attribute)
- **Genomic coordinates**: `start`, `end` (1-based), `strand` (+/-), `replicon` (seqid)
- **CDS length** = `end - start + 1`

### Parsing Module
`GenomicsProcessor.parse_gff()` in `preprocessing/genomics.py`

Parsing algorithm:
```
for each non-comment line in GFF:
    if type == 'CDS':
        parse attributes dict (semicolon-separated key=value pairs)
        extract locus_tag, protein_id, gene, product
        if no protein_id in attributes, try Dbxref=GenBank:{id}
        if no locus_tag, derive from Parent=gene-{tag}
        compute CDS_length = end - start + 1
        append record
filter rows with empty GeneID
sort by (replicon, start)
```

### Dataset Role
| Component | Contribution |
|-----------|-------------|
| Node creation | Gene node IDs, Protein node IDs |
| Edge construction | Genomic proximity edges via `compute_genomic_edges()` |
| Feature engineering | `start`, `CDS_length`, `strand` (one-hot) for gene node features |
| Cross-omics bridge | GeneID ↔ ProteinID mapping for integration |

### Dataset Limitation
The GFF contains the complete gene set (4,200+ CDS), but only genes with matched RNA-seq and proteomics data appear in the final graph (~2,283). Genes without expression data are excluded.

## 3.2 Transcriptomics — GEO GSE152966

| Property | Value |
|----------|-------|
| Dataset name | RNA-seq of *E. coli* B36 in RPMI vs pooled human serum |
| GEO accession | GSE152966 |
| Organism | *Escherichia coli* B36 |
| Experimental objective | Transcriptional response to bloodstream-like environment |
| Conditions | RPMI (baseline) vs 50% pooled human serum (stress) |
| Replicates | 6 per condition |
| Omics type | RNA-seq (gene-level counts) |

### Why Selected
The B36 strain was originally isolated from a urinary tract infection. Serum exposure mimics the transition from urine (RPMI-mimicking) to bloodstream — the key stress encountered during ascending UTI that leads to sepsis. The matched design (same strain, same conditions as proteomics) enables cross-omics concordance analysis.

### Files Used
| File | Format | Purpose | Key columns |
|------|--------|---------|-------------|
| `GSE152966_gene_counts.txt` | TSV | Gene-level raw counts (primary) | GeneID, sample columns |
| `GSE152966_gene_cpmlog2.txt.gz` | TSV.GZ | log2 CPM expression | Not currently used |

### Sample Mapping
| Sample IDs | Condition | Label |
|------------|-----------|-------|
| 50857, 50858, 51033, 51034, 51035, 51036 | RPMI (baseline) | `RNA_RPMI` |
| 50863, 50864, 51037, 51038, 51039, 51040 | Pooled human serum | `RNA_Sera` |

### Information Extracted
- **Gene ID**: Locus tag from first column
- **RNA_RPMI**: Mean of 6 RPMI replicate expression values
- **RNA_Sera**: Mean of 6 serum replicate expression values
- **RNA_logFC**: `RNA_Sera - RNA_RPMI`
- **RNA_Regulation**: `Up` / `Down` / `Stable` based on logFC threshold

### Transformation
The counts file is loaded via `TranscriptomicsProcessor.load_log2_cpm()`. Sample columns are read directly as-is (values are treated as expression values, no log2 conversion applied in the processor — column means are computed directly). `compute_condition_means()` averages replicates and computes logFC.

### Dataset Role
| Component | Contribution |
|-----------|-------------|
| Node features | `RNA_RPMI`, `RNA_Sera`, `RNA_logFC` as gene node attributes |
| Edge construction | Transcriptional correlation edges via Pearson r across replicate columns |
| Regulation analysis | Up/Down/Stable classification for summary dashboard |

## 3.3 Proteomics — SWATH-MS (DIA)

| Property | Value |
|----------|-------|
| Dataset name | SWATH-MS proteomics of *E. coli* B36 in RPMI vs serum |
| Organism | *Escherichia coli* B36 |
| Experimental objective | Proteomic response to bloodstream-like environment |
| Conditions | RPMI (baseline) vs 50% pooled human serum |
| Replicates | 6 per condition |
| Quantification type | Protein area (MS2 fragment ion intensity) |
| File | `Ecoli_B36_merged.xlsx`, sheet `Area - proteins` |

### Why Selected
Matched with the RNA-seq data (same strain, same growth conditions, same replicates design). This enables direct RNA↔protein logFC concordance analysis, which is a core biological feature captured in the `multi_omics_triple` hyperedge agreement score.

### Files Used
| File | Format | Purpose |
|------|--------|---------|
| `Ecoli_B36_merged.xlsx` | Excel | Primary SWATH-MS quantification |

### Column Naming Convention
SWATH-MS columns follow the pattern:
```
{6-digit_sample_id} (180525_P21506_SWATH_{sample_id}.wiff (sample 1))
```

### Information Extracted
- **ProteinID**: From `Protein` column, cleaned (`ref|WP_*` → `WP_*`)
- **Protein_RPMI**: `log2(mean(RPMI replicate areas) + 1)`
- **Protein_Sera**: `log2(mean(serum replicate areas) + 1)`
- **Protein_logFC**: `Protein_Sera - Protein_RPMI`
- **Protein_Regulation**: Up/Down/Stable based on logFC threshold

### Dataset Role
| Component | Contribution |
|-----------|-------------|
| Node features | `Protein_RPMI`, `Protein_Sera`, `Protein_logFC` as protein node attributes |
| Regulation analysis | Up/Down/Stable classification for summary dashboard |

## 3.4 External Knowledge Databases

### UniProt Knowledgebase
- **API endpoint**: `https://rest.uniprot.org/uniprotkb/search`
- **Data retrieved**: 31-column annotation table per protein
- **Key columns**: UniProtID, GeneName, ProteinName, Function, EC_number, GO_BP/MF/CC, InterPro, Pfam, KEGG, SubcellularLocation, AnnotationScore, Reviewed
- **Role**: Bridges RefSeq WP_ → UniProt accessions; provides functional metadata and KEGG cross-references

### STRING Database
- **API endpoint**: `https://string-db.org/api/json`
- **Data retrieved**: PPI edges with `combined_score`, `experimental_score`, `database_score`, `coexpression_score`, `cooccurrence_score`
- **Species fallback**: Primary = CFT073 (199310, UPEC), Fallback = K-12 (511145)
- **Score handling**: API returns scores as 0–1 floats; `required_score` parameter uses 0–1000 scale (400 = medium confidence)
- **Role**: Creates `(protein, ppi, protein)` edges and `ppi_cluster` hyperedges

### KEGG Database
- **API endpoint**: `https://rest.kegg.jp`
- **Organism code**: `eco` (*E. coli* K-12 — best-annotated reference)
- **Mapping chain**: Gene → KO (KEGG Orthology) → Pathway
- **Data retrieved**: `pathway_membership` dict `{pathway_id: [GeneIDs]}`
- **Role**: Creates `(gene, in_pathway, annotation)` edges and `kegg_pathway` hyperedges

### eggNOG Database
- **API status**: Real eggNOG-mapper 2.1.13 output (`emapper.py --dmnd_iterate no`) loaded via `load_precomputed()` for all five strains; `annotate_by_locus_tag()` remains only as an `S`-fallback for genes without an ortholog
- **Role**: Creates `(gene, cog_category, annotation)` edges and `cog_category` hyperedges. The emapper `KEGG_ko` / `KEGG_Pathway` columns are also carried through the TSV (via `scripts/convert_emapper_to_tsv.py`) and exposed by `EggNOGAnnotator.kegg_pathway_membership()`, so the graph gains a KEGG pathway layer for acquired/plasmid-borne AMR determinants that the `eco`-bridged UniProt KEGG cross-references cannot cover.

## 3.5 Dataset Limitations

The current datasets support multi-omics integration and AMR characterization, but **cannot support supervised AMR phenotype prediction** because:

1. **No AMR phenotype labels**: The datasets characterize the transcriptomic/proteomic response to serum stress, not antibiotic exposure. There are no resistance/susceptibility labels for supervised learning.
2. **No held-out test set**: The entire dataset is processed into one graph. No train/test/validation split is implemented.
3. **Small strain count**: Only four *E. coli* strain graphs carry AMR markers (B36, MS_14384, MS_14386, MS_14387) — insufficient for meaningful train/test evaluation.
4. **COG gaps**: A small set of aligned genes per strain has no eggNOG ortholog and falls back to COG 'S' (Function unknown).

---

# 4. Execution Flow

## 4.1 Sequence Diagram

```
main() CLI
  │
  ├── load_config(strain) → config dict
  │
  └── StrainPipeline(config)
        │
        ├── run_genomics()
        │   ├── GenomicsProcessor.parse_gff() → gene_table
        │   ├── gene_table.to_csv('genome_genes.csv')
        │   └── compute_genomic_edges() → genomic_edges CSV
        │
        ├── run_transcriptomics()
        │   ├── TranscriptomicsProcessor.load_log2_cpm()
        │   ├── compute_condition_means() → expression_table
        │   └── expression_table.to_csv('transcriptomics_expression.csv')
        │
        ├── run_proteomics()
        │   ├── Load SWATH-MS Excel
        │   ├── Compute Protein_RPMI, Protein_Sera, Protein_logFC
        │   ├── Assign Protein_Regulation (Up/Down/Stable)
        │   └── proto_agg.to_csv('proteomics_abundance.csv')
        │
        ├── run_integration()
        │   ├── Build ProteinID → GeneID map from genome
        │   ├── Join RNA + Protein on GeneID (inner)
        │   ├── Add genome features (start, end, strand, CDS_length)
        │   └── aligned.to_csv('aligned_multiomics.csv')
        │
        ├── run_annotation()
        │   ├── UniProtAnnotator.fetch_annotations() → uniprot_annotations.csv
        │   ├── STRINGClient (CFT073 + K-12 fallback) → string_ppi.csv
        │   ├── KEGGAnnotator.build_pathway_membership() → pathway_membership
        │   └── EggNOGAnnotator.annotate_by_locus_tag() → cog_annotations
        │
        ├── run_graph_construction()
        │   ├── HeterogeneousGraphBuilder.build()
        │   │   ├── add_gene_nodes()          # 7 features per gene
        │   │   ├── add_protein_nodes()       # 3 features per protein
        │   │   ├── add_genomic_proximity_edges()
        │   │   ├── add_transcriptional_correlation_edges()
        │   │   ├── add_encodes_edges()
        │   │   ├── add_abundance_correlation_edges()
        │   │   ├── add_protein_ppi_edges()
        │   │   ├── add_functional_edges(KEGG)  # in_pathway
        │   │   └── add_functional_edges(COG)   # cog_category
        │   ├── builder.save() → heterodata.pt
        │   └── Export edge CSVs
        │
        ├── run_hypergraph()
        │   ├── BiologicalHypergraphBuilder()
        │   ├── add_multi_omics_triples()
        │   ├── add_kegg_pathway_hyperedges()
        │   ├── add_cog_hyperedges()
        │   ├── add_ppi_clusters()           # Louvain communities
        │   ├── build_incidence_matrix() → H
        │   ├── build_hypergraph_laplacian()
        │   └── save() → hyperedges.csv, incidence_matrix.npy, hypergraph_laplacian.pt
        │
        ├── run_visualization()
        │   └── NetworkVisualizer → PNG figures
        │
        └── run_export()
            ├── aligned_multiomics.csv (re-export)
            ├── gene_features.csv
            ├── Edge list CSVs per relation type
            ├── hypergraph_incidence.npy
            ├── hypergraph_incidence_edges.csv
            ├── heterodata.pt
            └── hypergraph_laplacian.pt
```

## 4.2 CLI Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--strain` | str | `B36` | Strain identifier |
| `--all-strains` | flag | False | Process all configured strains |
| `--output` | str | `outputs` | Base output directory |
| `--correlation` | float | 0.7 | Pearson r threshold for co-expression edges |
| `--skip-annotation` | flag | False | Skip online API calls (UniProt, STRING, KEGG, eggNOG) |
| `--no-cache` | flag | False | Clear all cached API responses and `.pt` files before run |

## 4.3 Strain Configuration (`STRAIN_CONFIGS`)

```python
'B36': {
    'name': 'Escherichia coli B36',
    'ncbi_tax_id': 562,
    'string_species': 511145,                     # K-12 (backward compat)
    'string_species_primary': 199310,             # CFT073 (UPEC)
    'string_species_fallback': 511145,            # K-12 MG1655
    'regulation_std_multiplier': 0.0,             # 0 = absolute threshold; >0 = k * std(logFC)
    'kegg_code': 'eco',
    'genome_gff': BASE_DIR / 'genome/.../genomic.gff',
    'rna_counts_file': BASE_DIR / 'transcriptomic/GSE152966_gene_counts.txt',
    'prot_file': BASE_DIR / 'proteomics/Ecoli_B36_merged.xlsx',
    'prot_sheet': 'Area - proteins',
    'prot_protein_col': 'Protein',
    'rna_rpmi_samples': ['50857', '50858', ...],
    'rna_sera_samples': ['50863', '50864', ...],
    'prot_rpmi_ids': [50973, 50974, ...],
    'prot_sera_ids': [50979, 50980, ...],
}
```

---

# 5. Data Processing Pipeline

## 5.1 Stage 1: Genomics (`run_genomics()`)

**Input**: GFF3 file
**Processor**: `GenomicsProcessor`
**Output files**: `genome_genes.csv`, `genome_genomic_edges.csv`

### Gene Table Columns
| Column | Source | Type | Description |
|--------|--------|------|-------------|
| `GeneID` | `locus_tag` attribute | str | Primary gene identifier (e.g., `EW036_RS00010`) |
| `GeneName` | `gene` attribute | str | Gene symbol (e.g., `secY`) |
| `ProteinID` | `protein_id` attribute | str | RefSeq WP_ accession (e.g., `WP_000673464.1`) |
| `Product` | `product` attribute | str | Functional description |
| `replicon` | seqid column | str | Chromosome/plasmid identifier |
| `start` | start column | int | 1-based start coordinate |
| `end` | end column | int | End coordinate |
| `strand` | strand column | str | `+` or `-` |
| `CDS_length` | computed | int | `end - start + 1` |

### Genomic Proximity Edges
`compute_genomic_edges(max_gap=5000)` creates edges between consecutive genes on the same replicon if intergenic distance ≤ 5000 bp.

**Edge weight**: `1.0 / max(1.0, distance)` — inversely proportional to intergenic gap.

## 5.2 Stage 2: Transcriptomics (`run_transcriptomics()`)

**Input**: Gene counts file (raw counts)
**Processor**: `TranscriptomicsProcessor`
**Output file**: `transcriptomics_expression.csv`

### Computation
```
RNA_RPMI  = mean(counts of 6 RPMI sample columns)
RNA_Sera  = mean(counts of 6 serum sample columns)
RNA_logFC = RNA_Sera - RNA_RPMI
```

### Regulation Call
```
if std_multiplier > 0:
    threshold = std_multiplier * std(RNA_logFC)
else:
    threshold = regulation_threshold (default 1.0)

RNA_Regulation = 'Up'   if RNA_logFC > threshold
               = 'Down' if RNA_logFC < -threshold
               = 'Stable' otherwise
```

## 5.3 Stage 3: Proteomics (`run_proteomics()`)

**Input**: SWATH-MS Excel (`Area - proteins` sheet)
**Output file**: `proteomics_abundance.csv`

### SWATH-MS Path
1. Load Excel sheet
2. Clean Protein IDs: `ref|WP_*` → `WP_*`
3. Build RPMI/Sera column names from sample IDs + filename template
4. Validate columns exist in data
5. Compute:
   ```
   Protein_RPMI  = log2(mean(RPMI replicate areas) + 1)
   Protein_Sera  = log2(mean(serum replicate areas) + 1)
   Protein_logFC = Protein_Sera - Protein_RPMI
   ```

### Regulation Call
Same std-based or absolute threshold logic as transcriptomics.

### Protein Aggregation
Proteins appearing in multiple rows are aggregated by `ProteinID` using mean.

## 5.4 Stage 4: Integration (`run_integration()`)

**Input**: Gene table (genome), expression table (RNA), protein table
**Output file**: `aligned_multiomics.csv`

**Algorithm**:
```
1. Build ProteinID → GeneID map from genome gene_table
2. Map protein data to GeneID via ProteinID → GeneID bridge
3. Drop proteins without a GeneID mapping
4. Inner join: RNA[genes] + Protein[genes] on GeneID
5. Remove duplicate GeneID entries (keep first)
6. Add genome features: start, end, strand, CDS_length, replicon
7. Build gene_protein_map: {GeneID: ProteinID} for graph construction
```

**Aligned table columns**: `GeneID`, `GeneName`, `RNA_RPMI`, `RNA_Sera`, `RNA_logFC`, `ProteinID`, `Protein_RPMI`, `Protein_Sera`, `Protein_logFC`, `Protein_Regulation`, `start`, `end`, `strand`, `CDS_length`, `replicon`

## 5.5 Stage 5: Annotation (`run_annotation()`)

### UniProt
1. Extract RefSeq WP_ accessions from aligned protein IDs
2. `UniProtAnnotator.map_refseq_to_uniprot()` — batch search via REST API
3. `fetch_annotations()` — parse cached entries into 31-column DataFrame
4. Save to `uniprot_annotations.csv`

### STRING
1. Map RefSeq WP_ → UniProt accessions
2. Query CFT073 (199310) first, fallback to K-12 (511145) for additional coverage
3. Merge results, deduplicate pairs by highest `combined_score`
4. Map UniProt accessions back to original WP_ IDs (with version numbers)
5. Save to `string_ppi.csv`

### KEGG
1. Build `GeneID → KEGG gene ID` map from UniProt cross-references
2. `KEGGAnnotator.build_pathway_membership()`:
   - Map genes → KO numbers via KEGG link/ko API
   - Map KO → pathways via link/pathway API
   - Invert to `{pathway_id: [GeneIDs]}`
3. Store in `self.pathway_membership`

### eggNOG
1. `EggNOGAnnotator.annotate_by_locus_tag()` assigns COG 'S' to all genes
2. Store in `self.cog_annotations`

## 5.6 Stage 6: Graph Construction (`run_graph_construction()`)

**Input**: Aligned multi-omics data, annotations, PPIs
**Builder**: `HeterogeneousGraphBuilder`
**Output**: `heterodata.pt` (PyTorch Geometric `HeteroData`)

### Node Counts (B36)
| Node type | Count | Features |
|-----------|-------|----------|
| gene | 2,283 | [RNA_RPMI, RNA_Sera, RNA_logFC, start/1e6, CDS_length/1e3, strand_oh+] |
| protein | 2,283 | [Protein_RPMI, Protein_Sera, Protein_logFC] |
| annotation | 327 | name attribute (326 KEGG + 1 COG) |

### Edge Counts (B36)
| Edge type | Count | Edge attribute |
|-----------|-------|----------------|
| (gene, genomic_proximity, gene) | 1,165 | weight (1/distance) |
| (gene, transcriptional_correlation, gene) | varies | Pearson r |
| (gene, encodes, protein) | 2,283 | 1.0 (binary) |
| (protein, abundance_correlation, protein) | varies | Pearson r |
| (protein, ppi, protein) | 6,949 | combined_score (0–1) |
| (gene, in_pathway, annotation) | 2,100 | — |
| (gene, cog_category, annotation) | 2,283 | — |

See Section 10 for full graph construction documentation.

## 5.7 Stage 7: Hypergraph Construction (`run_hypergraph()`)

**Builder**: `BiologicalHypergraphBuilder`
**Output files**: `hyperedges.csv`, `incidence_matrix.npy`, `hypergraph_nodes.csv`, `hyperedge_type_counts.csv`, `hypergraph_laplacian.pt`

### Hyperedge Types (B36)
| Type | Count | Members |
|------|-------|---------|
| multi_omics_triple | 2,283 | (gene, RNA, protein) per aligned locus |
| kegg_pathway | 230 | Genes in same pathway |
| metabolic_pathway | 219 | Metabolites (+ genes) in same KEGG pathway |
| go_term | 779 | Proteins sharing a GO term (BP/MF/CC) |
| ppi_cluster | 10 | Louvain communities from PPI (min size 3) |
| cog_category | 20 | COG categories |

See Section 10 for full hypergraph documentation.

## 5.8 Stage 8: Visualization (`run_visualization()`)

Generates 6 PNG figures:
- `summary_dashboard.png` — 2x3 panel with layer sizes, regulation pies, graph stats
- `regulation_concordance.png` — RNA vs Protein logFC scatter + agreement histogram
- `cog_distribution.png` — COG category frequency bar chart
- `hyperedge_types.png` — Hyperedge type distribution bar chart
- `incidence_matrix.png` — Hypergraph incidence heatmap
- `heterogeneous_graph.png` — NetworkX layout of sampled subgraph

## 5.9 Stage 9: Export (`run_export()`)

Generates:
- `aligned_multiomics.csv` — re-export of aligned table
- `gene_features.csv` — gene table from genome
- `edges_{src}_{rel}_{dst}.csv` — one per edge type (node indices)
- `hypergraph_incidence.npy` — incidence matrix H
- `hypergraph_incidence_edges.csv` — sparse (node, hyperedge) pairs
- `heterodata.pt` — PyTorch Geometric HeteroData object
- `hypergraph_laplacian.pt` — Normalized hypergraph Laplacian matrix

---

# 6. Module Documentation

## 6.1 `preprocessing/genomics.py` — `GenomicsProcessor`

### `parse_gff(gff_path)`
- **Purpose**: Parse NCBI RefSeq GFF3, extract CDS entries
- **Parameters**: `gff_path` (str) — path to `.gff` file
- **Returns**: `pd.DataFrame` with CDS records
- **Exceptions**: `FileNotFoundError` if GFF missing
- **Algorithm**: Reads line by line, skips comments (`#`), parses only `type == 'CDS'`, parses semicolon-separated attributes into dict
- **Side effects**: Sets `self.gene_table`, calls `_build_protein_id_map()`
- **Time complexity**: O(n) where n = GFF lines

### `compute_genomic_edges(max_gap=5000)`
- **Purpose**: Create adjacency edges for genes within `max_gap` bp on same replicon
- **Parameters**: `max_gap` (int) — maximum intergenic distance
- **Returns**: `pd.DataFrame` with columns `GeneID_A`, `GeneID_B`, `replicon`, `distance`, `weight`
- **Algorithm**:
  ```
  for each replicon:
      sort genes by start position
      for each adjacent pair (g1, g2):
          dist = g2.start - g1.end
          if 0 < dist <= max_gap:
              weight = 1.0 / max(1.0, dist)
              add edge (g1, g2)
  ```

## 6.2 `preprocessing/transcriptomics.py` — `TranscriptomicsProcessor`

### `load_log2_cpm(filepath, sample_map, sheet_name)`
- Detects file format by extension (CSV, TSV, Excel, gzipped)
- Auto-detects GeneID and GeneName columns by header pattern matching
- Stores sample→condition mapping for replicate grouping

### `compute_condition_means(rpmi_samples, sera_samples)`
- **Purpose**: Compute per-condition means, logFC, and regulation calls
- **Returns**: `pd.DataFrame` with columns `GeneID`, `GeneName`, `RNA_RPMI`, `RNA_Sera`, `RNA_logFC`, `RNA_Regulation`
- **Algorithm**:
  ```
  RNA_RPMI  = mean(replicate columns)
  RNA_Sera  = mean(serum columns)
  RNA_logFC = RNA_Sera - RNA_RPMI
  threshold = std_multiplier * std(RNA_logFC)  (if std_multiplier > 0)
           or regulation_threshold              (otherwise)
  RNA_Regulation = np.select([logFC > threshold, logFC < -threshold], ['Up', 'Down'], 'Stable')
  ```

## 6.3 `preprocessing/proteomics.py` — `ProteomicsProcessor`

### `compute_condition_means(rpmi_cols, sera_cols)`
Same structure as transcriptomics, with added log2 transform:
```
Protein_RPMI  = log2(mean(RPMI replicate areas) + 1)
Protein_Sera  = log2(mean(serum replicate areas) + 1)
Protein_logFC = Protein_Sera - Protein_RPMI
```

## 6.4 `annotation/uniprot.py` — `UniProtAnnotator`

### `map_refseq_to_uniprot(wp_accessions, batch_size=30)`
- **Purpose**: Map RefSeq WP_ accessions to UniProt accessions
- **Algorithm**:
  1. Normalize IDs (strip version suffixes)
  2. Check cache → return if valid
  3. Batch search UniProt REST API (OR query per batch: `database:RefSeq {id}`)
  4. Build `{norm_wp: [(uniprot_id, entry)]}` index from search results
  5. Resolve multiple candidates by rank: organism match > reviewed > score
  6. Cache entry list and mapping dict
- **Returns**: `{WP_norm: UniProt_accession}` dict

### `fetch_annotations(wp_accessions)`
- **Purpose**: Fetch complete functional annotations for protein list
- **Algorithm**:
  1. Call `map_refseq_to_uniprot()` (also caches entry objects)
  2. Load cached entries, parse each via `_parse_entry()`
  3. Build 31-column DataFrame
  4. Map UniProt IDs back to WP_ accessions
  5. Validate, deduplicate, compute statistics
- **Returns**: `pd.DataFrame` with 31 columns (see Section 6.4.1)

### Parsed Annotation Columns
| Column | Description |
|--------|-------------|
| `UniProtID` | Primary UniProt accession |
| `GeneName` | Gene symbol |
| `ProteinName` | Recommended protein name |
| `AlternativeNames` | Alternative names (semicolon-separated) |
| `Function` | Functional description |
| `CatalyticActivity` | Catalytic activity comment |
| `EC_number` | Enzyme Commission numbers |
| `ProteinLength` | Amino acid length |
| `AnnotationScore` | UniProt annotation score (0–5) |
| `Reviewed` | Boolean: Swiss-Prot (True) or TrEMBL (False) |
| `EntryType` | Full entry type string |
| `ProteinExistence` | Evidence level |
| `Organism` | Scientific name |
| `SubcellularLocation` | Cellular location |
| `Keywords` | Semicolon-separated keywords |
| `KeywordCount` | Number of keywords |
| `GO_BP`, `GO_MF`, `GO_CC` | Semicolon-separated GO IDs by aspect |
| `GO_BP_count`, `GO_MF_count`, `GO_CC_count`, `GO_total` | GO term counts |
| `InterPro`, `Pfam` | Domain database cross-references |
| `KEGG` | KEGG gene IDs (used for pathway mapping) |
| `RefSeq` | RefSeq cross-references |
| `PDB`, `AlphaFold` | Structure cross-references |
| `ProteinID` | Original WP_ accession |

## 6.5 `annotation/string_api.py` — `STRINGClient`

### `get_string_ids(protein_ids)`
- **Purpose**: Map input IDs (UniProt accessions) to STRING internal IDs
- **API**: `POST /json/get_string_ids`
- **Caching**: Per species_id

### `fetch_interactions(protein_ids)`
- **Purpose**: Fetch PPI network for the mapped STRING IDs
- **Algorithm**:
  1. Check cache → return if exists
  2. Map UniProt IDs → STRING IDs via `get_string_ids()`
  3. Batch query `/json/network` (500 IDs per batch)
  4. Parse edges: `stringId_A`, `stringId_B`, `score`, `escore`, `dscore`, `tscore`, `pscore`
  5. Deduplicate pairs: group by `(sorted A, B)`, keep max `combined_score`
  6. Map STRING IDs back to original UniProt accessions
  7. Cache result
- **Returns**: `pd.DataFrame` with columns `ProteinID_A`, `ProteinID_B`, `combined_score`, `experimental_score`, `database_score`, `coexpression_score`, `cooccurrence_score`
- **Critical note**: The API's `required_score` parameter uses 0–1000 scale (pass `400` as integer). The returned scores are 0–1 floats (e.g., `0.404`). No scaling is applied to stored values.

### Species Fallback (in `main.py`)
```python
# Query primary species (CFT073)
all_ppi = []
for sp in [primary_sp, fallback_sp]:
    client = STRINGClient(species_id=sp)
    ppi_part = client.fetch_interactions(uniprot_ids)
    if not ppi_part.empty:
        all_ppi.append(ppi_part)

# Merge and deduplicate across species
ppi = pd.concat(all_ppi, ignore_index=True)
ppi['_pair_key'] = ppi.apply(lambda r: tuple(sorted([...])), axis=1)
ppi = ppi.loc[ppi.groupby('_pair_key')['combined_score'].idxmax()]
```

## 6.6 `annotation/kegg.py` — `KEGGAnnotator`

### `build_pathway_membership(gene_ids, gene_to_kegg_map=None)`
- **Purpose**: Build `{pathway_id: [GeneIDs]}` from KEGG KO mapping
- **Algorithm**:
  1. `get_ko_for_genes()`: Fetch `organism_code → KO` mappings via `/link/ko/{org}`
     - Direct match on GeneID
     - Via `gene_to_kegg_map` bridge if provided
     - Via short name fallback (last `_` component)
  2. `map_ko_to_pathways()`: Batch KO→pathway via `/link/pathway/{kos}`
  3. Invert to `pathway → [GeneIDs]`
  4. Also fetch pathway names via `list_pathways()`
- **Returns**: `{pathway_id: [GeneID, ...]}`

## 6.7 `annotation/eggnog.py` — `EggNOGAnnotator`

### `annotate_by_locus_tag(gene_ids)`
- **Current behavior**: Assigns COG 'S' (Function unknown) with empty description to every gene
- **Production path**: Replace with `load_precomputed()` using real eggNOG-mapper output

### `COG_CATEGORIES`
25 functional categories (J, K, L, D, V, T, M, N, Z, W, U, O, X, C, G, E, F, H, I, P, Q, R, S, A, B) with full descriptions.

---

# 7. Data Structures

## 7.1 `aligned` DataFrame (from `run_integration()`)
Shape: ~2283 × 15

| Column | Type | Source | Description |
|--------|------|--------|-------------|
| `GeneID` | str | Genome | Locus tag |
| `GeneName` | str | Genome | Gene symbol |
| `RNA_RPMI` | float | RNA | Mean RPMI expression |
| `RNA_Sera` | float | RNA | Mean serum expression |
| `RNA_logFC` | float | RNA | Serum − RPMI |
| `RNA_Regulation` | str | RNA | Up/Down/Stable |
| `ProteinID` | str | Proteomics | WP_ accession |
| `Protein_RPMI` | float | Proteomics | Mean RPMI abundance (log2) |
| `Protein_Sera` | float | Proteomics | Mean serum abundance (log2) |
| `Protein_logFC` | float | Proteomics | Serum − RPMI |
| `Protein_Regulation` | str | Proteomics | Up/Down/Stable |
| `start` | int | Genome | Gene start |
| `end` | int | Genome | Gene end |
| `strand` | str | Genome | + or − |
| `CDS_length` | int | Genome | CDS length |
| `replicon` | str | Genome | Chromosome/plasmid |

## 7.2 `gene_to_protein_map` dict
`{GeneID: ProteinID}` — built in `run_integration()` from genome gene_table. Used by `add_encodes_edges()`.

## 7.3 `pathway_membership` dict
`{pathway_id: [GeneID, ...]}` — built by `KEGGAnnotator.build_pathway_membership()`. Used for KEGG graph edges and hypergraph hyperedges.

## 7.4 `cog_annotations` dict
`{GeneID: {COG_category: str, COG_classes: [str], Description: str, Preferred_name: str}}` — built by `EggNOGAnnotator`. Used for COG graph edges and hypergraph hyperedges.

## 7.5 HeteroData Object (PyG)
```
{
    'gene': {
        'x': torch.Tensor [2283, 7]   # [RNA_RPMI, RNA_Sera, RNA_logFC, start/1e6, CDS_length/1e3, strand_oh[0], strand_oh[1]]
    },
    'protein': {
        'x': torch.Tensor [2283, 3]   # [Protein_RPMI, Protein_Sera, Protein_logFC]
    },
    'annotation': {
        'name': ['eco00010', 'eco00020', ..., 'S']  # 327 annotation names
    },
    ('gene', 'genomic_proximity', 'gene'): {
        'edge_index': torch.LongTensor [2, 1165],
        'edge_attr': torch.FloatTensor [1165, 1]
    },
    ('gene', 'transcriptional_correlation', 'gene'): {
        'edge_index': torch.LongTensor [2, N],
        'edge_attr': torch.FloatTensor [N, 1]
    },
    ('gene', 'encodes', 'protein'): {
        'edge_index': torch.LongTensor [2, 2283],
        'edge_attr': torch.FloatTensor [2283, 1]
    },
    ('protein', 'abundance_correlation', 'protein'): {
        'edge_index': torch.LongTensor [2, M],
        'edge_attr': torch.FloatTensor [M, 1]
    },
    ('protein', 'ppi', 'protein'): {
        'edge_index': torch.LongTensor [2, 6949],
        'edge_attr': torch.FloatTensor [6949, 1]
    },
    ('gene', 'in_pathway', 'annotation'): {
        'edge_index': torch.LongTensor [2, 2100]
    },
    ('gene', 'cog_category', 'annotation'): {
        'edge_index': torch.LongTensor [2, 2283]
    }
}
```

## 7.6 Hypergraph Data Structures

### `hyperedges` list
List of dicts, each:
```python
{
    'hyperedge_id': int,
    'type': str,              # 'multi_omics_triple', 'kegg_pathway', 'cog_category', 'ppi_cluster'
    'name': str,              # e.g., 'Triple_EW036_RS00010', 'COG_S: Function unknown'
    'nodes': [str, ...],      # node IDs
    'node_indices': [int, ...],  # indices into node_list
    'features': dict,
}
```

### `H` incidence matrix
`np.ndarray [n_nodes, n_hyperedges]` — binary: `H[i,j] = 1` if node i ∈ hyperedge j

### `hypergraph_laplacian` tensor
`torch.Tensor [n_nodes, n_nodes]` — computed as:
```
D_v = diag(H @ 1)          # node degrees
D_e = diag(1^T @ H)        # hyperedge degrees
L = D_v^{-1/2} @ H @ D_e^{-1} @ H^T @ D_v^{-1/2}
```

---

# 8. Feature Engineering

## 8.1 Gene Node Features (7 dimensions)

| Feature | Formula | Scaling | Source |
|---------|---------|---------|--------|
| `RNA_RPMI` | Mean of RPMI replicates | Raw | Transcriptomics |
| `RNA_Sera` | Mean of serum replicates | Raw | Transcriptomics |
| `RNA_logFC` | `RNA_Sera - RNA_RPMI` | Raw | Computed |
| `start_norm` | `start / 1e6` | Division by 1e6 | Genome |
| `length_norm` | `CDS_length / 1e3` | Division by 1e3 | Genome |
| `strand_plus` | 1.0 if `strand == '+'` else 0.0 | One-hot | Genome |
| `strand_minus` | 0.0 if `strand == '+'` else 1.0 | One-hot | Genome |

## 8.2 Protein Node Features (3 dimensions)

| Feature | Formula | Scaling | Source |
|---------|---------|---------|--------|
| `Protein_RPMI` | `log2(mean(RPMI areas) + 1)` | log2 | Proteomics |
| `Protein_Sera` | `log2(mean(serum areas) + 1)` | log2 | Proteomics |
| `Protein_logFC` | `Protein_Sera - Protein_RPMI` | Raw | Computed |

## 8.3 Missing Values
- Missing RNA features: default to 0.0
- Missing protein features: default to 0.0
- Missing genome features: default to 0.0
- NaN in correlation matrices: mean-imputed per column

---

# 9. Biological Relationship Modeling

## 9.1 Intra-Omics Relationships

### Genomics — Genomic Proximity
- **Relationship**: Consecutive genes on the same replicon within 5 kbp
- **Edge type**: `(gene, genomic_proximity, gene)`
- **Weight**: `1 / max(1, distance_in_bp)` — proximal genes get higher weight
- **Biological meaning**: Adjacent genes are often co-transcribed (operons), share regulatory regions, or are functionally related
- **Implementation**: `GenomicsProcessor.compute_genomic_edges(max_gap=5000)`

### Transcriptomics — Co-Expression
- **Relationship**: Genes with correlated expression across 12 biological replicates (6 RPMI + 6 serum)
- **Edge type**: `(gene, transcriptional_correlation, gene)`
- **Weight**: Pearson correlation coefficient `r`
- **Threshold**: `|r| >= correlation_threshold` (default 0.7)
- **Biological meaning**: Co-expressed genes are often co-regulated, share transcription factors, or participate in the same pathway
- **Implementation**: `HeterogeneousGraphBuilder.add_transcriptional_correlation_genes()`
  - Extracts replicate columns (all numeric columns except means/logFC/regulation)
  - Computes full Pearson correlation matrix
  - Applies upper-triangle mask, threshold filter
  - Mean-impures NaNs per column before `np.corrcoef()`

### Proteomics — Protein-Protein Interactions (STRING)
- **Relationship**: Physical or functional protein interactions
- **Edge type**: `(protein, ppi, protein)`
- **Weight**: `combined_score` (0–1 float from STRING API)
- **Biological meaning**: Direct physical binding, co-complex membership, or functional association
- **Implementation**: `STRINGClient.fetch_interactions()` + species fallback in `main.py`

### Proteomics — Abundance Correlation
- **Relationship**: Proteins with correlated abundance across conditions
- **Edge type**: `(protein, abundance_correlation, protein)`
- **Weight**: Pearson r across `Protein_RPMI` and `Protein_Sera` values
- **Implementation**: `HeterogeneousGraphBuilder.add_abundance_correlation_edges()`

## 9.2 Inter-Omics Relationships

### Gene → Protein (Central Dogma)
- **Relationship**: Each gene encodes a protein product
- **Edge type**: `(gene, encodes, protein)`
- **Weight**: 1.0 (binary — existence of translation relationship)
- **Bridge**: Built from genome GFF CDS entries mapping `locus_tag` → `protein_id`
- **Implementation**: `HeterogeneousGraphBuilder.add_encodes_edges()`

### Gene → Transcriptomics
- **Design choice**: RNA expression values are **attached as gene node features**, not as separate RNA nodes
- **Rationale**: In bacteria, transcription and translation are coupled; no splicing occurs. A separate RNA node type would add a redundant `gene → transcribes → RNA → translates → protein` chain. Instead, RNA features are attributes on the gene node, and the relationship is captured in the `multi_omics_triple` hyperedge

### Protein → Proteomics
- Abundance values are attached as direct protein node features

### Gene/Protein → Functional Annotation
- KEGG pathways and COG categories are represented as `annotation` nodes
- **Edge type**: `(gene, in_pathway, annotation)` and `(gene, cog_category, annotation)`
- These are bipartite edges linking biological entities to their functional groupings
- **Biological meaning**: A gene participates in a metabolic pathway or belongs to a functional ortholog group

## 9.3 Relationship Preservation During Graph Fusion

Graph fusion is implemented in `HeterogeneousGraphBuilder.build()`:

```
Input:
    gene_ids, protein_ids, rna_data, genome_data, protein_data,
    genomic_edges, gene_protein_map, ppi_edges,
    gene_to_ko, ko_to_pathway, gene_to_cog

Process:
    1. Create gene nodes with features from RNA + genome
    2. Create protein nodes with features from proteomics
    3. Add intra-omics edges (genomic, transcriptional, abundance, PPI)
    4. Add inter-omics edges (encodes, in_pathway, cog_category)

Output:
    HeteroData with all node types, features, and all edge types
    in a single unified data structure
```

Preservation guarantees:
- **Node identity preserved**: Each gene and protein has a unique index via `_map_nodes()`
- **Edge types preserved**: Each biological relationship is a separate edge type in HeteroData
- **Features preserved**: All omics-derived features attached to appropriate node type
- **Cross-layer mappings preserved**: `encodes` edges maintain gene↔protein correspondence
- **Annotation mappings preserved**: `in_pathway` and `cog_category` edges link entities to functional groups

---

# 10. Graph Construction

## 10.1 Heterogeneous Graph

### Purpose
Multi-relational graph for PyTorch Geometric heterogeneous GNNs (RGCN, HGT, etc.)

### Node Schema
```
┌──────────┐     ┌──────────┐     ┌──────────────┐
│   gene   │     │ protein  │     │  annotation  │
│  [7 feats]│     │ [3 feats]│     │  [name attr] │
└────┬─────┘     └────┬─────┘     └──────┬───────┘
     │                │                  │
     │   encodes      │                  │
     ├────────────────┤                  │
     │   (gene→prot)  │                  │
     │                │                  │
     │   in_pathway   │   cog_category   │
     └────────────────┴──────────────────┘
     (gene→annotation)  (gene→annotation)
```

### Edge Construction Algorithms

#### `add_genomic_proximity_edges(genomic_edges)`
**Input**: DataFrame with `GeneID_A`, `GeneID_B`, `weight`
**Output**: Edge index tensor [2, N] with edge_attr
**Complexity**: O(E) where E = number of consecutive gene pairs

#### `add_transcriptional_correlation_edges(gene_ids, rna_data)`
**Input**: Gene IDs, RNA data with replicate columns
**Algorithm**:
1. Identify replicate columns (all numeric excluding means/logFC/regulation)
2. Extract replicate matrix for aligned genes
3. Mean-impute NaN per column
4. `corr_mat = np.corrcoef(rep_mat)`
5. Upper-triangle mask: `|corr| >= threshold`
6. Build pairs and edge attributes
**Complexity**: O(n²) where n = number of genes — dominated by `np.corrcoef`

#### `add_encodes_edges(gene_protein_map)`
**Input**: `{GeneID: ProteinID}`
**Output**: Edge index with constant edge_attr = 1.0
**Complexity**: O(n) where n = number of mapped pairs

#### `add_protein_ppi_edges(ppi_edges)`
**Input**: PPI DataFrame with `ProteinID_A`, `ProteinID_B`, `combined_score`
**Output**: Edge index with score as edge_attr
**Complexity**: O(E) where E = PPI edges

#### `add_abundance_correlation_edges(protein_ids, protein_data)`
**Input**: Protein IDs, protein data with `Protein_RPMI`, `Protein_Sera`
**Algorithm**: Same as transcriptional correlation but only 2 replicates available

#### `add_functional_edges(gene_to_pathways, path_to_name, node_type, relation)`
**Input**: `{GeneID: [annotation_id]}`, annotation name map
**Algorithm**:
1. Collect unique annotation IDs → register as `annotation` node type
2. Build `(GeneID, annotation_id)` pairs
3. Add edge index
4. Append annotation names to `annotation.name`
**Complexity**: O(E) where E = number of gene-annotation assignments

## 10.2 Hypergraph

### Purpose
Higher-order biological modules for Hypergraph Neural Networks (HGNN).

### Incidence Matrix
`H ∈ {0,1}^(10552 × 3541)` — binary incidence of nodes in hyperedges (B36)

### Hyperedge Construction Algorithms

#### `add_multi_omics_triples(aligned_data)`
One hyperedge per aligned locus:
```
nodes = [GeneID, f"{GeneID}_RNA", ProteinID]
features = {'RNA_logFC': ..., 'Protein_logFC': ..., 'agreement_score': RNA_logFC × Protein_logFC}
```
**Total**: 2,283 hyperedges

#### `add_kegg_pathway_hyperedges(pathway_membership)`
One hyperedge per KEGG pathway with all participating genes as nodes.
**Total**: 206 hyperedges

#### `add_cog_hyperedges(cog_annotations, gene_to_protein)`
One hyperedge per COG category with all genes in that category (+ their protein IDs).
**Total**: 1 hyperedge (only COG 'S' is populated)

#### `add_ppi_clusters(ppi_edges, protein_list, min_cluster_size=3)`
1. Build NetworkX graph from PPI edges with `combined_score` as edge weight
2. Run Louvain community detection (`nx.community.louvain_communities` with seed=42, resolution=1.0)
3. Each community ≥ 3 nodes becomes a `ppi_cluster` hyperedge
**Total**: 10 hyperedges (sizes: 120, 116, 88, 69, 58, 49, 25, 11, 3, 3)

#### `add_go_hyperedges(protein_to_go, go_info)`
One hyperedge per GO term shared by ≥ 2 proteins. Built in `run_hypergraph()` from the UniProt annotation `GO_BP` / `GO_MF` / `GO_CC` columns (terms split on `;`, names/aspects from the annotator's `go_terms` registry).
**Total**: 779 hyperedges

### Laplacian
Standard normalized hypergraph Laplacian:
```
L = D_v^{-1/2} H D_e^{-1} H^T D_v^{-1/2}
```
Where:
- `D_v = diag(sum(H, axis=1))` — node degree diagonal
- `D_e = diag(sum(H, axis=0))` — hyperedge degree diagonal

Saved as `hypergraph_laplacian.pt` (6849 × 6849 symmetric tensor).

### Download/Export Format
- `incidence_matrix.npy`: Binary H matrix
- `hypergraph_incidence_edges.csv`: Sparse (node_idx, hyperedge_idx) pairs
- `hyperedges.csv`: Hyperedge metadata (ID, type, name, nodes, size)
- `hypergraph_nodes.csv`: Node IDs with inferred types (genome/proteome/transcriptome)
- `hyperedge_type_counts.csv`: Count per type

---

# 11. Configuration

## 11.1 Constants
| Constant | Value | Location | Purpose |
|----------|-------|----------|---------|
| `STRING_BASE` | `https://string-db.org/api` | `string_api.py:22` | STRING API root |
| `UNIPROT_BASE` | `https://rest.uniprot.org/uniprotkb` | `uniprot.py:24` | UniProt API root |
| `KEGG_BASE` | `https://rest.kegg.jp` | `kegg.py:26` | KEGG API root |
| `BATCH_SIZE` | 30 | `uniprot.py:27` | UniProt search batch size |
| `MAX_RETRIES` | 5 | `uniprot.py:28` | HTTP retry count |
| `B36_RPMI_SAMPLE_IDS` | `[50973-50978]` | `proteomics.py:26` | Proteomics sample IDs |
| `B36_SERA_SAMPLE_IDS` | `[50979-50984]` | `proteomics.py:27` | Proteomics sample IDs |

## 11.2 Thresholds
| Parameter | Default | Location | Purpose |
|-----------|---------|----------|---------|
| `regulation_threshold` | 1.0 | `transcriptomics.py:58`, `proteomics.py:57` | Absolute logFC cutoff |
| `std_multiplier` | 0.0 | `transcriptomics.py:59`, `proteomics.py:58` | If >0: k × std(logFC) replaces absolute threshold |
| `correlation_threshold` | 0.7 | `main.py:106`, `heterogeneous_graph.py:56` | Minimum |Pearson r| for co-expression/abundance edges |
| `confidence_threshold` | 400 | `string_api.py:47` | STRING minimum score (0-1000 API scale) |
| `min_cluster_size` | 3 | `hypergraph.py:299` | Minimum Louvain community for PPI hyperedge |
| `resolution` | 1.0 | `hypergraph.py:300` | Louvain resolution (higher = smaller communities) |
| `max_gap` | 5000 | `genomics.py:141` | Max intergenic distance (bp) for proximity edges |

## 11.3 Directories
| Path | Purpose |
|------|---------|
| `outputs/{strain}/` | Strain-specific output |
| `outputs/{strain}/figures/` | Generated figures |
| `outputs/{strain}/graph/` | Graph and hypergraph serialized files |
| `outputs/cache/` | API response cache (JSON) |

---

# 12. Dependencies

## 12.1 Python Packages
| Package | Version | Usage |
|---------|---------|-------|
| Python | ≥ 3.10 | Runtime |
| pandas | ≥ 1.5 | DataFrame operations |
| numpy | ≥ 1.24 | Matrix operations |
| torch | ≥ 2.0 | Tensor operations, HeteroData |
| networkx | ≥ 3.0 | Louvain community detection, graph visualization |
| requests | ≥ 2.28 | REST API calls |
| matplotlib | ≥ 3.7 | Visualization |
| (optional) torch_geometric | ≥ 2.3 | PyG HeteroData (falls back to dict) |
| (optional) openpyxl | — | Excel file reading |
| (optional) pytest | ≥ 7.0 | Unit tests |

## 12.2 External APIs
| API | Endpoint | Rate Limit | Dependency |
|-----|----------|------------|------------|
| UniProt | `rest.uniprot.org/uniprotkb/search` | ~3 req/s | Online annotation |
| STRING | `string-db.org/api/json` | None documented | Online annotation |
| KEGG | `rest.kegg.jp` | ~1 req/s (no concurrent) | Online annotation |

## 12.3 Data File Dependencies
| File | Source | Required |
|------|--------|----------|
| `genomic.gff` | NCBI RefSeq | Yes |
| `GSE152966_gene_counts.txt` | GEO | Yes |
| `Ecoli_B36_merged.xlsx` | SWATH-MS | Yes |

---

# 13. Output Artifacts

All files are written to `outputs/{strain}/` or subdirectories.

## 13.1 Core Data Files

| File | Format | Dimensions | Producer | Consumer |
|------|--------|------------|----------|----------|
| `genome_genes.csv` | CSV | ~4200 × 9 | `run_genomics()` | `run_integration()` |
| `genome_genomic_edges.csv` | CSV | ~4000 × 5 | `run_genomics()` | Graph construction |
| `transcriptomics_expression.csv` | CSV | ~4200 × 6 | `run_transcriptomics()` | `run_integration()` |
| `proteomics_abundance.csv` | CSV | ~2300 × 5 | `run_proteomics()` | `run_integration()` |
| `aligned_multiomics.csv` | CSV | ~2283 × 15 | `run_integration()` | Graph + Hypergraph |
| `uniprot_annotations.csv` | CSV | 1366 × 31 | `run_annotation()` | KEGG bridge |
| `string_ppi.csv` | CSV | 6949 × 7 | `run_annotation()` | Graph + Hypergraph |

## 13.2 Graph Files (`graph/`)

| File | Format | Dimensions | Description |
|------|--------|------------|-------------|
| `heterodata.pt` | PyTorch | See 7.5 | Full heterogeneous graph |
| `heterodata_metadata.txt` | Text | — | Human-readable graph schema |
| `edges_gene_genomic_proximity_gene.csv` | CSV | N × 3 | Edge list (indices) |
| `edges_gene_transcriptional_correlation_gene.csv` | CSV | N × 3 | Edge list |
| `edges_gene_encodes_protein.csv` | CSV | 2283 × 3 | Edge list |
| `edges_gene_in_pathway_annotation.csv` | CSV | 2100 × 3 | Edge list |
| `edges_gene_cog_category_annotation.csv` | CSV | 2283 × 3 | Edge list |
| `edges_protein_ppi_protein.csv` | CSV | 6949 × 3 | Edge list |
| `edges_protein_abundance_correlation_protein.csv` | CSV | N × 3 | Edge list |

## 13.3 Hypergraph Files (`graph/`)

| File | Format | Dimensions | Description |
|------|--------|------------|-------------|
| `hyperedges.csv` | CSV | 3541 × 5 | Hyperedge metadata |
| `incidence_matrix.npy` | NPY | 10552 × 3541 | Binary incidence (B36) |
| `hypergraph_nodes.csv` | CSV | 6849 × 2 | Node ID → type |
| `hyperedge_type_counts.csv` | CSV | 4 × 2 | Count per type |
| `hypergraph_laplacian.pt` | PyTorch | 6849 × 6849 | Normalized Laplacian |

## 13.4 Figures (`figures/`)

| File | Description |
|------|-------------|
| `summary_dashboard.png` | 6-panel overview dashboard |
| `regulation_concordance.png` | RNA vs Protein logFC scatter |
| `incidence_matrix.png` | Hypergraph incidence heatmap |
| `cog_distribution.png` | COG category bar chart |
| `hyperedge_types.png` | Hyperedge type distribution |
| `heterogeneous_graph.png` | Sampled subgraph network plot |

---

# 14. Current State of the Project

## 14.1 Completed Modules
- **All 9 pipeline stages**: Genomics through Export run end-to-end
- **50 unit tests**: All passing for UniProt annotation module
- **STRING PPI fallback**: Primary CFT073 + fallback K-12 with cross-species dedup
- **Louvain community detection**: Replaced connected components with weighted Louvain
- **Regulation threshold**: Data-driven std-based threshold option

## 14.2 Validated Modules
- Full B36 pipeline runs in ~7.5 minutes (~450 s, dominated by Metabolomics KEGG lookup and Visualization)
- All output files verified: correct dimensions, valid file formats
- HeteroData Laplacian verified symmetric
- PPI scores verified in 0–1 range
- Edge CSV columns verified unique
- All five strains (B36, MS_14384, MS_14385, MS_14386, MS_14387) produce complete outputs
- Real eggNOG COG annotations loaded for all strains (`load_precomputed`) and resolved onto locus-tag gene nodes
- AMR knowledge layer verified (manifest → hyperedges → highlight figures); `pytest` suite: 73 passing

## 14.3 Partially Implemented Modules
- None outstanding; MaxQuant fallback support was removed (SWATH-MS is the sole proteomics source).

## 14.4 Known Issues
1. **No separate RNA node type**: RNA features are embedded as gene node attributes. This is a design choice, not a bug, but users expecting explicit RNA nodes must add them.
2. **Annotation node redundancy**: KEGG pathways appear both as annotation nodes (hetero graph) and as hyperedges (hypergraph). Same information, two representations.
3. **Hypergraph flattens node types**: Gene and protein IDs share the same index space in the incidence matrix. A `kegg_pathway` hyperedge connecting only genes cannot be distinguished from a `ppi_cluster` hyperedge connecting only proteins at the matrix level.
4. **COG gap for a few genes**: ~3–15 aligned genes per strain have no eggNOG-mapper ortholog and fall back to COG `S` (Function unknown) via `annotate_by_locus_tag()`.

## 14.5 Current Limitations
1. **No GNN training**: The pipeline stops at graph export. No `train.py`, `model.py`, or `predict.py` modules exist.
2. **No train/test split**: The entire dataset is processed into one graph without any split for evaluation.
3. **No supervised AMR prediction**: With only four *E. coli* strain graphs, meaningful train/test evaluation is not feasible (`AMR_GRAPH_INTEGRATION.md` §8.4).
4. **Multi-species graphs deferred**: *K. pneumoniae*, *S. aureus*, *S. pyogenes* from the paper are not merged at this stage (§8.1).
5. **Data-driven threshold unused**: `regulation_std_multiplier` is set to 0.0 (off) in the B36 config — absolute threshold of 1.0 is still the default.

---

# 15. Future Work

Clearly implied by the current implementation:

1. **Implement GNN training module** — A `train.py` that loads `heterodata.pt` and runs RGCN/HGT/GCN for node classification or link prediction.

2. **Add UniProt features to protein nodes** — Enrich protein node features with ProteinLength, EC_number, KeywordCount, GO counts from the cached `uniprot_annotations.csv`.

3. **Multi-species ortholog representation** — eggNOG OG mapping across species to compare AMR mechanisms across *E. coli* and other species.

4. **Enable std-based threshold** — Set `regulation_std_multiplier: 1.5` in `STRAIN_CONFIGS` for data-driven regulation calls across strains.

> **Amendment (AMR phase):** This document predates the AMR + entity-resolution work.
> See **`AMR_GRAPH_INTEGRATION.md`** for the full documentation of Phase 1 (WP / COG
> entity resolution), Phase 2 (AMR knowledge hyperedges + visualization), and Phase 3
> (cross-strain AMR analysis), plus corrections made to the base pipeline. All five
> strains now carry real eggNOG-mapper 2.1.13 annotations; genome-only AMR markers
> remain represented as graph nodes and are drawn outlined in `amr_highlight.png`
> (see `scripts/amr_figures.py`).
