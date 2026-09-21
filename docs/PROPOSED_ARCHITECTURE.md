# Proposed Architecture — Multi-Omics Heterogeneous Graph & Hypergraph Fusion for Bacterial Systems Biology (with AMR Knowledge Layer)

**Status:** Proposed design document.
**Scope:** Full-stack proposal — data sources → per-layer processing → central-dogma alignment → external knowledge enrichment → heterogeneous graph → hypergraph → analysis scripts → export for GNN training. Describes the proposed architecture in detail and cross-references the existing implementation (`multiomics_graph/`) where the design is already realized.

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Design Goals & Non-Goals](#2-design-goals--non-goals)
3. [Conceptual Model](#3-conceptual-model)
4. [Repository Layout (Proposed)](#4-repository-layout-proposed)
5. [Data Sources](#5-data-sources)
6. [Pipeline Stage 1 — Per-Strain Configuration & Discovery](#6-pipeline-stage-1--per-strain-configuration--discovery)
7. [Pipeline Stage 2 — Omics Preprocessing](#7-pipeline-stage-2--omics-preprocessing)
8. [Pipeline Stage 3 — Cross-Omics Alignment & Integration](#8-pipeline-stage-3--cross-omics-alignment--integration)
9. [Pipeline Stage 4 — External Knowledge Enrichment](#9-pipeline-stage-4--external-knowledge-enrichment)
10. [Pipeline Stage 5 — Heterogeneous Graph Construction](#10-pipeline-stage-5--heterogeneous-graph-construction)
11. [Pipeline Stage 6 — Hypergraph Construction](#11-pipeline-stage-6--hypergraph-construction)
12. [Pipeline Stage 7 — AMR Knowledge Layer](#12-pipeline-stage-7--amr-knowledge-layer)
13. [Pipeline Stage 8 — Analysis Scripts](#13-pipeline-stage-8--analysis-scripts)
14. [Pipeline Stage 9 — Visualization](#14-pipeline-stage-9--visualization)
15. [Pipeline Stage 10 — Export & Downstream Interface](#15-pipeline-stage-10--export--downstream-interface)
16. [Data Model Summary](#16-data-model-summary)
17. [Design Decisions & Rationale](#17-design-decisions--rationale)
18. [Robustness & Correctness Measures](#18-robustness--correctness-measures)
19. [Implementation Status](#19-implementation-status)
20. [Open Questions & Future Work](#20-open-questions--future-work)
21. [Glossary](#21-glossary)

---

## 1. Executive Summary

The architecture transforms matched multi-omics bacterial measurements — **genomics, transcriptomics, proteomics, metabolomics** — into a unified computational representation composed of two complementary structures:

1. A **heterogeneous graph** (`HeteroData`) with typed biological nodes (gene, protein, metabolite, annotation) and typed, semantically distinct edges (genomic proximity, co-expression, encoding, PPI, pathway membership, COG, GO, AMR association, metabolic relations).
2. A **hypergraph** whose hyperedges encode higher-order biological modules (multi-omics central-dogma triples, KEGG pathways, COG categories, GO terms, PPI Louvain communities, metabolite–pathway joins, and curated **AMR mechanism** classes).

A curated **AMR knowledge layer** attaches known antimicrobial-resistance determinants to the graph as a biological annotation layer, enabling *characterization and analysis* of resistance mechanisms (not supervised phenotype prediction). The output is a file-based representation ready for downstream GNN / HGNN training (RGCN, HAN, HGNN) and mechanism-exploration analyses.

The design follows three principles:

- **Central-dogma as a computational graph** — a `gene --encodes--> protein` bridge connects layers so information flows across the genome → transcript → protein → metabolite axis.
- **Relational preservation over concatenation** — intra- and inter-omics relationships are first-class structures rather than flattened features.
- **Honest missing data** — genome-only markers (resistance genes without RNA/protein quantitation) remain represented; absent layers are reported as `unavailable` / `not_observed`, never fabricated zeros.

---

## 2. Design Goals & Non-Goals

### Goals
- Build a **layer-aligned, cross-omics graph** per bacterial strain from real data.
- Enrich with external knowledge (UniProt, STRING, KEGG, eggNOG/COG) so graph edges are biologically meaningful and cover genes invisible to a single bridging strategy.
- Provide **dual graph + hypergraph** representations for pairwise and set-level (module) reasoning.
- Attach **curated AMR knowledge** as an annotation layer; support cross-strain AMR presence/synergy reports and multi-hop molecular-context tracing.
- Produce **standardized export artifacts** (CSV, NPY, PT) consumable by GNN toolkits.
- Remain **dataset-agnostic**: adding a strain or species requires configuration, not code changes.

### Non-Goals (explicitly out of scope for the current design)
- **Supervised AMR phenotype prediction** — requires matched resistance-phenotype labels (MICs / binary), which the primary dataset does not provide. The current *E. coli* strain set (≤5 strain graphs) is also too small to supervise GNN training.
- Causal claims from path tracing — the tracing framework is **explanatory/characterization**, not inference of causation.
- Multi-species graphs (*K. pneumoniae*, *S. aureus*, *S. pyogenes*) in the first release — deferred (see §20).

---

## 3. Conceptual Model

Biological entities and their relationships are modeled with the **central dogma as the backbone**:

```
  GENOME          TRANSCRIPTOME          PROTEOME           METABOLOME
  ┌────────┐   expr  ┌────────┐  encodes ┌────────┐   catalyses / in  ┌────────────┐
  │  gene  │ ──────► │  RNA   │ ───────► │ protein│ ◄───────────────│ metabolite │
  └───┬────┘         └────────┘          └───┬────┘                 └────────────┘
      │                                      │
      ├── genomic proximity (gene–gene)      ├── PPI (protein–protein, STRING)
      ├── co-expression (gene–gene)          ├── abundance correlation
      └── pathway/COG/AMR membership         └── GO / KEGG / cluster membership
                                             └── AMR mechanism (curated)
```

- **Gene** = genomic locus (features: position, strand, CDS length + RNA state).
- **Protein** = product of a gene via the `encodes` edge (features: abundance, logFC).
- **Metabolite** = small-molecule layer (features: abundance, m/z, RT), bridged to proteins via EC-number (enzymes) and to annotations via KEGG pathways.
- **Annotation** = knowledge nodes (KEGG pathway, COG category, GO term, PPI cluster, AMR mechanism class) that group multiple biological entities.

The same entities reappear in the **hypergraph** as nodes, connected by module-level hyperedges.

---

## 4. Repository Layout (Proposed)

```
agent/
├── multiomics_graph/                 # Main Python package
│   ├── main.py                       # StrainPipeline orchestrator + CLI
│   ├── preprocessing/                # Parsers/normalizers per omics layer
│   │   ├── genomics.py               #   GenomicsProcessor
│   │   ├── transcriptomics.py        #   TranscriptomicsProcessor
│   │   ├── proteomics.py             #   ProteomicsProcessor
│   │   └── metabolomics.py           #   MetabolomicsProcessor
│   ├── annotation/                   # External knowledge retrievers
│   │   ├── uniprot.py                #   UniProtAnnotator (map + fetch + GO)
│   │   ├── string_api.py             #   STRINGClient (PPI, species fallback)
│   │   ├── kegg.py                   #   KEGGAnnotator (gene→KO→pathway)
│   │   ├── kegg_compound.py          #   KEGGCompoundResolver (ChEBI→compound→EC/path)
│   │   ├── eggnog.py                 #   EggNOGAnnotator (COG + KEGG recovery)
│   │   └── identifiers.py            #   Identifier normalization (WP versioning)
│   ├── amr/                          # AMR knowledge layer
│   │   ├── amr_manifest.json         #   Curated per-strain determinants
│   │   └── amr.py                    #   AMRManifest loader/query
│   ├── graph/
│   │   ├── heterogeneous_graph.py    #   HeterogeneousGraphBuilder (PyG HeteroData)
│   │   └── hypergraph.py             #   BiologicalHypergraphBuilder (incidence + Laplacian)
│   ├── visualization/network_plots.py#   NetworkVisualizer (incl. AMR highlight)
│   ├── scripts/                      # Analysis & maintenance tooling
│   │   ├── run_eggnog_mapper.py
│   │   ├── convert_emapper_to_tsv.py
│   │   ├── recover_amr_uniprot.py
│   │   ├── enrich_amr_graph.py
│   │   ├── amr_report.py
│   │   ├── amr_figures.py
│   │   └── amr_path_analysis.py
│   └── tests/                        # pytest suite (71 passing)
│
├── strains/
│   ├── TEMPLATE/                     # New-strain scaffold (config.json + subdirs)
│   ├── B36/                          # E. coli B36 (ST131) — primary validation
│   ├── MS_14384/ … MS_14387/         # Additional strains
│   └── (each contains config.json,
│         genome/, transcriptomic/,
│         proteomics/, eggnog_annotations.tsv)
│
├── metabolic/                        # MetaboLights MTBLS2015 MAF files
├── outputs/<strain>/                 # Generated artifacts (per strain)
│   ├── aligned_multiomics.csv        # Central-dogma aligned table
│   ├── genome_genes.csv / transcriptomics_expression.csv /
│   │   proteomics_abundance.csv / metabolomics_abundance.csv
│   ├── uniprot_annotations.csv / string_ppi.csv / gene_features.csv
│   ├── heterodata.pt                 # PyG HeteroData
│   ├── graph/                        # edge CSVs, metadata, hyperedges, incidence.npy
│   ├── figures/                      # visualization PNGs
│   └── cache/                        # annotation API cache
├── reports/                          # AMR analysis outputs (cross-strain)
└── docs/                             # Architecture documentation
```

**Conventions**
- Strains are self-contained: `strains/<name>/config.json` resolves all file paths; data files are auto-discovered by convention (`*.gff`, `*counts*`, `*.xlsx`, `*.faa`, `*maf.tsv`).
- Two config layouts are supported: per-strain folder (recommended) and a legacy centralized `strains_config.json`.
- All pipeline outputs are deterministic file artifacts under `outputs/<strain>/`.

---

## 5. Data Sources

| Layer | Source | Identifier / Format |
|---|---|---|
| Genome | NCBI Assembly (GFF3) | `GCF_900622635.1` (B36), `GCF_900622695.1` (MS_14384), `GCF_900622655.1` (MS_14385), `GCF_900622665.1` (MS_14386), `GCF_900622685.1` (MS_14387) |
| Transcriptome | GEO RNA-seq counts | GSE152966 (B36), GSE152967 (MS_14384); 6 RPMI + 6 sera replicates |
| Proteome | SWATH-MS | Excel `Area - proteins` sheet (bacterial-filtered) |
| Metabolome | MetaboLights MTBLS2015 | GC-MS / LC-MS MAF (`.maf.tsv`), sample-condition map |
| PPI | STRING API | *E. coli* CFT073 (tax 199310) primary, K-12 (tax 511145) fallback; `combined_score` 0–1 |
| Pathways | KEGG REST | organism code `eco`; gene → KO → pathway chain; pathway names |
| COG / orthology | eggNOG-mapper 2.1.13 | `strains/<name>/eggnog_annotations.tsv` (emapper output converted by `convert_emapper_to_tsv.py`) |
| Functional annotation | UniProtKB (REST) | RefSeq WP_ → UniProt mapping; GO, EC, names, cross-refs |
| AMR knowledge | Curated manifest (paper + GFF-verified) | `amr/amr_manifest.json` |

All online retrievals are **cached locally** (`outputs/<strain>/cache/`) and can be bypassed with `--skip-annotation` for offline runs.

---

## 6. Pipeline Stage 1 — Per-Strain Configuration & Discovery

**Entry point:** `python main.py --strain B36` (or `--all-strains`, `--list-strains`, `--discover`).

**Responsibilities of `load_config(strain)` (`multiomics_graph/main.py`):**
1. Resolve the strain definition from `strains/<name>/config.json` (primary) or the legacy centralized JSON (fallback).
2. Flatten nested sections (`transcriptomics`, `proteomics.swath`, `metabolomics`, `paths`) into a single pipeline config dict.
3. Resolve every data path relative to the strain folder, applying convention-based discovery when no explicit path is given:
   - `genome_gff` ← `genome/*.gff`
   - `protein_faa` ← `genome/*.faa`; `cds_fna` ← `genome/cds_from_genomic.fna`; `genome_fna` ← `genome/*.fna`
   - `rna_counts_file` ← `transcriptomic/*counts*`
   - `prot_file` ← `proteomics/*.xlsx`
   - `metabolomics_maf_files` ← `metabolic/*maf.tsv` or explicit list
4. Carry sample-ID wiring: RNA `rpmi_samples` / `sera_samples`, proteomics `rpmi_ids` / `sera_ids` plus column templates (with `{sid}` / `{strain}` placeholders), metabolomics sample lists.
5. Carry run parameters: KEGG code (`eco`), STRING primary/fallback tax IDs, `regulation_std_multiplier`, `correlation_threshold`.

**CLI flags**
- `--strain NAME`, `--all-strains`, `--output DIR`, `--correlation F`
- `--skip-annotation` (offline mode; fabricates minimal COG/pathway placeholders so the graph structure still builds)
- `--no-cache` (wipe annotation caches and recompute), `--list-strains`, `--discover`

---

## 7. Pipeline Stage 2 — Omics Preprocessing

Each layer is processed independently into a normalized, condition-aware table. All layers share the same two-condition experimental design (RPMI control vs. pooled sera treatment).

### 7.1 Genomics — `preprocessing/genomics.py::GenomicsProcessor`
- `parse_gff(path)`: parse the GFF3 annotation into a gene table with columns `GeneID` (locus tag), `GeneName`, `ProteinID` (RefSeq WP_ accessions, versioned), `start`, `end`, `strand`, `CDS_length`, `replicon`. Build bidirectional `GeneID ↔ ProteinID` maps.
- `compute_genomic_edges(max_gap=5000)`: undirected `gene–gene` **genomic proximity** edges between adjacent CDS on the same replicon whose intergenic gap is below the threshold; records distance and a weight.

**Outputs:** `genome_genes.csv`, `genome_genomic_edges.csv`.

### 7.2 Transcriptomics — `transcriptomics.py::TranscriptomicsProcessor`
- `load_log2_cpm(filepath)`: ingest gene-level counts, index by `GeneID`, expose the raw replicate matrix.
- `compute_condition_means(rpmi_samples, sera_samples)`: derive mean `RNA_RPMI`, `RNA_Sera`, `RNA_logFC`, and a regulation label (`Up`/`Down`/`Stable`). The threshold is either an absolute logFC cutoff or `std_multiplier × σ(logFC)` when `regulation_std_multiplier > 0`.
- `get_regulated_genes(direction)`, `summary()`.

**Output:** `transcriptomics_expression.csv`.

### 7.3 Proteomics — `proteomics.py::ProteomicsProcessor`
Loads the SWATH-MS Excel and aggregates by ProteinID:
- **SWATH-MS Excel** — sheet `Area - proteins`, protein column + template-generated RPMI/sera columns → log2(mean+1) intensity per condition.
- `filter_bacterial_proteins(valid_prefix="WP_")`, `compute_condition_means(...)`.

The regulation label uses the same std/absolute-threshold logic as RNA.

**Output:** `proteomics_abundance.csv` (aggregate `Protein_RPMI`, `Protein_Sera`, `Protein_logFC`, `Protein_Regulation`).

### 7.4 Metabolomics — `metabolomics.py::MetabolomicsProcessor`
- `parse_sample_file()`: map sample columns → condition from the sample sheet.
- `_load_maf` / `load_raw()`: parse MetaboLights MAF files per instrument; `_normalise_instruments` performs per-platform scale normalization so GC-MS and LC-MS are comparable.
- `_metabolite_id(row)`: canonicalize each metabolite to a stable ID (prefer `CHEBI:...`, else a name-derived key).
- `load(...)`: merge MAFs, assign samples to RPMI/sera, compute `Metabolite_RPMI`, `Metabolite_Sera`, `Metabolite_logFC`, and retain `MassToCharge` + `RetentionTime` as chemical features.
- `summary()`.

**Output:** `metabolomics_abundance.csv`.

---

## 8. Pipeline Stage 3 — Cross-Omics Alignment & Integration

**Method (`StrainPipeline.run_integration`)**
1. Build the `ProteinID → GeneID` bridge from the genome gene table.
2. Deduplicate genome rows per GeneID and index by GeneID.
3. Map proteomics onto GeneID via the bridge; drop unmapped proteins.
4. **3-layer merge (default):** inner join of RNA features (`RNA_RPMI`, `RNA_Sera`, `RNA_logFC`) and protein features (`Protein_*`, `Protein_Regulation`) on GeneID.
5. **2-layer fallback (no RNA):** align proteomics against the genome only; RNA columns are `NaN` (a true *absence*, not a fabricated zero) and flagged via `has_rna=False` downstream.
6. Attach genomic features (start, end, strand, CDS_length, replicon) to each aligned row.
7. Record `gene_protein_map` (versioned WP_) for graph construction.

**Output:** `aligned_multiomics.csv` — the central-dogma aligned table (gene | RNA | protein) that drives node features, `encodes` edges, and multi-omics-triple hyperedges. Metabolomics is *joined via annotation*, not row-aligned (metabolites have no GeneID).

**Design note:** the aligned genes define the "multi-omics-supported" set (drawn solid in AMR figures); the *gene node universe* is larger — genome-annotated ∪ aligned genes — so genome-only AMR determinants are never dropped (§12, §17).

---

## 9. Pipeline Stage 4 — External Knowledge Enrichment

Runs only when annotations are enabled (not `--skip-annotation`). Each retriever caches to `outputs/<strain>/cache/`.

### 9.1 UniProt — `annotation/uniprot.py::UniProtAnnotator`
- Batched RefSeq WP_ → UniProt mapping (`map_refseq_to_uniprot`) via UniProt REST; covers the **500-hits-per-batch** truncation for regular runs (recovered separately for AMR determinants, §13).
- `fetch_annotations`: per-protein functional fields — UniProt ID, gene/protein names, description, EC numbers, GO terms (BP/MF/CC) with a `go_terms` registry (name + aspect), subcellular location, cross-references (incl. KEGG), length, review status.
- `_validate` / `_compute_stats`: mapping coverage and annotation quality statistics.
- `_CacheManager`: local caching for offline reruns.

**Output:** `uniprot_annotations.csv`.

### 9.2 STRING PPI — `annotation/string_api.py::STRINGClient`
- `fetch_interactions(protein_ids, species_id)` over the STRING API with a `combined_score` confidence threshold.
- **Dual-species merge + fallback:** query the primary species (CFT073, tax 199310) *and* the fallback (K-12, tax 511145); concatenate and deduplicate on the unordered protein pair, keeping the max score.
- Map STRING (UniProt) IDs back to *versioned* RefSeq WP_ so graph node identity is consistent.
- `get_string_ids`, `summary()`.

**Output:** `string_ppi.csv`.

### 9.3 KEGG pathways — `annotation/kegg.py::KEGGAnnotator`
- Chain: **gene → KEGG gene ID → KO → pathway** via `get_ko_for_genes`, `map_ko_to_pathways`, `build_pathway_membership`, `get_pathway_name`.
- Gene→KEGG-Gene bridging for *E. coli* relies on UniProt cross-references (`eco:` KEGG IDs). Where the bridge is absent (acquired/plasmid AMR genes), coverage is recovered from **eggNOG KEGG_Pathway columns** (§9.5, §13).

**Output:** pathway-membership dicts (id → list of locus tags) + name map.

### 9.4 KEGG compounds — `annotation/kegg_compound.py::KEGGCompoundResolver`
- For each metabolite: **ChEBI → KEGG compound → enzymes (EC numbers) + pathways** via `resolve(chebi_ids)`; a name-based fallback `find_by_name(names)` covers ChEBI-less or unresolved metabolites.
- Results are re-keyed onto canonical `MetaboliteID`s for graph edges.

**Output:** `metabolite_to_pathway` and `metabolite_to_enzymes` maps (also exported as CSVs in §15).

### 9.5 eggNOG / COG — `annotation/eggnog.py::EggNOGAnnotator`
- `load_precomputed(filepath)`: read `strains/<name>/eggnog_annotations.tsv` (produced by eggNOG-mapper 2.1.13 `emapper.py --dmnd_iterate no`, converted by `scripts/convert_emapper_to_tsv.py`).
- Keyed by **unversioned** WP accessions; `as_locus_tag_map(gene_protein_map)` re-keys onto locus-tag gene nodes (the graph uses *versioned* WP_), so COG edges/hyperedges attach to real gene nodes.
- `kegg_pathway_membership(gene_protein_map)`: exposes the emapper `KEGG_ko` / `KEGG_Pathway` columns, closed gap for acquired/plasmid-borne determinants invisible to the `eco`-bridged UniProt route (§13).
- `annotate_by_locus_tag(gene_ids)`, `get_genes_by_cog(category)`, `summary()`.

### 9.6 Identifier normalization — `annotation/identifiers.py`
- `normalize_wp_id` and helpers to strip version suffixes for cross-database lookups and restore the *first versioned* Id when mapping back to graph nodes.

---

## 10. Pipeline Stage 5 — Heterogeneous Graph Construction

**Builder:** `graph/heterogeneous_graph.py::HeterogeneousGraphBuilder` → a PyTorch Geometric `HeteroData` object (dict fallback when PyG is absent). Node IDs map to consecutive indices via per-type maps.

### Node types and features
| Node type | Features | Notes |
|---|---|---|
| `gene` | `RNA_RPMI`, `RNA_Sera`, `RNA_logFC`, `start/1e6`, `CDS_length/1e3`, strand one-hot (`+`→`[1,0]`, `−`→`[0,1]`) | 7 features; gene universe = genome-annotated ∪ aligned (genome-only AMR markers included) |
| `protein` | `Protein_RPMI`, `Protein_Sera`, `Protein_logFC` | 3 features |
| `metabolite` | `Metabolite_RPMI`, `Metabolite_Sera`, `Metabolite_logFC`, `m/z/1e3`, `RT/1e6` | 5 features |
| `annotation` | name registry | shared type hosting KEGG pathways, COG categories, GO terms, PPI clusters, AMR mechanism nodes |

### Edge types (source --relation--> target)
| Relation | Semantics | Source |
|---|---|---|
| `(gene, genomic_proximity, gene)` | adjacent CDS on a replicon | GFF coordinates; edge attr = weight |
| `(gene, transcriptional_correlation, gene)` | |Pearson r| ≥ threshold across RNA replicates | RNA matrix |
| `(gene, encodes, protein)` | central dogma | GFF `protein_id` bridge |
| `(protein, abundance_correlation, protein)` | correlation of protein abundances | proteomics |
| `(protein, ppi, protein)` | STRING interaction (score attr) | STRING API |
| `(gene, in_pathway, annotation)` | KEGG pathway membership | KEGG (+ eggNOG recovery) |
| `(gene, cog_category, annotation)` | COG functional category | eggNOG TSV |
| `(protein, annotated_by, annotation)` | GO term | UniProt GO_BP/MF/CC (+ AMR recovery) |
| `(gene, associated_with, annotation)` | curated AMR mechanism class | `amr_manifest.json` |
| `(protein, member_of_cluster, annotation)` | Louvain PPI community | networkx on STRING |
| `(metabolite, in_pathway, annotation)` | metabolite → KEGG pathway | KEGG compound resolver |
| `(metabolite, metabolised_by, protein)` | metabolite → enzyme (EC) → protein | EC bridge |

Optional self-loops (`add_self_loops`) help GNN training. Edges of the same relation carry edge attributes where meaningful (correlation, PPI score, distance).

### PPI clusters
Weighted Louvain community detection (`resolution=1.0`, `seed=42`, min cluster ≥ 3) partitions the STRING graph into functional modules; each community becomes a `member_of_cluster` edge and a `ppi_cluster` hyperedge.

### Persistence (`builder.save`)
- `heterodata.pt` — the serialized `HeteroData`.
- `edges_<src>_<rel>_<dst>.csv` — per-relation edge lists (node IDs decoded from indices; same-type edges keep distinct `_target` columns).
- `heterodata_metadata.txt` — node/edge type summary.

---

## 11. Pipeline Stage 6 — Hypergraph Construction

**Builder:** `graph/hypergraph.py::BiologicalHypergraphBuilder`. A hyperedge groups ≥2 nodes that belong to a shared biological module. Representation: binary incidence matrix `H` (`n_nodes × n_hyperedges`, `H[i,j]=1` if node `i` ∈ hyperedge `j`).

### Hyperedge types
| Type | Contents |
|---|---|
| `multi_omics_triple` | per aligned locus: gene + `<gene>_RNA` + protein; features `RNA_logFC`, `Protein_logFC`, `has_rna`, `agreement_score` (product of logFCs; `None` when RNA absent — no pseudo-value) |
| `kegg_pathway` | all genes in a KEGG pathway |
| `cog_category` | genes (and their proteins) sharing a COG category; named `COG_<X>: <description>` |
| `go_term` | proteins sharing a GO term (from UniProt, incl. recovered AMR GO) |
| `ppi_cluster` | Louvain communities ≥ `min_cluster_size` from the STRING network |
| `amr_mechanism` | per curated mechanism class: its locus-tag genes + `AMR:<class>` mechanism node (keeps genome-only determinants represented) |
| `metabolic_pathway` | co-joins metabolites with pathway genes on shared KEGG metabolic pathways |

### Operators
- `build_incidence_matrix()` → `H`.
- `build_hypergraph_laplacian()` → `G = D_v^{-1/2} H D_e^{-1} H^T D_v^{-1/2}`, the propagation operator for hypergraph convolution (HGNN).

### Persistence
`hyperedges.csv`, `incidence_matrix.npy`, `hypergraph_nodes.csv` (with inferred node types: genome/transcriptome/proteome/amr_mechanism), `hyperedge_type_counts.csv`.

---

## 12. Pipeline Stage 7 — AMR Knowledge Layer

**Purpose.** Attach known resistance determinants as a biological *annotation layer* — not a supervised target. Resistance genes present in the genome but absent from RNA/protein measurements are real markers and remain represented (status `genome_only`, no invented values).

### Manifest — `amr/amr_manifest.json`
Curated per strain (B36, MS_14384, MS_14386, MS_14387; source = study paper + each strain's NCBI GFF). Each marker: `name` (e.g. `tet(A)`, `blaCTX-M-15`, `aac(6')-Ib-cr5`), `amr_class`, `locus_tag`, `protein_id`, `source`.

### Loader — `amr/amr.py::AMRManifest`
- `locus_to_amr(strain)` — locus → marker record.
- `class_to_loci(strain)` — mechanism class → locus tags.
- `amr_loci(strain)`, `markers`, `summary`, `class_name`.

Mechanism classes (`AMR_CLASS_NAMES`): β-lactam, aminoglycoside, sulfonamide, trimethoprim, tetracycline, macrolide, fluoroquinolone, phenicol, fosfomycin.

### Wiring into the graph
1. **Graph:** `gene_to_amr` inverted from `class_to_loci` → `(gene, associated_with, annotation)` edges + `AMR:<class>` annotation nodes (displayed as "&lt;Name&gt; Resistance").
2. **Hypergraph:** `add_amr_hyperedges` builds `amr_mechanism` hyperedges per class including the `AMR:<class>` node.
3. **Features/export:** `gene_features.csv` gains `is_amr_gene` and `amr_class` (semicolon-joined) columns.
4. **Figures:** `amr_highlight.png` — genome-only markers **outlined/hatched**, multi-omics-supported markers **solid**.

### Why the node universe must include genome-only markers
A determinant like `tet(A)` has no RNA/protein quantitation, yet it still carries `encodes`, COG, KEGG, GO, and genomic-proximity edges that connect it to the broader molecular state. Excluding it would sever those paths and bias the representation against exactly the markers AMR analysis cares about.

---

## 13. Pipeline Stage 8 — Analysis Scripts

### 13.1 eggNOG pipeline maintenance
- `scripts/run_eggnog_mapper.py` — regenerate per-strain eggNOG annotations (`emapper.py --dmnd_iterate no`, diamond).
- `scripts/convert_emapper_to_tsv.py` — convert emapper output to `strains/<name>/eggnog_annotations.tsv`, **preserving `KEGG_ko` / `KEGG_Pathway` columns** (previously dropped) so the graph gains a KEGG layer for acquired/plasmid-borne AMR genes.

### 13.2 AMR graph enrichment (offline, idempotent)
- `scripts/recover_amr_uniprot.py` — per-ID UniProt GO/function lookup for AMR determinants, bypassing the batch-mapping 500-hits truncation → `reports/amr_uniprot_recovery.csv`.
- `scripts/enrich_amr_graph.py` — merge eggNOG KEGG/COG, recovered GO, `encodes`, and genomic-proximity edges for the determinants into the *saved* per-strain graphs; updates kegg/cog/go hyperedges so BFS path tracing sees them.

### 13.3 AMR reporting
- `scripts/amr_report.py` — per-strain + cross-strain AMR presence/synergy reports from saved outputs.
- `scripts/amr_figures.py` — regenerate AMR highlight figures and verify **100%** of the paper's loci are drawn (B36 8/8, MS_14384 1/1, MS_14386 14/14, MS_14387 4/4).

### 13.4 8-layer determinant path tracing — `scripts/amr_path_analysis.py`
An **exploratory / mechanism-characterization** framework (no model, no prediction, no causal claims). For each curated determinant it describes, over **8 evidence layers** (`genome, transcriptome, proteome, metabolome, kegg, cog, go, ppi_amr`):
- which layers contain information about the marker,
- the **evidence taxonomy** — `directly_observed` (measured abundance/locus), `graph_derived` (proximity/PPI/cluster/hyperedge/bounded path), `annotation_mediated` (KEGG/COG/GO/AMR), `inferred` (metabolite via pathway/enzyme), `unavailable` (layer absent);
- bounded **≤4-hop** paths through the entity graph (BFS, max 400 rows per marker),
- per-strain differences and conserved context,
- **candidate AMR-associated elements** (non-AMR entities recurring in determinant context — hypotheses, not confirmed determinants), resolved to names (CHEBI → metabolite name; locus/protein → gene symbol/product) and de-biased: entities reaching ≥ `HUB_FRACTION` (0.75) of the maximum determinant count are flagged `hub_candidate` (ubiquitous pathway-shared metabolites) and ranked after specific candidates.

Missing data are never fabricated: a missing measurement is `not_observed`, a missing layer `unavailable`. RNA falls back to `transcriptomics_expression.csv` so genome-only markers report a *real* transcriptome layer rather than a fabricated zero.

**Outputs** (`reports/`): `amr_path_trace_<strain>.csv`, `amr_path_evidence_<strain>.csv`, `amr_path_hyperedge_<strain>.csv`, `amr_path_cross_strain.csv`, `amr_path_candidates.csv`, `amr_path_report.md`.

---

## 14. Pipeline Stage 9 — Visualization

**Module:** `visualization/network_plots.py::NetworkVisualizer` (writes to `outputs/<strain>/figures/`).

- **`plot_heterogeneous_graph`** — full network render from saved edge CSVs (typed node coloring), plus a sampled (≤500-node) variant when the graph is large.
- **`plot_amr_highlight`** — AMR determinants overlaid; mechanism nodes (`AMR:<class>`) linked to their loci; genome-only markers hatched vs. multi-omics-supported solid; verifies every manifest locus is drawn.
- **`plot_regulation_concordance`** — RNA vs. protein regulation concordance per aligned gene.
- **`plot_metabolite_regulation`** — metabolite up/down counts.
- **`plot_cog_distribution`** — COG category histogram (COG_CATEGORIES labels).
- **`plot_hyperedge_type_distribution`**, **`plot_incidence_matrix`** — hypergraph diagnostics.
- **`plot_summary_dashboard`** — stage-level summary panel.

---

## 15. Pipeline Stage 10 — Export & Downstream Interface

All artifacts under `outputs/<strain>/`:
- **Tables:** `aligned_multiomics.csv`, `genome_genes.csv`, `transcriptomics_expression.csv`, `proteomics_abundance.csv`, `metabolomics_abundance.csv`, `uniprot_annotations.csv`, `string_ppi.csv`, `gene_features.csv` (with `is_amr_gene`/`amr_class`), `metabolite_pathway_membership.csv`, `metabolite_enzymes.csv`.
- **Graph:** `heterodata.pt` (PyG `HeteroData`); per-relation `edges_*.csv`; `heterodata_metadata.txt`.
- **Hypergraph:** `hypergraph_incidence.npy`, `hypergraph_incidence_edges.csv` (sparse), `hypergraph_laplacian.pt`.
- **Consumed by:** RGCN / HAN (typed nodes/edges), HGNN (incidence + Laplacian), node classification, link prediction, and the AMR path-tracing scripts (§13.4).

---

## 16. Data Model Summary

| Component | Type | Cardinality (B36 current run) |
|---|---|---|
| Gene nodes | `gene` | ~2,283 aligned (+ genome-only markers) |
| Protein nodes | `protein` | ~2,283 |
| Metabolite nodes | `metabolite` | ~219 |
| Annotation nodes | `annotation` | KEGG pathways + COG + GO + clusters + AMR classes (hundreds) |
| Edge types | heterogeneous | 12 (genomic proximity, transcriptional correlation, encodes, abundance correlation, ppi, in_pathway, cog_category, annotated_by/GO, associated_with/AMR, member_of_cluster, metabolite in_pathway, metabolised_by) |
| Hyperedge types | hypergraph | 7 (multi_omics_triple, kegg_pathway, cog_category, go_term, ppi_cluster, amr_mechanism, metabolic_pathway) — ~3,500 hyperedges total |
| Incidence / Laplacian | `np.ndarray` / `torch.Tensor` | `n_nodes × n_hyperedges`; `n_nodes × n_nodes` |

Representative current graphs: B36/MS_14384/MS_14386/MS_14387 produce 2,246–2,348 aligned genes (MS_14385 runs the 2-layer genome+proteome alignment, no RNA). Exact counts live in `outputs/<strain>/graph/heterodata_metadata.txt`.

---

## 17. Design Decisions & Rationale

1. **Central dogma as an edge, not a join.** The `gene --encodes--> protein` edge preserves *who produces whom* as first-class structure; the alignment table only defines shared node features. This is what makes the fused graph relational rather than a concatenated feature matrix.
2. **Heterogeneous graph first, hypergraph as complement.** Pairwise edges support message passing; hyperedges capture module-level structure (pathways, categories, clusters, mechanisms) that pairwise edges cannot represent exactly. Both are exported.
3. **Genome-only AMR markers are never removed.** A determinant present in the genome but without RNA/protein quantitation is a real marker. It stays in the gene-node universe, gets `encodes`/COG/KEGG/GO/proximity edges, and is drawn hatched with status `genome_only`.
4. **Versioned WP_ as canonical protein identity.** eggNOG is keyed by unversioned accessions; UniProt/KEGG bridges and graph nodes use versioned IDs. A normalizer converts at the boundaries to avoid orphan nodes.
5. **Dual-bridge KEGG coverage.** Primary route: UniProt cross-reference → `eco:` KEGG gene → KO → pathway. Recovery route for acquired/plasmid determinants invisible to that bridge: emapper `KEGG_ko`/`KEGG_Pathway` columns carried through `convert_emapper_to_tsv.py`.
6. **Species fallback for STRING.** Primary UPEC species query + K-12 fallback, merged at max-confidence per pair — maximizes PPI coverage while staying within the *E. coli* clade.
7. **Config-over-code extensibility.** Adding a strain (or eventually a species) is a config exercise: new `strains/<name>/config.json` + data files + optional eggNOG TSV; no pipeline edits.
8. **Cache-everything, skip-annotation mode.** All external retrievals cache to `outputs/<strain>/cache/`; `--skip-annotation` permits fully offline reproduction with minimal placeholder annotations so the graph still builds.
9. **Honest missing-data semantics throughout.** `NaN` RNA in 2-layer strains, `has_rna=False`, `agreement_score=None`, `not_observed`/`unavailable` in path tracing — absence is represented, never imputed as zero.

---

## 18. Robustness & Correctness Measures

- **Unit tests:** pytest suite (71 passing) covering UniProt mapping/parsing, eggNOG/COG resolution, AMR manifest wiring, AMR path tracing, and corrections (`tests/`).
- **Offline reproducibility:** scripts in §13.2 read only saved outputs; `--skip-annotation` / `--no-cache` give explicit online vs. offline semantics.
- **Idempotent enrichment:** `enrich_amr_graph.py` merges deterministically so reruns do not duplicate edges.
- **Deterministic clustering:** Louvain with fixed `resolution` and `seed` (42); correlation edges use a fixed threshold (`--correlation`, default 0.7).
- **Coverage verification:** `amr_figures.py` asserts 100% of the paper's AMR loci are drawn per strain.
- **No fabricated data:** genome-only markers report `genome_only`; path tracing distinguishes `directly_observed` from `graph_derived` / `annotation_mediated` / `inferred` evidence.

---

## 19. Implementation Status

**Implemented and validated** on five *E. coli* strains (B36, MS_14384, MS_14385, MS_14386, MS_14387):
- All four omics preprocessing parsers and the 3-layer (or 2-layer) central-dogma alignment.
- UniProt, STRING (dual-species), KEGG, eggNOG/COG enrichment with caching.
- Heterogeneous graph (PyG `HeteroData`) + hypergraph (incidence, Laplacian).
- AMR knowledge layer (manifest, mechanism hyperedges, highlight figures, per-strain/cross-strain reports, 8-layer path tracing with candidates).
- Visualization and export artifacts as documented.

**Deferred / not yet implemented:**
- GNN training on the constructed graphs (RGCN/HAN/HGNN) — inputs are ready.
- Supervised AMR phenotype prediction (needs labeled phenotypes + larger strain panel).
- Multi-species graphs (*K. pneumoniae*, *S. aureus*, *S. pyogenes*) and species-specific PPI/KEGG mappings.

---

## 20. Open Questions & Future Work

- **GNN layer design.** Which model per task: RGCN (relation-aware message passing) for node/link tasks, HAN (metapath attention) for cross-omics aggregation, HGNN (incidence + Laplacian) for module-level learning? Metapath design (e.g., `gene → encodes → protein → ppi → protein`) needs empirical comparison.
- **Phenotype integration.** Adding MIC/binary resistance labels would unlock supervised AMR tasks; requires an expanded multi-omics strain panel.
- **Orthology / cross-species.** Representing AMR determinants by ortholog groups (eggNOG OGs) would allow multi-species graph fusion and knowledge transfer.
- **Community AMR databases.** Extending the manifest with CARD / ResFinder / ARG-ANNOT for broader, standardized coverage.
- **Candidate validation.** The candidate AMR-associated elements from path tracing are hypotheses; design experiments or orthogonal databases to test them.
- **Scale & memory.** Sparse incidence/Laplacian and metapath sampling for large panels; GPU-amenable formats for the incidence matrices.

---

## 21. Glossary

| Term | Meaning |
|---|---|
| `HeteroData` | PyTorch Geometric container for typed nodes and typed edge relations |
| Hyperedge | Set of ≥2 nodes forming one biological module (pathway, category, cluster, mechanism, triple) |
| Incidence matrix `H` | `H[i,j]=1` iff node `i` belongs to hyperedge `j` |
| Hypergraph Laplacian | `D_v^{-1/2} H D_e^{-1} H^T D_v^{-1/2}` — HGNN propagation operator |
| COG category | Broad functional class from eggNOG (J,K,L,D,…; `S` = function unknown) |
| KO (KEGG orthology) | Ortholog-group identifier linking genes to pathways |
| WP_ accession | RefSeq protein identifier; `.1` suffix is the version (graph uses versioned) |
| Locus tag | Per-strain gene identifier (e.g., `EW036_RS26095`) used as `GeneID` |
| Genome-only marker | AMR determinant in the genome without RNA/protein quantitation; kept with status `genome_only` |
| Evidence taxonomy | `directly_observed` / `graph_derived` / `annotation_mediated` / `inferred` / `unavailable` for path-trace evidence |
| Louvain community | Dense PPI subnetwork detected by weighted Louvain clustering; becomes a cluster edge + hyperedge |
| Multi-omics triple | (gene, RNA, protein) hyperedge per aligned locus — the central-dogma module |

---

*End of proposed architecture document. Cross-references: `multiomics_graph/TECHNICAL_DOCUMENTATION.md` (implementation), `multiomics_graph/AMR_GRAPH_INTEGRATION.md` (AMR layer phases), `multiomics_graph/METABOLOMICS_INTEGRATION.md` (metabolite layer), `docs/ARCHITECTURE_OVERVIEW.md` (high-level rationale).*
