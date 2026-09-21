# Multi-Omics Graph Fusion Architecture

## A Graph-Based Framework for Integrating Bacterial Genomic, Transcriptomic, and Proteomic Data

---

# 1. Introduction

Biological systems are governed by interactions that span multiple molecular levels — from genomic organization through transcriptional regulation to protein abundance and post-translational modification. No single omics layer captures this complexity in isolation. Genomics reveals the genetic potential of an organism but not which genes are actively expressed under a given condition. Transcriptomics captures gene activity but correlates imperfectly with protein abundance due to regulatory processes such as mRNA degradation, translational efficiency, and post-translational modifications. Proteomics measures the functional effectors of the cell but lacks the genomic context necessary to understand regulatory mechanisms.

The central dogma of molecular biology (DNA → RNA → Protein) provides a conceptual framework for relating these layers, but translating this linear model into a computational representation that preserves the inherently graph-like structure of biological relationships remains a challenge. Biological entities do not exist in isolation: genes are organized along chromosomes with shared regulatory contexts, transcripts form co-expression modules, and proteins participate in physical interaction networks. These relationships are naturally represented as graphs, where nodes correspond to biological entities and edges capture their functional, structural, or regulatory connections.

This document describes a multi-omics graph fusion architecture designed to integrate genomic, transcriptomic, and proteomic data from bacterial organisms into a unified graph representation. The architecture is motivated by the observation that biological systems are better understood through the relationships between molecular entities than through isolated measurements of individual components. By representing each omics layer as a graph and fusing these representations into a single heterogeneous graph, the architecture preserves both intra-omics relationships (e.g., genomic proximity, co-expression, protein–protein interactions) and inter-omics relationships (e.g., which gene encodes which protein, and how that gene's expression changes under different conditions).

---

# 2. Project Objective

The primary objective of the architecture is to construct a unified biological representation from matched multi-omics bacterial data. Rather than treating genomic, transcriptomic, and proteomic measurements as independent feature sets, the framework models them as interconnected layers of a single biological system.

The architecture is designed to support downstream graph neural network (GNN) learning for tasks such as gene function prediction, pathway activity inference, and antimicrobial resistance (AMR) analysis. The graph structure preserves the relational information that is typically lost when omics data are concatenated into flat feature vectors.

The current implementation uses a publicly available multi-omics dataset from *Escherichia coli* strains B36, MS_14384, MS_14385, MS_14386, and MS_14387. These provide matched genomic, transcriptomic, and proteomic measurements from the same organisms; B36, MS_14384, MS_14386, and MS_14387 carry the paper's AMR genotype descriptions (MS_14385 has none). The datasets are used to validate the alignment/integration pipeline and to **characterize AMR mechanisms**: curated resistance determinants are attached to the graph as an annotation layer (see §3.5), enabling cross-strain AMR analysis. They are not supervised AMR-phenotype datasets — direct AMR *prediction* would require matched resistance-phenotype labels.

---

# 3. External Knowledge Bases

The architecture integrates four external knowledge bases that enrich the graph with functional context beyond what the experimental measurements alone can provide. Each knowledge base contributes a distinct dimension of biological information and is accessed programmatically through its public API or local data files.

## 3.1 STRING — Protein–Protein Interaction Network

**Purpose.** STRING (Search Tool for the Retrieval of Interacting Genes/Proteins) provides known and predicted protein–protein interactions (PPIs) based on evidence from curated databases, experimental assays, co-expression analysis, genomic context, and automated text mining.

**Role in the architecture.** PPI edges are the primary mechanism for modeling relationships between protein nodes. Without STRING, the proteomic layer would rely solely on abundance correlation, which captures only co-regulation signals and misses physical and functional interactions that are invisible to correlation analysis.

**Implementation.** The architecture queries the STRING API for all protein identifiers in the dataset. Each query returns interaction pairs with a combined confidence score in the 0–1 range. A minimum confidence threshold is applied to retain only high-confidence interactions. Because the primary strain (*E. coli* B36, genomically closest to CFT073, tax ID 199310) may not have full STRING coverage, a fallback strategy is used: proteins unmapped in the primary species are re-queried against *E. coli* K-12 MG1655 (tax ID 511145). Results from both queries are merged and deduplicated, retaining the highest confidence score for each protein pair.

A community detection step (Louvain clustering) is applied to the PPI network to identify densely connected protein modules. These modules become hyperedges in the hypergraph, representing functional complexes or pathways that are not explicitly captured by pairwise interactions alone.

## 3.2 UniProt — Functional Protein Annotation

**Purpose.** UniProt provides comprehensive functional annotation for proteins, including gene ontology (GO) terms, protein names, EC numbers, transmembrane domains, and subcellular localization.

**Role in the architecture.** UniProt bridges the gap between protein identifiers (RefSeq WP_ accessions from the genome annotation) and functional descriptors. The architecture maps each WP_ accession to its UniProt entry and retrieves available functional annotations.

**Implementation.** RefSeq protein identifiers are mapped to UniProt accessions using the UniProt ID mapping service. For each mapped entry, functional annotations are retrieved from the UniProtKB REST API. The retrieved data are cached locally to avoid redundant API calls across pipeline runs.

The UniProt annotations serve as node features for protein nodes. This means each protein in the graph carries not only its abundance measurements from the proteomics experiment but also its known functional context from the curated literature.

## 3.3 KEGG — Pathway Membership

**Purpose.** KEGG (Kyoto Encyclopedia of Genes and Genomes) is a curated database of biological pathways, including metabolism, genetic information processing, environmental information processing, and cellular processes.

**Role in the architecture.** Pathway membership connects individual genes to the higher-order biological processes in which they participate. These set-level relationships are essential for the hypergraph, where each pathway becomes a hyperedge linking all genes that participate in it.

**Implementation.** The architecture follows the KEGG gene-to-pathway mapping chain: gene identifiers → KEGG orthologs (KO numbers) → pathway identifiers. The gene-to-KO mapping is retrieved from the KEGG API for the appropriate organism code (*eco* for *E. coli*). KO-to-pathway mappings are then used to assign each gene to its set of associated pathways.

In the heterogeneous graph, each gene–pathway pair produces an edge from the gene node to an annotation node representing the pathway. In the hypergraph, all genes belonging to the same pathway are grouped into a single hyperedge.

## 3.4 eggNOG — COG Functional Categories

**Purpose.** eggNOG (evolutionary genealogy of genes: Non-supervised Orthologous Groups) provides functional annotation through orthologous group assignments and the associated COG (Clusters of Orthologous Groups) functional categories — broad functional classes such as metabolism, information storage and processing, and cellular processes.

**Role in the architecture.** COG categories provide a coarse-grained functional classification for every gene in the genome. Unlike KEGG pathways, which are specific to known molecular processes, COG categories cover the entire gene catalogue, including genes of unknown function that are assigned to a broad category like "function unknown" (COG category S).

**Implementation.** eggNOG-mapper 2.1.13 (`emapper.py --dmnd_iterate no`, diamond search) is run locally per strain over the translated proteome (`protein.faa`). The resulting `.emapper.annotations` file is converted to `strains/<name>/eggnog_annotations.tsv` (via `scripts/convert_emapper_to_tsv.py`), keyed by **unversioned** RefSeq WP accessions. Because the graph/`genome_genes.csv` uses **versioned** WP accessions (`WP_000000542.1`), the annotation layer normalizes identifiers (`annotation/identifiers.py::normalize_wp_id`) and re-keys COG annotations onto locus-tag gene nodes (`annotation/eggnog.py::as_locus_tag_map`). COG edges and `cog_category` hyperedges therefore attach to real gene nodes. Each COG category becomes an annotation node in the heterogeneous graph and a hyperedge grouping all genes sharing that classification.

## 3.5 AMR Knowledge Layer (curated)

**Purpose.** Known antimicrobial-resistance determinants from the underlying study are curated per strain (`amr/amr_manifest.json`) and attached to the graph as a biological annotation layer — not as a supervised prediction target.

**Role in the architecture.** Each curated mechanism class (β-lactam, aminoglycoside, sulfonamide, trimethoprim, tetracycline, macrolide, fluoroquinolone, phenicol, fosfomycin) becomes an `amr_mechanism` hyperedge containing its locus-tag gene nodes plus an `AMR:<class>` mechanism node. Gene nodes are always included even for **genome-only** determinants with no RNA/protein quantitation — a resistance gene present in the genome is a real marker and must remain represented.

**Implementation.** `amr/amr.py` loads the manifest; `graph/hypergraph.py::add_amr_hyperedges` builds the class hyperedges; the AMR highlight figure distinguishes genome-only markers (outlined/hatched) from multi-omics-supported markers (solid) without inventing RNA/protein values; `gene_features.csv` flags `is_amr_gene`/`amr_class`. Cross-strain presence reports and figures are produced by `scripts/amr_report.py`, and `scripts/amr_figures.py` regenerates the highlights while verifying 100% of the paper's AMR loci are drawn. `scripts/amr_path_analysis.py` traces each determinant across the eight evidence layers and ranks candidate AMR-associated elements, resolving CHEBI/gene/protein ids to names and down-ranking generic pathway-shared "hub" candidates. Full details: `multiomics_graph/AMR_GRAPH_INTEGRATION.md`.

---

# 4. Primary Experimental Dataset

The architecture processes three aligned omics layers from the same bacterial strain, extended to four with metabolomics. The pipeline has been applied to five *E. coli* strains from the study (B36, MS_14384, MS_14385, MS_14386, MS_14387); this section documents the primary validation strain B36 and notes the multi-strain scope.

## 4.1 Genomics

**Source.** The genome of *Escherichia coli* B36 was obtained from the NCBI Assembly database (accession GCF_900622635.1). The assembly includes a complete genome sequence with annotation in GFF3 format. Each strain under analysis has its own NCBI genome (GCF_900622655.1 for MS_14385, correctly assigned to that strain; GCF_900622665.1 for MS_14386, etc.).

**Role.** The genome annotation provides the reference framework for the entire architecture. It defines the set of genes and their identifiers (locus tags), maps each gene to its encoded protein (via protein_id cross-references in CDS entries), and records genomic organization — replicon assignment, chromosomal coordinates, and strand orientation.

**Information extracted.** Gene identifiers, gene symbols, RefSeq protein identifiers (WP_ accessions), CDS boundaries, replicon assignment, and strand orientation. This set of 2,283 protein-coding genes forms the genomic node layer, with genomic proximity between adjacent genes encoding structural relationships.

## 4.2 Transcriptomics

**Source.** RNA sequencing data for *E. coli* B36 were obtained from the Gene Expression Omnibus (GEO accession GSE152966). The dataset contains gene-level count data from six biological replicates under RPMI laboratory growth conditions and six replicates under pooled human serum exposure.

**Role.** Transcriptomics captures which genes are actively transcribed and how expression changes between the two experimental conditions. The comparison between RPMI (control) and pooled sera (treatment) reveals the bacterial transcriptional response to host serum — a condition relevant to understanding infection biology and stress adaptation.

**Information extracted.** Log2-transformed counts per million per replicate, mean expression per condition, log2 fold change between conditions, and a categorical regulation label (up-regulated, down-regulated, or stable). Regulation thresholds are configurable and can be set as either absolute logFC cutoffs or as a multiple of the standard deviation of the observed logFC distribution.

## 4.3 Proteomics

**Source.** Proteomic data for *E. coli* B36 were acquired using SWATH-MS mass spectrometry. Six replicates per condition match the transcriptomic experimental design.

**Role.** Proteomics measures the abundance of proteins — the functional products of gene expression. Protein abundance frequently diverges from transcript abundance due to post-transcriptional regulation, differential translation efficiency, and protein degradation. Proteomics is therefore a necessary complement to transcriptomics for capturing the complete molecular response of the organism.

**Information extracted.** Log2-transformed protein intensity per replicate, mean abundance per condition, log2 fold change between conditions, and protein regulation status. Protein-level identifiers (WP_ accessions) enable direct cross-referencing with the genome annotation.

## 4.4 Dataset Scope and Limitations

The three omics layers together provide a comprehensive picture of the bacterial response to serum exposure at the transcriptional and translational levels. The dataset is well-suited for validating the multi-omics integration pipeline because all layers originate from the same strain under matched conditions, cross-omics alignment is unambiguous through shared identifiers, and the biological system shows measurable differential regulation across conditions.

However, this is not an antimicrobial resistance experiment. The dataset studies bacterial adaptation to human serum, not resistance to antibiotics. The architecture therefore cannot be directly validated for AMR prediction using this data — an AMR-specific task would require a dataset with matched multi-omics measurements from resistant and susceptible strains, phenotypic resistance labels (MIC values or binary classification), and coverage of known resistance genes and pathways.

The value of the current implementation lies in demonstrating that the graph construction and fusion pipeline produces a coherent, biologically interpretable representation from real multi-omics data. The architecture is designed to be dataset-agnostic: replacing the input data with an AMR-specific multi-omics dataset requires only configuration changes, not architectural modifications.

---

# 5. Overall Architecture

The architecture processes raw multi-omics data through a sequence of stages, each transforming the data toward the final unified graph representation.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                          INPUT DATA SOURCES                                 │
│                                                                             │
│  ┌─────────────────┐  ┌──────────────────┐  ┌──────────────────────┐       │
│  │    GENOMICS      │  │  TRANSCRIPTOMICS  │  │     PROTEOMICS       │       │
│  │  NCBI GFF/GCF    │  │  RNA-seq counts   │  │      SWATH-MS      │       │
│  │  (GCF_900622635) │  │  (GSE152966)      │  │  (6 RPMI + 6 sera)  │       │
│  └────────┬─────────┘  └────────┬──────────┘  └──────────┬───────────┘       │
│           │                     │                        │                   │
└───────────┼─────────────────────┼────────────────────────┼───────────────────┘
            │                     │                        │
            ▼                     ▼                        ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      PROCESSING & ALIGNMENT                                 │
│                                                                             │
│  ┌──────────────┐  ┌──────────────────┐  ┌────────────────────┐            │
│  │ GFF Parsing  │  │ log2 CPM +       │  │ log2 Intensity +   │            │
│  │ Gene/Protein │  │ Condition Means  │  │ Condition Means    │            │
│  │ ID Mapping   │  │ logFC +          │  │ logFC +            │            │
│  │ Proximity    │  │ Regulation Label │  │ Regulation Label   │            │
│  │ Edges        │  │ Co-expression    │  │ Abundance Corr.    │            │
│  └──────┬───────┘  └────────┬─────────┘  └────────┬──────────┘            │
│         │                   │                      │                        │
│         └───────────────────┼──────────────────────┘                        │
│                             │                                               │
│                             ▼                                               │
│              ┌──────────────────────────────┐                               │
│              │   CROSS-OMICS ALIGNMENT       │                               │
│              │  (Inner join on GeneID /      │                               │
│              │   ProteinID → GeneID bridge)  │                               │
│              └──────────────┬───────────────┘                               │
└─────────────────────────────┼───────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    EXTERNAL KNOWLEDGE ENRICHMENT                             │
│                                                                             │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐                   │
│  │  STRING  │  │  UniProt │  │   KEGG   │  │  eggNOG  │                   │
│  │ PPI      │  │ Function │  │ Pathways │  │ COG Cat. │                   │
│  │ (CFT073  │  │ Annot.   │  │ Gene→KO  │  │ Ortholog │                   │
│  │  + K-12  │  │ RefSeq→  │  │ →Pathway │  │ Groups   │                   │
│  │ fallback)│  │ UniProt  │  │ Mapping  │  │          │                   │
│  └────┬─────┘  └────┬─────┘  └────┬─────┘  └────┬─────┘                   │
│       │             │             │             │                           │
└───────┼─────────────┼─────────────┼─────────────┼───────────────────────────┘
        │             │             │             │
        ▼             ▼             ▼             ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                     HETEROGENEOUS GRAPH CONSTRUCTION                         │
│                                                                             │
│  ┌──────────────────────────────────────────────────────────────────┐      │
│  │                    NODE TYPES                                      │      │
│  │  ┌──────────┐  ┌──────────┐  ┌──────────────────────────┐        │      │
│  │  │   GENE   │  │ PROTEIN  │  │      ANNOTATION          │        │      │
│  │  │ (2,283)  │  │ (2,283)  │  │  (KEGG pathways + COG)   │        │      │
│  │  │ 7 feats  │  │ 3 feats  │  │       (327 nodes)        │        │      │
│  │  └────┬─────┘  └────┬─────┘  └───────────┬──────────────┘        │      │
│  └───────┼─────────────┼─────────────────────┼──────────────────────┘      │
│          │             │                     │                              │
│  ┌───────┴─────────────┴─────────────────────┴──────────────────────┐      │
│  │                    EDGE TYPES                                      │      │
│  │                                                                     │      │
│  │  ┌────────────────────┐  ┌──────────────┐  ┌──────────────────┐   │      │
│  │  │ Genomic Proximity  │  │   Encodes    │  │       PPI        │   │      │
│  │  │   gene ─── gene    │  │ gene ──►protein│  │ protein ── protein│   │      │
│  │  │    (1,165 edges)   │  │  (2,283 edges)│  │  (6,949 edges)   │   │      │
│  │  └────────────────────┘  └──────────────┘  └──────────────────┘   │      │
│  │                                                                     │      │
│  │  ┌────────────────────┐  ┌──────────────────────┐                  │      │
│  │  │   In Pathway       │  │    COG Category      │                  │      │
│  │  │ gene ───annotation │  │  gene ─── annotation │                  │      │
│  │  │  (2,100 edges)     │  │   (2,283 edges)      │                  │      │
│  │  └────────────────────┘  └──────────────────────┘                  │      │
│  └────────────────────────────────────────────────────────────────────┘      │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      HYPERGRAPH CONSTRUCTION                                 │
│                                                                             │
│   ┌─────────────────────┐  ┌──────────────────┐  ┌──────────────────┐      │
│   │ Multi-Omics Triples │  │ KEGG Pathway     │  │ COG Category     │      │
│   │ (gene + RNA + Prot) │  │ Hyperedges       │  │ Hyperedges       │      │
│   │   (2,283 hyperedges)│  │  (206 hyperedges) │  │  (1 hyperedge)   │      │
│   └─────────────────────┘  └──────────────────┘  └──────────────────┘      │
│                                                                             │
│   ┌──────────────────────────────┐                                          │
│   │   PPI Community Hyperedges   │                                          │
│   │   (Louvain clusters)         │                                          │
│   │   (10 hyperedges)            │                                          │
│   └──────────────────────────────┘                                          │
│                                                                             │
│   Total: 2,494 hyperedges | Incidence: 6,849 × 2,494 | Laplacian: 6,849²   │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                     UNIFIED BIOLOGICAL REPRESENTATION                        │
│                                                                             │
│   ┌─────────────────────────────────────┐  ┌──────────────────────────┐    │
│   │   Heterogeneous Graph (HeteroData)  │  │   Hypergraph             │    │
│   │   • Typed node types (3)            │  │   • Incidence matrix     │    │
│   │   • Typed edge types (5)            │  │   • Normalized Laplacian │    │
│   │   • Node feature tensors            │  │   • Hyperedge attributes │    │
│   │   • Edge index tensors              │  │                          │    │
│   └─────────────────────────────────────┘  └──────────────────────────┘    │
│                                                                             │
│   ┌─────────────────────────────────────────────────────────────────────┐  │
│   │                    EXPORT                                           │  │
│   │   CSV tables │ NPY matrices │ PT tensors │ Visualizations          │  │
│   └─────────────────────────────────────────────────────────────────────┘  │
│                                                                             │
│   ┌─────────────────────────────────────────────────────────────────────┐  │
│   │              READY FOR DOWNSTREAM GNN TRAINING                      │  │
│   │   • Node classification  • Link prediction  • Graph-level tasks    │  │
│   └─────────────────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────────────┘
```

**Stage 1 — Data Processing and Alignment.** Raw data from each omics layer are normalized independently. The genome annotation is parsed to establish the gene catalogue and gene–protein mappings. Transcriptomic counts are converted to log2 CPM, and mean expression is computed per condition. Proteomic intensities are log2-transformed and averaged per condition. The three data sources are then aligned experimentally: the transcriptome is indexed by gene identifier, the proteome by protein identifier, and the genome provides the bridge between them. Only entities with measurements in all three layers are retained in the final aligned table.

**Stage 2 — Individual Omics Representations.** Each omics layer is processed to extract layer-specific relationships. Genomic proximity between genes is computed based on chromosomal coordinates. Transcriptomic co-expression is inferred from Pearson correlation across replicates. Proteomic relationships combine abundance correlation with external protein–protein interaction data. These relationships are not yet embedded in a shared graph — they remain layer-specific at this stage.

**Stage 3 — Multi-Omics Graph Construction.** The aligned entities and their relationships are assembled into a heterogeneous graph with typed nodes and typed, directed edges. Gene nodes carry genomic features (position, strand, length), protein nodes carry expression-derived features (abundance, logFC, regulation status), and annotation nodes represent external knowledge (pathways, functional categories). Edge types include genomic proximity (gene–gene), co-expression (gene–gene), protein–protein interaction (protein–protein), encoding relationships (gene–protein), and functional associations (gene–annotation, protein–annotation).

**Stage 4 — Relationship Fusion.** The heterogeneous graph is complemented by a hypergraph that captures higher-order relationships — sets of nodes that together form a biological module. Each multi-omics triple (gene, its transcript, its protein) forms a hyperedge, as does each KEGG pathway, each COG functional category, and each PPI cluster. The hypergraph adds a grouping dimension that pairwise edges alone cannot represent.

**Stage 5 — Unified Biological Representation.** The final output consists of a PyTorch Geometric HeteroData object (the heterogeneous graph) and a hypergraph incidence matrix. These representations are ready for downstream GNN training, node classification, link prediction, or graph-level tasks.

---

# 6. Representation of Individual Omics Layers

Each omics layer is modeled as a distinct node type in the heterogeneous graph. The representation of each layer captures both the individual entities and their relationships.

## 6.1 Genomic Layer

The genomic layer represents genes as nodes. Each gene node is characterized by:
- Its genomic coordinates (replicon, start position, end position, strand)
- Protein-coding status and CDS length
- Its relationship to neighboring genes along the chromosome

The genomic layer does not model sequence content directly; rather, it captures genomic organization as a source of structural relationships. Genes that are adjacent on the chromosome often share regulatory elements, are co-transcribed as operons, or have functional coupling. This positional information is encoded as undirected edges between genes that are within a defined genomic distance threshold.

## 6.2 Transcriptomic Layer

The transcriptomic layer captures the expression state of each gene under the two experimental conditions. Each node in this layer corresponds to the same set of genes as the genomic layer, but its features are expression-derived rather than sequence-derived.

The transcriptomic representation includes:
- Mean expression level under control and treatment conditions
- Direction and magnitude of differential expression
- A categorical regulation label (up-regulated, down-regulated, or stable)

Relationships between transcriptomic nodes are established through expression pattern similarity. Genes with correlated expression profiles across replicates are likely to be co-regulated or functionally related. The architecture captures this through edges weighted by expression correlation.

## 6.3 Proteomic Layer

The proteomic layer models each protein as a node. Proteins are linked to their corresponding genes through the gene–protein mapping established during genome annotation parsing, providing the inter-omics bridge.

The proteomic representation includes:
- Mean protein abundance under control and treatment conditions
- Differential abundance between conditions
- Protein regulation status

Relationships between proteins derive from two sources: (1) abundance correlation across replicates, which can indicate shared regulatory mechanisms or complex membership, and (2) known protein–protein interactions retrieved from external databases (STRING), which capture physical and functional associations.

---

# 7. Capturing Intra-Omics Relationships

Intra-omics relationships model how entities within the same molecular layer relate to one another. These relationships carry specific biological meaning and are preserved as distinct edge types in the fused graph.

## 7.1 Genomic Proximity (Genomics)

Adjacent genes on a bacterial chromosome are often functionally linked. They may be organized into operons — co-transcribed units where a single promoter drives expression of multiple genes — or they may share regulatory regions. The architecture captures these relationships as genomic proximity edges connecting genes on the same replicon whose intergenic distance falls below a configurable threshold.

The biological rationale is that proximity implies potential regulatory or functional coupling. While not every adjacent gene pair is functionally related, the aggregate pattern of genomic organization provides meaningful structural context for the graph.

## 7.2 Expression Similarity (Transcriptomics)

Genes with correlated expression patterns across conditions and replicates are likely to be co-regulated or to participate in the same biological processes. The architecture computes pairwise Pearson correlation coefficients between all gene expression profiles and connects genes whose absolute correlation exceeds a threshold.

This approach captures relationships such as:
- Co-regulated stress response genes
- Members of the same metabolic pathway whose expression is coordinated
- Genes under the control of the same transcription factor

The resulting co-expression network is a standard systems biology tool; embedding it within the multi-omics graph adds the cross-omics context that a standalone co-expression analysis lacks.

## 7.3 Protein–Protein Interactions (Proteomics)

Physical and functional interactions between proteins are retrieved from the STRING database, which integrates evidence from experimental assays, curated databases, co-expression analysis, genomic context, and text mining. Each interaction carries a combined confidence score.

The architecture queries STRING for all proteins in the dataset, using the appropriate species identifier and applying a minimum confidence threshold. For the primary *E. coli* strain (B36, genomically closest to CFT073), proteins that are not mapped to the primary species are re-queried against a close relative (*E. coli* K-12 MG1655) as a fallback. Results from both queries are merged, deduplicated by retaining the highest confidence score for each pair.

PPI edges capture functional relationships that are not evident from genomic position or transcriptomic correlation alone, adding an orthogonal dimension to the graph.

---

# 8. Capturing Inter-Omics Relationships

Inter-omics relationships connect nodes across different molecular layers, creating cross-layer paths for information flow. These edges are the key innovation of the architecture: they enable the graph to represent the central dogma as a computational graph.

## 8.1 Gene–Protein Encoding Relationship

The most fundamental inter-omics relationship is the encoding relationship: each protein is produced from a specific gene. This relationship is established during genome annotation parsing, where CDS entries in the GFF file contain a protein_id attribute linking the genomic feature to its protein product.

In the heterogeneous graph, each encoding relationship creates a directed edge from a gene node to its corresponding protein node. This edge provides the bridge between the genomic and proteomic layers. Without it, the two layers would be disconnected subgraphs.

## 8.2 Transcript–Gene Correspondence

Transcriptomic measurements are inherently tied to genes — each expression value quantifies the transcription of a specific gene. In the architecture, this correspondence is used for alignment rather than as an explicit edge type. The same gene identifier indexes both the genomic node and the transcriptomic node features. The relationship is further captured through the multi-omics hyperedge (Section 9), which groups a gene, its expression profile, and its protein abundance into a single higher-order structure.

## 8.3 Cross-Layer Integration

The inter-omics edges ensure that information can propagate across molecular levels. In the fused graph, a protein node is connected to its encoding gene, which is connected to neighboring genes through genomic proximity, which in turn are connected to their expression profiles and protein products. This creates a connected representation where the graph structure mirrors the biological flow from genome to transcript to protein.

The integration is fundamentally different from simple feature concatenation. In a concatenated approach, each gene would be represented by a fixed vector (genomic features + expression values + protein abundance), and all relational information — which genes are adjacent, which proteins interact, which genes are co-expressed — would be discarded. The graph approach retains these relationships as first-class structures.

---

# 9. Graph Fusion Strategy

The architecture constructs separate representations of each omics layer before fusing them into a single heterogeneous graph. This staged approach preserves the distinct semantics of each relationship type while enabling cross-layer information flow.

## 9.1 Why Separate Graphs First

Each omics layer has its own measurement modality, noise characteristics, and relationship types. Genomic proximity is deterministic (based on coordinates), co-expression is statistical (based on correlation), and PPI is evidential (based on external database confidence scores). Processing these separately allows each to be treated with the appropriate method before integration.

## 9.2 The Fusion Mechanism

Fusion occurs by constructing a heterogeneous graph where:
- Node types correspond to biological entity types (gene, protein, annotation)
- Edge types correspond to relationship types (genomic_proximity, encodes, ppi, co_expression, in_pathway, cog_category)
- Each edge type carries its own semantic meaning and can be processed independently or jointly by downstream models

The fusion is biologically motivated: a heterogeneous graph with typed edges preserves the distinction between "these two genes are adjacent on the chromosome" and "these two proteins physically interact" and "this gene encodes this protein." A homogeneous graph that treated all edges identically would collapse these distinct biological relationships into a single undifferentiated connection.

## 9.3 Hypergraph Complement

The heterogeneous graph captures pairwise relationships well, but some biological structures are intrinsically higher-order. A KEGG pathway involves multiple genes; a COG category groups genes by functional annotation; a PPI community detected by clustering represents a functional module. These set-level relationships are captured by a complementary hypergraph, where each hyperedge connects the set of nodes that belong to a common biological module.

The hypergraph adds:
- Multi-omics triples (gene + its expression + its protein)
- Pathway hyperedges (all genes in a KEGG pathway)
- Functional category hyperedges (all genes in a COG category)
- PPI community hyperedges (dense interaction clusters detected by community detection)

The heterogeneous graph and hypergraph together provide a dual representation: pairwise relationships for message passing and set-level relationships for module-level reasoning.

---

# 10. Advantages of the Proposed Architecture

**Compared to single-omics analysis.** Single-omics studies (genomics alone, transcriptomics alone, proteomics alone) capture only one dimension of biological activity. Gene presence does not imply expression; transcript abundance does not guarantee protein abundance; protein abundance alone does not reveal regulatory mechanisms. The proposed architecture integrates all three views, allowing relationships across layers to inform the representation of each entity.

**Compared to feature concatenation.** The simplest multi-omics approach is to concatenate measurements from each layer into a fused feature vector per gene. This discards all relational information — both within layers (adjacency, co-expression, PPI) and between layers (who encodes whom). The graph representation preserves these relationships as explicit structures, enabling relational reasoning that concatenation cannot support.

**Compared to independent layer analysis.** Analyzing each omics layer independently and comparing results post hoc (e.g., overlapping lists of differentially expressed genes and differentially abundant proteins) treats the layers as separate experiments. The proposed architecture unifies them in a single computational framework where cross-omics information flow is built into the representation.

**Information complementarity.** Each omics layer contributes unique information that the others cannot provide:
- Genomics provides the invariant reference framework — the gene catalogue, structural organization, and functional potential encoded in the genome.
- Transcriptomics provides dynamic information about which genes are active under specific conditions.
- Proteomics provides information about the actual functional state of the cell, capturing post-transcriptional and post-translational regulation.
- External annotations (pathways, functional categories, PPIs) provide functional context not derivable from the omics measurements alone.

---

# 11. Current Status

The architecture is fully implemented and has been validated using real multi-omics data from five *E. coli* strains (B36, MS_14384, MS_14385, MS_14386, MS_14387) spanning genomics, transcriptomics, proteomics, and metabolomics (MetaboLights MTBLS2015).

**Completed components:**
- Multi-omics preprocessing pipeline (genomics, transcriptomics, proteomics, metabolomics parsers)
- Cross-omics alignment and integration (3-layer central-dogma merge, with a 2-layer genome+proteome fallback when RNA is unavailable)
- Cross-referencing with external knowledge bases (UniProt, STRING, KEGG, eggNOG-mapper)
- Heterogeneous graph construction (PyG HeteroData with typed nodes and edges)
- Hypergraph construction (multi-omics triples, pathway, COG, PPI cluster, GO, and AMR-mechanism hyperedges)
- AMR knowledge layer — curated resistance determinants per strain (`amr/amr_manifest.json`) attached as annotation hyperedges/nodes; genome-only markers remain represented
- Visualization modules (network layouts, feature distributions, incidence heatmaps, AMR highlight)
- Data export (CSV, NPY, PT formats for downstream use)
- Cross-strain AMR comparison (per-strain reports, presence matrix, synergies)

**Not yet implemented:**
- GNN-based learning on the constructed graph
- Training and evaluation on AMR phenotype prediction (only four *E. coli* strain graphs — insufficient for supervised train/test; see `AMR_GRAPH_INTEGRATION.md` §8.4)
- Multi-species graphs (*Klebsiella pneumoniae*, *Staphylococcus aureus*, *Streptococcus pyogenes* from the paper) — deferred (see `AMR_GRAPH_INTEGRATION.md` §8.1)

The architecture produces a file-based representation that is ready for GNN training per strain: a HeteroData object with gene/protein/metabolite/annotation node types and typed edges, complemented by a hypergraph of module-level hyperedges (see `outputs/<strain>/graph/`).

> **AMR update:** See `multiomics_graph/AMR_GRAPH_INTEGRATION.md` for the separate Phase 1–3 documentation of the AMR knowledge layer, WP/COG entity resolution, genome-only marker visualization, and cross-strain AMR analysis.

---

# 12. Future Work

**AMR phenotype prediction.** Curated resistance determinants are already attached as an annotation layer (§3.5), but the current *E. coli* strain set lacks matched resistance-phenotype labels and is too small for supervised GNN training. Applying the pipeline to a dataset with AMR phenotypic labels would enable the full prediction pipeline — from raw data through graph construction to AMR classification (see `AMR_GRAPH_INTEGRATION.md` §8.4).

**Expansion of curated resistance genes.** The curated manifest covers the study's known determinants. Incorporating community AMR gene databases (CARD, ResFinder, ARG-ANNOT) would broaden coverage and support orthology to multi-species contexts.

**Multi-species graphs.** The study also covers *Klebsiella pneumoniae*, *Staphylococcus aureus*, and *Streptococcus pyogenes*. Building separate strain-level graphs per species requires species-specific PPI, KEGG, and functional mappings (deferred — see `AMR_GRAPH_INTEGRATION.md` §8.1).

**Graph neural network learning.** The heterogeneous graph and hypergraph are structured inputs for heterogeneous GNNs (RGCN, HAN) and hypergraph neural networks (HGCN) with type-specific transformations.

**Explainability analysis.** GNN explainability methods (GNNExplainer, integrated gradients) could identify which nodes, edges, and features drive representations — potentially revealing biological mechanisms underlying resistance or identifying novel resistance determinants from the graph structure.

---

*Document generated for the Multi-Omics Graph Fusion Architecture project. This document describes methodology and rationale; for implementation details, refer to the technical documentation and source code.*
