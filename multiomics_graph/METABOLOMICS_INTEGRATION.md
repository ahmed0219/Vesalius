# Metabolomics Integration Layer

## Technical Documentation — Metabolite Layer

---

# 1. Overview

## 1.1 Objective

Add a **metabolomics layer** to the multi-omics heterogeneous graph + hypergraph. Metabolite abundances (RPMI vs Sera) are parsed from MetaboLights MAF files (MTBLS2015), and each metabolite is enriched with KEGG pathway and enzyme (EC number) membership via the ChEBI → KEGG compound resolver. The result is a fourth biological entity type — `metabolite` nodes — linked into the graph through:

- `(metabolite, in_pathway, annotation)` — metabolite → KEGG pathway edges (reusing the annotation node set shared with genes)
- `(metabolite, metabolised_by, protein)` — metabolite → enzyme (protein) edges bridged through EC numbers
- `(metabolite, abundance_correlation, metabolite)` — metabolite ↔ metabolite edges from Pearson correlation over the per-sample abundance columns, capped by the same mutual top-k rule as the RNA/protein correlation families (threshold `--correlation`, default 0.7; `--correlation-top-k`, default 10)
- `metabolic_pathway` hyperedges co-joining genes and metabolites on the same KEGG pathway
- `go_term` hyperedges co-joining proteins sharing the same Gene Ontology term (BP / MF / CC, from UniProt annotations)

## 1.2 Pipeline Position

The metabolomics stage runs after Proteomics and before Integration:

```
Genomics → Transcriptomics → Proteomics → Metabolomics →
Integration → Annotation → Graph Construction → Hypergraph Construction →
Visualization → Export
```

Implementation: `StrainPipeline.run_metabolomics()` in `main.py`.

## 1.3 Data Sources

| Item | Location | Description |
|------|----------|-------------|
| GC-MS MAF | `metabolic/m_MTBLS2015_GC-MS___metabolite_profiling_v2_maf.tsv` | 118 metabolite rows, 59 sample columns |
| LC-MS MAF | `metabolic/m_MTBLS2015_LC-MS_negative_hilic_metabolite_profiling_v2_maf.tsv` | 170 metabolite rows, 60 sample columns |
| Study design | `s_MTBLS2015.txt` | Sample → strain / growth-media mapping (B36, MS14384, MS14385, MS14386, MS14387) |
| KEGG REST | `https://rest.kegg.jp` | ChEBI cross-reference table + compound records |

**Sample assignment (from `s_MTBLS2015.txt`):**

| Strain | RPMI samples | Sera (Pooled sera) samples |
|--------|--------------|----------------------------|
| B36 | 50913–50918 | 50919–50924 |
| MS_14384 | 50925–50930 | 50931–50936 |

---

# 2. Module Documentation

## 2.1 `preprocessing/metabolomics.py` — `MetabolomicsProcessor`

Parses one or more MetaboLights MAF TSV files and produces a per-metabolite summary for one strain.

### Public Methods

| Method | Description |
|--------|-------------|
| `parse_sample_file()` | Read the `s_*.txt` sample sheet, map `{sample_key: growth_media}` |
| `load_raw()` | Load and concatenate all MAF files into a long DataFrame |
| `load(rpmi_samples, sera_samples, log2=True, min_detection=0.0, condition_map=None)` | Compute condition means and logFC; returns the summary table |
| `summary()` | `{metabolites, up_regulated, down_regulated, stable}` |

### Key Behaviour

- **Sample columns**: any MAF column matching `^\d+$` or `^102.\d+.\d+/\d+$` is treated as a sample abundance column.
- **Condition split**: explicit `rpmi_samples` / `sera_samples` lists (from the strain config) take precedence; otherwise inferred from the condition map (values starting with `RPMI`/`Sera`).
- **Metabolite ID**: `CHEBI:<accession>` when a `database_identifier` of form `CHEBI:NNNN` is present, else `M-<metabolite_identification>`.
- **Log transform**: means are `log2(mean_abundance + 1)`; `Metabolite_logFC = Metabolite_Sera - Metabolite_RPMI`.
- **Cross-instrument averaging**: identical MetaboliteIDs across instruments are merged by averaging the per-condition replicate means (and the logFC is recomputed from the averaged means), so evidence from both GC-MS and LC-MS is combined.
- **Abundance filtering**: `min_detection > 0` drops metabolites below the threshold in both groups.

### Output Table Columns

`MetaboliteID`, `MetaboliteName`, `Database`, `ChemicalFormula`, `MassToCharge`, `RetentionTime`, `Metabolite_RPMI`, `Metabolite_Sera`, `Metabolite_logFC`

## 2.2 `annotation/kegg_compound.py` — `KEGGCompoundResolver`

Maps ChEBI accessions → KEGG compound → enzymes (EC numbers) + pathways.

### Mapping Strategy

The ChEBI REST API is unreachable from the build environment, so the resolver uses the KEGG cross-reference table instead:

1. Download the full `rest.kegg.jp/conv/chebi/cpd` table (17 112 ChEBI → compound mappings), cached to `outputs/cache/kegg_chebi_cpd.json`.
2. Invert to `{chebi_accession: [kegg_compound_id]}`.
3. For each KEGG compound id, fetch `rest.kegg.jp/get/cpd:Cxxxxx` and parse the `ENZYME` (EC numbers) and `PATHWAY` blocks. Records cached to `outputs/cache/kegg_compounds/*.json`.

### Public Methods

| Method | Description |
|--------|-------------|
| `resolve(chebi_ids)` | Return `{chebi_id: {kegg_id, name, formula, enzymes, pathways}}` |
| `find_by_name(names)` | **Name-based fallback** — for metabolites without a ChEBI → KEGG mapping, query `rest.kegg.jp/find/compound/<name>` and resolve the first hit; returns the same `{name: {...}}` record shape |
| `fetch_compound(cpd_id)` | Fetch + cache a single KEGG compound record |
| `summary()` | Resolver statistics |

In `main.py`, `run_metabolomics()` first resolves by ChEBI, then applies `find_by_name` to the remaining (ChEBI-less / unresolved) metabolites and maps the name-keyed records back to `MetaboliteID` for edge construction.

### Parsing Notes

- Only valid EC numbers (`^\d+\.\d+\.\d+\.\d+$`) are retained; artifacts like `2.4.1.-` are dropped.
- Only KEGG pathway entries (`map*`) are retained.

## 2.3 `graph/heterogeneous_graph.py` — new methods

| Method | Description |
|--------|-------------|
| `add_metabolite_nodes(metabolite_ids, metabolite_data)` | Register `metabolite` nodes with 5-dim feature vectors |
| `add_metabolite_functional_edges(metabolite_to_func, func_to_name, relation='in_pathway')` | `(metabolite, in_pathway, annotation)` edges to shared annotation nodes |
| `add_metabolite_enzyme_edges(metabolite_to_enzymes, enzyme_to_protein)` | `(metabolite, metabolised_by, protein)` EC-bridged edges |
| `add_metabolite_correlation_edges(metabolite_ids, metabolite_data)` | `(metabolite, abundance_correlation, metabolite)` edges from Pearson correlation over the per-sample matrix (mutual top-k capped) |

`build()` gains optional parameters: `metabolite_ids`, `metabolite_data`, `metabolite_to_pathway`, `pathway_names`, `metabolite_to_enzymes`, `enzyme_to_protein`, `metabolite_replicates` (per-sample matrix retained by `MetabolomicsProcessor.raw_replicates` for the abundance-correlation family).

## 2.4 `graph/hypergraph.py` — new methods

`add_metabolite_hyperedges(metabolite_to_pathway, pathway_genes=None)` creates one `metabolic_pathway` hyperedge per KEGG pathway, containing its metabolites and (optionally) the genes already assigned to that pathway.

`add_go_hyperedges(protein_to_go, go_info)` creates one `go_term` hyperedge per Gene Ontology term (aspect-filtered), co-joining all proteins annotated with that term. In `run_hypergraph()` the protein → GO map is built from the UniProt annotation `GO_BP`, `GO_MF`, `GO_CC` columns and term names/aspects from the annotator's `go_terms` registry.

## 2.5 `visualization/network_plots.py` — new method

`plot_metabolite_regulation(metabolite_table, filename)` renders a logFC histogram plus a top-changed-metabolite bar chart. `metabolite` was added to the `COLORS` scheme.

---

# 3. Metabolite Node Features

Metabolite nodes carry a **5-dimension** feature vector in `HeteroData['metabolite'].x`:

| Index | Feature | Scaling |
|-------|---------|---------|
| 0 | `Metabolite_RPMI` | log2 mean abundance (RPMI replicates) |
| 1 | `Metabolite_Sera` | log2 mean abundance (Sera replicates) |
| 2 | `Metabolite_logFC` | Sera − RPMI |
| 3 | `MassToCharge` | / 1e3 (m/z ≈ 0.1-scale) |
| 4 | `RetentionTime` | / 1e6 (minutes) |

Node names are stored in `HeteroData['metabolite'].name`.

---

# 4. Edge Construction

## 4.1 `(metabolite, in_pathway, annotation)`

Each resolved KEGG pathway membership of a metabolite becomes an edge to the shared `annotation` node set — the same nodes that genes connect to via `(gene, in_pathway, annotation)`. This creates metabolite ↔ gene linkage through shared pathways without a dedicated pathway node type.

## 4.2 `(metabolite, metabolised_by, protein)`

The EC-number bridge:

1. `metabolite_to_enzymes` — `{MetaboliteID: [EC numbers]}` from the KEGG compound resolver.
2. `enzyme_to_protein` — `{EC: [ProteinID]}` built from the **UniProt `EC_number` column** of the loaded annotations.

Protein IDs are translated from UniProt's unversioned `WP_` form back to the versioned `WP_....1` ids used by the graph nodes before matching.

---

# 5. Configuration

Add a `metabolomics` block to `strains/<name>/config.json`:

```json
"metabolomics": {
  "enabled": true,
  "maf_files": [
    "../../metabolic/m_MTBLS2015_GC-MS___metabolite_profiling_v2_maf.tsv",
    "../../metabolic/m_MTBLS2015_LC-MS_negative_hilic_metabolite_profiling_v2_maf.tsv"
  ],
  "sample_file": "../../s_MTBLS2015.txt",
  "rpmi_samples": ["50913", "50914", "50915", "50916", "50917", "50918"],
  "sera_samples": ["50919", "50920", "50921", "50922", "50923", "50924"]
}
```

Config keys are resolved by `_load_strain_from_folder()` (paths relative to the strain directory). The legacy `strains_config.json` loader also supports the same keys.

---

# 6. Results (Current Run)

| Metric | B36 | MS_14384 |
|--------|-----|----------|
| Metabolites | 219 | 219 |
| Up-regulated (logFC > 0.5) | 46 | 53 |
| Down-regulated (logFC < −0.5) | 80 | 89 |
| KEGG-resolved (ChEBI + name fallback) | 178 | 178 |
| Metabolites with pathways | 168 | 168 |
| `in_pathway` edges | 1 674 | 1 674 |
| `metabolised_by` edges | 351 | 809 |
| `metabolic_pathway` hyperedges | 219 | 219 |
| `go_term` hyperedges | 779 | 779 |
| Proteins with EC annotation | 509 / 1299 | 509 / 1299 |

---

# 7. Output Artifacts

| File | Description |
|------|-------------|
| `outputs/{strain}/metabolomics_abundance.csv` | Per-metabolite summary table |
| `outputs/{strain}/metabolite_pathway_membership.csv` | `{MetaboliteID, KEGGPathway}` pairs (1 674 rows) |
| `outputs/{strain}/metabolite_enzymes.csv` | `{MetaboliteID, EC_number}` pairs (6 873 rows) |
| `outputs/{strain}/graph/edges_metabolite_in_pathway_annotation.csv` | Pathway edges (edge-list) |
| `outputs/{strain}/graph/edges_metabolite_metabolised_by_protein.csv` | Enzyme edges (edge-list) |
| `outputs/{strain}/figures/metabolite_regulation.png` | Regulation + top-changed figure |
| `outputs/cache/kegg_chebi_cpd.json` | ChEBI ↔ KEGG cross-reference (17 112 entries) |
| `outputs/cache/kegg_compounds/*.json` | Per-compound enzyme/pathway records |

---

# 8. Known Issues & Limitations

1. **ChEBI API unreachable**: the resolver relies on the KEGG cross-reference table plus a name-based `find/compound` fallback; compounds with neither a ChEBI → KEGG mapping nor a name match are not enriched (≈41 of 219 metabolites in the current run).
2. **Strain resolution**: the MAF files contain five strains; sample columns must be selected explicitly via `rpmi_samples` / `sera_samples`, otherwise all RPMI samples across all strains are averaged.
3. **Duplicate MAF metabolites**: GC-MS and LC-MS often detect the same compound; identical MetaboliteIDs are merged by averaging condition means rather than keeping the highest-detection row.
4. **EC match coverage**: enzyme edges depend on the overlap between KEGG compound EC lists and UniProt protein EC annotations; proteins without UniProt EC annotations are excluded.
5. **`_get_ec_numbers` fix**: UniProt's `CATALYTIC ACTIVITY` reaction is a single dict (not a list) in the API JSON; the parser now handles both shapes so `EC_number` populates correctly.
6. **Name lookup ambiguity**: `find/compound` returns the first (lexicographic) hit, which is not always the exact compound (e.g. `Succinic acid` resolves to oxaloacetate C00036); results are cached and recorded as-is.
