# AMR Graph Integration — Documentation

Multi-Omics Graph Fusion for Antimicrobial Resistance (AMR) Analysis
*Escherichia coli* — B36, MS_14384, MS_14385, MS_14386, MS_14387

---

## 1. Objective

Enrich and interpret the fused multi-omics graph using **known AMR determinants as
biological references**. The goal is *not* to build a direct AMR classifier; it is:

> Using multi-omics graph fusion to characterize and analyze antimicrobial
> resistance mechanisms.

AMR knowledge is attached to the graph as a **biological annotation layer** and used
for interpretation, visualization, and cross-strain comparison.

### Guiding principles
- **Entity resolution first.** Fix identifier integration (WP version-suffix mismatch)
  before any AMR analysis. Without it, COG/functional annotations disconnect from the
  graph.
- **AMR = knowledge, not a prediction target.** AMR mechanisms are added as annotation
  nodes (`AMR:<class>`) and hyperedges (`amr_mechanism`).
- **Genome-only resistance genes are valid.** A marker with no RNA/protein quantitation
  must remain represented (resistance can exist without detection).
- **No supervised GNN AMR prediction.** Only 4 *E. coli* strain graphs → insufficient
  samples for meaningful train/test evaluation.
- **Multi-species data is out of scope** for now (see §8.1).

---

## 2. Deliverables

| Artifact | Path |
|----------|------|
| Shared protein-ID normalizer | `annotation/identifiers.py` |
| eggNOG loader (dual-key + locus re-keying) | `annotation/eggnog.py` |
| AMR manifest (curated genotypes) | `amr/amr_manifest.json` |
| AMR manifest loader | `amr/amr.py` |
| AMR mechanism hyperedges | `graph/hypergraph.py::add_amr_hyperedges` |
| AMR highlight visualization | `visualization/network_plots.py::plot_amr_highlight` |
| AMR figures regen + 100% coverage verify | `scripts/amr_figures.py` |
| AMR reports + cross-strain comparison | `scripts/amr_report.py` |
| Validation tests | `tests/test_cog_resolution.py`, `tests/test_amr.py`, `tests/test_corrections.py` |

---

## 3. Phase 1 — Data Correctness & Entity Resolution

### 3.1 Problem

The genome/proteomics/graph nodes use **versioned** RefSeq protein IDs
(`WP_000000542.1`), while the eggNOG annotation file is keyed by **unversioned** IDs
(`WP_000000542`). Additionally, downstream COG code treated those unversioned WP keys
as if they were **locus_tag GeneIDs**. The consequence for every strain:

- COG edges attached to **orphan WP nodes** (no matching graph node).
- COG hyperedges **missed actual gene/protein nodes**.
- Functional annotations were **disconnected** from the multi-omics graph.

### 3.2 The biological identity hierarchy

```
locus_tag                        e.g. EW036_RS26095
    |
    |  gene_protein_mapping      {locus_tag: versioned_protein_id}
    v
protein_id (versioned)           e.g., WP_000000542.1
    |
    |  normalize_wp_id()
    v
protein_id (unversioned)         e.g., WP_000000542
    |
    v
eggNOG / COG annotation          keyed by unversioned WP accession
```

The fix **does not** globally replace identifiers. It maintains this hierarchy and simply
resolves COG annotations through the mapping at the join points.

### 3.3 Implementation

#### `annotation/identifiers.py`
```python
def normalize_wp_id(protein_id) -> str:
    """WP_000000542.1 -> WP_000000542 (strip numeric version suffix)."""
```
A single shared normalizer used by the eggNOG loader and downstream joins.

#### `annotation/eggnog.py` — `EggNOGAnnotator`
- `load_precomputed(filepath)`: for each row, builds a payload and registers it under
  **both** the normalized (unversioned) key and the original id, and stores
  `original_protein_id` for traceability. Compatibility with both id-spaces is preserved.
- New `as_locus_tag_map(gene_protein_map) -> {locus_tag: annotation}`: re-keys COG
  annotations onto real locus-tag gene nodes. For each `locus_tag → protein`, it looks
  up the annotation by the protein id (versioned) *or* its normalized form.

#### `main.py` — `run_annotation()`
- After loading eggNOG, resolves `self.cog_annotations` to **locus-tag keys** via
  `eggnog.as_locus_tag_map(self.gene_protein_map)`. Fallback to placeholder only if
  nothing resolves. This makes COG edges (gene→annotation) and COG hyperedges
  (`cog_category`) attach to real gene nodes.

### 3.4 Validation

`tests/test_cog_resolution.py`:
- `normalize_wp_id` behavior.
- eggNOG loader keys annotations under both normalized + original forms.
- `as_locus_tag_map` resolves via gene-protein mapping; unmatched proteins skipped.
- COG hyperedges contain real locus-tag gene nodes (plus their protein nodes).
- B36 `OXA-1` (a defense marker) falls in COG `V`.

### 3.5 Result

| Strain | Before | After |
|--------|--------|-------|
| B36 | 1 COG hyperedge | **21 COG category hyperedges**; 2283 annotations resolved |
| MS_14384 | 1 COG hyperedge | **21 COG category hyperedges**; 2348 annotations resolved |

B36's `blaOXA-1` gene (EW036_RS26250) now correctly maps to COG `V`
(defense mechanisms).

---

## 4. AMR verification across strains

All paper genotypes were traced into each strain's NCBI `genomic.gff`.

| Strain | Markers found | Notes |
|--------|---------------|-------|
| B36 (ST131) | tet(A), dfrA17, aadA5, sul1, mph(A), **blaCTX-M-15**, **blaOXA-1**, aac(6')-Ib-cr5 → **8/8** | plasmid-borne (NZ_LR130546.1) |
| MS_14384 (ST963) | blaCMY-2 → **1/1** | chromosomal ampC/BlaEC is distinct |
| MS_14386 (ST224) | blaTEM-1b, blaTEM-215, aac(3)-IId, aadA1*, aph(3')-Ia, strAB, sul2, sul3, dfrA12, dfrA14, floR, fosA4, mph(A) → **13/13** | *aadA1 = frameshifted pseudogene (null) |
| MS_14387 (ST69) | blaTEM-1b, strAB, sul2 → **3/3** | |

**Key caveats handled in the manifest:**
- `aadA1` is a frameshifted pseudogene with no protein accession → `protein_id: null`.
- `blaTEM-1b` vs `blaTEM-215` map to two distinct loci but share the GFF product name
  "TEM-1" → disambiguated by **locus_tag**, not name.
- `floR` (MS_14386), `strA` (MS_14387) present in genome but not intro aligned omics →
  treated as genome-only AMR nodes.

---

## 4. Phase 2 — AMR Knowledge Integration

### 4.1 `amr/amr_manifest.json` — Manifest (curated)

Structure (from paper + verified GFF):
```json
{
  "B36": {
    "blaCTX-M-15": {
      "name": "blaCTX-M-15",
      "amr_class": "beta_lactam",
      "locus_tag": "EW036_RS26225",
      "protein_id": "WP_000239590.1",
      "source": "paper (GFF: ESBL CTX-M-15)"
    }
  }
}
```
- `protein_id` may be `null` or a list (e.g. `strAB` covers two loci).
- Classes follow `AMR_CLASS_NAMES` (beta_lactam, aminoglycoside, sulfonamide,
  trimethoprim, tetracycline, macrolide, fluoroquinolone, phenicol, fosfomycin).

### 4.2 `amr/amr.py` — AMRManifest loader

API:
- `strains()` — strain keys.
- `markers(strain)` → normalized marker records.
- `locus_to_amr(strain)` → `{locus_tag: marker}`.
- `class_to_loci(strain)` → `{amr_class: [locus_tag, ...]}`.
- `amr_loci(strain)` → deduplicated AMR locus tags.
- `class_name(amr_class)` → human display name.
- `summary(strain)` → marker/locus/class counts.

### 4.3 AMR graph representation

The graph node-type set is extended from
`gene / protein / metabolite / annotation` to include:

```
gene  →  associated_with  →  AMR_mechanism   (e.g., beta_lactam_resistance)
```

`graph/hypergraph.py::add_amr_hyperedges(class_to_loci, amr_loci,
gene_to_protein=None, class_names=None)` creates one `amr_mechanism` hyperedge per
class, always attaching the **`AMR:<class>` mechanism node** so genome-only markers
remain represented:

```
beta_lactam_resistance: { blaCTX-M-15, blaOTH-1, associated genes, AMR:beta_lactam }
```

Node-type inference in `_infer_node_types` recognizes the `AMR:` prefix →
`amr_mechanism` node type.

### 4.4 Pipeline wiring (`main.py`)
- `__init__`: `self.amr_manifest = None`.
- `run_annotation`: loads `AMRManifest`, logs marker/locus/class counts.
- `run_hypergraph`: builds AMR mechanism hyperedges (step 3b).
- `run_export`: adds `is_amr_gene` (bool) and `amr_class` (semicolon-joined) to
  `gene_features.csv`.

### 4.5 Visualization
`visualization/network_plots.py`:
- `plot_amr_highlight(node_types, edge_types, amr_loci, amr_class_map, ...)` →
  `amr_highlight.png`. AMR genes enlarged and colored by mechanism class; background
  nodes faded.
- **Every manifest AMR locus is kept as a gene node**, even when it has no
  RNA/protein quantitation (genome-only determinants are real genomic markers
  and are never removed just because they are absent from
  `aligned_multiomics.csv`).
- **Visual distinction of omics support:**
  - solid node → AMR gene backed by multi-omics (RNA/protein) quantification;
  - outlined + hatched node → AMR gene detected only at the genome level
    (status stays explicitly `genome_only`; no RNA/protein values are invented).
- `AMR:<class>` mechanism nodes are drawn as **diamonds**, connected to their
  member genes via `amr_mechanism` hyperedge edges (β-lactam, aminoglycoside,
  sulfonamide, trimethoprim, tetracycline, macrolide, …).
- The caller merges the AMR layer into the draw graph (all loci as `gene_*`
  nodes + `AMR:<class>` nodes + hyperedge edges) before rendering; the draw
  graph therefore contains the full paper genotype set.
- Auto-verification: `plot_amr_highlight(..., verify=True)` asserts that every
  `amr_class_map` key is present in the drawn graph and returns the per-node
  state map (`multi_omics` | `genome_only`).

### 4.6 Verification
- B36: 8 markers, 8 loci, 7 classes → 7 AMR mechanism hyperedges; `amr_highlight.png`.
- MS_14386: 13 markers, 14 loci (aadA1 + strAB dual) → 7 AMR hyperedges.
- MS_14387: 3 markers, 4 loci → 3 AMR hyperedges.
- `gene_features.csv` correctly flags OXA-1: `is_amr_gene=True, amr_class=beta_lactam`.
- **100% coverage check**: `scripts/amr_figures.py` regenerates the figures and
  verifies that **all** paper AMR loci appear as drawn nodes (multi-omics +
  genome-only):
  - B36 8/8, MS_14384 1/1, MS_14386 14/14, MS_14387 4/4 — PASS.

---

## 5. Phase 3 — Biological Analysis

`scripts/amr_report.py` reads the finished `outputs/<strain>/` and the manifest
(read-only; it does not re-run the pipeline).

### 5.1 Per-strain report (`<strain>_amr_report.csv`)
Columns: `strain, marker, amr_class, locus_tag, protein_id, COG, n_kegg_pathways,
n_ppi_neighbors, in_genome, in_aligned_omics`.

### 5.2 Cross-strain comparison
- `amr_cross_strain_matrix.csv` — marker × strain binary presence.
- `amr_cross_strain_heatmap.png`.
- `amr_synergies.txt` — core / shared / strain-specific markers.

### 5.3 Biological takeaways
- **Beta-lactam** is the only core mechanism present in all four strains.
- Shared markers (≥2 strains): `blaTEM-1b`, `strAB`, `sul2`, `mph(A)`.
- 17 strain-specific markers reinforce the paper's multi-drug-resistant (MDR) strain
  definition (e.g. MS_14386's `floR`, `fosA4`, `sul3`, `dfrA14`).

### 5.4 Determinant path tracing (`scripts/amr_path_analysis.py`)
Reads the saved `outputs/<strain>/graph/` edge tables and traces each manifest
determinant across **eight evidence layers** (genome → transcriptome →
proteome → metabolome → KEGG → COG → GO → PPI/AMR mechanism), plus a bounded
multi-hop BFS context through gene/protein/metabolite/annotation nodes
(`reports/amr_path_*`). Layer statuses are *observed / not_observed /
unavailable* — missing measurements are never converted to zeros.

Coverage-recovery helpers for genome-only determinants (no network):
- `scripts/convert_emapper_to_tsv.py` now carries emapper `KEGG_ko` and
  `KEGG_Pathway` columns into `strains/<name>/eggnog_annotations.tsv`;
  `annotation/eggnog.py::kegg_pathway_membership()` exposes them re-keyed onto
  locus tags. `main.py` merges these into `pathway_membership`, so the graph's
  KEGG layer covers acquired/plasmid-borne genes the `eco`-bridged UniProt
  route misses.
- `main.py` builds gene nodes from `genome_genes ∪ aligned`, and merges
  recovered UniProt GO terms (below) into `protein_to_go`.
- `scripts/recover_amr_uniprot.py` — per-ID UniProt GO lookup for AMR
  determinants (bypasses the pipeline's 500-hits-per-batch truncation) →
  `reports/amr_uniprot_recovery.csv`.
- `scripts/enrich_amr_graph.py` — idempotently merges eggNOG KEGG/COG,
  recovered GO, encodes, and genomic-proximity edges for the determinants into
  the saved graphs, and updates the KEGG/COG/GO hyperedges so BFS path tracing
  connects through them.
- `scripts/amr_path_analysis.py` falls back to `transcriptomics_expression.csv`
  for RNA when a determinant has no aligned multi-omics row.

Candidate ranking (`reports/amr_path_candidates.csv`, report §G) resolves
CHEBI metabolite ids to names and gene/protein ids to gene symbols/products
from `genome_genes.csv`, and de-biases the hypothesis list: entities connected
to ≥ `HUB_FRACTION` (0.75) of the maximum determinant count reached by any
candidate are flagged `hub_candidate=True` (ubiquitous pathway-shared
metabolites such as central-metabolism compounds) and ordered after the
specific candidates, so informative elements such as `tolC`, `mrdA`, `cysK`
and `oppB` surface at the top.

Result: B36's genome-only determinants (`tet(A)`, `aadA5`, `blaCTX-M-15`)
move from `genome_only` (1/8 layers) to `multi_omics` (3–5/8 layers) with real
transcriptome and COG (and, for `blaCTX-M-15`, KEGG) evidence — no values
invented.

---

## 6. Tests

| File | Covers |
|------|--------|
| `tests/test_corrections.py` | 2-layer triple masking, 3-layer agreement, cross-instrument normalization, exact-name KEGG |
| `tests/test_cog_resolution.py` | WP-id normalization, dual-key loading, locus re-keying, COG `V` |
| `tests/test_amr.py` | manifest loading, class grouping, genome-only pseudogene, AMR hyperedges incl. mechanism node |
| `tests/test_amr_path_analysis.py` | eggNOG KEGG parsing, `_rna_row` transcriptomics fallback, layer status logic, enrichment idempotency, candidate hub de-biasing, hyperedge update bump |

Run: `python -m pytest tests/` (73 passing).

---

## 7. Usage

```bash
# run one strain through the pipeline (does Phase 1 + Phase 2 wiring)
python multiomics_graph/main.py --strain B36

# generate AMR per-strain + cross-strain reports from saved outputs
python scripts/amr_report.py --output reports

# regenerate AMR highlight figures + verify 100% paper-locus coverage
python scripts/amr_figures.py

# per-ID UniProt GO recovery for AMR determinants (offline, no network)
python scripts/recover_amr_uniprot.py

# merge eggNOG KEGG/COG + recovered GO + encodes/proximity edges into the
# saved graphs for the determinants (idempotent)
python scripts/enrich_amr_graph.py

# trace each determinant across the 8 evidence layers + path context
python scripts/amr_path_analysis.py

# visualise a single AMR determinant across the 8 omics layers + its
# propagation through the graph (two-panel figure per marker)
python scripts/amr_gene_trace_viz.py --strains B36 --marker blaOXA-1

# run all tests
python scripts (from multiomics_graph/): python -m pytest tests/
```

---

## 8. Limitations & Future Work

### 8.1 Multi-species data (deferred)
The paper also covers *Klebsiella pneumoniae*, *Staphylococcus aureus*, and
*Streptococcus pyogenes*. Do not merge at this stage. A future extension should build
**separate strain-level graphs per species**; cross-species comparison requires:
species-specific PPI, KEGG mapping, functional annotation, and orthology/conserved
representation.

### 8.2 Multi-species ortholog representation (future)
To compare mechanisms across E. coli and other species, orthologous (conserved) gene
sets must be added (e.g. eggNOG OG mapping), allowing a common representation across
species.

### 8.3 eggNOG annotations
All five strains (B36, MS_14384, MS_14385, MS_14386, MS_14387) now have real
eggNOG-mapper annotations (`strains/<name>/eggnog_annotations.tsv`, 20–22 COG
categories, 4426–4820 proteins each), generated with fresh `emapper.py`
`--dmnd_iterate no` runs and converted via
`scripts/convert_emapper_to_tsv.py`. MS_14385's genome was corrected to its
proper assembly GCF_900622655.1 (`EW040_RS` locus prefix; 5322 CDS; distinct
proteome), replacing the previously mis-assigned GCF_900622665.1 (which
belonged to MS_14386). Reproduce the conversion with:

```bash
python scripts/convert_emapper_to_tsv.py MS_14385
```

The converter now also carries the emapper `KEGG_ko` and `KEGG_Pathway`
columns through to the TSV. `annotation/eggnog.py` parses these, and
`EggNOGAnnotator.kegg_pathway_membership()` exposes
`{pathway_id: [locus_tag, ...]}` re-keyed onto graph gene nodes, so
acquired/plasmid-borne AMR determinants (invisible to the `eco`-bridged
UniProt KEGG cross-references) still get a KEGG pathway layer. The saved
per-strain graphs are enriched in-place with these memberships plus eggNOG
COG, recovered UniProt GO, encodes, and genomic-proximity edges for AMR
determinants by:

```bash
python scripts/enrich_amr_graph.py
```

### 8.4 Supervised GNN AMR prediction (out of scope)
With only 4 strain graphs, meaningful train/test evaluation is not feasible. The
appropriate objective remains *graph fusion to characterize and analyze AMR
mechanisms*.