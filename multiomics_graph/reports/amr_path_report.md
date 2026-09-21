# AMR Multi-Omics Path Analysis

> **Scope.** This is an **exploratory, mechanism/context characterization** framework. It is **not an AMR classifier** and **not a causal inference system**. Graph connectivity does not imply causality; annotation-mediated relationships are not equivalent to experimental measurements; candidate elements are hypotheses, not confirmed resistance determinants.

## A. Dataset / analysis summary

- Strains analyzed: B36, MS_14384, MS_14386, MS_14387
- The analysis reads only saved pipeline outputs (`outputs/<strain>/` and `outputs/<strain>/graph/`).
- Bounded molecular-context paths: <= 4 hops (shortest first, max 400 path records per determinant).
- Eight layers are traced per determinant; layer coverage is a **descriptive integration measure** and is **not** an evidence strength or biological importance score.

## B. Per-strain AMR determinant table

| Strain | Determinants | Multi-omics | Genome-only | Mean layer coverage |
|---|---|---|---|---|
| B36 | 8 | 8 | 0 | 6.12/8 |
| MS_14384 | 1 | 1 | 0 | 8.00/8 |
| MS_14386 | 13 | 13 | 0 | 6.54/8 |
| MS_14387 | 3 | 3 | 0 | 7.33/8 |

## C. Layer coverage summary

Layer coverage counts how many of the 8 layers contain information about a determinant. The layers are heterogeneous and not equivalent: a determinant with 6/8 layers is **not** necessarily biologically more important than one with 4/8.

| Strain | Determinants | Mean layer coverage |
|---|---|---|
| B36 | 8 | 6.12/8 |
| MS_14384 | 1 | 8.00/8 |
| MS_14386 | 13 | 6.54/8 |
| MS_14387 | 3 | 7.33/8 |

## D. Evidence-level summary

Evidence is classified as **directly observed** (measured abundance / represented locus), **graph-derived** (edges, hyperedge membership, bounded paths), **annotation-mediated** (KEGG/COG/GO/curated mechanism), or **inferred** (metabolite relationships via pathway/enzyme information). Missing layers are reported as **unavailable**, never fabricated as zeros.

| Strain | Direct | Graph-derived | Annotation-mediated | Inferred | Missing layers (marker-level) |
|---|---|---|---|---|---|
| B36 | 340 | 3284 | 59 | 327 | 0 |
| MS_14384 | 13 | 414 | 10 | 11 | 0 |
| MS_14386 | 678 | 5346 | 100 | 654 | 0 |
| MS_14387 | 171 | 1239 | 28 | 165 | 0 |

**Important:** annotation-mediated evidence (KEGG membership, GO annotation, COG assignment) indicates association in a reference database, not observed activity in the strain under study. STRING PPI edges are database/network evidence, not strain-specific experimental validation.

## E. RNA/protein concordance summary

Concordance compares the sign of `RNA log2(Sera/RPMI)` (derived) with the sign of `Protein_logFC`; both are log scales. It is a descriptive observation, not proof of post-transcriptional regulation.

| Strain | Concordant | Discordant | Neutral | Unavailable |
|---|---|---|---|---|
| B36 | 4 | 1 | 0 | 3 |
| MS_14384 | 0 | 1 | 0 | 0 |
| MS_14386 | 9 | 2 | 0 | 2 |
| MS_14387 | 0 | 3 | 0 | 0 |

## F. Hyperedge summary

| Hyperedge type | Determinant-memberships across strains |
|---|---|
| amr_mechanism | 25 |
| cog_category | 24 |
| multi_omics_triple | 20 |
| go_term | 18 |
| kegg_pathway | 13 |
| metabolic_pathway | 11 |

Hyperedge membership is graph-derived context, not evidence of functional activity.

## G. Top candidate AMR-associated elements

Candidates are **non-AMR** genes, proteins, or metabolites that repeatedly occur in the molecular context of known AMR determinants. They are ranked by transparent descriptive criteria (number of distinct AMR determinants connected, number of strains, minimum path length, direct observations) and are **hypotheses for follow-up, not confirmed resistance determinants**. Candidates connected to >= 75% of the maximum determinant count reached by any candidate are labelled **hub candidates** (e.g. ubiquitous pathway-shared metabolites) and are ordered after the specific candidates.

Candidates also carry a one-sided **over-representation test**: the probability that another entity of the same type would connect to >= that many determinants purely by chance (binomial tail, `p_value`, corrected across all candidates with **Benjamini-Hochberg FDR** into `q_value`). A `q_value <= 0.05` marks a candidate as **significantly enriched** among the entities of its type that touch AMR context. Significance is reported alongside the descriptive ranking and does not change the ordering; it is distinct from the `hub` flag, which measures non-specificity.

| Rank | Candidate | Name | Type | #AMR determinants | #Strains | Min path | Hub ? | p-value | q-value (FDR) | Significant | Evidence types |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | CHEBI:30797 | L-Malic acid | metabolite | 9 | 4 | 2 |  | 0.366 | 0.855 |  | directly_observed; inferred |
| 2 | CHEBI:15741 | Succinic acid | metabolite | 9 | 4 | 2 |  | 0.366 | 0.855 |  | directly_observed; inferred |
| 3 | CHEBI:4170 | Glucose 6-phosphate | metabolite | 9 | 4 | 2 |  | 0.366 | 0.855 |  | directly_observed; inferred |
| 4 | CHEBI:17053 | L-Aspartic acid | metabolite | 9 | 4 | 2 |  | 0.366 | 0.855 |  | directly_observed; inferred |
| 5 | CHEBI:18012 | Fumaric acid | metabolite | 9 | 4 | 2 |  | 0.366 | 0.855 |  | directly_observed; inferred |
| 6 | CHEBI:15344 | Acetoacetic acid | metabolite | 9 | 4 | 2 |  | 0.366 | 0.855 |  | directly_observed; inferred |
| 7 | CHEBI:30769 | Citric acid | metabolite | 9 | 4 | 2 |  | 0.366 | 0.855 |  | directly_observed; inferred |
| 8 | CHEBI:17634 | D-Glucose | metabolite | 9 | 4 | 2 |  | 0.366 | 0.855 |  | directly_observed; inferred |
| 9 | CHEBI:16027 | Adenosine monophosphate | metabolite | 9 | 4 | 2 |  | 0.366 | 0.855 |  | directly_observed; inferred |
| 10 | CHEBI:17659 | Uridine 5'-diphosphate | metabolite | 9 | 4 | 2 |  | 0.366 | 0.855 |  | directly_observed; inferred |
| 11 | WP_000074527.1 | alaA | protein | 9 | 4 | 2 |  | 4.1e-05 | 0.024 | yes | annotation_mediated |
| 12 | CHEBI:16467 | L-Arginine | metabolite | 9 | 3 | 2 |  | 0.366 | 0.855 |  | directly_observed; inferred |
| 13 | EW035_RS00170 | ilvN | gene | 9 | 1 | 2 |  | 0.000183 | 0.024 | yes | annotation_mediated |
| 14 | EW035_RS00460 | gpmM | gene | 9 | 1 | 2 |  | 0.000183 | 0.024 | yes | annotation_mediated |
| 15 | EW035_RS00680 | avtA | gene | 9 | 1 | 2 |  | 0.000183 | 0.024 | yes | annotation_mediated |
| 16 | EW035_RS01515 | gntK | gene | 9 | 1 | 2 |  | 0.000183 | 0.024 | yes | annotation_mediated |
| 17 | EW035_RS01530 | asd | gene | 9 | 1 | 2 |  | 0.000183 | 0.024 | yes | annotation_mediated |
| 18 | EW035_RS01700 | pckA | gene | 9 | 1 | 2 |  | 0.000183 | 0.024 | yes | annotation_mediated |
| 19 | EW035_RS01770 | aroB | gene | 9 | 1 | 2 |  | 0.000183 | 0.024 | yes | annotation_mediated |
| 20 | EW035_RS01785 | rpe | gene | 9 | 1 | 2 |  | 0.000183 | 0.024 | yes | annotation_mediated |

## H. Cross-strain observations

Determinants present in >= 2 strains (differences are descriptive, not evidence of biological causation):

| Marker | Class | B36 coverage | MS_14384 coverage | MS_14386 coverage | MS_14387 coverage |
|---|---|---|---|---|---|
| blaTEM-1b | beta_lactam |  |  | 8 | 8 |
| mph(A) | macrolide | 6 |  | 6 |  |
| strAB | aminoglycoside |  |  | 6 | 6 |
| sul2 | sulfonamide |  |  | 8 | 8 |

## I. Missing-data and data-quality notes

No layer is unavailable for any strain; all 8 layers exist for every strain in this analysis.

Data-correctness checks:

- **B36** KEGG consistency: healthy (472 distinct pathway annotations)
- B36: WP versioning -- all protein ids are versioned (WP_....N); no normalization needed
- **MS_14384** KEGG consistency: healthy (476 distinct pathway annotations)
- MS_14384: WP versioning -- all protein ids are versioned (WP_....N); no normalization needed
- **MS_14386** KEGG consistency: healthy (480 distinct pathway annotations)
- MS_14386: WP versioning -- all protein ids are versioned (WP_....N); no normalization needed
- **MS_14387** KEGG consistency: healthy (472 distinct pathway annotations)
- MS_14387: WP versioning -- all protein ids are versioned (WP_....N); no normalization needed

Where a saved edge table has duplicate columns, it is flagged and not silently trusted. Suspicious KEGG results are reported, not silently normalized away.

Recovered UniProt annotations:

- **20** determinant rows have GO terms recovered by per-ID UniProt lookup (`reports/amr_uniprot_recovery.csv`). The upstream pipeline maps RefSeq -> UniProt in ~1000-ID batches but only keeps the first 500 matching UniProt entries per request; later IDs in a batch are silently dropped, which is why some multi-omics-quantified AMR determinants previously had no GO annotation. Re-querying each WP accession individually (same search, same candidate resolution and parsing) recovers those annotations. The recovered terms are real UniProt GO annotations and are labelled annotation-mediated with their own evidence source.

Recovered eggNOG KEGG pathways:

- AMR determinants now carry KEGG pathway memberships merged from the eggNOG-mapper `KEGG_Pathway` column (regenerated TSVs via `scripts/convert_emapper_to_tsv.py`, applied to the saved graphs by `scripts/enrich_amr_graph.py`). This closes the previous coverage gap where `eco`-bridged UniProt cross-references skipped plasmid-borne genes.

Known reference-coverage limitations (KEGG / PPI):
- **KEGG for acquired AMR determinants.** The `eco`-bridged UniProt KEGG cross-references do not cover plasmid-borne genes, but the eggNOG-mapper output now carries a `KEGG_Pathway` column (see `scripts/convert_emapper_to_tsv.py`), and those memberships are merged into the graph (`scripts/enrich_amr_graph.py`). KEGG is still genuinely absent where eggNOG assigns no pathway (e.g. many efflux/modification determinants with only a `-` KEGG_Pathway).
- **PPI = 0 for AMR determinants.** STRING is queried against E. coli K-12/UPEC reference proteins; the clinical-strain WP accessions are not hosted by STRING, and the acquired resistance genes have no K-12 orthologs. Zero STRING edges for a determinant reflects reference coverage, not an absence of physical interactions.

## J. Interpretation guidelines and limitations

- This is exploratory **molecular-context characterization**, not AMR prediction and not causal inference.
- Graph connectivity does not imply causality. A path between a determinant and an element means a *graph-derived association*, nothing more.
- KEGG membership does **not** mean the pathway is active; GO annotation does **not** mean functional activity was observed; STRING PPI is **not** strain-specific experimental validation.
- A metabolite sharing a pathway does **not** mean it directly interacts with the AMR gene; such links are labeled *inferred*.
- A determinant with higher layer coverage is **not** automatically more important; coverage is a descriptive integration measure.
- Candidate AMR-associated elements are **hypotheses**, not confirmed resistance determinants.
- Discordant RNA/protein signs are descriptive observations and are **not** interpreted as proof of post-transcriptional regulation.
- Missing measurements are never converted to zeros; they are reported as *not observed* (layer exists, no value) or *unavailable* (layer data do not exist).
