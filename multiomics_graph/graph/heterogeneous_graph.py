"""
heterogeneous_graph.py
Build a heterogeneous multi-omics graph using PyTorch Geometric HeteroData.

Node types:
  - gene:    Genomic loci with transcriptomic features
  - protein: Protein products with proteomic features

Edge types (source → relation → target):
  - gene → genomic_proximity → gene
  - gene → transcriptional_correlation → gene
  - protein → abundance_correlation → protein
  - protein → ppi → protein            (STRING)
  - gene → encodes → protein           (central dogma)
  - protein → annotated_by → kegg_pathway
  - protein → annotated_by → go_term
  - gene → annotated_by → cog_category

The HeteroData object is directly compatible with
PyTorch Geometric's heterogeneous GNN operators.
"""

import pandas as pd
import numpy as np
import torch
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def _isnan(v) -> bool:
    """True for NaN values without raising on non-numeric input."""
    try:
        return bool(pd.isna(v))
    except Exception:
        return False

try:
    from torch_geometric.data import HeteroData
except ImportError:
    # Fallback if PyG not installed: use dict-based simulation
    HeteroData = dict
    print("  [PyG] torch_geometric not available — using dict fallback")


class HeterogeneousGraphBuilder:
    """
    Build a heterogeneous multi-omics graph.

    Central dogma architecture:
        gene ──encodes──> protein

    Each node type has omics-derived features.
    Each edge type represents a biological relationship.

    Parameters
    ----------
    add_self_loops : bool
        Add self-loops for all node types (helps GNN training)
    correlation_threshold : float
        Minimum |Pearson r| for correlation edges
    """

    def __init__(self, add_self_loops: bool = True,
                 correlation_threshold: float = 0.7,
                 correlation_top_k: int = 10):
        self.add_self_loops = add_self_loops
        self.correlation_threshold = correlation_threshold
        self.correlation_top_k = correlation_top_k
        self.data = None
        self._node_maps = {}  # {node_type: {node_id: index}}

    @staticmethod
    def _correlation_edges(corr_mat: np.ndarray, threshold: float,
                           top_k: int):
        """Upper-triangle pairs passing |r| >= threshold (+ mutual top-k cap).

        The plain |r| >= threshold rule is far too dense on real RNA-seq
        (roughly a third of all gene pairs exceed 0.7 with 12 replicates), so
        the mutual top-k criterion additionally keeps an edge only when the
        pair is among the ``top_k`` most strongly correlated partners of both
        endpoints. This yields a sparse, symmetric co-expression graph of the
        same order as the other edge families.

        Returns
        -------
        (i_idx, j_idx, values) : tuple of arrays
        """
        abs_c = np.abs(corr_mat)
        n = corr_mat.shape[0]
        mask = (abs_c >= threshold) & np.triu(np.ones((n, n), dtype=bool), k=1)
        if top_k and top_k > 0:
            order = np.argsort(-abs_c, axis=1)
            mutual = np.zeros((n, n), dtype=bool)
            for i in range(n):
                for j in order[i, :top_k]:
                    if j != i:
                        mutual[i, j] = True
                        mutual[j, i] = True
            mask &= mutual
        i_idx, j_idx = np.where(mask)
        values = [[float(corr_mat[i, j])] for i, j in zip(i_idx, j_idx)]
        return i_idx, j_idx, values

    def _map_nodes(self, ids: List[str], node_type: str) -> torch.Tensor:
        """Create and cache a mapping from node IDs to consecutive indices."""
        if node_type not in self._node_maps:
            self._node_maps[node_type] = {
                node_id: idx for idx, node_id in enumerate(ids)
            }
        else:
            existing = self._node_maps[node_type]
            next_idx = len(existing)
            for nid in ids:
                if nid not in existing:
                    existing[nid] = next_idx
                    next_idx += 1
        return torch.arange(len(self._node_maps[node_type]))

    def _attach_mask(self, node_type: str, mask_list: List[List[bool]]):
        """Store a boolean feature mask alongside a node type's ``x`` tensor.

        A ``True`` entry marks a feature slot that was NOT measured for that
        node (e.g. the three RNA slots of a genome-only gene) and was therefore
        zero-filled in ``x``. Downstream models can use this mask to avoid
        treating such zeros as observed values — preserving the "honest missing
        data" contract end to end without fabricating measurements.
        """
        mask_t = torch.tensor(mask_list, dtype=torch.bool)
        if HeteroData is not dict:
            self.data[node_type].x_mask = mask_t
        else:
            self.data.setdefault(node_type, {})['x_mask'] = mask_t


    def _build_edge_index(self, pairs: List[Tuple[str, str]],
                           src_type: str, dst_type: str) -> torch.Tensor:
        """
        Convert node-pair list to edge index tensor [2, n_edges].

        Parameters
        ----------
        pairs : list of (src_id, dst_id)
        src_type : str
        dst_type : str

        Returns
        -------
        torch.LongTensor of shape [2, n_edges]
        """
        src_map = self._node_maps.get(src_type, {})
        dst_map = self._node_maps.get(dst_type, {})

        edges = []
        for s, d in pairs:
            if s in src_map and d in dst_map:
                edges.append([src_map[s], dst_map[d]])

        if not edges:
            return torch.zeros((2, 0), dtype=torch.long)

        return torch.tensor(edges, dtype=torch.long).t()

    def add_gene_nodes(
        self, gene_ids: List[str],
        rna_data: pd.DataFrame,
        genome_data: pd.DataFrame,
    ) -> 'HeterogeneousGraphBuilder':
        """
        Add gene nodes with transcriptomic and genomic features.

        Features:
          - RNA_RPMI, RNA_Sera, RNA_logFC (transcriptomic)
          - start, end, CDS_length (genomic)
          - strand (one-hot: + → [1,0], - → [0,1])

        Parameters
        ----------
        gene_ids : list of str
            Gene identifiers to include as nodes
        rna_data : pd.DataFrame
            RNA expression data (indexed by GeneID)
        genome_data : pd.DataFrame
            Genome annotation data (indexed by GeneID)
        """
        self._map_nodes(gene_ids, 'gene')

        feat_list = []
        mask_list = []
        for gid in gene_ids:
            rna = rna_data.loc[gid] if gid in rna_data.index else pd.Series()
            if isinstance(rna, pd.DataFrame):
                rna = rna.iloc[0]
            gen = genome_data.loc[gid] if gid in genome_data.index else pd.Series()
            if isinstance(gen, pd.DataFrame):
                gen = gen.iloc[0]

            rna_missing = (_isnan(rna.get('RNA_RPMI', 0))
                           or _isnan(rna.get('RNA_Sera', 0))
                           or _isnan(rna.get('RNA_logFC', 0)))
            rpmi = float(rna.get('RNA_RPMI', 0)) if not _isnan(rna.get('RNA_RPMI', 0)) else 0.0
            sera = float(rna.get('RNA_Sera', 0)) if not _isnan(rna.get('RNA_Sera', 0)) else 0.0
            logfc = float(rna.get('RNA_logFC', 0)) if not _isnan(rna.get('RNA_logFC', 0)) else 0.0
            start = float(gen.get('start', 0)) / 1e6
            length = float(gen.get('CDS_length', 0)) / 1e3
            strand = gen.get('strand', '+')
            strand_oh = [1.0, 0.0] if strand == '+' else [0.0, 1.0]

            feat_list.append([rpmi, sera, logfc, start, length] + strand_oh)
            # Mask the three transcriptomic slots when RNA is absent for this
            # gene (genome-only genes); the genomic slots are always present.
            mask_list.append([rna_missing, rna_missing, rna_missing,
                              False, False, False, False])

        if HeteroData is not dict:
            self.data['gene'].x = torch.tensor(feat_list, dtype=torch.float)
            self._attach_mask('gene', mask_list)
        else:
            self.data = self.data or {}
            self.data.setdefault('gene', {})['x'] = torch.tensor(
                feat_list, dtype=torch.float
            )
            self._attach_mask('gene', mask_list)

        return self

    def add_protein_nodes(
        self, protein_ids: List[str],
        protein_data: pd.DataFrame,
    ) -> 'HeterogeneousGraphBuilder':
        """
        Add protein nodes with proteomic features.

        Features:
          - Protein_RPMI, Protein_Sera, Protein_logFC

        Parameters
        ----------
        protein_ids : list of str
        protein_data : pd.DataFrame
            Protein abundance data (indexed by ProteinID)
        """
        self._map_nodes(protein_ids, 'protein')

        feat_list = []
        mask_list = []
        for pid in protein_ids:
            prot = protein_data.loc[pid] if pid in protein_data.index else pd.Series()
            if isinstance(prot, pd.DataFrame):
                prot = prot.iloc[0]
            missing = (_isnan(prot.get('Protein_RPMI', 0))
                       or _isnan(prot.get('Protein_Sera', 0))
                       or _isnan(prot.get('Protein_logFC', 0)))
            feat_list.append([
                float(prot.get('Protein_RPMI', 0)),
                float(prot.get('Protein_Sera', 0)),
                float(prot.get('Protein_logFC', 0)),
            ])
            mask_list.append([missing, missing, missing])

        if HeteroData is not dict:
            self.data['protein'].x = torch.tensor(feat_list, dtype=torch.float)
            self._attach_mask('protein', mask_list)
        else:
            self.data.setdefault('protein', {})['x'] = torch.tensor(
                feat_list, dtype=torch.float
            )
            self._attach_mask('protein', mask_list)

        return self

    def add_genomic_proximity_edges(
        self, genomic_edges: pd.DataFrame,
    ) -> 'HeterogeneousGraphBuilder':
        """
        Add gene ↔ gene edges based on genomic proximity.

        Edge features: genomic_distance (inverse normalised)

        Edge type: (gene, genomic_proximity, gene)

        Parameters
        ----------
        genomic_edges : pd.DataFrame
            Columns: GeneID_A, GeneID_B, distance, weight
        """
        pairs = []
        edge_feats = []
        for _, row in genomic_edges.iterrows():
            pairs.append((row['GeneID_A'], row['GeneID_B']))
            edge_feats.append([float(row.get('weight', 1.0))])

        edge_index = self._build_edge_index(pairs, 'gene', 'gene')
        edge_attr = torch.tensor(edge_feats, dtype=torch.float) if edge_feats else None

        key = ('gene', 'genomic_proximity', 'gene')
        if HeteroData is not dict:
            self.data[key].edge_index = edge_index
            if edge_attr is not None:
                self.data[key].edge_attr = edge_attr
        else:
            self.data[key] = {'edge_index': edge_index}
            if edge_attr is not None:
                self.data[key]['edge_attr'] = edge_attr

        return self

    def add_transcriptional_correlation_edges(
        self, gene_ids: List[str],
        rna_data: pd.DataFrame,
    ) -> 'HeterogeneousGraphBuilder':
        """
        Add gene ↔ gene edges based on co-expression correlation.

        Edge created when |Pearson r| > threshold across RNA replicates.

        Edge type: (gene, transcriptional_correlation, gene)

        Parameters
        ----------
        gene_ids : list of str
        rna_data : pd.DataFrame
            Must contain per-replicate columns (not just means)
        """
        # Find replicate columns (numeric, not GeneID/GeneName)
        exclude = {'GeneID', 'GeneName', 'RNA_RPMI', 'RNA_Sera',
                    'RNA_logFC', 'RNA_Regulation'}
        rep_cols = [
            c for c in rna_data.columns
            if c not in exclude and pd.api.types.is_numeric_dtype(rna_data[c])
        ]

        if len(rep_cols) < 3:
            print("  [HeteroGraph] Too few replicate columns for correlation")
            return self

        # Subset replicate matrix for aligned genes (vectorised correlation)
        gene_set = [g for g in gene_ids if g in rna_data.index]
        rep_mat = rna_data.loc[gene_set, rep_cols].values.astype(float)

        # Mean-impute remaining NaNs so corrcoef doesn't fail
        col_mean = np.nanmean(rep_mat, axis=0)
        inds = np.where(np.isnan(rep_mat))
        rep_mat[inds] = np.take(col_mean, inds[1])

        corr_mat = np.corrcoef(rep_mat)
        i_idx, j_idx, corr_vals = self._correlation_edges(
            corr_mat, self.correlation_threshold, self.correlation_top_k)
        pairs = [(gene_set[i], gene_set[j]) for i, j in zip(i_idx, j_idx)]

        edge_index = self._build_edge_index(pairs, 'gene', 'gene')
        edge_attr = torch.tensor(corr_vals, dtype=torch.float) if corr_vals else None

        key = ('gene', 'transcriptional_correlation', 'gene')
        if HeteroData is not dict:
            self.data[key].edge_index = edge_index
            if edge_attr is not None:
                self.data[key].edge_attr = edge_attr
        else:
            self.data[key] = {'edge_index': edge_index}
            if edge_attr is not None:
                self.data[key]['edge_attr'] = edge_attr

        return self

    def add_encodes_edges(
        self, gene_protein_map: Dict[str, str]
    ) -> 'HeterogeneousGraphBuilder':
        """
        Add gene → protein edges representing the central dogma.

        Edge type: (gene, encodes, protein)

        Parameters
        ----------
        gene_protein_map : dict
            {GeneID: ProteinID}
        """
        pairs = []
        for gene_id, prot_id in gene_protein_map.items():
            if prot_id and not pd.isna(prot_id):
                pairs.append((gene_id, prot_id))

        edge_index = self._build_edge_index(pairs, 'gene', 'protein')
        # Edge feature: 1.0 for all (binary relationship)
        edge_attr = torch.ones((edge_index.size(1), 1), dtype=torch.float) \
            if edge_index.size(1) > 0 else None

        key = ('gene', 'encodes', 'protein')
        if HeteroData is not dict:
            self.data[key].edge_index = edge_index
            if edge_attr is not None:
                self.data[key].edge_attr = edge_attr
        else:
            self.data[key] = {'edge_index': edge_index}
            if edge_attr is not None:
                self.data[key]['edge_attr'] = edge_attr

        return self

    def add_protein_ppi_edges(
        self, ppi_edges: pd.DataFrame,
        protein_col_a: str = 'ProteinID_A',
        protein_col_b: str = 'ProteinID_B',
        score_col: str = 'combined_score',
    ) -> 'HeterogeneousGraphBuilder':
        """
        Add protein ↔ protein PPI edges from STRING.

        Edge type: (protein, ppi, protein)

        Parameters
        ----------
        ppi_edges : pd.DataFrame
        protein_col_a, protein_col_b : str
            Column names for protein identifiers
        score_col : str
            Column with interaction confidence score
        """
        pairs = []
        scores = []
        for _, row in ppi_edges.iterrows():
            pairs.append((row[protein_col_a], row[protein_col_b]))
            scores.append([float(row.get(score_col, 0))])

        edge_index = self._build_edge_index(pairs, 'protein', 'protein')
        edge_attr = torch.tensor(scores, dtype=torch.float) if scores else None

        key = ('protein', 'ppi', 'protein')
        if HeteroData is not dict:
            self.data[key].edge_index = edge_index
            if edge_attr is not None:
                self.data[key].edge_attr = edge_attr
        else:
            self.data[key] = {'edge_index': edge_index}
            if edge_attr is not None:
                self.data[key]['edge_attr'] = edge_attr

        return self

    def add_abundance_correlation_edges(
        self, protein_ids: List[str],
        protein_data: pd.DataFrame,
    ) -> 'HeterogeneousGraphBuilder':
        """
        Add protein ↔ protein edges based on abundance correlation.

        Edge created when |Pearson r| > threshold across the per-sample
        replicate columns (per-sample intensities), mirroring the RNA
        co-expression family.

        Edge type: (protein, abundance_correlation, protein)
        """
        # Find replicate columns (numeric, not the summarised means/labels)
        exclude = {'ProteinID', 'Protein', 'GeneSymbol', 'GeneName',
                    'Protein_RPMI', 'Protein_Sera', 'Protein_logFC',
                    'Protein_Regulation'}
        rep_cols = [
            c for c in protein_data.columns
            if c not in exclude and pd.api.types.is_numeric_dtype(protein_data[c])
        ]

        if len(rep_cols) < 3:
            return self

        prot_set = [p for p in protein_ids if p in protein_data.index]
        prot_mat = protein_data.loc[prot_set, rep_cols].values.astype(float)

        col_mean = np.nanmean(prot_mat, axis=0)
        inds = np.where(np.isnan(prot_mat))
        prot_mat[inds] = np.take(col_mean, inds[1])

        corr_mat = np.corrcoef(prot_mat)
        i_idx, j_idx, corr_vals = self._correlation_edges(
            corr_mat, self.correlation_threshold, self.correlation_top_k)
        pairs = [(prot_set[i], prot_set[j]) for i, j in zip(i_idx, j_idx)]

        edge_index = self._build_edge_index(pairs, 'protein', 'protein')
        edge_attr = torch.tensor(corr_vals, dtype=torch.float) if corr_vals else None

        key = ('protein', 'abundance_correlation', 'protein')
        if HeteroData is not dict:
            self.data[key].edge_index = edge_index
            if edge_attr is not None:
                self.data[key].edge_attr = edge_attr
        else:
            self.data[key] = {'edge_index': edge_index}
            if edge_attr is not None:
                self.data[key]['edge_attr'] = edge_attr

        return self

    def add_metabolite_correlation_edges(
        self, metabolite_ids: List[str],
        metabolite_data: pd.DataFrame,
    ) -> 'HeterogeneousGraphBuilder':
        """
        Add metabolite ↔ metabolite edges based on abundance correlation.

        Same rule as the RNA/protein correlation families: |Pearson r| above
        the threshold, capped by the mutual top-k criterion, computed over the
        per-sample abundance columns (MetaboliteID x replicates).

        Edge type: (metabolite, abundance_correlation, metabolite)
        """
        exclude = {'MetaboliteID', 'MetaboliteName', 'Metabolite_RPMI',
                    'Metabolite_Sera', 'Metabolite_logFC',
                    'MassToCharge', 'RetentionTime', 'ChemicalFormula',
                    'Database', 'Instruments', 'n_instruments'}
        rep_cols = [
            c for c in metabolite_data.columns
            if c not in exclude and pd.api.types.is_numeric_dtype(
                metabolite_data[c])
        ]

        if len(rep_cols) < 3:
            return self

        metab_set = [m for m in metabolite_ids if m in metabolite_data.index]
        metab_mat = metabolite_data.loc[metab_set, rep_cols].values.astype(float)

        col_mean = np.nanmean(metab_mat, axis=0)
        inds = np.where(np.isnan(metab_mat))
        metab_mat[inds] = np.take(col_mean, inds[1])

        corr_mat = np.corrcoef(metab_mat)
        i_idx, j_idx, corr_vals = self._correlation_edges(
            corr_mat, self.correlation_threshold, self.correlation_top_k)
        pairs = [(metab_set[i], metab_set[j]) for i, j in zip(i_idx, j_idx)]

        edge_index = self._build_edge_index(pairs, 'metabolite', 'metabolite')
        edge_attr = torch.tensor(corr_vals, dtype=torch.float) if corr_vals else None

        key = ('metabolite', 'abundance_correlation', 'metabolite')
        if HeteroData is not dict:
            self.data[key].edge_index = edge_index
            if edge_attr is not None:
                self.data[key].edge_attr = edge_attr
        else:
            self.data[key] = {'edge_index': edge_index}
            if edge_attr is not None:
                self.data[key]['edge_attr'] = edge_attr

        return self

    def add_functional_edges(
        self, gene_to_pathways: Dict[str, List[str]],
        path_to_name: Dict[str, str],
        node_type: str = 'gene',
        relation: str = 'annotated_by',
    ) -> 'HeterogeneousGraphBuilder':
        """
        Add edges from genes/proteins to functional annotations.

        Edge type: (node_type, relation, annotation_type)

        Where annotation_type is 'kegg_pathway', 'go_term', or 'cog_category'.
        For simplicity, annotation nodes are stored as a single 'annotation'
        node type with a mapping.

        Parameters
        ----------
        gene_to_pathways : dict
            {GeneID: [annotation_id, ...]}
        path_to_name : dict
            {annotation_id: display_name}
        node_type : str
            Source node type ('gene' or 'protein')
        relation : str
            Edge relation name
        """
        # Collect all unique annotation IDs
        all_annotations = set()
        for anns in gene_to_pathways.values():
            all_annotations.update(anns)
        ann_list = sorted(all_annotations)

        if not ann_list:
            return self

        self._map_nodes(ann_list, 'annotation')

        pairs = []
        for gene_id, anns in gene_to_pathways.items():
            for ann in anns:
                if ann in self._node_maps.get('annotation', {}):
                    pairs.append((gene_id, ann))

        edge_index = self._build_edge_index(pairs, node_type, 'annotation')
        if edge_index.size(1) == 0:
            return self

        key = (node_type, relation, 'annotation')
        ann_names = [path_to_name.get(a, a) for a in ann_list]
        if HeteroData is not dict:
            self.data[key].edge_index = edge_index
            existing = getattr(self.data['annotation'], 'name', [])
            self.data['annotation'].name = list(existing) + ann_names
        else:
            self.data[key] = {'edge_index': edge_index}
            existing = self.data.setdefault('annotation', {}).get('name', [])
            self.data['annotation']['name'] = list(existing) + ann_names

        return self

    def add_metabolite_nodes(
        self, metabolite_ids: List[str],
        metabolite_data: pd.DataFrame,
    ) -> 'HeterogeneousGraphBuilder':
        """
        Add metabolite nodes with abundance + chemical features.

        Features:
          - Metabolite_RPMI, Metabolite_Sera, Metabolite_logFC
          - MassToCharge (scaled), RetentionTime (scaled)

        Parameters
        ----------
        metabolite_ids : list of str
        metabolite_data : pd.DataFrame
            Indexed by MetaboliteID
        """
        self._map_nodes(metabolite_ids, 'metabolite')

        feat_list = []
        mask_list = []
        for mid in metabolite_ids:
            m = metabolite_data.loc[mid] if mid in metabolite_data.index else pd.Series()
            if isinstance(m, pd.DataFrame):
                m = m.iloc[0]
            mz = float(m.get('MassToCharge', 0))
            rt = float(m.get('RetentionTime', 0))
            missing = (_isnan(m.get('Metabolite_RPMI', 0))
                       or _isnan(m.get('Metabolite_Sera', 0))
                       or _isnan(m.get('Metabolite_logFC', 0)))
            feat_list.append([
                float(m.get('Metabolite_RPMI', 0)),
                float(m.get('Metabolite_Sera', 0)),
                float(m.get('Metabolite_logFC', 0)),
                float(mz) / 1e3,          # m/z ~ hundreds, scale to ~0.1
                float(rt) / 1e6,          # retention time ~ minutes
            ])
            mask_list.append([missing, missing, missing, False, False])

        if HeteroData is not dict:
            self.data['metabolite'].x = torch.tensor(feat_list, dtype=torch.float)
            self.data['metabolite'].name = list(metabolite_ids)
            self._attach_mask('metabolite', mask_list)
        else:
            self.data.setdefault('metabolite', {})['x'] = torch.tensor(
                feat_list, dtype=torch.float
            )
            self.data['metabolite']['name'] = list(metabolite_ids)
            self._attach_mask('metabolite', mask_list)

        return self

    def add_metabolite_functional_edges(
        self, metabolite_to_func: Dict[str, List[str]],
        func_to_name: Dict[str, str],
        relation: str = 'in_pathway',
    ) -> 'HeterogeneousGraphBuilder':
        """
        Add edges from metabolites to functional/pathway annotations.

        Edge type: (metabolite, relation, annotation)

        Parameters
        ----------
        metabolite_to_func : dict
            {MetaboliteID: [annotation_id, ...]}
        func_to_name : dict
            {annotation_id: display_name}
        relation : str
            Edge relation name (default 'in_pathway')
        """
        all_anns = set()
        for anns in metabolite_to_func.values():
            all_anns.update(anns)
        ann_list = sorted(all_anns)
        if not ann_list:
            return self

        self._map_nodes(ann_list, 'annotation')

        pairs = []
        for mid, anns in metabolite_to_func.items():
            for ann in anns:
                if ann in self._node_maps.get('annotation', {}):
                    pairs.append((mid, ann))

        edge_index = self._build_edge_index(pairs, 'metabolite', 'annotation')
        if edge_index.size(1) == 0:
            return self

        key = ('metabolite', relation, 'annotation')
        ann_names = [func_to_name.get(a, a) for a in ann_list]
        if HeteroData is not dict:
            self.data[key].edge_index = edge_index
            existing = getattr(self.data['annotation'], 'name', [])
            self.data['annotation'].name = list(existing) + ann_names
        else:
            self.data[key] = {'edge_index': edge_index}
            existing = self.data.setdefault('annotation', {}).get('name', [])
            self.data['annotation']['name'] = list(existing) + ann_names

        return self

    def add_metabolite_enzyme_edges(
        self, metabolite_to_enzymes: Dict[str, List[str]],
        enzyme_to_protein: Dict[str, List[str]],
    ) -> 'HeterogeneousGraphBuilder':
        """
        Add edges from metabolites to proteins that catalyse them.

        Edge type: (metabolite, metabolised_by, protein)

        The EC-number bridge: metabolite -> enzymes (EC) -> proteins
        annotated with those EC numbers.

        Parameters
        ----------
        metabolite_to_enzymes : dict
            {MetaboliteID: [EC, ...]}
        enzyme_to_protein : dict
            {EC: [ProteinID, ...]}
        """
        pairs = []
        prots_used = set()
        for mid, ecs in metabolite_to_enzymes.items():
            for ec in ecs:
                for pid in enzyme_to_protein.get(ec, []):
                    pairs.append((mid, pid))
                    prots_used.add(pid)

        if not pairs:
            return self

        edge_index = self._build_edge_index(pairs, 'metabolite', 'protein')
        if edge_index.size(1) == 0:
            return self

        key = ('metabolite', 'metabolised_by', 'protein')
        edge_attr = torch.ones((edge_index.size(1), 1), dtype=torch.float)
        if HeteroData is not dict:
            self.data[key].edge_index = edge_index
            self.data[key].edge_attr = edge_attr
        else:
            self.data[key] = {'edge_index': edge_index, 'edge_attr': edge_attr}

        return self

    def build(
        self, gene_ids: List[str], protein_ids: List[str],
        rna_data: pd.DataFrame, genome_data: pd.DataFrame,
        protein_data: pd.DataFrame,
        genomic_edges: pd.DataFrame,
        gene_protein_map: Dict[str, str],
        ppi_edges: Optional[pd.DataFrame] = None,
        gene_to_ko: Optional[Dict[str, List[str]]] = None,
        ko_to_pathway: Optional[Dict[str, str]] = None,
        gene_to_cog: Optional[Dict[str, List[str]]] = None,
        metabolite_ids: Optional[List[str]] = None,
        metabolite_data: Optional[pd.DataFrame] = None,
        metabolite_replicates: Optional[pd.DataFrame] = None,
        metabolite_to_pathway: Optional[Dict[str, List[str]]] = None,
        pathway_names: Optional[Dict[str, str]] = None,
        metabolite_to_enzymes: Optional[Dict[str, List[str]]] = None,
        enzyme_to_protein: Optional[Dict[str, List[str]]] = None,
        protein_to_go: Optional[Dict[str, List[str]]] = None,
        go_info: Optional[Dict[str, dict]] = None,
        class_to_loci: Optional[Dict[str, List[str]]] = None,
        class_names: Optional[Dict[str, str]] = None,
        ppi_clusters: Optional[List[List[str]]] = None,
        rna_replicates: Optional[pd.DataFrame] = None,
        protein_replicates: Optional[pd.DataFrame] = None,
    ) -> 'HeteroData':
        """
        Build the complete heterogeneous graph.

        Parameters
        ----------
        gene_ids, protein_ids : list of str
        rna_data : pd.DataFrame
        genome_data : pd.DataFrame
        protein_data : pd.DataFrame
        genomic_edges : pd.DataFrame
        gene_protein_map : dict
        ppi_edges : pd.DataFrame, optional
        gene_to_ko : dict, optional
        ko_to_pathway : dict, optional
        gene_to_cog : dict, optional
        rna_replicates : pd.DataFrame, optional
            Gene x replicate-column matrix (raw per-sample counts/log2CPM).
            When provided, the RNA co-expression edges are computed from it;
            otherwise ``rna_data`` (condition means only) is used and the
            co-expression family degrades to no edges.
        protein_replicates : pd.DataFrame, optional
            ProteinID x replicate-column matrix (per-sample intensities).
            Same semantics as ``rna_replicates`` for the protein abundance
            correlation family.

        Returns
        -------
        HeteroData (or dict fallback)
        """
        self.data = HeteroData() if HeteroData is not dict else {}

        # Add nodes with features
        self.add_gene_nodes(gene_ids, rna_data, genome_data)
        self.add_protein_nodes(protein_ids, protein_data)

        # Add core biological edges
        self.add_genomic_proximity_edges(genomic_edges)
        rna_for_corr = rna_replicates if rna_replicates is not None else rna_data
        self.add_transcriptional_correlation_edges(gene_ids, rna_for_corr)
        self.add_encodes_edges(gene_protein_map)
        prot_for_corr = protein_replicates if protein_replicates is not None else protein_data
        self.add_abundance_correlation_edges(protein_ids, prot_for_corr)

        if ppi_edges is not None and not ppi_edges.empty:
            self.add_protein_ppi_edges(ppi_edges)

        if gene_to_ko:
            self.add_functional_edges(
                gene_to_ko, ko_to_pathway or {},
                node_type='gene', relation='in_pathway',
            )

        if gene_to_cog:
            self.add_functional_edges(
                gene_to_cog, {}, node_type='gene', relation='cog_category',
            )

        # GO Terms (from UniProt)
        if protein_to_go:
            go_names = {}
            if go_info:
                for go_id, info in go_info.items():
                    go_names[go_id] = info.get('name', go_id)
            self.add_functional_edges(
                protein_to_go, go_names,
                node_type='protein', relation='annotated_by',
            )

        # AMR mechanisms (curated resistance knowledge)
        if class_to_loci:
            # Invert {class: [locus_tag, ...]} to {locus_tag: [class, ...]}
            gene_to_amr = {}
            for amr_class, loci in class_to_loci.items():
                for locus in loci:
                    gene_to_amr.setdefault(locus, []).append(amr_class)
            
            # Map using display names if available
            amr_names = {}
            if class_names:
                for amr_class, name in class_names.items():
                    amr_names[amr_class] = f"{name} Resistance"
            self.add_functional_edges(
                gene_to_amr, amr_names,
                node_type='gene', relation='associated_with',
            )

        # PPI clusters (Louvain communities)
        if ppi_clusters:
            protein_to_cluster = {}
            cluster_names = {}
            for i, cluster in enumerate(ppi_clusters):
                cluster_id = f"Cluster_{i}"
                cluster_names[cluster_id] = f"PPI Louvain Cluster {i} (size {len(cluster)})"
                for p in cluster:
                    protein_to_cluster.setdefault(p, []).append(cluster_id)
            self.add_functional_edges(
                protein_to_cluster, cluster_names,
                node_type='protein', relation='member_of_cluster',
            )

        # Metabolomics layer
        if metabolite_ids and metabolite_data is not None:
            self.add_metabolite_nodes(metabolite_ids, metabolite_data)
            metab_for_corr = (metabolite_replicates
                              if metabolite_replicates is not None
                              else metabolite_data)
            self.add_metabolite_correlation_edges(metabolite_ids, metab_for_corr)
            if metabolite_to_pathway:
                self.add_metabolite_functional_edges(
                    metabolite_to_pathway, pathway_names or {},
                    relation='in_pathway',
                )
            if metabolite_to_enzymes and enzyme_to_protein:
                self.add_metabolite_enzyme_edges(
                    metabolite_to_enzymes, enzyme_to_protein,
                )

        return self.data

    def save(self, output_dir: str):
        """Save HeteroData, edge CSVs, and metadata to disk."""
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        if hasattr(self.data, 'store'):
            # PyG ≥ 2.0
            torch.save(self.data, out_dir / 'heterodata.pt')
        else:
            torch.save(self.data, out_dir / 'heterodata.pt')

        # Save edge CSVs for visualization by extracting from HeteroData
        idx_to_id = {
            ntype: {idx: nid for nid, idx in mapping.items()}
            for ntype, mapping in self._node_maps.items()
        }

        # Handle both dict fallback and PyG HeteroData
        if isinstance(self.data, dict):
            data_keys = [k for k in self.data if isinstance(k, tuple) and len(k) == 3]
        else:
            data_keys = list(self.data.edge_types) if hasattr(self.data, 'edge_types') else []

        for key in data_keys:
            if not isinstance(key, tuple) or len(key) != 3:
                continue
            src_type, rel, dst_type = key
            store = self.data[key]
            ei = store.get('edge_index') if isinstance(store, dict) else \
                 getattr(store, 'edge_index', None)
            if ei is None:
                continue
            n_edges = ei.size(1) if ei.dim() > 1 else 0
            if n_edges == 0:
                continue
            src_ids = [idx_to_id.get(src_type, {}).get(int(i), f"n{i}") for i in ei[0]]
            dst_ids = [idx_to_id.get(dst_type, {}).get(int(i), f"n{i}") for i in ei[1]]
            src_name = f'{src_type}_id'
            dst_name = f'{dst_type}_id'
            if src_name == dst_name:
                # Same-type edges (e.g., protein--ppi-->protein): keep both
                # columns distinct so neither side is silently dropped.
                dst_name = f'{dst_type}_target'
            csv_path = out_dir / f'edges_{src_type}_{rel}_{dst_type}.csv'
            pd.DataFrame({src_name: src_ids, dst_name: dst_ids, 'relation': rel}).to_csv(csv_path, index=False)
            print(f"  [Saved] {csv_path.name} ({n_edges} edges)")

        # Save node and edge type metadata
        meta_file = out_dir / 'heterodata_metadata.txt'
        with open(meta_file, 'w') as f:
            f.write(self._metadata_str())

        print(f"  [Saved] HeteroData -> {out_dir / 'heterodata.pt'}")

    def _metadata_str(self) -> str:
        """Generate human-readable metadata string."""
        lines = ["=== Heterogeneous Graph Metadata ===\n"]

        # Node types
        if hasattr(self.data, 'node_types'):
            lines.append(f"\nNode types ({len(self.data.node_types)}):")
            for nt in self.data.node_types:
                store = self.data[nt]
                x = store.get('x', None)
                if hasattr(x, 'size'):
                    n = x.size(0)
                    d = x.size(1) if x.dim() > 1 else 0
                else:
                    # Annotation nodes carry no numeric feature tensor; their
                    # identity is the semantic `name` list. Count them from it.
                    names = store.get('name', None)
                    n = len(names) if names is not None else 0
                    d = 0
                lines.append(f"  {nt}: {n} nodes, {d} features")
                # Report missing-feature masks (honest missing-data contract).
                x_mask = store.x_mask if hasattr(store, 'x_mask') else store.get('x_mask') if isinstance(store, dict) else None
                if x_mask is not None:
                    n_masked = int(x_mask.any(dim=1).sum()) if hasattr(x_mask, 'dim') else 0
                    lines.append(f"       {n_masked} nodes with masked (unobserved) features")
        elif isinstance(self.data, dict):
            lines.append(f"\nNode types ({len(self.data)}):")
            for nt in self.data:
                if isinstance(self.data[nt], dict) and 'x' in self.data[nt]:
                    x = self.data[nt]['x']
                    n = x.size(0)
                    d = x.size(1) if x.dim() > 1 else 0
                    lines.append(f"  {nt}: {n} nodes, {d} features")

        # Edge types
        if hasattr(self.data, 'edge_types'):
            edge_types = list(self.data.edge_types)
        elif isinstance(self.data, dict):
            edge_types = [k for k in self.data if isinstance(k, tuple) and len(k) == 3]
        else:
            edge_types = []

        lines.append(f"\nEdge types ({len(edge_types)}):")
        for et in edge_types:
            store = self.data[et]
            ei = store.get('edge_index', torch.zeros((2, 0))) if isinstance(store, dict) else getattr(store, 'edge_index', torch.zeros((2, 0)))
            n = ei.size(1) if hasattr(ei, 'size') else 0
            lines.append(f"  {et[0]} --{et[1]}--> {et[2]}: {n} edges")

        return '\n'.join(lines)

    def summary(self) -> dict:
        """Return graph summary statistics."""
        if self.data is None:
            return {'total_nodes': 0, 'total_edges': 0}
        n_nodes = 0; n_edges = 0
        if hasattr(self.data, 'node_types'):
            for nt in self.data.node_types:
                store = self.data[nt]
                x = store.get('x') if isinstance(store, dict) else getattr(store, 'x', None)
                if x is not None:
                    n_nodes += x.size(0)
        if hasattr(self.data, 'edge_types'):
            for et in self.data.edge_types:
                store = self.data[et]
                ei = store.get('edge_index') if isinstance(store, dict) else getattr(store, 'edge_index', None)
                if ei is not None:
                    n_edges += ei.size(1)
        return {
            'total_nodes': n_nodes,
            'total_edges': n_edges,
        }
