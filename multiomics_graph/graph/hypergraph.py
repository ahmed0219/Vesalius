"""
hypergraph.py
Build biological-module hypergraphs from meta-dimensional bacterial data.

A hyperedge represents a biological module that connects multiple
biological entities simultaneously:

  - KEGG pathways: all genes in glycolysis share a hyperedge
  - GO biological processes: genes involved in DNA replication
  - Protein complexes: physically interacting proteins
  - COG categories: genes with shared evolutionary function
  - Multi-omics triples: gene + RNA + protein per locus

The hypergraph is represented as an incidence matrix H where:
  H[i, j] = 1 if node i belongs to hyperedge j

This representation is compatible with Hypergraph Neural Networks (HGNN)
and can be converted to a bipartite graph for message passing.
"""

import pandas as pd
import numpy as np
import torch
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict


class BiologicalHypergraphBuilder:
    """
    Build hypergraphs from biological modules.

    Hyperedge types:
      1. multi_omics_triple — (Gene, RNA, Protein) per aligned locus
      2. kegg_pathway      — genes in same KEGG pathway
      3. cog_category      — genes with same COG functional category
      4. go_term           — proteins sharing a GO term
      5. ppi_cluster       — Louvain communities from STRING PPI network

    Attributes
    ----------
    node_list : list
        All nodes (gene + protein identifiers) across all hyperedges
    hyperedge_list : list of dict
        Each: {hyperedge_id, type, name, nodes, features}
    H : np.ndarray
        Incidence matrix (n_nodes × n_hyperedges)
    node_to_idx : dict
        {node_id: row_index in H}
    """

    def __init__(self):
        self.node_list = []
        self.node_to_idx = {}
        self.hyperedges = []
        self.H = None
        self.hyperedge_type_map = defaultdict(list)

    def _register_node(self, node_id: str) -> int:
        """Register a node and return its index."""
        if node_id not in self.node_to_idx:
            self.node_to_idx[node_id] = len(self.node_list)
            self.node_list.append(node_id)
        return self.node_to_idx[node_id]

    def _add_hyperedge(self, nodes: List[str],
                        hyperedge_type: str,
                        name: str = '',
                        features: Optional[dict] = None):
        """
        Add a hyperedge connecting the given nodes.

        Parameters
        ----------
        nodes : list of str
            Node identifiers in this hyperedge
        hyperedge_type : str
            Type identifier (e.g., 'kegg_pathway', 'multi_omics_triple')
        name : str
            Human-readable name
        features : dict, optional
            Hyperedge-level features
        """
        node_indices = [self._register_node(n) for n in nodes if n]
        if len(node_indices) < 2:
            return  # Need at least 2 nodes for a meaningful hyperedge

        he_id = len(self.hyperedges)
        self.hyperedges.append({
            'hyperedge_id': he_id,
            'type': hyperedge_type,
            'name': name,
            'nodes': nodes,
            'node_indices': node_indices,
            'features': features or {},
        })
        self.hyperedge_type_map[hyperedge_type].append(he_id)

    def build_incidence_matrix(self) -> np.ndarray:
        """
        Build the binary incidence matrix H.

        H[i, j] = 1 if node i is in hyperedge j

        Returns
        -------
        np.ndarray of shape (n_nodes, n_hyperedges)
        """
        n_nodes = len(self.node_list)
        n_hyperedges = len(self.hyperedges)
        H = np.zeros((n_nodes, n_hyperedges), dtype=np.float32)

        for j, he in enumerate(self.hyperedges):
            for i in he['node_indices']:
                H[i, j] = 1.0

        self.H = H
        return H

    def build_hypergraph_laplacian(self) -> torch.Tensor:
        """
        Build the hypergraph Laplacian: D_v^{-1/2} H D_e^{-1} H^T D_v^{-1/2}

        This is the core propagation operator for hypergraph convolution.

        Returns
        -------
        torch.Tensor of shape (n_nodes, n_nodes)
        """
        if self.H is None:
            self.build_incidence_matrix()

        n_nodes, n_hyperedges = self.H.shape
        d_v = np.maximum(np.sum(self.H, axis=1), 1e-10)
        d_e = np.maximum(np.sum(self.H, axis=0), 1e-10)

        D_v_inv_sqrt = np.diag(1.0 / np.sqrt(d_v))
        D_e_inv = np.diag(1.0 / d_e)

        G = D_v_inv_sqrt @ self.H @ D_e_inv @ self.H.T @ D_v_inv_sqrt
        return torch.FloatTensor(G)

    # ------------------------------------------------------------------
    # Hyperedge builders
    # ------------------------------------------------------------------

    def add_multi_omics_triples(
        self, aligned_data: pd.DataFrame,
        gene_col: str = 'GeneID',
        protein_col: str = 'ProteinID',
    ):
        """
        Create hyperedges from aligned multi-omics data.

        Each hyperedge connects:
          gene_node -- rna_node -- protein_node

        representing the central dogma flow for one genomic locus.

        Hyperedge features: RNA_logFC, Protein_logFC, agreement_score

        Parameters
        ----------
        aligned_data : pd.DataFrame
            Multi-omics aligned table
        gene_col, protein_col : str
        """
        for _, row in aligned_data.iterrows():
            gene_id = row[gene_col]
            rna_id = f"{gene_id}_RNA"
            prot_id = row.get(protein_col, '')

            nodes = [gene_id, rna_id]
            if prot_id and not pd.isna(prot_id):
                nodes.append(prot_id)

            # RNA may be absent in 2-layer (genome + proteome) strains. Detect
            # missing RNA and mask the concordance signal accordingly so NaNs
            # do not propagate into downstream feature/agreement consumers.
            rna_val = row.get('RNA_logFC', 0)
            if pd.isna(rna_val):
                rna_val = 0.0
                has_rna = False
            else:
                rna_val = float(rna_val)
                has_rna = True

            features = {
                'RNA_logFC': rna_val,
                'Protein_logFC': float(row.get('Protein_logFC', 0)),
                'has_rna': has_rna,
            }
            # Agreement is only meaningful when both layers carry signal; leave
            # it as None (not a pseudo-value) for RNA-less 2-layer strains.
            features['agreement_score'] = (
                features['RNA_logFC'] * features['Protein_logFC']
                if has_rna else None
            )

            self._add_hyperedge(
                nodes=nodes,
                hyperedge_type='multi_omics_triple',
                name=f"Triple_{gene_id}",
                features=features,
            )

    def add_kegg_pathway_hyperedges(
        self, pathway_membership: Dict[str, List[str]],
        pathway_names: Optional[Dict[str, str]] = None,
    ):
        """
        Create hyperedges from KEGG pathway membership.

        Each KEGG pathway becomes a hyperedge containing all
        bacterial genes that participate in that pathway.

        Parameters
        ----------
        pathway_membership : dict
            {pathway_id: [list of GeneIDs]}
        pathway_names : dict, optional
            {pathway_id: display_name}
        """
        for path_id, gene_ids in pathway_membership.items():
            name = (pathway_names or {}).get(path_id, path_id)
            # Add both gene and protein nodes if available
            nodes = list(gene_ids)
            self._add_hyperedge(
                nodes=nodes,
                hyperedge_type='kegg_pathway',
                name=name,
                features={'pathway_id': path_id},
            )

    def add_cog_hyperedges(
        self, cog_annotations: Dict[str, dict],
        gene_to_protein: Optional[Dict[str, str]] = None,
        cog_categories: Optional[Dict[str, str]] = None,
    ):
        """
        Create hyperedges from COG/eggNOG functional categories.

        Each COG category (e.g., 'J' = Translation) becomes a
        hyperedge containing all genes with that annotation.

        Parameters
        ----------
        cog_annotations : dict
            {GeneID: {COG_classes: [...], ...}}
        gene_to_protein : dict, optional
            {GeneID: ProteinID} to also include protein nodes
        """
        cog_to_genes = defaultdict(list)
        for gene, ann in cog_annotations.items():
            for cog in ann.get('COG_classes', []):
                cog_to_genes[cog].append(gene)

        cat = cog_categories or {}

        for cog, genes in cog_to_genes.items():
            if not genes:
                continue
            nodes = list(genes)
            if gene_to_protein:
                for g in genes:
                    pid = gene_to_protein.get(g, '')
                    if pid:
                        nodes.append(pid)

            name = cat.get(cog, f"COG_{cog}")
            self._add_hyperedge(
                nodes=list(set(nodes)),
                hyperedge_type='cog_category',
                name=f"COG_{cog}: {name}",
                features={'cog': cog, 'description': name},
            )

    def add_amr_hyperedges(self, class_to_loci: Dict[str, List[str]],
                            amr_loci: List[str],
                            gene_to_protein: Optional[Dict[str, str]] = None,
                            class_names: Optional[Dict[str, str]] = None):
        """
        Create hyperedges from AMR mechanism classes.

        Each AMR mechanism (e.g. beta_lactam) becomes a hyperedge containing
        all genes carrying a resistance determinant of that class. Gene nodes
        (locus tags) are always included so a genome-only resistance gene
        remains represented even without RNA/protein quantification.

        Parameters
        ----------
        class_to_loci : dict
            {amr_class: [locus_tag, ...]} from the AMR manifest.
        amr_loci : list
            All AMR locus tags (used to trap/highlight membership).
        gene_to_protein : dict, optional
            {locus_tag: protein_id} to also attach protein nodes.
        class_names : dict, optional
            {amr_class: display name}.
        """
        for amr_class, loci in class_to_loci.items():
            if not loci:
                continue
            nodes = list(loci)
            if gene_to_protein:
                for g in loci:
                    pid = gene_to_protein.get(g, '')
                    if pid:
                        nodes.append(pid)
            name = (class_names or {}).get(amr_class, amr_class)
            # Always attach the AMR_mechanism node so the class is represented
            # even when a genome-only resistance gene has no quantified
            # protein/RNA (biologically valid) and the hyperedge is not empty.
            mechanism_node = f"AMR:{amr_class}"
            nodes.append(mechanism_node)
            self._add_hyperedge(
                nodes=list(set(nodes)),
                hyperedge_type='amr_mechanism',
                name=f"{amr_class}_resistance",
                features={
                    'amr_class': amr_class,
                    'description': name,
                    'n_loci': len(loci),
                },
            )

    def add_go_hyperedges(
        self, protein_to_go: Dict[str, List[str]],
        go_info: Optional[Dict[str, dict]] = None,
    ):
        """
        Create hyperedges from GO term annotations.

        Each GO term becomes a hyperedge with all proteins
        annotated with that term.

        Parameters
        ----------
        protein_to_go : dict
            {ProteinID: [GO:xxxxxx, ...]}
        go_info : dict, optional
            {GO:xxxxxx: {name, aspect}}
        """
        go_to_proteins = defaultdict(list)
        for prot, go_terms in protein_to_go.items():
            for go in go_terms:
                go_to_proteins[go].append(prot)

        for go, proteins in go_to_proteins.items():
            if len(proteins) < 2:
                continue
            name = (go_info or {}).get(go, {}).get('name', go)
            aspect = (go_info or {}).get(go, {}).get('aspect', '')
            self._add_hyperedge(
                nodes=list(set(proteins)),
                hyperedge_type='go_term',
                name=f"{go}: {name}",
                features={'go_id': go, 'aspect': aspect},
            )

    def add_ppi_clusters(self, ppi_edges: pd.DataFrame,
                          protein_list: List[str],
                          min_cluster_size: int = 3,
                          resolution: float = 1.0,
                          seed: int = 42):
        """
        Create hyperedges from Louvain communities in the PPI network.

        Uses weighted Louvain community detection (combined_score as edge
        weight) to partition the PPI graph into biologically coherent
        protein modules. Each community >= min_cluster_size becomes a
        hyperedge.

        Parameters
        ----------
        ppi_edges : pd.DataFrame
            PPI edges with ProteinID_A, ProteinID_B, combined_score
        protein_list : list of str
            All protein identifiers in the dataset
        min_cluster_size : int
            Minimum community size to create a hyperedge
        resolution : float
            Louvain resolution parameter (default 1.0).
            Higher values yield more, smaller communities.
        seed : int
            Random seed for reproducible community assignments
        """
        import networkx as nx

        G = nx.Graph()
        for _, row in ppi_edges.iterrows():
            a, b = row.get('ProteinID_A', ''), row.get('ProteinID_B', '')
            w = row.get('combined_score', 1.0)
            if a and b:
                G.add_edge(a, b, weight=w)

        if G.number_of_nodes() == 0:
            return

        communities = nx.community.louvain_communities(
            G, weight='weight', resolution=resolution, seed=seed
        )
        for community in communities:
            if len(community) >= min_cluster_size:
                self._add_hyperedge(
                    nodes=list(community),
                    hyperedge_type='ppi_cluster',
                    name=f"PPI_cluster_{len(self.hyperedges)}",
                    features={'size': len(community)},
                )

    def add_metabolite_hyperedges(
        self, metabolite_to_pathway: Dict[str, List[str]],
        pathway_genes: Optional[Dict[str, List[str]]] = None,
    ):
        """
        Create hyperedges from metabolite-to-pathway membership.

        Each KEGG metabolic pathway becomes a hyperedge containing all
        metabolites that participate in that pathway, plus (optionally)
        the genes already known to belong to it.

        Parameters
        ----------
        metabolite_to_pathway : dict
            {MetaboliteID: [pathway_id, ...]}
        pathway_genes : dict, optional
            {pathway_id: [GeneID, ...]} to co-join genes on the same path
        """
        path_to_metabs = defaultdict(list)
        for mid, paths in metabolite_to_pathway.items():
            for p in paths:
                path_to_metabs[p].append(mid)

        for path_id, metabods in path_to_metabs.items():
            nodes = list(set(metabods))
            if pathway_genes:
                nodes += list(pathway_genes.get(path_id, []))
            nodes = list(set(nodes))
            if len(nodes) < 2:
                continue
            self._add_hyperedge(
                nodes=nodes,
                hyperedge_type='metabolic_pathway',
                name=f"Metabolic_{path_id}",
                features={'pathway_id': path_id, 'n_metabolites': len(metabods)},
            )

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------

    def save(self, output_dir: str):
        """Save hypergraph to CSV files."""
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        # Hyperedges
        records = []
        for he in self.hyperedges:
            records.append({
                'hyperedge_id': he['hyperedge_id'],
                'type': he['type'],
                'name': he['name'],
                'nodes': ';'.join(he['nodes']),
                'n_nodes': len(he['nodes']),
            })
        he_df = pd.DataFrame(records)
        he_df.to_csv(out_dir / 'hyperedges.csv', index=False)
        print(f"  [Saved] {out_dir / 'hyperedges.csv'} ({len(he_df)} hyperedges)")

        # Incidence matrix
        if self.H is not None:
            np.save(out_dir / 'incidence_matrix.npy', self.H)

        # Node order
        node_df = pd.DataFrame({
            'node_id': self.node_list,
            'node_type': self._infer_node_types(),
        })
        node_df.to_csv(out_dir / 'hypergraph_nodes.csv', index=False)
        print(f"  [Saved] {out_dir / 'hypergraph_nodes.csv'} ({len(node_df)} nodes)")

        # Hyperedge type counts
        type_counts = pd.Series({
            k: len(v) for k, v in self.hyperedge_type_map.items()
        }).sort_values(ascending=False)
        type_counts.to_csv(out_dir / 'hyperedge_type_counts.csv', header=['count'])

    def _infer_node_types(self) -> List[str]:
        """Infer node types from identifier patterns."""
        types = []
        for node in self.node_list:
            if node.endswith('_RNA'):
                types.append('transcriptome')
            elif node.startswith('AMR:'):
                types.append('amr_mechanism')
            elif node.startswith('WP_') or '.' in node:
                types.append('proteome')
            else:
                types.append('genome')
        return types

    def summary(self) -> dict:
        """Return hypergraph summary."""
        if self.H is not None:
            density = np.sum(self.H > 0) / max(self.H.size, 1)
        else:
            density = 0.0

        type_counts = {
            k: len(v) for k, v in self.hyperedge_type_map.items()
        }

        return {
            'n_nodes': len(self.node_list),
            'n_hyperedges': len(self.hyperedges),
            'incidence_density': float(f"{density:.6f}"),
            'hyperedge_types': type_counts,
            'mean_hyperedge_size': float(
                np.mean([len(he['nodes']) for he in self.hyperedges])
            ) if self.hyperedges else 0,
        }
