"""
build_aligned_omics.py
Process real multi-omics data files into the aligned table format.

Identifier bridging strategy:
  - Proteomics identifies proteins by WP_* accession (Protein IDs)
  - Genome GFF links WP_* (via CDS protein_id) to locus_tag (GeneID) and gene symbol (GeneName)
  - RNA counts use locus_tag (GeneID) and gene symbol (GeneName)

Merge pipeline:
  1. Parse genome GFF CDS entries to build: ProteinID -> GeneID, GeneName, replicon, coords
  2. Parse proteomics: ProteinID + Intensity values
  3. Merge proteomics -> genome on ProteinID  (inner join)
  4. Merge result -> RNA on GeneID           (inner join)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def parse_genome_protein_map(gff_path: str) -> pd.DataFrame:
    """
    Parse the GFF and build a ProteinID -> gene feature mapping.

    Extracts from both 'gene' and 'CDS' entries:
      - ProteinID (WP_* from CDS protein_id)
      - GeneID (locus_tag from gene entry)
      - GeneName (gene symbol from Name/gene attribute)
      - Replicon, coordinates, strand, CDS length
    """
    gene_info = {}      # locus_tag -> info
    cds_info = {}       # locus_tag or protein_id -> info

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            ftype = parts[2]
            attr = dict(a.split('=', 1) for a in parts[8].split(';') if '=' in a)

            if ftype == 'gene':
                lt = attr.get('locus_tag', '')
                if lt:
                    symbol = attr.get('gene', attr.get('Name', lt))
                    gene_info[lt] = {
                        'GeneID': lt,
                        'GeneName': symbol,
                        'replicon': parts[0],
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'strand': parts[6],
                    }

            elif ftype == 'CDS':
                pid = attr.get('protein_id', '').strip()
                lt = attr.get('locus_tag', '')
                if pid:
                    symbol = attr.get('gene', attr.get('Name', ''))
                    parent = attr.get('Parent', '').replace('gene-', '')
                    key = lt if lt else parent
                    cds_info[key] = {
                        'ProteinID': pid,
                        'CDS_locus_tag': lt,
                        'CDS_GeneName': symbol,
                    }
                    # Also index by ProteinID directly
                    cds_info[pid] = {
                        'ProteinID': pid,
                        'CDS_locus_tag': lt,
                        'CDS_GeneName': symbol,
                    }

    # Build final table: for each gene, attach ProteinID from CDS
    rows = []
    for lt, g in gene_info.items():
        cds = cds_info.get(lt, {})
        pid = cds.get('ProteinID', '')
        cds_symbol = cds.get('CDS_GeneName', '')
        # Use CDS gene symbol when available (it's often more specific)
        gene_symbol = cds_symbol if cds_symbol else g['GeneName']

        rows.append({
            'GeneID': g['GeneID'],
            'GeneName': gene_symbol,
            'ProteinID': pid,
            'replicon': g['replicon'],
            'start': g['start'],
            'end': g['end'],
            'strand': g['strand'],
            'CDS_length': g['end'] - g['start'] + 1,
        })

    df = pd.DataFrame(rows)

    # Also build a direct ProteinID -> GeneID lookup for proteins mapped by WP_* only
    protein_to_gene = {}
    for r in rows:
        if r['ProteinID']:
            protein_to_gene[r['ProteinID']] = r['GeneID']

    df = df.sort_values(['replicon', 'start']).reset_index(drop=True)
    return df


def process_transcriptomics(counts_path: str) -> pd.DataFrame:
    """
    Process RNA-seq count data.
    RPMI: 50857, 50858, 51033, 51034, 51035, 51036
    Sera: 50863, 50864, 51037, 51038, 51039, 51040
    """
    df = pd.read_csv(counts_path, sep='\t')

    rpmi_cols = ['50857', '50858', '51033', '51034', '51035', '51036']
    sera_cols = ['50863', '50864', '51037', '51038', '51039', '51040']

    for col in rpmi_cols + sera_cols:
        total = df[col].sum()
        cpm = df[col] / total * 1e6 if total > 0 else 0.0
        df[col + '_log2'] = np.log2(cpm + 1)

    df['RNA_RPMI'] = df[[c + '_log2' for c in rpmi_cols]].mean(axis=1)
    df['RNA_Sera'] = df[[c + '_log2' for c in sera_cols]].mean(axis=1)
    df['RNA_logFC'] = df['RNA_Sera'] - df['RNA_RPMI']

    result = df[['GeneID', 'GeneName', 'RNA_RPMI', 'RNA_Sera', 'RNA_logFC']].copy()
    result = result.sort_values('GeneID').reset_index(drop=True)
    return result


def process_proteomics(proteomics_path: str) -> pd.DataFrame:
    """
    Process MaxQuant proteinGroups.txt.

    Intensity columns:
      RPMI: 50973..50978
      Sera: 50979..50984

    Returns DataFrame indexed by ProteinID with mean log2 intensity per condition.
    """
    df = pd.read_csv(proteomics_path, sep='\t', low_memory=False)

    rpmi_ids = [50973, 50974, 50975, 50976, 50977, 50978]
    sera_ids = [50979, 50980, 50981, 50982, 50983, 50984]

    # Try LFQ first, fall back to raw Intensity
    rpmi_cols = [f'LFQ intensity {sid}_B36_RPMI' for sid in rpmi_ids]
    sera_cols = [f'LFQ intensity {sid}_B36_Pooled sera' for sid in sera_ids]
    missing = [c for c in rpmi_cols + sera_cols if c not in df.columns]
    if missing:
        rpmi_cols = [f'Intensity {sid}_B36_RPMI' for sid in rpmi_ids]
        sera_cols = [f'Intensity {sid}_B36_Pooled sera' for sid in sera_ids]

    for col in rpmi_cols + sera_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
        df[col + '_log2'] = np.log2(df[col] + 1)

    # Extract primary Protein ID (strip REV__ or other prefixes)
    df['ProteinID_raw'] = df['Protein IDs'].str.split(';').str[0]
    df['ProteinID'] = df['ProteinID_raw'].str.replace('^REV__', '', regex=True)

    # Keep only E. coli proteins (filter by FASTA header)
    is_ecoli = df['Fasta headers'].str.contains('Escherichia|E. coli|B36',
                                                  case=False, na=False, regex=True)
    if is_ecoli.sum() == 0:
        print("  [Warning] No E. coli proteins identified by FASTA header — using all")
        is_ecoli = pd.Series([True] * len(df))

    df_ecoli = df[is_ecoli].copy()
    print(f"    Total protein groups: {len(df)}")
    print(f"    E. coli protein groups: {len(df_ecoli)}")

    # Compute per-replicate log2 intensity columns and mean per condition
    rpmi_log2 = [c + '_log2' for c in rpmi_cols]
    sera_log2 = [c + '_log2' for c in sera_cols]

    df_ecoli['Protein_RPMI'] = df_ecoli[rpmi_log2].mean(axis=1)
    df_ecoli['Protein_Sera'] = df_ecoli[sera_log2].mean(axis=1)
    df_ecoli['Protein_logFC'] = df_ecoli['Protein_Sera'] - df_ecoli['Protein_RPMI']

    # Aggregate: mean per ProteinID (some WP_* appear in multiple rows/contaminants)
    df_agg = df_ecoli.groupby('ProteinID', as_index=False).agg({
        'Protein_RPMI': 'mean',
        'Protein_Sera': 'mean',
        'Protein_logFC': 'mean',
    })

    return df_agg


def build_aligned_table(genome_df, rna_df, protein_df, gff_path=None) -> pd.DataFrame:
    """
    Build aligned multi-omics table.

    Strategy:
      1. proteomics (ProteinID=WP_*) -> genome (ProteinID) on ProteinID
         -> yields GeneID + GeneName + protein values
      2. Merged -> RNA on GeneID -> yields RNA values

    This ensures every aligned row has measurements in all 3 omics.
    """
    print(f"\n    Merging proteomics -> genome on ProteinID...")
    # Note: genome_df already has ProteinID column from parse_genome_protein_map
    prot_genome = protein_df.merge(
        genome_df[['GeneID', 'GeneName', 'ProteinID', 'replicon', 'start', 'end', 'strand', 'CDS_length']],
        on='ProteinID',
        how='inner',
    )
    print(f"    Proteomics-genome merge: {len(prot_genome)} aligned "
          f"(from {len(protein_df)} protein groups)")

    if len(prot_genome) == 0:
        print("    [ERROR] No proteins matched genome. Debugging...")
        print(f"    Sample protein IDs: {protein_df['ProteinID'].head(5).tolist()}")
        print(f"    Genome ProteinID samples: {genome_df['ProteinID'].dropna().head(5).tolist()}")
        return pd.DataFrame()

    # Drop rows without GeneID
    prot_genome = prot_genome[prot_genome['GeneID'].notna() & (prot_genome['GeneID'] != '')]

    # Now merge with RNA on GeneID
    print(f"    Merging -> RNA on GeneID...")
    aligned = prot_genome.merge(
        rna_df[['GeneID', 'RNA_RPMI', 'RNA_Sera', 'RNA_logFC']],
        on='GeneID',
        how='inner',
    )
    
    print(f"    RNA merge: {len(aligned)} fully aligned triples")

    # Build final table
    result = pd.DataFrame()
    result['GeneID'] = aligned['GeneID']
    result['GeneName'] = aligned['GeneName']
    result['RNA_RPMI'] = aligned['RNA_RPMI']
    result['RNA_Sera'] = aligned['RNA_Sera']
    result['RNA_logFC'] = aligned['RNA_logFC']
    result['ProteinID'] = aligned['ProteinID']
    result['Protein_RPMI'] = aligned['Protein_RPMI']
    result['Protein_Sera'] = aligned['Protein_Sera']
    result['Protein_logFC'] = aligned['Protein_logFC']

    for col in ['replicon', 'start', 'end', 'strand', 'CDS_length']:
        if col in aligned.columns:
            result[col] = aligned[col]

    return result.drop_duplicates(subset='GeneID').reset_index(drop=True)


def main():
    base_dir = Path(__file__).parent.absolute()

    gff_path = (
        base_dir / 'genome' / '900622635.1' / 'ncbi_dataset' / 'data'
        / 'GCF_900622635.1' / 'genomic.gff'
    )
    counts_path = base_dir / 'transcriptomic' / 'GSE152966_gene_counts.txt'
    proteomics_path = (
        base_dir / 'proteomics' / 'Ecoli_B36_MQtxt' / 'txt' / 'proteinGroups.txt'
    )
    output_path = base_dir / 'project' / 'data' / 'aligned_omics_real.csv'

    print("=" * 60)
    print("BUILDING ALIGNED MULTI-OMICS TABLE FROM REAL DATA")
    print("=" * 60)

    print("\n[1] Parsing genome GFF (ProteinID -> GeneID mapping)...")
    genome_df = parse_genome_protein_map(str(gff_path))
    print(f"    {len(genome_df)} genes with {genome_df['ProteinID'].astype(bool).sum()} annotated ProteinIDs")

    print("\n[2] Processing transcriptomics counts...")
    rna_df = process_transcriptomics(str(counts_path))
    print(f"    {len(rna_df)} RNA profiles (GeneID: {rna_df['GeneID'].iloc[0]} .. {rna_df['GeneID'].iloc[-1]})")

    print("\n[3] Processing proteomics protein groups...")
    protein_df = process_proteomics(str(proteomics_path))
    print(f"    {len(protein_df)} E. coli protein profiles aggregated")

    print("\n[4] Building aligned multi-omics table...")
    aligned = build_aligned_table(genome_df, rna_df, protein_df, gff_path)

    if len(aligned) == 0:
        print("\n[ERROR] Alignment failed — no overlapping genes across all 3 omics.")
        return

    aligned.to_csv(output_path, index=False)
    print(f"\n[5] Saved to: {output_path}")
    print(f"    Shape: {aligned.shape}")
    print(f"    Columns: {list(aligned.columns)}")

    print("\n" + "-" * 50)
    print("SUMMARY STATISTICS")
    print("-" * 50)
    print(f"  Total aligned:         {len(aligned)}")
    print(f"  Unique replicons:      {aligned['replicon'].nunique()}")
    print(f"  RNA_RPMI range:        {aligned['RNA_RPMI'].min():.2f} - {aligned['RNA_RPMI'].max():.2f}")
    print(f"  RNA_Sera range:        {aligned['RNA_Sera'].min():.2f} - {aligned['RNA_Sera'].max():.2f}")
    print(f"  RNA_logFC range:       {aligned['RNA_logFC'].min():.3f} - {aligned['RNA_logFC'].max():.3f}")
    print(f"  Protein_RPMI range:    {aligned['Protein_RPMI'].min():.2f} - {aligned['Protein_RPMI'].max():.2f}")
    print(f"  Protein_Sera range:    {aligned['Protein_Sera'].min():.2f} - {aligned['Protein_Sera'].max():.2f}")
    print(f"  Protein_logFC range:   {aligned['Protein_logFC'].min():.3f} - {aligned['Protein_logFC'].max():.3f}")

    n_missing = aligned.isna().sum().sum()
    if n_missing > 0:
        print(f"\n  [Warning] {n_missing} missing values")
        for col in aligned.columns:
            nm = aligned[col].isna().sum()
            if nm > 0:
                print(f"    {col}: {nm}")
    else:
        print("\n  No missing values")

    print("\n  Sample rows:")
    print(aligned[['GeneID', 'GeneName', 'RNA_logFC', 'Protein_logFC']].head(10).to_string(index=False))

    return aligned


if __name__ == '__main__':
    main()
