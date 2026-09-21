"""
identifiers.py
Shared RefSeq protein identifier normalisation utilities.

The pipeline must reconcile several identifier forms of the same protein:

    locus_tag                      (e.g. EW036_RS26095)
      |
      |  gene_protein_mapping
      v
    protein_id (versioned)         (e.g. WP_000000542.1)
      |
      |  normalize_wp_id()
      v
    protein_id (unversioned)       (e.g. WP_000000542)
      |
      v
    eggNOG / COG annotation        (keyed by unversioned WP accession)

These helpers centralise the version-suffix stripping so that joins between
the genome / proteomics (versioned) and eggNOG / UniProt (unversioned) layers
align everywhere.
"""

__all__ = ['normalize_wp_id']


def normalize_wp_id(protein_id) -> str:
    """Strip the version suffix from a RefSeq protein accession.

    ``WP_000000542.1`` -> ``WP_000000542``

    Parameters
    ----------
    accession : str
        A RefSeq protein accession, possibly version-suffixed.

    Returns
    -------
    str
        The normalized (unversioned) accession, or the original string if the
        input carries no version suffix / is not a WP accession.
    """
    if not protein_id:
        return ''
    s = str(protein_id).strip()
    if not s:
        return ''
    if s.startswith('WP_'):
        # Only strip the numeric version tail (e.g. `.1`).
        return s.split('.')[0].strip()
    return s