"""Fetch AlphaMissense predictions."""

from typing import Iterator
import requests

# AlphaMissense data is available from:
# https://console.cloud.google.com/storage/browser/dm_alphamissense
# Or via EBI: https://alphafold.ebi.ac.uk/download

ALPHAMISSENSE_URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz"


def fetch_alphamissense(gene: str) -> Iterator[dict]:
    """
    Fetch AlphaMissense scores for a gene.
    
    Yields dicts with:
        - hgvs_p: e.g., "p.Ala123Val"
        - alphamissense_score: float 0-1
        - alphamissense_class: "likely_benign" | "ambiguous" | "likely_pathogenic"
    """
    # TODO: Implement - need to download and parse the TSV
    # The full file is ~4GB, so we may want to:
    # 1. Download once and cache locally
    # 2. Filter by gene during parsing
    # 3. Or use a pre-filtered per-gene approach
    raise NotImplementedError("AlphaMissense fetcher not yet implemented")
