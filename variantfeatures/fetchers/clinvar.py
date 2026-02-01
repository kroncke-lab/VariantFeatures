"""Fetch ClinVar annotations."""

from typing import Iterator
import requests

# ClinVar FTP: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
# Useful files:
#   variant_summary.txt.gz - all variants with clinical significance
#   submission_summary.txt.gz - individual submissions

CLINVAR_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"


def fetch_clinvar(gene: str) -> Iterator[dict]:
    """
    Fetch ClinVar annotations for a gene.
    
    Yields dicts with:
        - hgvs_p: e.g., "p.Arg528His"
        - clinvar_id: int
        - clinvar_significance: e.g., "Pathogenic", "Likely benign"
        - clinvar_review_status: e.g., "criteria provided, single submitter"
        - clinvar_last_evaluated: date
    """
    # TODO: Implement
    # Options:
    # 1. Download variant_summary.txt.gz, filter by gene
    # 2. Use NCBI E-utilities API (slower but no bulk download)
    # 3. Use ClinVar's VCF files
    raise NotImplementedError("ClinVar fetcher not yet implemented")
