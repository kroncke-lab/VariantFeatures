"""Small UniProt resolver used by protein-level annotators.

The tool is gene-agnostic: no gene list is hardcoded. Any human gene symbol is
resolved to a reviewed UniProtKB accession via UniProt's REST search. Optional
runtime overrides may be supplied via the `VARIANTFEATURES_UNIPROT_MAP` /
`ALPHAMISSENSE_UNIPROT` env vars (comma-separated `GENE:ACCESSION`) or the
`extra` argument -- e.g. to pin a specific isoform accession or to run offline.
"""

from __future__ import annotations

import os
from typing import Optional

import requests


UNIPROT_REST = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_TIMEOUT = 30


class UniProtError(Exception):
    pass


def gene_uniprot_map(extra: Optional[dict[str, str]] = None) -> dict[str, str]:
    """Return configured GENE -> UniProt accession overrides (no hardcoded genes).

    Overrides come from the `VARIANTFEATURES_UNIPROT_MAP` / `ALPHAMISSENSE_UNIPROT`
    env vars (comma-separated `GENE:ACCESSION`) and the `extra` argument. Genes
    not listed here are resolved dynamically via UniProt REST -- the tool does
    not carry a built-in gene list.
    """
    out: dict[str, str] = {}
    env_pairs = os.environ.get("VARIANTFEATURES_UNIPROT_MAP") or os.environ.get("ALPHAMISSENSE_UNIPROT")
    if env_pairs:
        for pair in env_pairs.split(","):
            if ":" not in pair:
                continue
            gene, accession = pair.split(":", 1)
            gene = gene.strip().upper()
            accession = accession.strip()
            if gene and accession:
                out[gene] = accession
    if extra:
        for gene, accession in extra.items():
            out[gene.upper()] = accession
    return out


def resolve_uniprot_accession(
    gene_symbol: str,
    *,
    timeout: int = DEFAULT_TIMEOUT,
    extra: Optional[dict[str, str]] = None,
) -> str:
    """Resolve a human gene symbol to a reviewed UniProtKB accession."""
    gene = gene_symbol.strip().upper()
    configured = gene_uniprot_map(extra)
    if gene in configured:
        return configured[gene]

    params = {
        "query": f"(gene_exact:{gene}) AND (organism_id:9606) AND (reviewed:true)",
        "fields": "accession,id,gene_primary,reviewed,organism_id,length",
        "format": "json",
        "size": "1",
    }
    resp = requests.get(UNIPROT_REST, params=params, timeout=timeout)
    if resp.status_code == 404:
        raise UniProtError(f"No UniProt entry found for gene {gene!r}")
    resp.raise_for_status()
    data = resp.json()
    results = data.get("results") if isinstance(data, dict) else None
    if not results:
        raise UniProtError(
            f"No reviewed human UniProt entry found for gene {gene!r}; "
            "set VARIANTFEATURES_UNIPROT_MAP=GENE:ACCESSION to override."
        )
    accession = results[0].get("primaryAccession")
    if not accession:
        raise UniProtError(f"UniProt response for {gene!r} did not contain primaryAccession")
    return accession
