"""Small UniProt resolver used by protein-level annotators.

The built-in map covers the current cardiac/metabolic/amyloid grant targets and
common ACMG genes. For any other human gene, fall back to UniProt's REST search.
"""

from __future__ import annotations

import os
from typing import Optional

import requests


UNIPROT_REST = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_TIMEOUT = 30


BUILTIN_GENE_TO_UNIPROT: dict[str, str] = {
    "KCNH2": "Q12809",
    "KCNQ1": "P51787",
    "SCN5A": "Q14524",
    "RYR2": "Q92736",
    "CACNA1C": "Q13936",
    "KCNJ2": "P63252",
    "KCNE1": "P15382",
    "KCNE2": "Q9Y6J6",
    "BRCA1": "P38398",
    "BRCA2": "P51587",
    "TP53": "P04637",
    "MYH7": "P12883",
    "MYBPC3": "Q14896",
    "TTN": "Q8WZ42",
    "LDLR": "P01130",
    "APOB": "P04114",
    "PCSK9": "Q8NBP7",
    "TTR": "P02766",
    "ALPL": "P05186",
    "RB1": "P06400",
    "VHL": "P40337",
    "ATP7B": "P35670",
    "RYR1": "P21817",
    "CACNA1S": "Q13698",
    "SDHB": "P21912",
    "SDHC": "Q99643",
    "SDHD": "O14521",
    "MEN1": "O00255",
    "RET": "P07949",
    "PTEN": "P60484",
    "STK11": "Q15831",
    "MUTYH": "Q9UIF7",
    "MLH1": "P40692",
    "MSH2": "P43246",
    "MSH6": "P52701",
    "PMS2": "P54278",
    "APC": "P25054",
}


class UniProtError(Exception):
    pass


def gene_uniprot_map(extra: Optional[dict[str, str]] = None) -> dict[str, str]:
    """Return the configured GENE -> UniProt accession map.

    `VARIANTFEATURES_UNIPROT_MAP` may add or override pairs as
    comma-separated `GENE:ACCESSION` values.
    """
    out = dict(BUILTIN_GENE_TO_UNIPROT)
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
