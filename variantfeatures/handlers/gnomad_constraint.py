"""gnomAD gene-level constraint handler.

Fetches the per-gene constraint metrics that gnomAD computes from sequencing
data: pLI, LOEUF (oe_lof_upper), missense Z (mis_z), observed/expected ratios,
etc. These are the single most predictive features for novel LoF variants —
a nonsense in a gene with pLI > 0.9 is far more likely pathogenic than the
same nonsense in an LoF-tolerant gene.

Stored per gene in `gene_constraint`, joinable via `variant_consequences.gene_symbol`.
"""

from __future__ import annotations

import time
from typing import Optional

import requests


SOURCE = "gnomad_api"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
DEFAULT_TIMEOUT = 30
DEFAULT_DATASET = "gnomad_v4"
_MAX_RETRIES_429 = 3


class HandlerError(Exception):
    pass


_QUERY = """
query GeneConstraint($symbol: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    gene_id
    canonical_transcript_id
    gnomad_constraint {
      pli
      lof_z
      mis_z
      syn_z
      oe_lof
      oe_lof_lower
      oe_lof_upper
      oe_mis
      oe_mis_upper
      exp_lof
      obs_lof
    }
  }
}
"""


def fetch_constraint(gene_symbol: str, *, timeout: int = DEFAULT_TIMEOUT) -> Optional[dict]:
    """Return gnomAD's gene record for `gene_symbol`, or None if not present."""
    body_payload = {"query": _QUERY, "variables": {"symbol": gene_symbol}}

    backoff = 2.0
    for attempt in range(_MAX_RETRIES_429 + 1):
        resp = requests.post(
            GNOMAD_API,
            json=body_payload,
            headers={"Content-Type": "application/json"},
            timeout=timeout,
        )
        if resp.status_code == 429 and attempt < _MAX_RETRIES_429:
            retry_after = resp.headers.get("Retry-After")
            try:
                wait = float(retry_after) if retry_after else backoff
            except ValueError:
                wait = backoff
            time.sleep(min(wait, 30.0))
            backoff *= 2
            continue
        resp.raise_for_status()
        body = resp.json()
        break

    if "errors" in body:
        msgs = "; ".join(e.get("message", "?") for e in body["errors"])
        # Variant-not-found-style errors -> silent absent
        if "not found" in msgs.lower():
            return None
        raise HandlerError(f"gnomAD API error for gene {gene_symbol!r}: {msgs}")
    gene = (body.get("data") or {}).get("gene")
    return gene


def annotate_gene(db, gene_symbol: str, *, dataset: str = DEFAULT_DATASET, timeout: int = DEFAULT_TIMEOUT) -> bool:
    """Fetch + persist constraint metrics for one gene. Returns True if a row was written."""
    record = fetch_constraint(gene_symbol, timeout=timeout)
    if not record:
        return False
    constraint = record.get("gnomad_constraint")
    if not constraint:
        return False

    db.upsert_gene_constraint(
        gene_symbol,
        dataset=dataset,
        pli=constraint.get("pli"),
        lof_z=constraint.get("lof_z"),
        mis_z=constraint.get("mis_z"),
        syn_z=constraint.get("syn_z"),
        oe_lof=constraint.get("oe_lof"),
        oe_lof_lower=constraint.get("oe_lof_lower"),
        oe_lof_upper=constraint.get("oe_lof_upper"),
        oe_mis=constraint.get("oe_mis"),
        oe_mis_upper=constraint.get("oe_mis_upper"),
        exp_lof=constraint.get("exp_lof"),
        obs_lof=constraint.get("obs_lof"),
        source=SOURCE,
    )
    return True
