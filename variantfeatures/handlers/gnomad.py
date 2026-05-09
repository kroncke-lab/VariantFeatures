"""gnomAD GraphQL handler.

Per-variant lookup via the gnomAD GraphQL API. Pulls the full population
breakdown (afr/amr/asj/eas/fin/mid/nfe/sas/remaining) for both exomes and
genomes plus the popmax stats, and writes one annotations_population row
per (dataset, pop).

Fills the gnomAD coverage gap that MyVariant.info has for many variants.

API: https://gnomad.broadinstitute.org/api
"""

from __future__ import annotations

import time
from typing import Iterable, Optional

import requests

from .clingen_ar import normalize_chromosome


SOURCE = "gnomad"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
DEFAULT_RATE_LIMIT_SEC = 0.6  # gnomAD prefers ~2 req/sec or less
DEFAULT_TIMEOUT = 60
DEFAULT_DATASET = "gnomad_r4"

# GraphQL error messages that mean "this variant simply isn't in gnomAD" — treat as
# a silent absent, not a failure. (gnomAD returns these as `errors` not `data.variant: null`.)
_NOT_FOUND_PHRASES = ("variant not found", "no variant found")
_MAX_RETRIES_429 = 3

# gnomAD v4 population codes. 'mid' (Middle Eastern) is new in v4.
GNOMAD_POPS = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "remaining"]


class HandlerError(Exception):
    pass


_QUERY = """
query VariantQuery($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    rsids
    exome {
      ac
      an
      af
      homozygote_count
      filters
      populations {
        id
        ac
        an
        homozygote_count
      }
      faf95 { popmax popmax_population }
    }
    genome {
      ac
      an
      af
      homozygote_count
      filters
      populations {
        id
        ac
        an
        homozygote_count
      }
      faf95 { popmax popmax_population }
    }
  }
}
"""


def fetch_gnomad(chromosome: str, position: int, ref: str, alt: str,
                 *, dataset: str = DEFAULT_DATASET, timeout: int = DEFAULT_TIMEOUT) -> Optional[dict]:
    """Fetch one variant's gnomAD record. Returns None if not in gnomAD.

    Retries up to 3 times on HTTP 429 with exponential backoff, then re-raises.
    "Variant not found" GraphQL errors are returned as None (silent absent).
    """
    chrom = normalize_chromosome(chromosome)
    variant_id = f"{chrom}-{position}-{ref}-{alt}"
    body_payload = {"query": _QUERY, "variables": {"variantId": variant_id, "dataset": dataset}}

    backoff = 2.0
    for attempt in range(_MAX_RETRIES_429 + 1):
        resp = requests.post(
            GNOMAD_API,
            json=body_payload,
            headers={"Content-Type": "application/json"},
            timeout=timeout,
        )
        if resp.status_code == 429 and attempt < _MAX_RETRIES_429:
            # Respect Retry-After header if present, else exponential backoff.
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
        msgs = [e.get("message", "?") for e in body["errors"]]
        joined_lc = "; ".join(m.lower() for m in msgs)
        if any(p in joined_lc for p in _NOT_FOUND_PHRASES):
            return None
        raise HandlerError(f"gnomAD API error: {'; '.join(msgs)}")
    return (body.get("data") or {}).get("variant")


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def _per_pop_rows(payload: dict, dataset_label: str) -> Iterable[dict]:
    """Yield {dataset, pop, af, ac, an, n_homozygotes, filter_status} rows for one assay (exome / genome)."""
    if not payload:
        return
    filters = payload.get("filters") or []
    filter_status = ",".join(filters) if filters else "PASS"

    ac = payload.get("ac")
    an = payload.get("an")
    af = payload.get("af")
    if af is None and ac is not None and an:
        af = ac / an if an else None
    yield {
        "dataset": dataset_label,
        "pop": "all",
        "af": af,
        "ac": ac,
        "an": an,
        "n_homozygotes": payload.get("homozygote_count"),
        "filter_status": filter_status,
    }

    # popmax (filtering allele frequency 95%)
    faf = payload.get("faf95") or {}
    popmax_af = faf.get("popmax")
    popmax_pop = faf.get("popmax_population")
    if popmax_af is not None:
        yield {
            "dataset": dataset_label,
            "pop": "popmax",
            "af": popmax_af,
            "ac": None,
            "an": None,
            "n_homozygotes": None,
            "filter_status": popmax_pop,  # store which pop was max in the filter slot
        }

    # Per-population breakdown
    pops = payload.get("populations") or []
    for entry in pops:
        pop = (entry.get("id") or "").lower()
        if not pop:
            continue
        # gnomAD v4 returns sex breakdowns as standalone "XX"/"XY" plus per-pop
        # combos like "afr_XX". Skip both forms.
        if pop in {"xx", "xy"} or pop.endswith("_xx") or pop.endswith("_xy"):
            continue
        ac = entry.get("ac")
        an = entry.get("an")
        af = (ac / an) if (ac is not None and an) else None
        yield {
            "dataset": dataset_label,
            "pop": pop,
            "af": af,
            "ac": ac,
            "an": an,
            "n_homozygotes": entry.get("homozygote_count"),
            "filter_status": filter_status,
        }


def parse_gnomad_variant(variant: dict, *, dataset_version: str = "v4") -> Iterable[dict]:
    """Yield population rows from a gnomAD variant payload (exome + genome combined)."""
    if not variant:
        return
    exome = variant.get("exome")
    if exome:
        yield from _per_pop_rows(exome, f"gnomad_exomes_{dataset_version}")
    genome = variant.get("genome")
    if genome:
        yield from _per_pop_rows(genome, f"gnomad_genomes_{dataset_version}")


def parse_aliases(variant: dict) -> Iterable[dict]:
    """Yield alias rows from a gnomAD variant: gnomad_id and any rsIDs."""
    if not variant:
        return
    vid = variant.get("variant_id")
    if vid:
        yield {"alias_type": "gnomad_id", "alias_value": vid}
    for rs in variant.get("rsids") or []:
        if rs:
            rs_str = rs if str(rs).startswith("rs") else f"rs{rs}"
            yield {"alias_type": "rsid", "alias_value": rs_str}


# ---------------------------------------------------------------------------
# Persistence + worker entry
# ---------------------------------------------------------------------------

def persist(db, variant_id: int, variant: dict, *, dataset_version: str = "v4") -> dict:
    """Write all population rows + aliases for one gnomAD payload. Returns counts."""
    counts = {"population": 0, "aliases": 0}
    for row in parse_gnomad_variant(variant, dataset_version=dataset_version):
        db.upsert_population(variant_id, source=SOURCE, **row)
        counts["population"] += 1
    aliases = list(parse_aliases(variant))
    if aliases:
        counts["aliases"] = db.add_aliases(variant_id, aliases, source=SOURCE)
    return counts


def handle(db, variant_id: int, payload: Optional[str] = None,
           *, dataset: str = DEFAULT_DATASET, timeout: int = DEFAULT_TIMEOUT) -> None:
    """Worker entry point: fetch one variant from gnomAD and persist population rows."""
    cur = db.conn.execute(
        "SELECT chromosome, position, ref, alt FROM variants WHERE id = ?",
        [variant_id],
    )
    row = cur.fetchone()
    if row is None:
        raise HandlerError(f"variant_id {variant_id} not found in variants table")

    variant = fetch_gnomad(row["chromosome"], row["position"], row["ref"], row["alt"],
                           dataset=dataset, timeout=timeout)
    # If absent in gnomAD, mark job done with no rows (job is "succeeded — variant not present").
    if variant is None:
        return
    persist(db, variant_id, variant)
