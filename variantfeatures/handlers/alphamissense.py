"""AlphaMissense local-file handler.

Reads the published `AlphaMissense_aa_substitutions.tsv.gz` (~1.2GB) which is
keyed by `(uniprot_id, protein_variant)` rather than genomic coordinates. We
match pending jobs by joining `annotation_jobs` -> `variants` ->
`variant_consequences` to recover `(gene_symbol, aa_ref, aa_pos, aa_alt)`,
then convert to the AlphaMissense short form (e.g. `KCNH2 + A614V`).

Per-gene UniProt resolution is gene-agnostic: each gene symbol is resolved to a
reviewed UniProtKB accession via `variantfeatures.uniprot.resolve_uniprot_accession`
(UniProt REST), with optional overrides via the `gene_uniprot` argument or the
`ALPHAMISSENSE_UNIPROT` env var (comma-separated `GENE:UNIPROT` pairs). No gene
list is hardcoded.

File format:
    uniprot_id  protein_variant  am_pathogenicity  am_class

Configure file path via `ALPHAMISSENSE_FILE` env var or `file_path=` kwarg.
Default: `data/alphamissense/AlphaMissense_aa_substitutions.tsv.gz`.

Download: https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path
from typing import Optional

import requests

from ..uniprot import UniProtError, resolve_uniprot_accession


SOURCE = "alphamissense"
DEFAULT_TIMEOUT = 30

DEFAULT_REL_PATH = Path("data") / "alphamissense" / "AlphaMissense_aa_substitutions.tsv.gz"
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


class HandlerError(Exception):
    pass


# ---------------------------------------------------------------------------
# Path / config helpers
# ---------------------------------------------------------------------------

def resolve_file_path(explicit: Optional[str] = None) -> Path:
    if explicit:
        return Path(explicit)
    env_path = os.environ.get("ALPHAMISSENSE_FILE")
    if env_path:
        return Path(env_path)
    return PROJECT_ROOT / DEFAULT_REL_PATH


def file_present(file_path: Optional[str] = None) -> bool:
    return resolve_file_path(file_path).exists()


def gene_uniprot_map(extra: Optional[dict[str, str]] = None) -> dict[str, str]:
    """Return GENE -> UniProt accession *overrides* (no hardcoded gene list).

    Overrides come from the `ALPHAMISSENSE_UNIPROT` env var and the `extra`
    argument. Genes not listed here are resolved dynamically via UniProt REST
    (see `run_batch`), so the handler works for any gene without configuration.
    """
    out: dict[str, str] = {}
    env_pairs = os.environ.get("ALPHAMISSENSE_UNIPROT")
    if env_pairs:
        for pair in env_pairs.split(","):
            if ":" in pair:
                g, u = pair.split(":", 1)
                out[g.strip().upper()] = u.strip()
    if extra:
        for g, u in extra.items():
            out[g.upper()] = u
    return out


# ---------------------------------------------------------------------------
# Batch runner
# ---------------------------------------------------------------------------

def run_batch(
    db,
    *,
    file_path: Optional[str] = None,
    limit: Optional[int] = None,
    gene_filter: Optional[str] = None,
    progress_every: int = 5_000_000,
    progress_callback=None,
    gene_uniprot: Optional[dict[str, str]] = None,
    timeout: int = DEFAULT_TIMEOUT,
) -> dict:
    """Drain pending alphamissense jobs by single-pass over the AM TSV.

    Gene symbols are resolved to UniProt accessions dynamically (UniProt REST),
    honouring `gene_uniprot` / `ALPHAMISSENSE_UNIPROT` overrides first.

    Returns: {"claimed": int, "matched": int, "failed": int, "lines_scanned": int}
    """
    path = resolve_file_path(file_path)
    if not path.exists():
        raise HandlerError(
            f"AlphaMissense file not found at {path}. "
            f"Set ALPHAMISSENSE_FILE or download from "
            f"https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz"
        )

    pending = _load_pending_index(
        db,
        limit=limit,
        gene_filter=gene_filter,
        overrides=gene_uniprot_map(gene_uniprot),
        timeout=timeout,
    )
    if not pending["by_key"]:
        # Either no pending jobs, or none of them matched a known gene -> mark them failed.
        failed = _mark_unmappable(db, pending["unmappable"])
        return {"claimed": failed, "matched": 0, "failed": failed, "lines_scanned": 0}

    matched = 0
    lines_scanned = 0

    by_key = pending["by_key"]

    with _open_maybe_gzip(path) as fh:
        header = None
        for line in fh:
            lines_scanned += 1
            if progress_callback and lines_scanned % progress_every == 0:
                progress_callback({"lines": lines_scanned, "matched": matched, "remaining": len(by_key)})

            # Skip license/copyright comment lines.
            if line.startswith("#"):
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                # Sanity check: expect (uniprot_id, protein_variant, am_pathogenicity, am_class)
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            uniprot_id, protein_variant, score_str, am_class = cols[0], cols[1], cols[2], cols[3]
            key = (uniprot_id, protein_variant)
            jobs = by_key.get(key)
            if not jobs:
                continue
            try:
                score = float(score_str)
            except ValueError:
                continue
            category = _AM_CLASS_TO_CODE.get(am_class, am_class)
            for job in jobs:
                _persist_match(db, job=job, score=score, category=category, am_class=am_class)
                matched += 1
            del by_key[key]
            if not by_key:
                break

    # Anything still in by_key was not in the AM file (extremely rare). Mark failed.
    failed = 0
    for jobs in by_key.values():
        for job in jobs:
            db.mark_job_failed(job["job_id"], "no AlphaMissense entry for this (uniprot, AA-substitution)")
            failed += 1

    # Plus jobs we couldn't even map to a uniprot to begin with (unknown gene).
    failed += _mark_unmappable(db, pending["unmappable"])

    return {"claimed": matched + failed, "matched": matched, "failed": failed, "lines_scanned": lines_scanned}


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

_AM_CLASS_TO_CODE = {
    # Newer AlphaMissense releases use bare labels; older releases prefixed with `likely_`.
    "benign": "B",
    "likely_benign": "B",
    "ambiguous": "A",
    "pathogenic": "P",
    "likely_pathogenic": "P",
}


def _open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _load_pending_index(
    db,
    *,
    limit: Optional[int],
    gene_filter: Optional[str],
    overrides: dict[str, str],
    timeout: int,
) -> dict:
    """Claim pending alphamissense jobs and index them by (uniprot, AM-style variant).

    Each gene symbol is resolved to a UniProt accession -- using `overrides`
    first, then UniProt REST -- so the handler is gene-agnostic. Jobs whose gene
    can't be resolved (or that aren't missense) are returned in the 'unmappable'
    list so the caller can fail them without scanning the big TSV.
    """
    jobs = db.claim_pending_jobs(
        source=SOURCE,
        limit=limit if limit is not None else 1_000_000,
        gene_filter=gene_filter,
    )
    if not jobs:
        return {"by_key": {}, "unmappable": []}

    ids = [j["id"] for j in jobs]
    placeholders = ",".join("?" * len(ids))
    cur = db.conn.execute(
        f"""
        SELECT j.id AS job_id, j.variant_id AS variant_id,
               c.gene_symbol AS gene_symbol,
               c.aa_ref AS aa_ref, c.aa_pos AS aa_pos, c.aa_alt AS aa_alt,
               c.consequence AS consequence
        FROM annotation_jobs j
        JOIN variant_consequences c ON c.variant_id = j.variant_id
                                       AND c.source = 'enumerated'
        WHERE j.id IN ({placeholders})
        """,
        ids,
    )

    # Resolve each distinct gene once; cache accession-or-failure-reason.
    resolution_cache: dict[str, tuple[Optional[str], Optional[str]]] = {}

    def _resolve(gene: str) -> tuple[Optional[str], Optional[str]]:
        if gene in resolution_cache:
            return resolution_cache[gene]
        if gene in overrides:
            result: tuple[Optional[str], Optional[str]] = (overrides[gene], None)
        else:
            try:
                result = (resolve_uniprot_accession(gene, timeout=timeout), None)
            except (UniProtError, requests.RequestException) as e:
                result = (None, f"could not resolve UniProt accession for gene {gene!r}: {e}")
        resolution_cache[gene] = result
        return result

    by_key: dict[tuple, list[dict]] = {}
    unmappable: list[dict] = []
    for r in cur.fetchall():
        d = dict(r)
        consequence = d.get("consequence")
        # AM only covers missense.
        if consequence != "missense_variant":
            d["_reason"] = f"alphamissense covers missense only (consequence={consequence!r})"
            unmappable.append(d)
            continue
        gene = (d.get("gene_symbol") or "").upper()
        uniprot, reason = _resolve(gene)
        if not uniprot:
            d["_reason"] = reason or f"no UniProt accession for gene {gene!r}"
            unmappable.append(d)
            continue
        if not (d.get("aa_ref") and d.get("aa_pos") and d.get("aa_alt")):
            d["_reason"] = "consequence missing aa_ref / aa_pos / aa_alt"
            unmappable.append(d)
            continue
        key = (uniprot, f"{d['aa_ref']}{d['aa_pos']}{d['aa_alt']}")
        by_key.setdefault(key, []).append(d)

    return {"by_key": by_key, "unmappable": unmappable}


def _persist_match(db, *, job: dict, score: float, category: Optional[str], am_class: str) -> None:
    db.upsert_pathogenicity(
        job["variant_id"],
        "alphamissense",
        score=score,
        category=category,
        source=SOURCE,
    )
    db.mark_job_done(job["job_id"])


def _mark_unmappable(db, unmappable: list[dict]) -> int:
    for j in unmappable:
        db.mark_job_failed(j["job_id"], j.get("_reason", "unmappable"))
    return len(unmappable)
