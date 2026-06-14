"""AlphaFold DB residue-confidence annotator.

The AlphaFold DB API returns one or more model fragments for a UniProt accession.
Each model links to a confidence JSON with per-residue pLDDT scores. This handler
stores the pLDDT at the affected amino-acid / truncation position in
`annotations_structure`.
"""

from __future__ import annotations

from typing import Optional

import requests

from ..uniprot import UniProtError, resolve_uniprot_accession


SOURCE = "alphafold"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
DEFAULT_TIMEOUT = 30
DEFAULT_RATE_LIMIT_SEC = 0.0


class HandlerError(Exception):
    pass


def fetch_predictions(uniprot_accession: str, *, timeout: int = DEFAULT_TIMEOUT) -> list[dict]:
    """Return AlphaFold DB prediction metadata records for one UniProt accession."""
    resp = requests.get(f"{ALPHAFOLD_API}/{uniprot_accession}", timeout=timeout)
    if resp.status_code == 404:
        return []
    resp.raise_for_status()
    data = resp.json()
    if not isinstance(data, list):
        raise HandlerError(f"Unexpected AlphaFold response for {uniprot_accession!r}")
    return data


def fetch_confidence(plddt_doc_url: str, *, timeout: int = DEFAULT_TIMEOUT) -> dict:
    """Fetch an AlphaFold confidence JSON document."""
    resp = requests.get(plddt_doc_url, timeout=timeout)
    resp.raise_for_status()
    data = resp.json()
    if not isinstance(data, dict):
        raise HandlerError(f"Unexpected AlphaFold confidence JSON from {plddt_doc_url}")
    return data


def run_batch(
    db,
    *,
    limit: Optional[int] = None,
    gene_filter: Optional[str] = None,
    timeout: int = DEFAULT_TIMEOUT,
    uniprot_overrides: Optional[dict[str, str]] = None,
) -> dict:
    """Drain pending AlphaFold jobs.

    Returns: {"claimed": int, "annotated": int, "failed": int}
    """
    jobs = db.claim_pending_jobs(
        source=SOURCE,
        limit=limit if limit is not None else 1_000_000,
        gene_filter=gene_filter,
    )
    if not jobs:
        return {"claimed": 0, "annotated": 0, "failed": 0}

    variant_context = _load_variant_context(db, jobs)
    predictions_by_uniprot: dict[str, list[dict]] = {}
    confidence_by_url: dict[str, dict] = {}

    annotated = 0
    failed = 0

    for job in jobs:
        ctx = variant_context.get(job["variant_id"])
        if not ctx:
            db.mark_job_failed(job["id"], "no consequence row with gene_symbol and aa_pos")
            failed += 1
            continue

        gene = ctx["gene_symbol"]
        aa_pos = int(ctx["aa_pos"])
        try:
            uniprot = resolve_uniprot_accession(gene, timeout=timeout, extra=uniprot_overrides)
            if uniprot not in predictions_by_uniprot:
                predictions_by_uniprot[uniprot] = fetch_predictions(uniprot, timeout=timeout)
            prediction = _pick_prediction(predictions_by_uniprot[uniprot], aa_pos)
            if not prediction:
                raise HandlerError(f"No AlphaFold model covers {uniprot}:{aa_pos}")
            plddt_url = prediction.get("plddtDocUrl")
            if not plddt_url:
                raise HandlerError(f"AlphaFold model {prediction.get('entryId') or uniprot} lacks plddtDocUrl")
            if plddt_url not in confidence_by_url:
                confidence_by_url[plddt_url] = fetch_confidence(plddt_url, timeout=timeout)

            plddt = _plddt_at_residue(confidence_by_url[plddt_url], aa_pos, prediction)
            feature_version = _feature_version(prediction)
            db.upsert_structure(
                job["variant_id"],
                "alphafold_plddt",
                feature_version=feature_version,
                protein_accession=uniprot,
                residue_number=aa_pos,
                score=plddt,
                category=_plddt_category(plddt),
                source=SOURCE,
            )
            db.mark_job_done(job["id"])
            annotated += 1
        except (requests.RequestException, UniProtError, HandlerError, ValueError) as e:
            db.mark_job_failed(job["id"], str(e))
            failed += 1

    return {"claimed": len(jobs), "annotated": annotated, "failed": failed}


def _load_variant_context(db, jobs: list[dict]) -> dict[int, dict]:
    """Choose one protein-position context per variant, preferring enumerated rows."""
    ids = [j["id"] for j in jobs]
    placeholders = ",".join("?" * len(ids))
    cur = db.conn.execute(
        f"""
        SELECT j.variant_id AS variant_id,
               c.gene_symbol AS gene_symbol,
               c.aa_pos AS aa_pos,
               c.source AS source,
               c.consequence AS consequence
        FROM annotation_jobs j
        JOIN variant_consequences c ON c.variant_id = j.variant_id
        WHERE j.id IN ({placeholders})
          AND c.gene_symbol IS NOT NULL
          AND c.aa_pos IS NOT NULL
        ORDER BY
          CASE c.source
            WHEN 'enumerated' THEN 0
            WHEN 'vep' THEN 1
            WHEN 'annovar' THEN 2
            ELSE 3
          END
        """,
        ids,
    )
    out: dict[int, dict] = {}
    for row in cur.fetchall():
        d = dict(row)
        out.setdefault(d["variant_id"], d)
    return out


def _pick_prediction(predictions: list[dict], residue_number: int) -> Optional[dict]:
    """Pick the latest AlphaFold model fragment that covers a residue."""
    covering = []
    for p in predictions:
        start = p.get("sequenceStart") or p.get("uniprotStart") or 1
        end = p.get("sequenceEnd") or p.get("uniprotEnd") or 0
        try:
            if int(start) <= residue_number <= int(end):
                covering.append(p)
        except (TypeError, ValueError):
            continue
    if not covering:
        return None
    return max(covering, key=lambda p: int(p.get("latestVersion") or 0))


def _plddt_at_residue(confidence: dict, residue_number: int, prediction: dict) -> float:
    residues = confidence.get("residueNumber")
    scores = confidence.get("confidenceScore")
    if not isinstance(residues, list) or not isinstance(scores, list) or len(residues) != len(scores):
        raise HandlerError("AlphaFold confidence JSON lacks residueNumber/confidenceScore arrays")

    for residue, score in zip(residues, scores):
        try:
            if int(residue) == residue_number:
                return float(score)
        except (TypeError, ValueError):
            continue

    # Some fragment files number confidence arrays from 1 even when the model
    # covers a later UniProt interval. Fall back to the fragment offset.
    start = int(prediction.get("sequenceStart") or prediction.get("uniprotStart") or 1)
    offset = residue_number - start
    if 0 <= offset < len(scores):
        return float(scores[offset])
    raise HandlerError(f"No pLDDT score found for residue {residue_number}")


def _feature_version(prediction: dict) -> str:
    version = prediction.get("latestVersion")
    return f"AFDB_v{version}" if version is not None else ""


def _plddt_category(score: float) -> str:
    if score >= 90:
        return "very_high_confidence"
    if score >= 70:
        return "confident"
    if score >= 50:
        return "low_confidence"
    return "very_low_confidence"
