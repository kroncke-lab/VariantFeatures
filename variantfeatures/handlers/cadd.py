"""CADD REST API handler for normalized variants."""

from __future__ import annotations

from typing import Optional

from ..fetchers.cadd import CADD_VERSION, fetch_cadd_single


SOURCE = "cadd"
SOURCE_LABEL = "cadd_api"
DEFAULT_RATE_LIMIT_SEC = 0.2


class HandlerError(Exception):
    pass


def handle(db, variant_id: int, payload: Optional[str] = None) -> None:
    """Fetch CADD PHRED/raw for one canonical variant and persist it."""
    cur = db.conn.execute(
        "SELECT chromosome, position, ref, alt FROM variants WHERE id = ?",
        [variant_id],
    )
    row = cur.fetchone()
    if row is None:
        raise HandlerError(f"variant_id {variant_id} not found in variants table")

    result = fetch_cadd_single(
        row["chromosome"],
        int(row["position"]),
        ref=row["ref"],
        alt=row["alt"],
        version=CADD_VERSION,
    )
    if result is None:
        return

    if result.get("cadd_phred") is not None:
        db.upsert_pathogenicity(
            variant_id,
            "cadd_phred",
            predictor_version=CADD_VERSION,
            score=result["cadd_phred"],
            source=SOURCE_LABEL,
        )
    if result.get("cadd_raw") is not None:
        db.upsert_pathogenicity(
            variant_id,
            "cadd_raw",
            predictor_version=CADD_VERSION,
            score=result["cadd_raw"],
            source=SOURCE_LABEL,
        )
